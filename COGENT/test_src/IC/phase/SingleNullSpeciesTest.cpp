#include "SingleNullSpeciesTest.H"

#include "CONSTANTS.H"

#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

#if 1  // warning, OS dependencies, will not work on all platforms
#include <sys/stat.h>
#endif

#include "Box.H"
#include "BoxIterator.H"
#include "HDF5Portable.H"
#include "LoadBalance.H"
#include "MayDay.H"
#include "newMappedGridIO.H"
#include "Vector.H"

#include "Directions.H"
#include "KineticSpeciesICFactory.H"
// these will go away
#include "IBCInterface.H";
#include "CLSInterface.H";
#include "VirtualCLS.H";

#undef CH_SPACEDIM
#define CH_SPACEDIM VEL_DIM
#include "Box.H"
#include "ProblemDomain.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "Box.H"
#include "BoxIterator.H"
#include "LoadBalance.H"
#include "MagCoordSys.H"
#include "MagBlockCoordSys.H"
#include "SingleNullCoordSys.H"
#include "SingleNullBlockCoordSys.H"
#include "newMappedGridIO.H"
#include "MultiBlockLevelExchangeAverage.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#define NORADDECOMP
#undef NORADDECOMP

#include "NamespaceHeader.H"


inline
bool notZeroVector( const std::vector<int>& a_x )
{
   const std::vector<int> z( a_x.size(), 0 );
   return std::lexicographical_compare( z.begin(), z.end(), a_x.begin(), a_x.end() );
}


SingleNullSpeciesTest::SingleNullSpeciesTest( ParmParse& a_pp )
{
   parseParameters( a_pp );

   createConfigurationSpace();

   createVelocitySpace();

   createPhaseSpace( a_pp );

   m_solution.resize( NUM_IC );
   for (int n(0); n<NUM_IC; n++) {
      std::string name( "hydrogen" );
      Real mass(1.0);
      Real charge(1.0);
      PhaseGeom* species_geom = new PhaseGeom( *m_phase_geom, mass, charge );
      IBCInterface* ics = NULL;
      CLSInterface* cls = new VirtualCLS();
      m_solution[n] = new KineticSpecies( name,
                                          mass,
                                          charge,
                                          *species_geom,
                                          *ics,
                                          *ics);
      m_solution[n]->distributionFunction().define( m_phase_geom->grids(),
                                                    1,
                                                    IntVect::Zero );
   }

   m_ic_name.resize( NUM_IC );
   // To add a new IC to test, add a name to m_ic_name that corresponds to
   // name in ICFactory.  Also, add a new section "IC.name..." in the input
   // file.
   m_ic_name[ZERO]      = "zero";
   m_ic_name[CONSTANT]  = "constant";
   m_ic_name[LOCALIZED] = "localized";
   m_ic_name[MAXWELLIAN] = "maxwellian";

   KineticSpeciesICFactory ic_factory;
   m_ic.resize( NUM_IC );
   for (int n(0); n<NUM_IC; n++) {
      std::string prefix( "IC." + m_ic_name[n] );
      ParmParse ppic( prefix.c_str() );
      m_ic[n] = ic_factory.create( m_ic_name[n], ppic, 1 );
   }

   m_fixed_plotindices.resize( 5 );
   for (int i(0); i<5; i++) {
      m_fixed_plotindices[i] = 0;
   }
   m_fixed_plotindices[POLOIDAL_DIR] = 64;
//   m_fixed_plotindices[RADIAL_DIR] = m_num_cells[RADIAL_DIR] / 2;
}



void
SingleNullSpeciesTest::initialize()
{
   if (procID()==0) cout << "Entering initialize" << endl;

   setInitialConditions();
   MPI_Barrier(MPI_COMM_WORLD);

   if (procID()==0) cout << "Exiting initialize" << endl;
}



SingleNullSpeciesTest::~SingleNullSpeciesTest()
{
   for (int n(0); n<NUM_IC; n++) {
      delete m_solution[n];
      delete m_ic[n];
   }
   delete m_mag_geom;
   delete m_mag_geom_coords;
   delete m_velocity_coords;
   delete m_phase_geom;
   delete m_phase_coords;
}



void SingleNullSpeciesTest::createConfigurationSpace()
{
   ParmParse pp_grid( "singlenull" );
   std::string prefix( "magnetic_geometry_mapping."
                       + CFG::SingleNullBlockCoordSys::pp_name );
   ParmParse pp_geom( prefix.c_str() );
   m_mag_geom_coords = new CFG::SingleNullCoordSys( pp_grid, pp_geom );

   CFG::DisjointBoxLayout grids;
   createConfigurationSpaceDisjointBoxLayout( grids );

   if (procID()==0) {
      std::cout << "Constructing magnetic geometry" << std::endl;
   }

   int space_order(4); // = 2*L with L as in Peter's notes
   int ghosts(0);
   m_mag_geom = new CFG::MultiBlockLevelGeom();
   m_mag_geom->define( m_mag_geom_coords, grids, ghosts, space_order );

   if (procID()==0) {
      std::cout << "Done constructing magnetic geometry"
                << std::endl << std::endl;
   }
}


void SingleNullSpeciesTest::createConfigurationSpaceDisjointBoxLayout(
   CFG::DisjointBoxLayout& grids )
{
   Vector<CFG::Box> boxes;
   for (int block(0); block<m_mag_geom_coords->numBlocks(); ++block) {

      const CFG::MagBlockCoordSys* mag_block_coords
         = (CFG::MagBlockCoordSys *)m_mag_geom_coords->getCoordSys( block );
      const CFG::ProblemDomain& domain( mag_block_coords->domain() );
      const CFG::IntVect decomp( m_mag_geom_coords->getDecomposition( block ) );
      const CFG::Box& domain_box( domain.domainBox() );

      /*
        Chop up the configuration space domain box over the number of processors specified
        for this block.  At this point, we insist that the box decomposes uniformly, or an
        error is thrown.
      */
      int nproc( decomp.product() );
      CFG::IntVect n_loc( D_DECL6(0,0,0,0,0,0) );

#if CFG_DIM==2

      for (int dir(0); dir<CFG_DIM; ++dir) {
         int decomp_dir( decomp[dir] );
         if (domain_box.size(dir)%decomp_dir!=0) {
            stringstream msg( "Decomposition in configuration direction ", ios_base::out|ios_base::ate );
            msg << dir << " does not evenly divide domain dimension";
            MayDay::Error( msg.str().c_str() );
         }
         else {
            n_loc[dir] = domain_box.size(dir) / decomp_dir;
         }
      }

      if (n_loc[RADIAL_DIR]!=0 && n_loc[POLOIDAL_DIR]!=0) {

         CFG::IntVect box_size( n_loc[RADIAL_DIR], n_loc[POLOIDAL_DIR] );
         CFG::Box patch( domain_box.smallEnd(), domain_box.smallEnd() + box_size-1 );
         CFG::Box skeleton( CFG::IntVect::Zero, CFG::IntVect( domain_box.size() / n_loc - 1 ) );
         for (CFG::BoxIterator bit( skeleton ); bit.ok();++bit) {
            CFG::Box thisBox( patch + bit() * box_size );
            boxes.push_back( thisBox );
         }
      }
      else {
         MayDay::Error( "Configuration domain box cannot be load balanced" );
      }
#else
      MayDay::Error( "SingleNullSpeciesTest::createConfigurationSpace has only been implemented for 2D" );
#endif
   }

   // Make the layout.  This is where boxes are assigned to processes.
   CFG::Vector<int> procMap;
   CFG::LoadBalance( procMap, boxes );
   CFG::ProblemDomain prob_domain;

   CFG::Box bounding_box;
   for (int n(0); n<boxes.size(); n++) {
      bounding_box = minBox( bounding_box, boxes[n] );
   }

   m_domain = CFG::ProblemDomain( bounding_box );
   grids.define( boxes, procMap, m_domain );
   grids.close();
}


void
SingleNullSpeciesTest::createVelocitySpace()
{
   const VEL::ProblemDomain& domain( getVelocityDomain() );

   VEL::DisjointBoxLayout grids;
   createVelocitySpaceDisjointBoxLayout( grids, domain );

   ParmParse pppsm( "phase_space_mapping" );

   Real v_parallel_max;
   pppsm.get( "v_parallel_max", v_parallel_max );
   // v_parallel runs from -v_parallel_max to + v_parallel_max
   Real dv_parallel( 2 * v_parallel_max / domain.size(0) );

   Real mu_max;
   pppsm.get( "mu_max", mu_max );
   Real dmu( mu_max / domain.size(1) );

   const VEL::RealVect dv( dv_parallel, dmu );

   std::string prefix( string("gksystem.") + string(VEL::VelCoordSys::pp_name) );
   ParmParse pp_vel( prefix.c_str() );

   m_velocity_coords = new VEL::VelCoordSys( pp_vel, grids, domain, dv );
}


void
SingleNullSpeciesTest::createVelocitySpaceDisjointBoxLayout(
   VEL::DisjointBoxLayout& grids,
   const VEL::ProblemDomain& domain )
{
   VEL::Box domain_box( domain.domainBox() );

   /*
     Chop up the velocity space domain box over the current number of processors.
     At this point, we insist that the box decomposes uniformly, or an error is thrown.
   */
   int nproc( numProc() );
   VEL::IntVect n_loc( VEL::IntVect::Zero );

   int nproc_vel(1);
   for (int dir(0); dir<VEL_DIM; ++dir) {
      nproc_vel *= m_velocity_decomposition[dir];
   }

   if (nproc_vel>nproc) {
      MayDay::Error( "Requested velocity space decomposition exceeds available processors" );
   }
   for (int dir(0); dir<VEL_DIM; ++dir) {
      if (domain_box.size(dir)%m_velocity_decomposition[dir] != 0) {
         stringstream msg( "Decomposition in velocity direction ",
                           ios_base::out|ios_base::ate );
         msg << dir << " does not evenly divide domain dimension";
         MayDay::Error( msg.str().c_str() );
      }
      else {
         n_loc[dir] = domain_box.size(dir) / m_velocity_decomposition[dir];
      }
   }

   Vector<VEL::Box> boxes;
   if ( (n_loc[0]!=0) && (n_loc[1]!=0) ) {
      VEL::Box patch( domain_box.smallEnd(),
                      domain_box.smallEnd() + n_loc - 1 );
      VEL::IntVect upper( domain_box.size() / n_loc - 1 );
      VEL::Box skeleton( VEL::IntVect::Zero, upper );
      for (VEL::BoxIterator bit( skeleton ); bit.ok(); ++bit) {
         VEL::Box thisBox( patch + bit() * n_loc );
         boxes.push_back( thisBox );
      }
   }
   else {
      MayDay::Error( "Velocity domain box cannot be load balanced" );
   }

   // Make the layout.  This is where boxes are assigned to processes.
   Vector<int> procMap;
   VEL::LoadBalance( procMap, boxes );
   grids.define( boxes, procMap, domain );
   grids.close();
}


void
SingleNullSpeciesTest::createPhaseSpace( ParmParse& a_pp )
{
  // Build a vector of ProblemDomains, one per configuration space block

   const VEL::ProblemDomain& vel_domain( getVelocityDomain() );

   Vector<ProblemDomain> domains;
   for (int block(0); block<m_mag_geom_coords->numBlocks(); ++block) {

      // Construct the phase space block domain from the configuration and velocity domains

      const CFG::MagBlockCoordSys* mag_block_coords
         = (CFG::MagBlockCoordSys *)m_mag_geom_coords->getCoordSys( block );

      const CFG::ProblemDomain& cfg_domain( mag_block_coords->domain() );
      const CFG::Box& cfg_box( cfg_domain.domainBox() );

      IntVect loVect, hiVect;
      bool isPeriodic[PDIM];

      for (int n(0); n<CFG_DIM; n++) {
         loVect[n] = cfg_box.smallEnd( n );
         hiVect[n] = cfg_box.bigEnd( n );
         isPeriodic[n] = cfg_domain.isPeriodic( n );
      }
      const VEL::Box& vel_box( vel_domain.domainBox() );
      for (int n(CFG_DIM); n<PDIM; n++) {
         loVect[n] = vel_box.smallEnd( n - CFG_DIM );
         hiVect[n] = vel_box.bigEnd( n - CFG_DIM );
         isPeriodic[n] = vel_domain.isPeriodic( n - CFG_DIM );
      }

      Box domain_box( loVect, hiVect );
      ProblemDomain domain( domain_box, isPeriodic );
      domains.push_back( domain );
   }

   ParmParse pp( "singlenull.decomp" );
   m_phase_coords
      = new SingleNullPhaseCoordSys( pp,
                                     *(CFG::SingleNullCoordSys*)m_mag_geom_coords,
                                     *m_velocity_coords,
                                     domains );

    DisjointBoxLayout phase_space_grid;
    createPhaseSpaceDisjointBoxLayout( phase_space_grid, domains );

    int ghosts(0);
    m_phase_geom = new PhaseGeom( a_pp,
                                  m_phase_coords,
                                  phase_space_grid,
                                  *m_mag_geom,
                                  *m_velocity_coords,
                                  ghosts,
                                  1 );
//                                m_units->larmorNumber() );

   MPI_Barrier( MPI_COMM_WORLD );

}



void
SingleNullSpeciesTest::createPhaseSpaceDisjointBoxLayout(
   DisjointBoxLayout&           a_grids,
   const Vector<ProblemDomain>& a_domains )
{
   Vector<Box> boxes;
   for (int block(0); block<a_domains.size(); ++block) {

      const Box& domain_box( a_domains[block].domainBox() );

      IntVect decomp( m_phase_coords->getDecomposition( block ) );

      /*
        Chop up the domain box in the configuration space directions over the current
        number of processors.
      */
      int nproc( decomp.product() );
      IntVect n_loc( IntVect::Zero );

      for (int dir(0); dir<PDIM; ++dir) {
         int decomp_dir( decomp[dir] );
         if (domain_box.size(dir)%decomp_dir != 0) {
            stringstream msg("Decomposition in phase direction ",
                             ios_base::out|ios_base::ate);
            msg << dir << " does not evenly divide domain dimension";
            MayDay::Error( msg.str().c_str() );
         }
         else {
            n_loc[dir] = domain_box.size(dir) / decomp_dir;
         }
      }

      if ( n_loc!=IntVect::Zero ) {
         for (int dir(0); dir<PDIM; ++dir) {
            if (n_loc[dir]<4) {
               MayDay::Error( "Phase space box is less than 4 cells wide" );
            }
         }

         Box patch( domain_box.smallEnd(),
                    domain_box.smallEnd() + n_loc - IntVect::Unit );
         IntVect upper( domain_box.size() / n_loc - 1 );
         Box skeleton( IntVect::Zero, upper );

         BoxIterator bit( skeleton );
         for (BoxIterator bit( skeleton ); bit.ok(); ++bit) {
            Box thisBox( patch + bit() * n_loc );
            boxes.push_back( thisBox );
         }
      }
      else {
         MayDay::Error( "Phase space domain box cannot be load balanced" );
      }
   }

   // Make the layout.  This is where boxes are assigned to processes.
   Vector<int> procMap;
   LoadBalance( procMap, boxes );

   Box bounding_box;
   for (int n(0); n<boxes.size(); n++) {
      bounding_box = minBox( bounding_box, boxes[n] );
   }

   ProblemDomain prob_domain( bounding_box );
   a_grids.define( boxes, procMap, prob_domain );
   a_grids.close();
}


void SingleNullSpeciesTest::setInitialConditions()
{
   for (int n(0); n<NUM_IC; n++) {
      m_ic[n]->initialize( *(m_solution[n]), 0.0 );
   }
}


void SingleNullSpeciesTest::writePlotFile( const std::string& prefix )
{
   if (procID()==0) cout << "Writing plotfile" << endl;

   for (int n(0); n<NUM_IC; n++) {

      // Get solution distribution function for the current species
      KineticSpecies& soln_species
         = dynamic_cast<KineticSpecies&>( *( m_solution[n] ) );
      LevelData<FArrayBox>& soln_dfn( soln_species.distributionFunction() );
      const PhaseGeom& geom( soln_species.phaseSpaceGeometry() );

      // Need to divide the soln, which is Bstarpar * f, by Bstarpar.
      LevelData<FArrayBox> dfn( soln_dfn.getBoxes(),
                                soln_dfn.nComp(),
                                soln_dfn.ghostVect() );
      for (DataIterator dit( dfn.dataIterator() ); dit.ok(); ++dit) {
         dfn[dit].copy( soln_dfn[dit] );
      }
      geom.divideBStarParallel( dfn );
      geom.divideJonValid( dfn ); // technically, lose order, but just for plotting

      // Plot f versus vparallel and theta at fixed psi and mu.
      std::string dir_prefix( prefix );
      std::string iter_str( prefix + "_vpar_poloidal_plots" );
      if (procID() == 0) {
         mkdir( iter_str.c_str(), 0777 );
      }
      dir_prefix = iter_str + "/" + prefix;

      std::string filename( dir_prefix + "." + m_ic_name[n] + ".0000.2d.hdf5" );

      geom.plotVParPoloidalData( filename.c_str(),
                                 m_fixed_plotindices[RADIAL_DIR],
#if (CFG_DIM>2)
                                 m_fixed_plotindices[TOROIDAL_DIR],
#else
                                 0,
#endif
                                 m_fixed_plotindices[MU_DIR],
                                 dfn );

      // Plot f versus psi and theta at fixed mu and vparallel.
      dir_prefix = prefix;
      iter_str = prefix + "_r_theta_plots";
      if (procID() == 0) {
         mkdir( iter_str.c_str(), 0777 );
      }
      dir_prefix = iter_str + "/" + prefix;

      filename = ( dir_prefix + "." + m_ic_name[n] + ".0000." );

      VEL::IntVect vspace_index( m_fixed_plotindices[VPARALLEL_DIR],
                                 m_fixed_plotindices[MU_DIR] );
      geom.plotAtVelocityIndex( filename.c_str(), vspace_index, dfn );


      // Plot f versus vparallel and mu at fixed configuration coordinate
      dir_prefix = prefix;
      iter_str = prefix + "_vpar_mu_plots";
      if (procID() == 0) {
         mkdir( iter_str.c_str(), 0777 );
      }
      dir_prefix = iter_str + "/" + prefix;

      filename = ( dir_prefix + "." + m_ic_name[n] + ".0000." );

      CFG::IntVect cspace_index( m_fixed_plotindices[RADIAL_DIR],
                                 m_fixed_plotindices[POLOIDAL_DIR] );
      geom.plotAtConfigurationIndex( filename.c_str(), cspace_index, soln_dfn);
   }

   if (procID()==0) cout << "Done writing plotfile" << endl;
}


void SingleNullSpeciesTest::parseParameters( ParmParse& a_pp )
{
   m_num_velocity_cells.resize( VEL_DIM );
   m_num_velocity_cells.assign( VEL_DIM, 0 );
   a_pp.getarr( "num_velocity_cells", m_num_velocity_cells, 0, m_num_velocity_cells.size() );
   CH_assert( notZeroVector( m_num_velocity_cells ) );

   if (a_pp.contains( "velocity_decomp" )) {
      m_velocity_decomposition.resize( VEL_DIM );
      m_velocity_decomposition.assign( VEL_DIM, 0 );
      a_pp.queryarr( "velocity_decomp",
                     m_velocity_decomposition,
                     0, m_velocity_decomposition.size() );
      CH_assert( notZeroVector( m_velocity_decomposition ) );
   }
}


VEL::ProblemDomain SingleNullSpeciesTest::getVelocityDomain() const
{
   /*
     For now, the number of cells in the vparallel direction must be even.
     This is due to the assumed convention used to compute the physical
     velocity directly from the global index and mesh size.  We intend to
     remove this assumption, and the associated limitation tested for here,
     in the future.
   */
   if ( m_num_velocity_cells[0]%2 != 0 ) {
      MayDay::Warning( "vparallel dimension must be even" );
   }

   int half_vp_dim( m_num_velocity_cells[0] / 2 );
   VEL::IntVect lo( -half_vp_dim, 0 );
   VEL::IntVect hi( half_vp_dim-1, m_num_velocity_cells[1]-1 );

   bool isPeriodic[VEL_DIM];
   for (int n(0); n<VEL_DIM; n++) {
      isPeriodic[n] = false;
   }

   return VEL::ProblemDomain( VEL::Box( lo, hi ), isPeriodic );
}


#include "NamespaceFooter.H"
