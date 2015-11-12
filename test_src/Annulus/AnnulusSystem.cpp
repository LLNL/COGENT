#include "AnnulusSystem.H"
#include "AnnulusSystemF_F.H"
#include "CONSTANTS.H"

#include <fstream>
#include <sstream>
#if 1  // warning, OS dependencies, will not work on all platforms
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "CH_HDF5.H"
#include "HDF5Portable.H"
#include "BoxIterator.H"
#include "LoadBalance.H"
#include "MagBlockCoordSys.H"
#include "MillerBlockCoordSys.H"
#include "MillerCoordSys.H"
#include "newMappedGridIO.H"
#include "FourthOrderUtil.H"
#include "Directions.H"

#define NORADDECOMP
#undef NORADDECOMP


#include "NamespaceHeader.H"


AnnulusSystem::AnnulusSystem( ParmParse& a_pp )
   : m_verbosity(0),
     m_hdf_density(true),
     m_hdf_pressure(false),
     m_hdf_vpartheta(false),
     m_hdf_frtheta(false),
     m_hdf_vparmu(false),
     m_ghostVect(4*IntVect::Unit),
     m_ppgksys("gksystem"),
     m_enforce_stage_positivity(false),
     m_enforce_step_positivity(false)
{
   // Read the input
   parseParameters( m_ppgksys );

   /*
     Create the configuration space (magnetic geometry) domains,
     data layout and coordinate system.
   */

   createConfigurationSpace();

   m_operator = new Advect( ParmParse(Advect::pp_name) );

   m_solution.define( m_mag_geom->grids(), 1, IntVect::Zero);

   ParmParse pp_ibc("IBC");

   m_ibc = new LocSpatialIBC( pp_ibc, *m_mag_geom, 0);
}



void
AnnulusSystem::initialize(const int cur_step)
{
   if (procID()==0) cout << "Entering initialize" << endl;

   // Set the initial conditions if this is the first step
   if (cur_step == 0) {
      setInitialConditions( m_solution );
   }

   // Set the boundary data.  m_solution is not modified by this call (since
   // it now contains the initial condition or the solution obtain from a restart),
   // but simply provides a template for the boundary condition object.
   setBoundaryData( m_solution );

   // Make the solution prototype to hold old mapped and real data.
   m_solution_mapped.define(m_solution.disjointBoxLayout(), m_solution.nComp(), m_ghostVect);
   m_solution_mapped_old.define(m_solution.disjointBoxLayout(), m_solution.nComp(), m_ghostVect);

   DataIterator dit = m_solution.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      m_solution_mapped[dit].copy(m_solution[dit]);
      m_solution_mapped_old[dit].copy(m_solution[dit]);
   }

   divideJ(m_solution, m_solution);

   // Create the integrator
   m_integrator = new RK4Integrator<LevelData<FArrayBox> >(*this, m_solution_mapped);

   MPI_Barrier(MPI_COMM_WORLD);
   if (procID()==0) cout << "Exiting initialize" << endl;
}



AnnulusSystem::~AnnulusSystem()
{
   delete m_mag_geom;
   delete m_mag_geom_coords;
   delete m_operator;
   delete m_integrator;
}



void AnnulusSystem::createConfigurationSpace()
{
   string prefix = string("gksystem.magnetic_geometry_mapping.")
     + string(MillerBlockCoordSys::pp_name);

   ParmParse pp( prefix.c_str() );

   // Construct the multiblock magnetic coordinate system
   m_mag_geom_coords
     = new MillerCoordSys(pp,
                          m_num_cells,
                          m_is_periodic,
                          m_configuration_decomposition);

   DisjointBoxLayout grids;
   createDisjointBoxLayout(grids);

   // Construct the magnetic geometry

   int ghosts = 0;
   for (int n=0; n<CFG_DIM; ++n) {
      if (m_ghostVect[n] > ghosts) ghosts = m_ghostVect[n];
   }

   m_mag_geom = new MultiBlockLevelGeom(m_mag_geom_coords, grids, ghosts);
}



void AnnulusSystem::createDisjointBoxLayout(DisjointBoxLayout& grids)
{
   Vector<Box> boxes;
   const MagBlockCoordSys* mag_block_coords
     = (MagBlockCoordSys *)m_mag_geom_coords->getCoordSys(0);

   const ProblemDomain& domain = mag_block_coords->domain();

   const IntVect decomp = m_mag_geom_coords->getDecomposition(0);

   const Box& domain_box = domain.domainBox();

   /*
     Chop up the configuration space domain box over the number of processors specified
     for this block.  At this point, we insist that the box decomposes uniformly, or an
     error is thrown.
   */
   int nproc = 1;
   for (int dir=0; dir<CFG_DIM; ++dir) {
     nproc *= decomp[dir];
   }

   int n_loc[CFG_DIM];
   for (int dir=0; dir<CFG_DIM; ++dir) {
     n_loc[dir] = 0;
   }

#if CFG_DIM==2

   for (int dir=0; dir<CFG_DIM; ++dir) {
     int decomp_dir = decomp[dir];
     if (domain_box.size(dir)%decomp_dir != 0) {
       stringstream msg("Decomposition in configuration direction ", ios_base::out|ios_base::ate);
       msg << dir << " does not evenly divide domain dimension";
       MayDay::Error( msg.str().c_str() );
     }
     else {
       n_loc[dir] = domain_box.size(dir) / decomp_dir;
     }
   }

   if (n_loc[RADIAL_DIR] != 0 && n_loc[POLOIDAL_DIR] != 0) {

     IntVect box_size(n_loc[RADIAL_DIR],n_loc[POLOIDAL_DIR]);
     Box patch(domain_box.smallEnd(), domain_box.smallEnd() + box_size-1);
     Box skeleton(IntVect::Zero, IntVect(domain_box.size(0)/n_loc[RADIAL_DIR]-1,domain_box.size(1)/n_loc[POLOIDAL_DIR]-1));
     BoxIterator bit(skeleton);
     for (bit.begin();bit.ok();++bit) {
       Box thisBox = patch + bit()*box_size;
       boxes.push_back(thisBox);
     }
   }
   else {
     MayDay::Error( "Configuration domain box cannot be load balanced" );
   }
#else
   MayDay::Error( "AnnulusSystem::createConfigurationSpace has only been implemented for 2D" );
#endif

   // Make the layout.  This is where boxes are assigned to processes.
   Vector<int> procMap;
   LoadBalance( procMap, boxes );

   m_domain = domain;

   grids.define( boxes, procMap, m_domain );
   grids.close();

#ifdef CH_MPI
   if (procID()==0) {
#endif
     if (m_verbosity>0) {
       for (int n=0; n<boxes.size(); n++) {
         const Box& local_box = boxes[n];
           cout << "   Configuration space box " << local_box << " is assigned to process " << n << endl;
       }
     }
#ifdef CH_MPI
   }
#endif
}



void AnnulusSystem::defineRHSData( LevelData<FArrayBox>&       a_rhs,
                                   const LevelData<FArrayBox>& a_prototype )
{
   if ( a_rhs.isDefined() ) {
      MayDay::Error("AnnulusSystem::defineRHSData(): a_rhs is already defined");
   }

   a_rhs.define( a_prototype.disjointBoxLayout(),
                 a_prototype.nComp(),
                 IntVect::Zero );

   DataIterator dit = a_prototype.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_rhs[dit].copy(a_prototype[dit]);
   }
}



void AnnulusSystem::defineSolnData( LevelData<FArrayBox>&       a_soln,
                                    const LevelData<FArrayBox>& a_prototype )
{
   if ( a_soln.isDefined() ) {
      MayDay::Error("AnnulusSystem::defineSolnData(): a_soln is already defined");
   }

   a_soln.define( a_prototype.disjointBoxLayout(),
                  a_prototype.nComp(),
                  IntVect::Zero );

   DataIterator dit = a_prototype.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_soln[dit].copy(a_prototype[dit]);
   }
}



bool AnnulusSystem::validSolnData( const LevelData<FArrayBox>& a_soln,
                                   const LevelData<FArrayBox>& a_protoSoln )
{
   return (a_protoSoln.disjointBoxLayout() == a_soln.disjointBoxLayout()) &&
          (a_protoSoln.nComp()             == a_soln.nComp())             &&
          (a_protoSoln.ghostVect()         == a_soln.ghostVect());
}



bool AnnulusSystem::validRHSData( const LevelData<FArrayBox>& a_rhs,
                                  const LevelData<FArrayBox>& a_protoRHS )
{
   return (a_protoRHS.disjointBoxLayout() == a_rhs.disjointBoxLayout()) &&
          (a_protoRHS.nComp()             == a_rhs.nComp());
}



void AnnulusSystem::zeroSolnData( LevelData<FArrayBox>& a_soln )
{
   DataIterator dit = a_soln.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_soln[dit].setVal(0.);
   }
}



void AnnulusSystem::addSolnData( LevelData<FArrayBox>&       a_soln,
                                 const LevelData<FArrayBox>& a_increment,
                                 const Real                  a_scale )
{
   DataIterator dit = a_soln.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_soln[dit].plus(a_increment[dit], a_scale);
   }
}



void AnnulusSystem::copySolnData( LevelData<FArrayBox>&       a_dstSoln,
                                  const LevelData<FArrayBox>& a_srcSoln )
{
   CH_assert( validSolnData(a_dstSoln, a_srcSoln) );

   DataIterator dit = a_srcSoln.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_dstSoln[dit].copy(a_srcSoln[dit]);
   }
}



void AnnulusSystem::evalRHS( LevelData<FArrayBox>&       a_rhs,
                             const LevelData<FArrayBox>& a_soln,
                             const int                   a_step_number,
                             const Real                  a_time,
                             const int                   a_stage)
{
   LevelData<FArrayBox> grown_soln(a_soln.disjointBoxLayout(),
                                   a_soln.nComp(),
                                   m_ghostVect);
   DataIterator dit = a_soln.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      grown_soln[dit].copy(a_soln[dit]);
   }
   grown_soln.exchange();
   m_ibc->ghostCellBC(grown_soln, *m_mag_geom, 0., a_time);

   divideJ(grown_soln, grown_soln);

   m_operator->evalRHS( *m_mag_geom, a_rhs, grown_soln, a_time );
}



Real AnnulusSystem::stableDt(const int a_step_number)
{
   return m_operator->computeDt(*m_mag_geom, m_solution_mapped);
}



void AnnulusSystem::enforcePositivity( LevelData<FArrayBox>& a_soln )
{
   m_positivity_post_processor.enforce( a_soln, 1.0 );
}



Real AnnulusSystem::advance( const Real a_cur_time,
                        const Real a_dt,
                        const int  a_step_number)
{
   double new_time = m_integrator->advance( m_solution_mapped,
                                            m_solution_mapped_old,
                                            a_step_number,
                                            a_cur_time,
                                            a_dt );

   if (m_enforce_step_positivity) {
      enforcePositivity( m_solution_mapped );
   }

   copySolnData( m_solution_mapped_old, m_solution_mapped );

   divideJ(m_solution_mapped, m_solution);

   return new_time;
}



void AnnulusSystem::setInitialConditions( LevelData<FArrayBox>& a_soln )
{
   m_ibc->initialize(a_soln, *m_mag_geom, 0., 0.);
}



void AnnulusSystem::setBoundaryData( LevelData<FArrayBox>& a_soln )
{
   m_ibc->setBoundaryData(a_soln, *m_mag_geom, 0., 0.);
}



void AnnulusSystem::multJ( const LevelData<FArrayBox>& a_soln_physical,
                           LevelData<FArrayBox>&       a_soln_mapped)
{
   // Need to grow the soln factor since we need transverse gradients for the fourth-order product
   const DisjointBoxLayout& grids = a_soln_physical.disjointBoxLayout();
   IntVect grown_vect = a_soln_mapped.ghostVect() + IntVect::Unit;
   LevelData<FArrayBox> grown_soln(grids, 1, grown_vect);

   LevelData<FArrayBox> J(grids, 1, grown_vect);
   m_mag_geom_coords->getJ(J);

   DataIterator dit = grown_soln.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox& this_grown_soln = grown_soln[dit];
      Box this_box = a_soln_mapped[dit].box();

      this_grown_soln.copy(a_soln_physical[dit]);

      for (int dir=0; dir<SpaceDim; ++dir) {
         secondOrderTransExtrap(this_grown_soln, dir, this_box);
      }
   }
   grown_soln.exchange();

   // Compute fourth-order product
   fourthOrderCellProd(a_soln_mapped, grown_soln, J);
}



void AnnulusSystem::divideJ( const LevelData<FArrayBox>& a_soln_mapped,
                             LevelData<FArrayBox>&       a_soln_physical)
{
   CH_assert( a_soln_mapped.ghostVect() >= a_soln_physical.ghostVect() );

   const DisjointBoxLayout& grids = a_soln_mapped.disjointBoxLayout();
   IntVect grown_vect = a_soln_physical.ghostVect() + IntVect::Unit;
   LevelData<FArrayBox> grown_soln(grids, 1, grown_vect);

   LevelData<FArrayBox> J(grids, 1, grown_vect);
   m_mag_geom_coords->getJ(J);

   DataIterator dit = a_soln_physical.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
     FArrayBox& this_grown_soln = grown_soln[dit];
     Box this_box = a_soln_physical[dit].box();

     this_grown_soln.copy(a_soln_mapped[dit]);

     for (int dir=0; dir<SpaceDim; ++dir) {
       secondOrderTransExtrap(this_grown_soln, dir, this_box);
     }
   }
   grown_soln.exchange();

   for (dit.begin(); dit.ok(); ++dit) {
     cellFGToCellF(a_soln_physical[dit], grown_soln[dit], J[dit], a_soln_physical[dit].box());
   }
}



void AnnulusSystem::writePlotFile(const char *prefix, const int cur_step )
{
   if (procID()==0) cout << "Writing plotfile" << endl;
   char iter_str[100], dir_prefix[100];

   stringstream filename;

   sprintf(dir_prefix, "%s", prefix);

#if 1  // warning, OS dependencies, will not work on all platforms
   // place density plots in their own directory;
   sprintf(iter_str, "%s_density_plots", prefix);
#ifdef CH_MPI
   if (procID() == 0) {
#endif
      mkdir(iter_str, 0777); // only works the first time, subsequent failure is normal and expected
#ifdef CH_MPI
      }
#endif
   sprintf(dir_prefix, "%s//%s", iter_str, prefix);
#endif

   stringstream filename1;
   sprintf(iter_str, "%s.density%04d.%d.%s.", dir_prefix, cur_step, 0, "hydrogen");
   filename1 << iter_str;

   WriteMappedUGHDF5(filename1.str().c_str(), m_solution.disjointBoxLayout(), m_solution, *m_mag_geom_coords, m_domain.domainBox());
   if (procID()==0) cout << "Done writing plotfile" << endl;
}



void AnnulusSystem::writeCheckpointFile( HDF5Handle&  handle,
                                    const int    cur_step,
                                    const double cur_time,
                                    const double cur_dt )
{
}



void AnnulusSystem::readCheckpointFile( HDF5Handle& handle, int* cur_step, double* cur_time, double* cur_dt )
{
}


inline void AnnulusSystem::printParameters()
{
   if (procID() == 0 && m_verbosity) {
      cout << "enforce_positivity = " << (m_enforce_step_positivity||m_enforce_stage_positivity) << endl;
      std::string ptype("stage");
      if (m_enforce_step_positivity)
        ptype = "step";
      cout << "enforce_positivity_type = " << ptype << endl;
   }
}



void AnnulusSystem::writeFieldHistory(int cur_step, double cur_time, bool startup_flag)
{
}




void AnnulusSystem::parseParameters( ParmParse& a_ppgksys)
{
   // This determines the amount of diagnositic output generated
   a_ppgksys.query( "verbosity", m_verbosity );
   CH_assert( m_verbosity >= 0 );

   // Set the grid size
   m_num_cells.resize( CFG_DIM );
   for (int i=0; i<CFG_DIM; ++i) m_num_cells[i] = 0;
   a_ppgksys.getarr( "num_cells", m_num_cells, 0, CFG_DIM );
   for (int i=0; i<CFG_DIM; ++i) CH_assert( m_num_cells[i]>0 );

   // Determine which spatial directions are periodic
   m_is_periodic.resize(CFG_DIM);
   vector<int> isPeriodic( CFG_DIM ); // why should I have to do this?
   a_ppgksys.getarr( "is_periodic", isPeriodic, 0, CFG_DIM );
   for (int dim=0; dim<SpaceDim; dim++)  {
     m_is_periodic[dim] = (isPeriodic[dim] == 1);
   }

   // Get the domain decomposition parameters
   if (a_ppgksys.contains("configuration_decomp")) {
     m_configuration_decomposition.resize( CFG_DIM );
     for (int i=0; i<CFG_DIM; ++i) m_configuration_decomposition[i] = 0;
     a_ppgksys.getarr( "configuration_decomp", m_configuration_decomposition, 0, CFG_DIM );
     for (int i=0; i<CFG_DIM; ++i) CH_assert( m_configuration_decomposition[i]>0 );
   }

   // Should we make hdf files for density?
   a_ppgksys.query("hdf_density",m_hdf_density);
   // Should we make hdf files for f versus vparallel and poloidal angle?
   a_ppgksys.query("hdf_vpartheta",m_hdf_vpartheta);
   // Should we make hdf files for f versus radius and poloidal angle at a specified vpar, mu?
   a_ppgksys.query("hdf_frtheta",m_hdf_frtheta);
   // Should we make hdf files for pressure?
   a_ppgksys.query("hdf_pressure",m_hdf_pressure);
   // Should we make hdf files for vparallel-mu?
   a_ppgksys.query("hdf_vparmu",m_hdf_vparmu);
   // At what fixed phase space indices should I plot?  (Indices plotted against in a given plot
   //   are ignored.  Specify in 5D; toroidal index ignored in 4D and set to zero in arguments
   //   of hdf write methods.
#if 0
   m_fixed_plotindices.resize( 5 );
   for (int i=0; i<5; ++i) m_fixed_plotindices[i]=0;
   a_ppgksys.queryarr("fixed_plot_indices",m_fixed_plotindices,0,5);
   for (int block=0; block<a_config_blocks.size(); ++block) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         CH_assert( a_config_blocks[block].loMappedIndex()[dir] <= m_fixed_plotindices[dir] );
         CH_assert( a_config_blocks[block].hiMappedIndex()[dir] >  m_fixed_plotindices[dir] );
      }
   }
#endif

   bool enforce_positivity(false);
   if (a_ppgksys.contains("enforce_positivity")) {
     a_ppgksys.get("enforce_positivity", enforce_positivity);
   }

   if (enforce_positivity) {

     std::string ptype("stage");
     m_enforce_stage_positivity = true;
     if (a_ppgksys.contains("enforce_positivity_type")) {
       a_ppgksys.get("enforce_positivity_type", ptype);
       if (ptype=="step") {
         m_enforce_step_positivity = true;
         m_enforce_stage_positivity = false;
       }
       else if (ptype!="stage") {
         MayDay::Error("Invalid positivity enforcement type");
       }
     }

     int n_iter(5);
     if (a_ppgksys.contains("max_positivity_iter")) {
       a_ppgksys.get("max_positivity_iter", n_iter);
     }

     bool verbose(false);
     if (a_ppgksys.contains("positivity_verbose_output")) {
       a_ppgksys.get("positivity_verbose_output", verbose);
     }

     int width(2);
     if (m_enforce_step_positivity) width++;
     IntVect halo( width*IntVect::Unit );
     m_positivity_post_processor.define( halo, n_iter, verbose );
   }

   if (m_verbosity) {
      printParameters();
   }
}



void AnnulusSystem::printDiagnostics()
{
}



void AnnulusSystem::postStageAdvance( LevelData<FArrayBox>& a_soln )
{
   if (m_enforce_stage_positivity) {
      enforcePositivity( a_soln );
   }
}



#include "NamespaceFooter.H"
