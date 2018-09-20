#include "AnnulusTest.H"

#include "CONSTANTS.H"

#include <iostream>
#include <sstream>

#if 1  // warning, OS dependencies, will not work on all platforms
#include <sys/stat.h>
#endif

#include "Box.H"
#include "BoxIterator.H"
#include "IntVect.H"
#include "HDF5Portable.H"
#include "LoadBalance.H"
#include "MayDay.H"
#include "newMappedGridIO.H"
#include "Vector.H"

#include "Directions.H"
#include "MagBlockCoordSys.H"
#include "MillerBlockCoordSys.H"
#include "MillerCoordSys.H"
#include "ICFactory.H"

#define NORADDECOMP
#undef NORADDECOMP

#include "NamespaceHeader.H"


AnnulusTest::AnnulusTest( ParmParse& a_pp )
{
   parseParameters( a_pp );

   createConfigurationSpace();

   m_solution.resize( NUM_IC );
   for (int n(0); n<NUM_IC; n++) {
      m_solution[n] = new LevelData<FArrayBox>( m_mag_geom->grids(), 1, IntVect::Zero );
   }

   m_ic_name.resize( NUM_IC );
   // To add a new IC to test, add a name to m_ic_name that corresponds to
   // name in ICFactory.  Also, add a new section "IC.name..." in the input
   // file.
   m_ic_name[ZERO]      = "zero";
   m_ic_name[CONSTANT]  = "constant";
   m_ic_name[COSINE]    = "cosine";
   m_ic_name[LOCALIZED] = "localized";
   m_ic_name[TANH]      = "tanh";

   ICFactory ic_factory;

   m_ic.resize( NUM_IC );
   for (int n(0); n<NUM_IC; n++) {
      std::string prefix( "IC." + m_ic_name[n] );
      ParmParse ppic( prefix.c_str() );
      m_ic[n] = ic_factory.create( m_ic_name[n], ppic, 1 );
   }
}



void
AnnulusTest::initialize()
{
   if (procID()==0) cout << "Entering initialize" << endl;

   setInitialConditions();
   MPI_Barrier(MPI_COMM_WORLD);

   if (procID()==0) cout << "Exiting initialize" << endl;
}



AnnulusTest::~AnnulusTest()
{
   for (int n(0); n<NUM_IC; n++) {
      delete m_solution[n];
      delete m_ic[n];
   }
   delete m_mag_geom;
   delete m_mag_geom_coords;
}



void AnnulusTest::createConfigurationSpace()
{
   string prefix( "magnetic_geometry_mapping." + MillerBlockCoordSys::pp_name );

   ParmParse pp( prefix.c_str() );

   m_mag_geom_coords
     = new MillerCoordSys( pp, m_num_cells, m_is_periodic, m_decomposition );

   DisjointBoxLayout grids;
   createDisjointBoxLayout( grids );

   const int ghosts(1);
   m_mag_geom = new MultiBlockLevelGeom( m_mag_geom_coords, grids, ghosts );

   std::cout << std::endl;
}



void AnnulusTest::createDisjointBoxLayout(DisjointBoxLayout& grids)
{
   Vector<Box> boxes;
   const MagBlockCoordSys* mag_block_coords
     = (MagBlockCoordSys *)m_mag_geom_coords->getCoordSys(0);
   const ProblemDomain& domain( mag_block_coords->domain() );
   const IntVect decomp( m_mag_geom_coords->getDecomposition(0) );
   const Box& domain_box( domain.domainBox() );

   /*
     Chop up the configuration space domain box over the number of processors specified
     for this block.  At this point, we insist that the box decomposes uniformly, or an
     error is thrown.
   */
   int nproc(1);
   int n_loc[CH_SPACEDIM];
   for (int dir(0); dir<CH_SPACEDIM; ++dir) {
     nproc *= decomp[dir];
     n_loc[dir] = 0;
   }

#if CH_SPACEDIM==2
   for (int dir(0); dir<CH_SPACEDIM; ++dir) {
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

     IntVect box_size( n_loc[RADIAL_DIR], n_loc[POLOIDAL_DIR] );
     Box patch( domain_box.smallEnd(), domain_box.smallEnd() + box_size-1 );
     Box skeleton( IntVect::Zero, IntVect(domain_box.size(0)/n_loc[RADIAL_DIR]-1,domain_box.size(1)/n_loc[POLOIDAL_DIR]-1) );
//     Box skeleton( IntVect::Zero, IntVect( domain_box.size() / n_loc - 1 ) );
     for (BoxIterator bit(skeleton); bit.ok(); ++bit) {
        Box thisBox( patch + bit() * box_size );
       boxes.push_back( thisBox );
     }
   }
   else {
     MayDay::Error( "Configuration domain box cannot be load balanced" );
   }
#else
   MayDay::Error( "AnnulusTest::createConfigurationSpace has only been implemented for 2D" );
#endif

   // Make the layout.  This is where boxes are assigned to processes.
   Vector<int> procMap;
   LoadBalance( procMap, boxes );
   m_domain = domain;
   grids.define( boxes, procMap, m_domain );
   grids.close();
}


void AnnulusTest::setInitialConditions()
{
   for (int n(0); n<NUM_IC; n++) {
      m_ic[n]->initialize( *(m_solution[n]), *m_mag_geom, 0.0 );
   }
}


void AnnulusTest::writePlotFile( const std::string& prefix )
{
   if (procID()==0) cout << "Writing plotfile" << endl;

   std::string dir_prefix( prefix );

#if 1
   std::string iter_str( prefix + "_plots" );
   if (procID() == 0) {
      mkdir( iter_str.c_str(), 0777 );
   }
   dir_prefix = iter_str + "/" + prefix;
#endif

   for (int n(0); n<NUM_IC; n++) {
      std::string filename( dir_prefix + "." + m_ic_name[n] + ".0000." );
      WriteMappedUGHDF5( filename.c_str(),
                         m_solution[n]->disjointBoxLayout(),
                         *(m_solution[n]),
                         *m_mag_geom_coords,
                         m_domain.domainBox() );
   }

   if (procID()==0) cout << "Done writing plotfile" << endl;
}



void AnnulusTest::parseParameters( ParmParse& a_pp )
{
   m_num_cells.resize( CH_SPACEDIM );
   for (int i(0); i<CH_SPACEDIM; ++i) {
      m_num_cells[i] = 0;
   }

   a_pp.getarr( "num_cells", m_num_cells, 0, m_num_cells.size() );
   for (int i(0); i<CH_SPACEDIM; ++i) {
      CH_assert( m_num_cells[i]>0 );
   }

   m_is_periodic.resize( CH_SPACEDIM );
   a_pp.getarr( "is_periodic", m_is_periodic, 0, m_is_periodic.size() );

   if (a_pp.contains("configuration_decomp")) {
      m_decomposition.resize( CH_SPACEDIM );
      for (int i(0); i<CH_SPACEDIM; ++i) {
         m_decomposition[i] = 0;
      }
      a_pp.getarr( "configuration_decomp", m_decomposition, 0, m_decomposition.size() );
      for (int i(0); i<CH_SPACEDIM; ++i) {
         CH_assert( m_decomposition[i]>0 );
      }
   }
}




#include "NamespaceFooter.H"
