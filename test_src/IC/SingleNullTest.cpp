#include "SingleNullTest.H"

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
#include "SingleNullBlockCoordSys.H"
#include "SingleNullCoordSys.H"
#include "ICFactory.H"

#define NORADDECOMP
#undef NORADDECOMP

#include "NamespaceHeader.H"


SingleNullTest::SingleNullTest( ParmParse& a_pp )
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
   m_ic_name[LOCALIZED] = "localized";

   ICFactory ic_factory;

   m_ic.resize( NUM_IC );
   for (int n(0); n<NUM_IC; n++) {
      std::string prefix( "IC." + m_ic_name[n] );
      ParmParse ppic( prefix.c_str() );
      m_ic[n] = ic_factory.create( m_ic_name[n], ppic, 1 );
   }
}



void
SingleNullTest::initialize()
{
   if (procID()==0) cout << "Entering initialize" << endl;

   setInitialConditions();
   MPI_Barrier(MPI_COMM_WORLD);

   if (procID()==0) cout << "Exiting initialize" << endl;
}



SingleNullTest::~SingleNullTest()
{
   for (int n(0); n<NUM_IC; n++) {
      delete m_solution[n];
      delete m_ic[n];
   }
   delete m_mag_geom;
   delete m_mag_geom_coords;
}



void SingleNullTest::createConfigurationSpace()
{
   ParmParse pp_grid( "singlenull" );
   std::string prefix( "magnetic_geometry_mapping."
                       + SingleNullBlockCoordSys::pp_name );
   ParmParse pp_geom( prefix.c_str() );
   m_mag_geom_coords = new SingleNullCoordSys( pp_grid, pp_geom );

   DisjointBoxLayout grids;
   createDisjointBoxLayout( grids );

   if (procID()==0) {
      std::cout << "Constructing magnetic geometry" << std::endl;
   }

   int space_order(4); // = 2*L with L as in Peter's notes
   int ghosts(0);
   m_mag_geom = new MultiBlockLevelGeom();
   m_mag_geom->define( m_mag_geom_coords, grids, ghosts, space_order );

   if (procID()==0) {
      std::cout << "Done constructing magnetic geometry"
                << std::endl << std::endl;
   }
}


void SingleNullTest::createDisjointBoxLayout( DisjointBoxLayout& grids )
{
   Vector<Box> boxes;
   for (int block(0); block<m_mag_geom_coords->numBlocks(); ++block) {

      const MagBlockCoordSys* mag_block_coords
         = (MagBlockCoordSys *)m_mag_geom_coords->getCoordSys( block );
      const ProblemDomain& domain( mag_block_coords->domain() );
      const IntVect decomp( m_mag_geom_coords->getDecomposition( block ) );
      const Box& domain_box( domain.domainBox() );

      /*
        Chop up the configuration space domain box over the number of processors specified
        for this block.  At this point, we insist that the box decomposes uniformly, or an
        error is thrown.
      */
      int nproc(1);
      for (int dir(0); dir<CH_SPACEDIM; ++dir) {
         nproc *= decomp[dir];
      }
      IntVect n_loc( D_DECL6(0,0,0,0,0,0) );

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
         Box skeleton( IntVect::Zero, IntVect( domain_box.size() / n_loc - 1 ) );
         for (BoxIterator bit( skeleton ); bit.ok();++bit) {
            Box thisBox( patch + bit() * box_size );
            boxes.push_back( thisBox );
         }
      }
      else {
         MayDay::Error( "Configuration domain box cannot be load balanced" );
      }
#else
      MayDay::Error( "GKSystem::createConfigurationSpace has only been implemented for 2D" );
#endif
   }

   // Make the layout.  This is where boxes are assigned to processes.
   Vector<int> procMap;
   LoadBalance( procMap, boxes );
   ProblemDomain prob_domain;

   Box bounding_box;
   for (int n(0); n<boxes.size(); n++) {
      bounding_box = minBox( bounding_box, boxes[n] );
   }
   prob_domain = ProblemDomain( bounding_box );

   grids.define( boxes, procMap, prob_domain );
   grids.close();
}


void SingleNullTest::setInitialConditions()
{
   for (int n(0); n<NUM_IC; n++) {
      m_ic[n]->initialize( *(m_solution[n]), *m_mag_geom, 0.0 );
   }
}


void SingleNullTest::writePlotFile( const std::string& prefix )
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


void SingleNullTest::parseParameters( ParmParse& a_pp )
{
}

#include "NamespaceFooter.H"
