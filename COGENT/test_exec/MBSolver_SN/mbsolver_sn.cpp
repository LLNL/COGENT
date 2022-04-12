#define CH_SPACEDIM 2

#include "ParmParse.H"
#include "LoadBalance.H"

#include "Diffusion.H"
#include "SingleNullCoordSys.H"

#include "MagGeom.H"

#include "parstream.H"
#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#include "MBHypreSolver.H"
#include "MBSolverUtil.H"

#include "UsingNamespace.H"


inline int checkCommandLineArgs( int a_argc, char* a_argv[] )
{
   // Check for an input file
   if (a_argc<=1) {
#ifdef CH_MPI
      if (procID()==0) {
#endif
      cout << "Usage:  cogent...ex <inputfile>" << endl;
      cout << "No input file specified" << endl;
#ifdef CH_MPI
      }
#endif
      return -1;
   }
   return 0;
}


void getDisjointBoxLayout( const SingleNullCoordSys*  a_mag_geom_coords,
                           DisjointBoxLayout&         a_grids )
{
  Vector<Box> boxes;
  for (int block=0; block<a_mag_geom_coords->numBlocks(); ++block) {

    const MagBlockCoordSys* mag_block_coords
       = (MagBlockCoordSys *)(a_mag_geom_coords->getCoordSys(block));

    const ProblemDomain& domain = mag_block_coords->domain();
    const Box& domain_box = domain.domainBox();

    const IntVect decomp = a_mag_geom_coords->getDecomposition(block);

    /*
      Chop up the configuration space domain box over the number of processors specified
      for this block.  At this point, we insist that the box decomposes uniformly, or an
      error is thrown.
    */

    IntVect n_loc = IntVect::Zero;

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

    int box_cell_num(1);
    for (int dir=0; dir<CFG_DIM; ++dir) {
      box_cell_num *= n_loc[dir];
    }
    
    if (box_cell_num > 0) {
      IntVect box_size(n_loc);
      Box patch(domain_box.smallEnd(), domain_box.smallEnd() + box_size-1);
      IntVect lo = IntVect::Zero;
      IntVect hi;
      for (int dir=0; dir<CFG_DIM; ++dir) {
        hi[dir] = domain_box.size(dir)/n_loc[dir]-1;
      }
      Box skeleton(lo, hi);
      BoxIterator bit(skeleton);
      for (bit.begin();bit.ok();++bit) {
        Box thisBox = patch + bit()*box_size;
        boxes.push_back(thisBox);
      }
    }
    else {
      MayDay::Error( "Configuration domain box cannot be load balanced" );
    }
  }

  // Make the layout.  This is where boxes are assigned to processes.
  Vector<int> procMap;
  LoadBalance( procMap, boxes );

  Box bounding_box;
  for (int n=0; n<boxes.size(); n++) {
     bounding_box = minBox(bounding_box, boxes[n]);
  }
  ProblemDomain prob_domain(bounding_box);

  a_grids.define( boxes, procMap, prob_domain );
  a_grids.close();

#ifdef CH_MPI
  if (procID()==0) {
#endif
    if (true) {
      for (int n=0; n<boxes.size(); n++) {
        const Box& local_box = boxes[n];
        cout << "   ";
        cout  << "Configuration space box " 
              << local_box << " is assigned to process " 
              << procMap[n] << endl;
      }
    }
#ifdef CH_MPI
  }
#endif
}


int main( int a_argc, char* a_argv[] )
{
#ifdef CH_MPI
   // Start MPI
   MPI_Init( &a_argc, &a_argv );
   setChomboMPIErrorHandler();
#endif

   int status = checkCommandLineArgs( a_argc, a_argv );

   if (status==0) {
      ParmParse pp( a_argc-2, a_argv+2, NULL, a_argv[1] );

      ParmParse pp_grid("grid");
      ParmParse pp_geom("geom");
      SingleNullCoordSys* mb_coord_sys = new SingleNullCoordSys(pp_grid, pp_geom);

      DisjointBoxLayout grids;
      getDisjointBoxLayout(mb_coord_sys, grids);

      int ghosts = 4;
      MagGeom geom(pp, RefCountedPtr<MagCoordSys>(mb_coord_sys), grids, ghosts);

      int order = 4;
#ifdef CH_MPI
      if (procID()==0) {
#endif
         cout << "Order = " << order << endl;
#ifdef CH_MPI
      }
#endif

      MBHypreSolver solver(geom, 1, order);

      ParmParse ppls("linear_solver");
      solver.setMethodParams(ppls);
      solver.setConvergenceParams(ppls);

      LevelData<FArrayBox> rhs(grids, 1, IntVect::Zero);
      LevelData<FArrayBox> phi(grids, 1, IntVect::Zero);

      double rhs_val;
      pp.get("rhs", rhs_val);

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         rhs[dit].setVal(rhs_val, 0);
         phi[dit].setVal(0., 0);
      }

      LevelData<FluxBox> coef(grids, SpaceDim*SpaceDim, IntVect::Unit);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            coef[dit][dir].setVal(1.,0);
            coef[dit][dir].setVal(0.,1);
            coef[dit][dir].setVal(0.,2);
            coef[dit][dir].setVal(1.,3);
         }
      }

      ParmParse ppsp( "BC" );
      int nblocks = mb_coord_sys->numBlocks();
      EllipticOpBC* bc = new SingleNullEllipticOpBC( "potential", ppsp, nblocks, nblocks, true );

      MBSolverUtil solver_util(solver);
      LevelData<FluxBox> mapped_coef;
      solver_util.computeMappedCoefficients(coef, mapped_coef);

      Diffusion op(pp, geom);
      op.setNeumannNaturalFactor(coef, *bc);
      op.modifyForNeumannAndNaturalBCs(*bc, coef, mapped_coef );

      geom.plotFaceData( "coef", coef, 0.);
      geom.plotFaceData( "mapped_coef", mapped_coef, 0.);

      LevelData<FArrayBox> alpha(grids, 1, IntVect::Zero);
      solver_util.getCellVolumes(alpha);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         alpha[dit].invert(1.);
      }

      solver.initializeMatrix();
      solver.constructMatrixBlock(0, 0, alpha, mapped_coef, *bc);
      solver.finalizeMatrix();

      bool homogeneous_bcs = false;
      solver.solve(rhs, phi, homogeneous_bcs);

      geom.plotCellData( "solution", phi, 0.);
   }

#ifdef CH_MPI
   MPI_Finalize();
#endif

   return status;
}
