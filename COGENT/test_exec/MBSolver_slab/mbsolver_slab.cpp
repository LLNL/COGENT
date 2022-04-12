#define CH_SPACEDIM 2

#include "ParmParse.H"
#include "LoadBalance.H"

#include "Diffusion.H"
#include "LogRectCoordSys.H"
#include "LogRectEllipticOpBC.H"


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


void
getDisjointBoxLayout( const MultiBlockCoordSys& a_mb_coords,
                      Vector<IntVect>&          a_decomposition,
                      DisjointBoxLayout&        a_grids )
{
  Vector<Box> boxes;
  const Vector<Box>& mapping_blocks = a_mb_coords.mappingBlocks();
  for (int block=0; block<a_mb_coords.numBlocks(); ++block) {
    const Box& domain_box = mapping_blocks[block];
    const IntVect decomp = a_decomposition[block];

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

    if (n_loc[0] != 0 && n_loc[1] != 0) {

      IntVect box_size(n_loc[0],n_loc[1]);
      Box patch(domain_box.smallEnd(), domain_box.smallEnd() + box_size-1);
      Box skeleton(IntVect::Zero, IntVect(domain_box.size(0)/n_loc[0]-1,domain_box.size(1)/n_loc[1]-1));
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

  ProblemDomain prob_domain;

  Box bounding_box;
  for (int n=0; n<boxes.size(); n++) {
     bounding_box = minBox(bounding_box, boxes[n]);
  }
  prob_domain = ProblemDomain(bounding_box);

  a_grids.define( boxes, procMap, prob_domain );
  a_grids.close();

#ifdef CH_MPI
  if (procID()==0) {
#endif
     for (int n=0; n<boxes.size(); n++) {
        const Box& local_box = boxes[n];
        cout << "   Box " << local_box << " is assigned to process " << procMap[n] << endl;
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

      ParmParse pp_slab("slab");

      std::vector<int> num_cells(SpaceDim);
      pp_slab.getarr("num_cells", num_cells, 0, SpaceDim);
      std::vector<bool> is_periodic(SpaceDim);
      std::vector<int> decomp(SpaceDim);
      for (int n=0; n<SpaceDim; ++n) {
         is_periodic[n] = false;
         decomp[n] = 1;
      }

      LogRectCoordSys* mb_coord_sys = new LogRectCoordSys(pp_slab, "slab", num_cells, is_periodic, decomp);

      // Create a DisjointBoxLayout
      Vector<IntVect> decomposition(mb_coord_sys->numBlocks());
      for (int n=0; n<decomposition.size(); ++n) {
         decomposition[n] = IntVect::Unit;
      }
      DisjointBoxLayout grids;
      getDisjointBoxLayout( *mb_coord_sys, decomposition, grids );

      int ghosts = 4;
      MagGeom geom(pp, RefCountedPtr<MagCoordSys>(mb_coord_sys), grids, ghosts);

      int order = 2;
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

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         rhs[dit].setVal(1.e-2);
         phi[dit].setVal(0.);
      }

      LevelData<FluxBox> coef(grids, SpaceDim*SpaceDim, IntVect::Unit);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            coef[dit][dir].setVal(3.,0);
            coef[dit][dir].setVal(1.,1);
            coef[dit][dir].setVal(1.,2);
            coef[dit][dir].setVal(1.,3);
         }
      }

      bool periodic[CH_SPACEDIM];
      for (int i=0; i<SpaceDim; ++i) periodic[i] = is_periodic[i];
      ParmParse ppsp( "BC.potential" );
      EllipticOpBC* bc = new LogRectEllipticOpBC( "potential", ppsp, periodic, true );

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

      geom.plotCellData("solution", phi, 0.);
   }

#ifdef CH_MPI
   MPI_Finalize();
#endif

   return status;
}
