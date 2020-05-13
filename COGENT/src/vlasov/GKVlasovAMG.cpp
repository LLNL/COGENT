#include "GKVlasovAMG.H"

// Hypre headers
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_sstruct_mv.h"

#include "NamespaceHeader.H"


GKVlasovAMG::GKVlasovAMG( const ParmParse& a_pp )
   : m_hypre_allocated(false),
     m_A(NULL),
     m_AMG_solver_allocated(false),
     m_discretization_order(2)
{
   if ( a_pp.contains("restrict_type") ) {
      a_pp.get("restrict_type", m_restrict_type);
   }
   else {
      m_restrict_type = 2;
   }

   double tol;
   if ( a_pp.contains("tol") ) {
      a_pp.get("tol", tol);
   }
   else {
      tol = 0.;
   }

   int max_iter;
   if ( a_pp.contains("max_iter") ) {
      a_pp.get("max_iter", max_iter);
   }
   else {
      max_iter = 2;
   }

   bool verbose;
   if ( a_pp.contains("verbose") ) {
      a_pp.get("verbose", verbose);
   }
   else {
      verbose = false;
   }

   setConvergenceParams(tol, max_iter, verbose);
}

void GKVlasovAMG::createHypreData( const DisjointBoxLayout&   a_grids,
                                   const MultiBlockCoordSys*  a_coord_sys,
                                   const string&              a_face_avg_type )
{
   if (m_hypre_allocated) {
      MayDay::Error( "GKVlasovAMG:createHypreData:: data already allocated; deallocate first!" );
   }

   if ( a_face_avg_type == "uw1" ) {
      m_stencil_radius = 1;
   }
   else if ( a_face_avg_type == "uw3" ) {
      m_stencil_radius = 2;
   }
   else {
      MayDay::Error("GKVlasovAMG::createHypreData(): Unrecognized face_avg_type");
   }

   m_hypre_object_type = HYPRE_PARCSR;

   Box stencil_box(IntVect::Zero, IntVect::Zero);
   stencil_box.grow(m_stencil_radius);
   stencil_box.shift(m_stencil_radius*IntVect::Unit);

   int num_blocks = a_coord_sys->numBlocks();

   int nvar = 1;
   int var = 0;

   // Create the grid

   {
      HYPRE_SStructGridCreate(MPI_COMM_WORLD, SpaceDim, num_blocks, &m_grid);

      for (DataIterator dit(a_grids); dit.ok(); ++dit) {
         const Box & box = a_grids[dit];
         int block_number = a_coord_sys->whichBlock(box);
         IntVect lower(box.loVect());
         IntVect upper(box.hiVect());

         HYPRE_SStructGridSetExtents(m_grid, block_number, lower.dataPtr(), upper.dataPtr());
      }

      HYPRE_SStructVariable vartypes[1] = {HYPRE_SSTRUCT_VARIABLE_CELL};

      for (int block=0; block<num_blocks; ++block) {
         HYPRE_SStructGridSetVariables(m_grid, block, nvar, vartypes);
      }

      // Specify periodicity.  This is only applicable for single block geometry,
      // since there doesn't seem to be any support for individually periodic blocks
      // in MultiBlockLevelGeom.

      if (num_blocks == 1) {

         const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = a_coord_sys->boundaries();

         const ProblemDomain& domain = a_grids.physDomain();
         IntVect periodic = domain.size();
         int block_number = 0;

         for (int dir=0; dir<SpaceDim; dir++) {
            bool this_dir_periodic = domain.isPeriodic(dir);
#if 0
            if ( !this_dir_periodic || 
                 (block_boundaries[block_number][dir].isInterface() && a_mblex_Ptr != NULL) ||
                 (block_boundaries[block_number][dir + SpaceDim].isInterface() && a_mblex_Ptr != NULL) ) {
                 (block_boundaries[block_number][dir].isInterface() ) ||
                 (block_boundaries[block_number][dir + SpaceDim].isInterface() ) ) {
               periodic[dir] = 0;
            }
#endif
            if ( !this_dir_periodic ) {
               periodic[dir] = 0;
            }
            else {
               // This direction is periodic.  Make sure that the number of
               // cells in this direction is a power of 2, which the Hypre
               // multigrid solver requires.
               int ncells = periodic[dir];
               while (ncells%2 == 0) {
                  ncells /= 2;
               }
               if (ncells != 1) {
                  if (procID()==0) {
                     MayDay::Error( "GKVlasovAMG::createHypreData(): Periodic direction must be a power of 2!" );
                  }
               }
            }
         }

         HYPRE_SStructGridSetPeriodic(m_grid, block_number, periodic.dataPtr());
      }

      int num_ghost[CH_SPACEDIM*2];
      for (int n=0; n<2*SpaceDim; ++n) {
         num_ghost[n] = m_stencil_radius;
      }

      HYPRE_SStructGridSetNumGhost(m_grid, num_ghost);

      HYPRE_SStructGridAssemble(m_grid);
   }

   // Define the discretization stencil for the system matrix A

   {
      /* Define the geometry of the stencil. Each represents a
         relative offset (in the index space). */

      int stencil_size = stencil_box.numPts();

      int** offsets = new int*[stencil_size];
      int n = -1;
      for (BoxIterator bit(stencil_box - m_stencil_radius*IntVect::Unit); bit.ok(); ++bit) {
         n++;
         offsets[n] = new int[SpaceDim];
         IntVect iv = bit();
         for (int m=0; m<SpaceDim; ++m) {
            offsets[n][m] = iv[m];
         }
      }

      m_A_diagonal_offset = (stencil_size - 1) / 2;

      /* Create an empty stencil object */
      HYPRE_SStructStencilCreate(SpaceDim, stencil_size, &m_A_stencil);

      /* Assign each of the stencil entries */
      for (int entry = 0; entry < stencil_size; entry++)
        HYPRE_SStructStencilSetEntry(m_A_stencil, entry, offsets[entry], var);

      for (int n=0; n<stencil_size; ++n) {
         delete [] offsets[n];
      }
      delete [] offsets;
   }

   // Set up the graph

   {
      HYPRE_SStructGraphCreate(MPI_COMM_WORLD, m_grid, &m_A_graph);

      HYPRE_SStructGraphSetObjectType(m_A_graph, m_hypre_object_type);

      for (int block=0; block<num_blocks; ++block) {
         HYPRE_SStructGraphSetStencil(m_A_graph, block, var, m_A_stencil);
      }

      HYPRE_SStructGraphAssemble(m_A_graph);
   }

   // Set up b and x vectors

   {
      /* Create an empty vector object */
      HYPRE_SStructVectorCreate(MPI_COMM_WORLD, m_grid, &m_b);
      HYPRE_SStructVectorCreate(MPI_COMM_WORLD, m_grid, &m_x);

      HYPRE_SStructVectorSetObjectType(m_b, m_hypre_object_type);
      HYPRE_SStructVectorSetObjectType(m_x, m_hypre_object_type);

      HYPRE_SStructVectorInitialize(m_b);
      HYPRE_SStructVectorInitialize(m_x);
   }

   m_hypre_allocated = true;
}


void GKVlasovAMG::AMGSetup( const HYPRE_SStructMatrix& a_matrix )
{
   CH_TIME("GKVlasovAMG::AMGSetup");

   if ( m_AMG_solver_allocated ) {
      HYPRE_BoomerAMGDestroy(m_par_AMG_solver);
   }
   
   HYPRE_BoomerAMGCreate(&m_par_AMG_solver);

   /*
     Instructions from Ruipeng Li (11/6/2019):

     To use the AIR (approximate ideal restriction) algorithm in hypre, you need to set
     HYPRE_BoomerAMGSetRestriction(amg_solver, restri_type) with restri_type = 1, or 2 for distance-1
     and -2 AIR respectively. Distance-2 is more expensive to build and leads to a denser restriction
     operator than the distance-1 but should approximate the ideal restriction better. Do not use aggressive
     coarsening with AIR ( i.e., HYPRE_ParCSRHybridSetAggNumLevels(amg, 0) ). 

     For hyperbolic problems like advection dominant, highly nonsymmetric, using very simple interpolation
     like 1-pt interpolation was often found enough, set HYPRE_BoomerAMGSetInterpType(amg, 100) for that.
     Also,  C/F relaxations can be helpful, in the paper, we used FFC, which can be set by
     HYPRE_BoomerAMGSetGridRelaxPoints(amg_solver, grid_relax_points); with          
         grid_relax_points[2][0] = -1; // F
         grid_relax_points[2][1] = -1; // F
         grid_relax_points[2][2] =  1; // C

     These options can be found in hypre, src/test/ij.c, search for air or AIR to see the default settings.
   */

   HYPRE_BoomerAMGSetRestriction(m_par_AMG_solver, m_restrict_type);

   HYPRE_Int **grid_relax_points = NULL;
   HYPRE_Int ns_down = 1;
   HYPRE_Int ns_up = 3;
   HYPRE_Int ns_coarse = 1;
   /* this is a 2-D 4-by-k array using Double pointers */
   grid_relax_points = hypre_CTAlloc(HYPRE_Int*, 4, HYPRE_MEMORY_HOST);
   grid_relax_points[0] = NULL;
   grid_relax_points[1] = hypre_CTAlloc(HYPRE_Int, ns_down, HYPRE_MEMORY_HOST);
   grid_relax_points[2] = hypre_CTAlloc(HYPRE_Int, ns_up, HYPRE_MEMORY_HOST);
   grid_relax_points[3] = hypre_CTAlloc(HYPRE_Int, ns_coarse, HYPRE_MEMORY_HOST);
   /* down cycle: C */
   for (int i=0; i<ns_down; i++) {
      grid_relax_points[1][i] = 0;//1;
   }
   /* up cycle: F */
   //for (int i=0; i<ns_up; i++)
   //{
   if (ns_up == 3) {
      grid_relax_points[2][0] = -1; // F
      grid_relax_points[2][1] = -1; // F
      grid_relax_points[2][2] =  1; // C
   }
   else if (ns_up == 2) {
      grid_relax_points[2][0] = -1; // F
      grid_relax_points[2][1] = -1; // F
   }
   //}
   /* coarse: all */
   for (int i=0; i<ns_coarse; i++) {
      grid_relax_points[3][i] = 0;
   }
   HYPRE_BoomerAMGSetGridRelaxPoints(m_par_AMG_solver, grid_relax_points);

   HYPRE_BoomerAMGSetRelaxType(m_par_AMG_solver, 0);  // Jacobi

   HYPRE_Int interp_type = 100;  // 1-pt
   HYPRE_BoomerAMGSetInterpType(m_par_AMG_solver, interp_type);
   
   HYPRE_ParCSRHybridSetAggNumLevels(m_par_AMG_solver, 0);

   HYPRE_Int coarse_threshold = 0;
   HYPRE_ParCSRHybridSetMaxCoarseSize(m_par_AMG_solver, coarse_threshold);

   //   HYPRE_BoomerAMGSetStrongThreshold(m_par_AMG_solver, 0.25);
   //   HYPRE_BoomerAMGSetStrongThresholdR(m_par_AMG_solver, 0.25);
   //   HYPRE_BoomerAMGSetCoarsenType(m_par_AMG_solver, 6);  // Falgout coarsening
   //   HYPRE_BoomerAMGSetCycleRelaxType(m_par_AMG_solver, 3, 3);  // hybrid Gauss-Seidel on coarsest level
   //   HYPRE_BoomerAMGSetCycleRelaxType(m_par_AMG_solver, 9, 3);  // Gaussian elimination on coarsest level
   HYPRE_BoomerAMGSetCycleType(m_par_AMG_solver, 2);  // V = 1, W = 2
   //   HYPRE_BoomerAMGSetCycleNumSweeps(m_par_AMG_solver, 1, 1);  // 1 sweep on down cycle
   //   HYPRE_BoomerAMGSetCycleNumSweeps(m_par_AMG_solver, 1, 2);  // 1 sweep on up cycle
   //   HYPRE_BoomerAMGSetCycleNumSweeps(m_par_AMG_solver, 1, 3);  // 1 sweeps on coarsest level

   m_AMG_solver_allocated = true;

   HYPRE_ParCSRMatrix par_A;
   HYPRE_SStructMatrixGetObject(a_matrix, (void **) &par_A);

   HYPRE_ParVector par_b;  // not used, but needs to be passed to HYPRE_BoomerAMGSetup
   HYPRE_ParVector par_x;  // not used, but needs to be passed to HYPRE_BoomerAMGSetup

   HYPRE_BoomerAMGSetup(m_par_AMG_solver, par_A, par_b, par_x);
}


void GKVlasovAMG::constructMatrix( const MultiBlockCoordSys*               a_coord_sys,
                                   const string&                           a_face_avg_type,
                                   LevelData<BaseFab<Vector<IntVect> > >&  a_couplings,
                                   LevelData<BaseFab<Vector<Real> > >&     a_weights )
{
   const DisjointBoxLayout& grids = a_couplings.disjointBoxLayout();
   CH_assert(grids == a_weights.disjointBoxLayout());

   if ( !m_hypre_allocated ) {
      m_coord_sys_ptr = a_coord_sys;

      createHypreData(grids, m_coord_sys_ptr, a_face_avg_type);
   }

   HYPRE_SStructMatrixDestroy(m_A);
   HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, m_A_graph, &m_A);
   HYPRE_SStructMatrixSetObjectType(m_A, m_hypre_object_type);

   HYPRE_SStructMatrixInitialize(m_A);

   int var = 0;

   const ProblemDomain& problem_domain = grids.physDomain();
   const Box& domain_box = problem_domain.domainBox();

   bool at_lo_boundary[CH_SPACEDIM];
   bool at_hi_boundary[CH_SPACEDIM];

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box& box = grids[dit];
      int block_number = m_coord_sys_ptr->whichBlock(box);
      const BaseFab<Vector<IntVect> >& this_couplings = a_couplings[dit];
      const BaseFab<Vector<Real> >& this_weights = a_weights[dit];

      // Set the intra-block matrix entries

      for (BoxIterator bit(box); bit.ok(); ++bit) {
         IntVect iv = bit();

         Box stencil_box(iv,iv);
         stencil_box.grow(m_stencil_radius);

         for (int dir=0; dir<SpaceDim; ++dir) {
            at_lo_boundary[dir] = iv[dir] == domain_box.smallEnd(dir) && !problem_domain.isPeriodic(dir);
            at_hi_boundary[dir] = iv[dir] == domain_box.bigEnd(dir) && !problem_domain.isPeriodic(dir);
         }

         int* index = iv.dataPtr();
         int nentries = 1;

         for (int n=0; n<this_couplings(iv,0).size(); ++n) {
            IntVect coupling = this_couplings(iv,0)[n];

            bool add_entry = true;
#if 0
            for (int dir=0; dir<SpaceDim; ++dir) {
               if ( (at_lo_boundary[dir] && coupling[dir] < iv[dir]) ||
                    (at_hi_boundary[dir] && coupling[dir] > iv[dir]) ) {
                  add_entry = false;
                  break;
               }
            }
#else
            for (int dir=0; dir<SpaceDim; ++dir) {
               if ( (coupling[dir] < domain_box.smallEnd(dir) && !problem_domain.isPeriodic(dir)) ||
                    (coupling[dir] > domain_box.bigEnd(dir)   && !problem_domain.isPeriodic(dir) ) ) {
                  add_entry = false;
                  break;
               }
            }
#endif
            
            if ( add_entry ) {
               int entries = findHypreEntry(stencil_box, coupling);
               double values = this_weights(iv,0)[n];
               HYPRE_SStructMatrixSetValues(m_A, block_number, index, var, nentries, &entries, &values);
            }
         }
      }
   }

   /* This is a collective call finalizing the matrix assembly.
      The matrix is now ``ready to be used'' */
   HYPRE_SStructMatrixAssemble(m_A);

   //   HYPRE_SStructMatrixPrint("HYPRE_MATRIX.A", m_A, 0);

   AMGSetup(m_A);
}


int GKVlasovAMG::findHypreEntry(const Box&         a_stencil_box,
                                const IntVect&     a_iv) const
{
   /*
     The Hypre matrix coefficients corresponding to the equation centered at cell a_iv are
     identified by an entry number, which is returned by this function.  Regular "stencil"
     couplings to cells in a_stencil_box are enumerated first (in the order with which a 
     BoxIterator built on a_stencil_box traverses cells), followed by the unstructured couplings
     to cells in other blocks contained in a_unstructured_ivs.
   */

   bool found_entry = false;

   int entry = 0;
   BoxIterator bit(a_stencil_box);
   for (bit.begin(); bit.ok(); ++bit) {
      if (bit() == a_iv) {
         found_entry = true;
         break;
      }
      entry++;
   }

   CH_assert(found_entry);

   return entry;
}


void GKVlasovAMG::multiplyMatrix( const LevelData<FArrayBox>&  a_in,
                                  LevelData<FArrayBox>&        a_out ) const
{
   if (!m_hypre_allocated) {
      MayDay::Error( "GKVlasovAMG::multiplyMatrix(): Operator has not yet been initialized!" );
   }

   /* Create temporary Hypre vector objects */

   HYPRE_SStructVector in_vector;
   HYPRE_SStructVector out_vector;

   HYPRE_SStructVectorCreate(MPI_COMM_WORLD, m_grid, &in_vector);
   HYPRE_SStructVectorCreate(MPI_COMM_WORLD, m_grid, &out_vector);

   HYPRE_SStructVectorSetObjectType(in_vector, m_hypre_object_type);
   HYPRE_SStructVectorSetObjectType(out_vector, m_hypre_object_type);

   HYPRE_SStructVectorInitialize(in_vector);
   HYPRE_SStructVectorInitialize(out_vector);

   // Copy the input LevelData to its Hypre vector

   int var = 0;

   const DisjointBoxLayout& grids = a_in.disjointBoxLayout();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box & box = grids[dit];
      int block_number = m_coord_sys_ptr->whichBlock(grids[dit]);

      IntVect lower(box.loVect());
      IntVect upper(box.hiVect());

      FArrayBox tmp(box,1);
      tmp.copy(a_in[dit]);

      HYPRE_SStructVectorSetBoxValues(in_vector, block_number, lower.dataPtr(), upper.dataPtr(),
                                      var, tmp.dataPtr());
   }

   HYPRE_SStructVectorAssemble(in_vector);
   HYPRE_SStructVectorAssemble(out_vector);

   // Do the matrix-vector multiply

   void *matvec_vdata;
   hypre_SStructMatvecCreate( &matvec_vdata );

   hypre_SStructMatvecSetup(matvec_vdata, m_A, in_vector);

   hypre_SStructMatvecCompute(matvec_vdata, 1., m_A, in_vector, 1., out_vector);

   hypre_SStructMatvecDestroy (matvec_vdata);

   if (m_hypre_object_type == HYPRE_PARCSR) {
     HYPRE_SStructVectorGather(out_vector);
   }

   // Copy the output Hypre vector to its LevelData

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box & box = grids[dit];
      int block_number = m_coord_sys_ptr->whichBlock(grids[dit]);

      IntVect lower(box.loVect());
      IntVect upper(box.hiVect());

      FArrayBox tmp(box,1);

      HYPRE_SStructVectorGetBoxValues(out_vector, block_number,
                                      lower.dataPtr(), upper.dataPtr(), var, tmp.dataPtr());

      a_out[dit].copy(tmp);
   }

   HYPRE_SStructVectorDestroy(out_vector);
   HYPRE_SStructVectorDestroy(in_vector);
}


void GKVlasovAMG::setConvergenceParams( const double  a_tol,
                                        const int     a_max_iter,
                                        const bool    a_verbose )
{
   m_tol       = a_tol;
   m_max_iter  = a_max_iter;
   m_verbose   = a_verbose;

   m_convergence_params_set = true;
}


void GKVlasovAMG::solve( const LevelData<FArrayBox>&  a_rhs,
                         LevelData<FArrayBox>&        a_solution,
                         bool                         a_homogeneous_bcs )
{
   if ( !m_convergence_params_set ) {
      MayDay::Error("GKVlasovAMG::solve(): solver convergence parameters have not been set");
   }

   if ( a_homogeneous_bcs ) {
      copyToHypreVector(a_rhs, m_b);
   }
   else {
      LevelData<FArrayBox> rhs;
      rhs.define(a_rhs);

#if 0
      for (DataIterator dit(rhs.dataIterator()); dit.ok(); ++dit) {
         rhs[dit] += m_rhs_from_bc[dit];
      }
#endif

      copyToHypreVector(rhs, m_b);
   }

   copyToHypreVector(a_solution, m_x);

   /* This is a collective call finalizing the vector assembly.
      The vectors are now ``ready to be used'' */
   HYPRE_SStructVectorAssemble(m_b);
   HYPRE_SStructVectorAssemble(m_x);

   AMG( m_A, m_b, m_tol, m_max_iter, m_verbose, m_x );

   copyFromHypreVector(m_x, a_solution);

   a_solution.exchange();
}


void GKVlasovAMG::AMG( const HYPRE_SStructMatrix&  a_matrix,
                       const HYPRE_SStructVector&  a_b,
                       const double                a_tol,
                       const int                   a_max_iter,
                       const bool                  a_verbose,
                       const HYPRE_SStructVector&  a_x )
{
   CH_TIMERS("MBHypreSolver::AMG");
   CH_TIMER("AMG_solve",t_AMG_solve);

   HYPRE_ParCSRMatrix    par_A;
   HYPRE_ParVector       par_b;
   HYPRE_ParVector       par_x;

   HYPRE_SStructMatrixGetObject(a_matrix, (void **) &par_A);
   HYPRE_SStructVectorGetObject(a_b, (void **) &par_b);
   HYPRE_SStructVectorGetObject(a_x, (void **) &par_x);

   HYPRE_BoomerAMGSetTol(m_par_AMG_solver, a_tol);
   HYPRE_BoomerAMGSetMaxIter(m_par_AMG_solver, a_max_iter);
   HYPRE_BoomerAMGSetPrintLevel(m_par_AMG_solver, 1);

   CH_START(t_AMG_solve);
   HYPRE_BoomerAMGSolve(m_par_AMG_solver, par_A, par_b, par_x);
   CH_STOP(t_AMG_solve);

   int num_iterations;
   HYPRE_BoomerAMGGetNumIterations(m_par_AMG_solver, &num_iterations);
   double final_res_norm;
   HYPRE_BoomerAMGGetFinalRelativeResidualNorm(m_par_AMG_solver, &final_res_norm);
   if (a_verbose && procID()==0) {
      cout << "        --> GKVlasovAMG solver residual = " << final_res_norm << " after " << num_iterations << " iterations" << endl;
   }

   if (m_hypre_object_type == HYPRE_PARCSR) {
      HYPRE_SStructVectorGather(a_x);
   }
}

   
void GKVlasovAMG::copyToHypreVector( const LevelData<FArrayBox>&  a_in,
                                     HYPRE_SStructVector&         a_out ) const
{
   int var = 0;

   const DisjointBoxLayout & grids = a_in.disjointBoxLayout();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box & box = grids[dit];
      int block_number = m_coord_sys_ptr->whichBlock(grids[dit]);

      IntVect lower(box.loVect());
      IntVect upper(box.hiVect());

      FArrayBox tmp(box,1);
      tmp.copy(a_in[dit]);

      HYPRE_SStructVectorSetBoxValues(a_out, block_number, lower.dataPtr(), upper.dataPtr(), var, tmp.dataPtr());
   }
} 



void GKVlasovAMG::copyFromHypreVector( const HYPRE_SStructVector&  a_in,
                                       LevelData<FArrayBox>&       a_out ) const
{
   int var = 0;

   const DisjointBoxLayout & grids = a_out.disjointBoxLayout();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box & box = grids[dit];
      int block_number = m_coord_sys_ptr->whichBlock(grids[dit]);

      FArrayBox tmp(box,1);

      IntVect lower(box.loVect());
      IntVect upper(box.hiVect());

      HYPRE_SStructVectorGetBoxValues(a_in, block_number, lower.dataPtr(), upper.dataPtr(), var, tmp.dataPtr());

      a_out[dit].copy(tmp);
   }
} 


#include "NamespaceFooter.H"

