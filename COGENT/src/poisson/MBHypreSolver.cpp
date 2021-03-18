#include "MBHypreSolver.H"
#include "MBStencilIterator.H"
#include "BlockRegister.H"
#include "SparseCoupling.H"
#include "MBSolverF_F.H"

#include "_hypre_sstruct_mv.h"


#include "NamespaceHeader.H"



const std::string MBHypreSolver::pp_name = "MBHypreSolver";


MBHypreSolver::MBHypreSolver( const MultiBlockLevelGeom&      a_geom,
                              const int                       a_nvar,  
                              const int                       a_discretization_order,
                              MultiBlockLevelExchangeCenter*  a_mblex_ptr)
   : MBSolver(a_geom, a_discretization_order, a_mblex_ptr),
     m_hypre_allocated(false),
     m_A(NULL),
     m_nvar(a_nvar),
     m_AMG_solver_allocated(false),
     m_ILU_solver_allocated(false),
     m_MGR_solver_allocated(false),
     m_matrix_initialized(false),
     m_matrix_finalized(false)
{
   createHypreData();
}
      

MBHypreSolver::~MBHypreSolver()
{
   destroyHypreData();
}


void
MBHypreSolver::multiplyMatrix( const LevelData<FArrayBox>&  a_in,
                               LevelData<FArrayBox>&        a_out ) const
{
   if (!m_hypre_allocated) {
      MayDay::Error( "MBHypreSolver::multiplyMatrix(): Operator has not yet been initialized!" );
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

   const DisjointBoxLayout& grids = a_in.disjointBoxLayout();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box & box = grids[dit];
      int block_number = m_coord_sys_ptr->whichBlock(grids[dit]);

      IntVect lower(box.loVect());
      IntVect upper(box.hiVect());

      FArrayBox tmp(box,1);

      for ( int var=0; var<m_nvar; ++var ) {
         tmp.copy(a_in[dit],var,0,1);

         HYPRE_SStructVectorSetBoxValues(in_vector, block_number, lower.dataPtr(),
                                         upper.dataPtr(), var, tmp.dataPtr());
      }
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

      for (int var=0; var<m_nvar; ++var) {
         HYPRE_SStructVectorGetBoxValues(out_vector, block_number,
                                         lower.dataPtr(), upper.dataPtr(), var, tmp.dataPtr());

         a_out[dit].copy(tmp,0,var,1);
      }
   }

   HYPRE_SStructVectorDestroy(out_vector);
   HYPRE_SStructVectorDestroy(in_vector);
}


void
MBHypreSolver::setMethodParams( const ParmParse&  a_pp )
{
   if ( a_pp.query("method", m_method) == 0 ) m_method = "AMG";
   if ( a_pp.query("precond_method", m_precond_method) == 0 ) m_precond_method = "";

   if ( m_method == "AMG" ) {
      if ( a_pp.query("coarsen_type", m_AMG_coarsen_type) == 0 ) m_AMG_coarsen_type = 6;
      if ( a_pp.query("strong_threshold", m_AMG_strong_threshold) == 0 ) m_AMG_strong_threshold = 0.25;
      if ( a_pp.query("amg_old_setup", m_AMG_old_setup) == 0 ) m_AMG_old_setup = false;
   }
   else if ( m_method == "MGR" ) {
      if ( a_pp.query("cpoint", m_MGR_cpoint) == 0 ) m_MGR_cpoint = 1;
      if ( a_pp.query("print_level", m_MGR_print_level) == 0 ) m_MGR_print_level = 0;
      if ( a_pp.query("amg_print_level", m_MGR_amg_print_level) == 0 ) m_MGR_amg_print_level = 0;
      if ( a_pp.query("amg_max_iter", m_MGR_amg_max_iter) == 0 ) m_MGR_amg_max_iter = 1;
   }
   else if ( m_method == "ILU" ) {
      if ( a_pp.query("type", m_ILU_type) == 0 ) m_ILU_type = 1;
      if ( a_pp.query("fill_level", m_ILU_fill_level) == 0 ) m_ILU_fill_level = 1;
      if ( a_pp.query("drop_threshold", m_ILU_drop_threshold) == 0 ) m_ILU_drop_threshold = 1.e-2;
      if ( a_pp.contains("drop_threshold_array") ) {
         Vector<double> tmp;
         a_pp.getarr( "drop_threshold_array", tmp, 0, CFG_DIM);
         for (int i=0; i<CFG_DIM; ++i) {
            m_ILU_drop_threshold_array[i] = tmp[i];
         }
      }
      else {
         for (int i=0; i<CFG_DIM; ++i) {
            m_ILU_drop_threshold_array[i] = 1.e-2;
         }
      }
      if ( a_pp.query("max_nnz_per_row", m_ILU_max_nnz_per_row) ) m_ILU_max_nnz_per_row = 1000;
      if ( a_pp.query("max_iter", m_ILU_max_iter) ) m_ILU_max_iter = 20;
      if ( a_pp.query("max_schur_iter", m_ILU_max_schur_iter) ) m_ILU_max_schur_iter = 5;
      if ( a_pp.query("tol", m_ILU_tol) ) m_ILU_tol = 1.e-7;
      if ( a_pp.query("print_level", m_ILU_print_level) ) m_ILU_print_level = 0;
   }
   else {
      MayDay::Error("MBHypreSolver::setMethodParams(): unknown method");
   }

   m_method_params_set = true;
}


void
MBHypreSolver::solve( const LevelData<FArrayBox>&  a_rhs,
                      LevelData<FArrayBox>&        a_solution,
                      bool                         a_homogeneous_bcs )
{
   CH_TIMERS("MBHypreSolver::solve");
   CH_TIMER("HYPRE_SStructVectorAssemble",t_HYPRE_SStructVectorAssemble);
   CH_assert(a_rhs.nComp() == m_nvar);
   CH_assert(a_solution.nComp() == m_nvar);
   
   if ( !m_matrix_finalized ) {
      MayDay::Error("MBHypreSolver::solve(): Hypre matrix has not been finalized");
   }

   if ( !m_convergence_params_set ) {
      MayDay::Error("MBHypreSolver::solve(): solver convergence parameters have not been set");
   }

   if ( a_homogeneous_bcs ) {
      copyToHypreVector(a_rhs, m_b);
   }
   else {
      LevelData<FArrayBox> rhs;
      rhs.define(a_rhs);

      for (DataIterator dit(rhs.dataIterator()); dit.ok(); ++dit) {
         rhs[dit] += m_rhs_from_bc[dit];
      }

      copyToHypreVector(rhs, m_b);
   }

   copyToHypreVector(a_solution, m_x);

   /* This is a collective call finalizing the vector assembly.
      The vectors are now ``ready to be used'' */
   CH_START(t_HYPRE_SStructVectorAssemble);
   HYPRE_SStructVectorAssemble(m_b);
   HYPRE_SStructVectorAssemble(m_x);
   CH_STOP(t_HYPRE_SStructVectorAssemble);

   if ( m_method == "GMRES" ) {
      AMG_preconditioned_GMRES( m_A, m_A, m_b, m_method_tol, m_method_max_iter,
                                m_precond_tol, m_precond_max_iter, m_method_verbose, m_x );
   }
   else if ( m_method == "AMG" ) {
      AMG( m_A, m_b, m_method_tol, m_method_max_iter, m_method_verbose, m_x );
   }
   else if ( m_method == "ILU" ) {
      ILU( m_A, m_b, m_method_tol, m_method_max_iter, m_method_verbose, m_x );
   }
   else if ( m_method == "MGR" ) {
      MGR( m_A, m_b, m_method_tol, m_method_max_iter, m_method_verbose, m_x );
   }
   else {
      MayDay::Error("MBHypreSolver::solve(): Unknown method, only GMRES, AMG, ILU or MGR recognized");
   }

   copyFromHypreVector(m_x, a_solution);

   a_solution.exchange();
}



void
MBHypreSolver::dumpMatrix(const string&  a_file_name) const
{
   HYPRE_SStructMatrixPrint(a_file_name.c_str(), m_A, 0);
}



void
MBHypreSolver::initializeMatrix()
{
   HYPRE_SStructMatrixDestroy(m_A);  // Just in case
   HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, m_A_graph, &m_A);
   HYPRE_SStructMatrixSetObjectType(m_A, m_hypre_object_type);

   HYPRE_SStructMatrixInitialize(m_A);

   m_matrix_initialized = true;
}



void
MBHypreSolver::finalizeMatrix()
{
   HYPRE_SStructMatrixAssemble(m_A);
   m_matrix_finalized = true;

   if ( m_method == "AMG" ) {
      AMGSetup(m_A);
   }
   else if ( m_method == "ILU" ) {
      ILUSetup(m_A);
   }
   else if ( m_method == "MGR" ) {
      MGRSetup(m_A);
   }

   //   dumpMatrix("HYPRE_MATRIX.A");
}



void
MBHypreSolver::constructMatrixGeneral( LevelData<FArrayBox>&  a_alpha_coefficient,
                                       LevelData<FluxBox>&    a_tensor_coefficient,
                                       LevelData<FArrayBox>&  a_beta_coefficient,
                                       const EllipticOpBC&    a_bc )
{
   bool fourth_order = (m_discretization_order == 4);

   if ( fourth_order && !(a_tensor_coefficient.ghostVect() >= IntVect::Unit) ) {
      MayDay::Error("MBHypreSolver::constructMatrixGeneral(): Fourth-order solve requires tensor coefficient with one transverse ghost cell");
   } 

   constructHypreMatrix(a_alpha_coefficient, a_tensor_coefficient, a_beta_coefficient, a_bc,
                        m_A_stencil_values, m_A_diagonal_offset, m_A_unstructured_coupling,
                        fourth_order, m_rhs_from_bc);
}


void
MBHypreSolver::constructMatrixBlockGeneral( const int              a_block_row,
                                            const int              a_block_column,
                                            LevelData<FArrayBox>&  a_alpha_coefficient,
                                            LevelData<FluxBox>&    a_tensor_coefficient,
                                            LevelData<FArrayBox>&  a_beta_coefficient,
                                            const EllipticOpBC&    a_bc )
{
   bool fourth_order = (m_discretization_order == 4);

   if ( fourth_order && !(a_tensor_coefficient.ghostVect() >= IntVect::Unit) ) {
      MayDay::Error("MBHypreSolver::constructMatrixBlockGeneral(): Fourth-order solve requires tensor coefficient with one transverse ghost cell");
   } 

   constructHypreMatrixBlock(a_block_row, a_block_column, a_alpha_coefficient, a_tensor_coefficient,
                             a_beta_coefficient, a_bc, m_A_graph, m_A_stencil_values, m_A_diagonal_offset,
                             m_A_unstructured_coupling, fourth_order, m_A, m_rhs_from_bc);
}


void
MBHypreSolver::createHypreData()
{
   if (m_hypre_allocated) {
      MayDay::Error( "MBHypreData:createHypreData:: data already allocated; deallocate first!" );
   }

   m_hypre_object_type = HYPRE_PARCSR;

   IntVect stencil_box_lo(IntVect::Zero);
   IntVect stencil_box_hi;
   int radius;
   if (m_discretization_order == 4) {
      stencil_box_hi = 4*IntVect::Unit;
      radius = 2;
   }
   else {
      stencil_box_hi = 2*IntVect::Unit;
      radius = 1;
   }
   Box stencil_box(stencil_box_lo, stencil_box_hi);
          
   // This is a scratch container for the stencil couplings at a grid cell
   m_A_stencil_values.define(stencil_box, 1);

   const DisjointBoxLayout& grids = m_geometry.grids();

   int num_blocks = m_coord_sys_ptr->numBlocks();

   // Create the grid

   {
      HYPRE_SStructGridCreate(MPI_COMM_WORLD, SpaceDim, num_blocks, &m_grid);

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& box = grids[dit];
         int block_number = m_coord_sys_ptr->whichBlock(box);
         IntVect lower(box.loVect());
         IntVect upper(box.hiVect());

         HYPRE_SStructGridSetExtents(m_grid, block_number, lower.dataPtr(), upper.dataPtr());
      }

      HYPRE_SStructVariable* vartypes = new HYPRE_SStructVariable[m_nvar];
      for (int var=0; var<m_nvar; ++var) {
         vartypes[var] = HYPRE_SSTRUCT_VARIABLE_CELL;
      }
      
      for (int block=0; block<num_blocks; ++block) {
         HYPRE_SStructGridSetVariables(m_grid, block, m_nvar, vartypes);
      }

      delete [] vartypes;

      // Specify periodicity.  This is only applicable for single block geometry,
      // since there doesn't seem to be any support for individually periodic blocks
      // in MultiBlockLevelGeom.

      if (num_blocks == 1) {

         const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries 
            = m_coord_sys_ptr->boundaries();

         const ProblemDomain& domain = m_geometry.gridsFull().physDomain();
         IntVect periodic = domain.size();
         int block_number = 0;

         for (int dir=0; dir<SpaceDim; dir++) {
            bool this_dir_periodic = domain.isPeriodic(dir);
            if ( !this_dir_periodic || 
                 (block_boundaries[block_number][dir].isInterface() && m_mblex_potential_Ptr != NULL) ||
                 (block_boundaries[block_number][dir + SpaceDim].isInterface() && m_mblex_potential_Ptr != NULL) ) {
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
                     MayDay::Error( "MBHypreSolver::createHypreData(): Periodic direction must be a power of 2!" );
                  }
               }
            }
         }

         HYPRE_SStructGridSetPeriodic(m_grid, block_number, periodic.dataPtr());
      }

      int num_ghost[] = {radius, radius, radius, radius, radius, radius};

      HYPRE_SStructGridSetNumGhost(m_grid, num_ghost);

      HYPRE_SStructGridAssemble(m_grid);
   }

   // Create the graph
   createGraph(stencil_box, radius);

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


void
MBHypreSolver::createGraph( const Box& a_stencil_box,
                            const int  a_radius )
{
   // Define the discretization stencil for the system matrix A

   {
      /* Define the geometry of the stencil. Each represents a
         relative offset (in the index space). */

      int stencil_size = a_stencil_box.numPts();

      int** offsets = new int*[stencil_size];
      int n = -1;
      for (BoxIterator bit(a_stencil_box - a_radius*IntVect::Unit); bit.ok(); ++bit) {
         n++;
         offsets[n] = new int[SpaceDim];
         IntVect iv = bit();
         for (int m=0; m<SpaceDim; ++m) {
            offsets[n][m] = iv[m];
         }
      }

      m_A_diagonal_offset = (stencil_size - 1) / 2;

      m_A_stencil = new HYPRE_SStructStencil[m_nvar];

      for (int row_var=0; row_var<m_nvar; ++row_var) {

         /* Create an empty stencil object */
         HYPRE_SStructStencilCreate(SpaceDim, stencil_size*m_nvar, &(m_A_stencil[row_var]));

         int block_entry_start = 0;
         for (int column_var=0; column_var<m_nvar; ++column_var) {
            /* Assign each of the stencil entries */
            for (int entry = 0; entry < stencil_size; entry++) {
               HYPRE_SStructStencilSetEntry(m_A_stencil[row_var], entry+block_entry_start,
                                            offsets[entry], column_var);
            }
            block_entry_start += stencil_size;
         }
      }

      for (int n=0; n<stencil_size; ++n) {
         delete [] offsets[n];
      }
      delete [] offsets;
   }

   // Set up the graph

   {
      HYPRE_SStructGraphCreate(MPI_COMM_WORLD, m_grid, &m_A_graph);

      HYPRE_SStructGraphSetObjectType(m_A_graph, m_hypre_object_type);

      int num_blocks = m_coord_sys_ptr->numBlocks();

      for (int var=0; var<m_nvar; ++var) {
         for (int block=0; block<num_blocks; ++block) {
            HYPRE_SStructGraphSetStencil(m_A_graph, block, var, m_A_stencil[var]);
         }
      }

      getUnstructuredCouplings(a_radius, m_A_unstructured_coupling);

      for (int row_var=0; row_var<m_nvar; ++row_var) {
         for (int column_var=0; column_var<m_nvar; ++column_var) {
            addUnstructuredGraphEntries(a_radius, m_A_unstructured_coupling, 
                                        row_var, column_var, m_A_graph);
         }
      }

      HYPRE_SStructGraphAssemble(m_A_graph);
   }
}

void
MBHypreSolver::destroyHypreData()
{
   if (m_hypre_allocated) {
      if (m_AMG_solver_allocated) {
         HYPRE_BoomerAMGDestroy(m_par_AMG_solver);
         m_AMG_solver_allocated = false;
      }
      if (m_ILU_solver_allocated) {
         HYPRE_ILUDestroy(m_par_ILU_solver);
         m_ILU_solver_allocated = false;
      }
      if (m_MGR_solver_allocated) { 
         if(m_mgr_num_cindexes)
            hypre_TFree(m_mgr_num_cindexes, HYPRE_MEMORY_HOST);
         m_mgr_num_cindexes = NULL;

         if(m_mgr_cindexes) {
            for (int i=0; i<m_mgr_nlevels; i++) {
               if(m_mgr_cindexes[i])
                  hypre_TFree(m_mgr_cindexes[i], HYPRE_MEMORY_HOST);
            }
            hypre_TFree(m_mgr_cindexes, HYPRE_MEMORY_HOST);
            m_mgr_cindexes = NULL;
         }

         if (m_mgr_amg_solver) HYPRE_BoomerAMGDestroy(m_mgr_amg_solver);
         HYPRE_MGRDestroy(m_par_MGR_solver);
         m_MGR_solver_allocated = false;
      }

      if (m_A) {
         HYPRE_SStructMatrixDestroy(m_A);
         m_A = NULL;
      }
      HYPRE_SStructVectorDestroy(m_x);
      HYPRE_SStructVectorDestroy(m_b);
      HYPRE_SStructGraphDestroy(m_A_graph);
      for (int var=0; var<m_nvar; ++var) {
         HYPRE_SStructStencilDestroy(m_A_stencil[var]);
      }
      delete [] m_A_stencil;
      HYPRE_SStructGridDestroy(m_grid);

      m_hypre_allocated = false;
   }
}


int
MBHypreSolver::findHypreEntry( const int          a_block_column,
                               const Box&         a_stencil_box,
                               const IntVectSet&  a_unstructured_ivs,
                               const IntVect&     a_iv ) const
{
   /*
     The Hypre matrix coefficients corresponding to the equation centered at cell a_iv are
     identified by an entry number, which is returned by this function.  Regular "stencil"
     couplings to cells in a_stencil_box are enumerated first (in the order with which a 
     BoxIterator built on a_stencil_box traverses cells), followed by the unstructured couplings
     to cells in other blocks contained in a_unstructured_ivs.
   */

   bool found_entry = false;

   // Search the structured entries
   
   int entry = a_block_column * a_stencil_box.numPts();

   for (BoxIterator bit(a_stencil_box); bit.ok(); ++bit) {
      if (bit() == a_iv) {
         found_entry = true;
         break;
      }
      entry++;
   }

   // Search the unstructured entries

   if (!found_entry) {

      entry = m_nvar * a_stencil_box.numPts() + a_block_column * a_unstructured_ivs.numPts();

      for (IVSIterator ivit(a_unstructured_ivs); ivit.ok(); ++ivit) {
         if (a_iv == ivit()) {
            found_entry = true;
            break;
         }
         entry++;
      }

      CH_assert(found_entry);
   }

   return entry;
}



void
MBHypreSolver::addUnstructuredGraphEntries( const int                                 a_radius,
                                            const LayoutData< BaseFab<IntVectSet> >&  a_unstructured_coupling,
                                            const int                                 a_var,
                                            const int                                 a_to_var,
                                            HYPRE_SStructGraph&                       a_graph ) const
{
   if (m_mblex_potential_Ptr) {

      const Box stencil_offsets(IntVect(D_DECL6(-a_radius,-a_radius,-a_radius,-a_radius,-a_radius,-a_radius)),
                                IntVect(D_DECL6(a_radius,a_radius,a_radius,a_radius,a_radius,a_radius)));

      const DisjointBoxLayout & grids = m_geometry.grids();

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         int dst_block_number = m_coord_sys_ptr->whichBlock(grids[dit]);

         // Find all cells that may be coupled through a block interface
         IntVectSet dst_ivs = getInterBlockCoupledCells(dst_block_number, a_radius, grids[dit]);

         for (IVSIterator it(dst_ivs); it.ok(); ++it) {
            IntVect iv = it();
            Box stencil_box = stencil_offsets + iv;

            const IntVectSet& unstructured_couplings = a_unstructured_coupling[dit](iv);
            for (IVSIterator sivsit(unstructured_couplings); sivsit.ok(); ++sivsit) {
               IntVect siv = sivsit();
               int src_block_number = m_coord_sys_ptr->whichBlock(siv);
               CH_assert(src_block_number >= 0 && src_block_number < m_geometry.coordSysPtr()->numBlocks());
               HYPRE_SStructGraphAddEntries(a_graph, dst_block_number, iv.dataPtr(), a_var,
                                            src_block_number, siv.dataPtr(), a_to_var);
            }
         }
      }
   }
}



void
MBHypreSolver::addUnstructuredMatrixEntries( const int                                 a_block_row,
                                             const int                                 a_block_column,
                                             const LevelData<FArrayBox>&               a_alpha_coefficient,
                                             const LevelData<FluxBox>&                 a_tensor_coefficient,
                                             const EllipticOpBC&                       a_bc,
                                             FArrayBox&                                a_stencil_values,
                                             const bool                                a_fourth_order,
                                             const LayoutData< BaseFab<IntVectSet> >&  a_unstructured_coupling,
                                             Vector< Vector<CoDim1Stencil> >&          a_codim1_stencils,
                                             Vector< Vector<CoDim2Stencil> >&          a_codim2_stencils,
                                             HYPRE_SStructMatrix&                      a_matrix ) const
{
   if (m_mblex_potential_Ptr) {

      const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = m_coord_sys_ptr->boundaries();

      Box stencil_offsets = a_stencil_values.box();
      int radius = (stencil_offsets.size(0)-1)/2;
      stencil_offsets -= radius*IntVect::Unit;

      const DisjointBoxLayout & grids = m_geometry.grids();

      const LayoutData< IntVectSet >& ghosts = m_mblex_potential_Ptr->ghostCells();
      const LayoutData< RefCountedPtr< IVSFAB<MBStencil> > >& stencil = m_mblex_potential_Ptr->stencils();

      BlockBaseRegister<BaseFab<SparseCoupling> > blockRegister(m_coord_sys_ptr, grids, 0);

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box & box = grids[dit];

         int block_number = m_coord_sys_ptr->whichBlock(box);

         const NewCoordSys* block_coord_sys = m_coord_sys_ptr->getCoordSys(block_number);
         RealVect dx = block_coord_sys->dx();

         const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];
         bool lo_block_interface[SpaceDim];
         bool hi_block_interface[SpaceDim];
         for (int dir=0; dir<SpaceDim; ++dir) {
            lo_block_interface[dir] = this_block_boundaries[dir].isInterface();
            hi_block_interface[dir] = this_block_boundaries[dir + SpaceDim].isInterface();
         }

         const FArrayBox* alphaPtr = NULL;
         if ( a_alpha_coefficient.isDefined() ) {
            alphaPtr = &a_alpha_coefficient[dit];
         }

         const Box& domainBox = (m_coord_sys_ptr->mappingBlocks())[block_number];
         const IntVectSet& this_ghosts_fab = ghosts[dit];
         if ( stencil[dit] != NULL ) {
            const IVSFAB<MBStencil>& this_stencil = *stencil[dit];

            BaseFab<SparseCoupling> data[SpaceDim];
            for( int dir = 0; dir < SpaceDim; ++dir ) {
               data[dir].define(surroundingNodes(box, dir), 1);
            }

            // Here, we loop over all cells whose stencil is affected by a block boundary
            IntVectSet interblock_ivs = getInterBlockCoupledCells(block_number, radius, box);
            for (IVSIterator ivsit(interblock_ivs); ivsit.ok(); ++ivsit) {
               IntVect iv = ivsit();
               Box stencil_box = stencil_offsets + iv;
               Box stencil_box_valid = stencil_box & domainBox;
               IntVectSet extra_block_ghosts = stencil_box & this_ghosts_fab;

               double alpha = alphaPtr? alphaPtr->operator()(iv): 1.;

               // Get all of the indices that are sparsely coupled to this cell
               const IntVectSet& unstructured_couplings = a_unstructured_coupling[dit](iv);

               for (int dir=0; dir<SpaceDim; ++dir) {
                  bool at_lo_block_interface = lo_block_interface[dir] && (iv[dir] == domainBox.smallEnd(dir));
                  bool at_hi_block_interface = hi_block_interface[dir] && (iv[dir] == domainBox.bigEnd(dir));
                  bool force_codim2_condense = lo_block_interface[dir] || hi_block_interface[dir];
               
                  for (SideIterator sit; sit.ok(); ++sit) {
                     Side::LoHiSide side = sit();

                     IntVect iv_face;

                     // Halve the contributions on face contained in block interfaces to effect averaging
                     double fac;
                     if ( (at_lo_block_interface && side == Side::LoHiSide::Lo) ||
                          (at_hi_block_interface && side == Side::LoHiSide::Hi) ) {
                        fac = m_flux_average[dir]? 0.5: 1.0;
                        iv_face = iv;
                        iv_face.shift(dir,side);
                     }
                     else {
                        fac = 1.;
                     }

                     for (int dir2=0; dir2<SpaceDim; ++dir2) {

                        a_stencil_values.setVal(0.);

                        accumStencilMatrixEntries(iv, dir, side, dir2, a_tensor_coefficient[dit],
                                                  dx, a_fourth_order, a_stencil_values);

                        FArrayBox dummy;
                        modifyStencilForBCs( a_codim1_stencils[block_number], a_codim2_stencils[block_number],
                                             iv, a_stencil_values, dummy, false, force_codim2_condense );

                        BoxIterator bit(stencil_box);
                        for (bit.begin(); bit.ok(); ++bit) {
                           IntVect stencil_box_iv = bit();

                           double s = fac * a_stencil_values( stencil_box_iv - iv + radius*IntVect::Unit );

                           if ( s != 0. ) {

                              if ( stencil_box_valid.contains(stencil_box_iv) ) {

                                 // If the current side of cell iv is on a block interface, add the stencil
                                 // contribution to the BlockRegister data.  The minus sign accounts for the
                                 // normal component, which is negated on the adacent cell.  Also, since we
                                 // have not yet accounted for the regular stencil contributions, do that now.
                                 if ( (at_lo_block_interface && side == Side::LoHiSide::Lo) ||
                                      (at_hi_block_interface && side == Side::LoHiSide::Hi) ) {
                                    data[dir](iv_face).add(stencil_box_iv, -s);

                                    int entry = findHypreEntry(a_block_column, stencil_box, 
                                                               unstructured_couplings, stencil_box_iv);
                                    double value = s * alpha;

                                    int num_entries = 1;
                                    HYPRE_SStructMatrixAddToValues(a_matrix, block_number, iv.dataPtr(),
                                                                   a_block_row, num_entries, &entry, &value);
                                 }
                              }
                              else if ( extra_block_ghosts.contains(stencil_box_iv) ) {

                                 const MBStencil& elements = this_stencil(stencil_box_iv, 0);
                                 int num_unstructured_entries = elements.size();
                                 double * values = new double[num_unstructured_entries];
                                 int * entries = new int[num_unstructured_entries];

                                 MBStencilIterator stit(elements);
                                 int k = 0;
                                 for (stit.begin(); stit.ok(); ++stit) {
                                    const MBStencilElement& stencilElement = stit();
                                    IntVect interp_cell = stencilElement.cell();

                                    entries[k] = findHypreEntry(a_block_column, stencil_box,
                                                                unstructured_couplings, interp_cell);
                                    values[k] = s * stencilElement.weight();

                                    // If we are on a block interface, add the stencil contribution to the
                                    // BlockRegister data.  The minus sign accounts for the normal component,
                                    // which is negated on the adacent cell.
                                    if ( (at_lo_block_interface && side == Side::LoHiSide::Lo) ||
                                         (at_hi_block_interface && side == Side::LoHiSide::Hi) ) {
                                       data[dir](iv_face).add(interp_cell, -values[k]);
                                    }

                                    values[k] *= alpha;
                                    k++;
                                 }

                                 HYPRE_SStructMatrixAddToValues(a_matrix, block_number, iv.dataPtr(),
                                                                a_block_row, num_unstructured_entries,
                                                                entries, values);

                                 delete [] entries;
                                 delete [] values;
                              }
                              else {

                                 for (IVSIterator givsit(extra_block_ghosts); givsit.ok(); ++givsit) {
                                    cout << iv << " " << stencil_box_iv << " " << givsit() << endl;
                                 }

                                 // If we've landed here, then stencil_box_iv is neither a valid cell nor an
                                 // extrablock ghost cell (relative to the center cell iv).  The only remaining
                                 // possibility is that it's a ghost cell at a codim = 1 or codim = 2 physical
                                 // boundary, in which case the stencil modification to account for boundary
                                 // conditions should have zeroed out the corresponding stencil entry (i.e.,
                                 // s = 0).  Something has therefore gone wrong and needs to be investigated.

                                 MayDay::Error("MBHypreSolver::addUnstructuredMatrixEntries(): Invalid stencil_box_iv");
                              }
                           }
                        }
                     } // dir2
                  }  // sit
               }  // dir
            }  // ivsit

            // Store the contributions at block interfaces
            for (int dir = 0; dir < SpaceDim; dir++) {
               if ( m_flux_average[dir] ) {
                  for (SideIterator sit; sit.ok(); ++sit) {
                     Side::LoHiSide side = sit();
                     if (blockRegister.hasInterface(dit(), dir, side)) {
                        blockRegister.store(data[dir], dit(), dir, side);
                     }
                  }
               }
            }
         }
      }

      // Exchange the block interface contributions and accumulate them to the matrix

      blockRegister.exchange();

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         int block_number = m_coord_sys_ptr->whichBlock(grids[dit]);

         const FArrayBox* alphaPtr = NULL;
         if ( a_alpha_coefficient.isDefined() ) {
            alphaPtr = &a_alpha_coefficient[dit];
         }

         for (SideIterator sit; sit.ok(); ++sit) {
            Side::LoHiSide side = sit();
            for (int dir = 0; dir < SpaceDim; dir++) {
               if ( m_flux_average[dir] && blockRegister.hasInterface(dit(), dir, side) ) {
                  Box faceBox = adjCellBox(grids[dit], dir, side, 1);
                  // if Lo, then shift +1; if Hi, then shift -1
                  faceBox.shiftHalf(dir, -sign(side));
                  BaseFab<SparseCoupling> otherFab(faceBox, 1);
                  blockRegister.fill(otherFab, dit(), dir, side, side);
                  BoxIterator bit(faceBox);
                  for (bit.begin(); bit.ok(); ++bit) {
                     IntVect iv_face = bit();
                     IntVect iv(iv_face);
                     iv.shift(dir,-side);
                     Box stencil_box = stencil_offsets + iv;
                     const IntVectSet& unstructured_couplings = a_unstructured_coupling[dit](iv);

                     double alpha = alphaPtr? alphaPtr->operator()(iv): 1.;

                     SparseCoupling& coupling = otherFab(iv_face,0);
                     int num_entries = coupling.size();
                     double * values = new double[num_entries];
                     int * entries = new int[num_entries];

                     int k = 0;
                     for (SparseCouplingIterator it(coupling); it.ok(); ++it) {
                        IntVect index = coupling[it()];
                        entries[k] = findHypreEntry(a_block_column, stencil_box, unstructured_couplings, index);
                        values[k] = coupling.weight(index) * alpha;
                        k++;
                     }

                     HYPRE_SStructMatrixAddToValues(a_matrix, block_number, iv.dataPtr(),
                                                    a_block_row, num_entries, entries, values);

                     delete [] entries;
                     delete [] values;

                  }  // bit
               }  // if interface
            }  // dir
         }  // sit
      }  // dit
   }
}


void
MBHypreSolver::constructHypreMatrix( LevelData<FArrayBox>&               a_alpha_coefficient, 
                                     LevelData<FluxBox>&                 a_tensor_coefficient,
                                     LevelData<FArrayBox>&               a_beta_coefficient,
                                     const EllipticOpBC&                 a_bc,
                                     FArrayBox&                          a_stencil_values,
                                     const int                           a_diagonal_offset,
                                     LayoutData< BaseFab<IntVectSet> >&  a_unstructured_coupling,
                                     const bool                          a_fourth_order,
                                     LevelData<FArrayBox>&               a_rhs_from_bc )
{
   initializeMatrix();

   constructHypreMatrixBlock(0, 0, a_alpha_coefficient, a_tensor_coefficient, a_beta_coefficient,
                             a_bc, m_A_graph, a_stencil_values, a_diagonal_offset,
                             a_unstructured_coupling, a_fourth_order, m_A, a_rhs_from_bc );

   finalizeMatrix();
   
#if 0

   testMatrixConstruct( a_alpha_coefficient, 
                        a_tensor_coefficient,
                        a_beta_coefficient,
                        block_boundaries,
                        codim1_stencils,
                        codim2_stencils,
                        a_stencil_values,
                        a_fourth_order,
                        a_rhs_from_bc );
#endif
}
   
void
MBHypreSolver::constructHypreMatrixBlock( const int                           a_block_row,
                                          const int                           a_block_column,
                                          LevelData<FArrayBox>&               a_alpha_coefficient, 
                                          LevelData<FluxBox>&                 a_tensor_coefficient,
                                          LevelData<FArrayBox>&               a_beta_coefficient,
                                          const EllipticOpBC&                 a_bc,
                                          HYPRE_SStructGraph&                 a_graph,
                                          FArrayBox&                          a_stencil_values,
                                          const int                           a_diagonal_offset,
                                          LayoutData< BaseFab<IntVectSet> >&  a_unstructured_coupling,
                                          const bool                          a_fourth_order,
                                          HYPRE_SStructMatrix&                a_matrix,
                                          LevelData<FArrayBox>&               a_rhs_from_bc )
{
   /* Set the matrix coefficients for the (a_block_row,a_block_column) block.  Each processor
      assigns coefficients for the boxes in the grid that it owns. Note that the coefficients
      associated with each stencil entry may vary from grid point to grid point if desired.
      Here, we first set the same stencil entries for each grid point.  Then we make modifications
      to grid points near the boundary. */

   // Make sure the coefficient ghost cell values are set.
   // The ghost cell values at physical boundaries aren't used.
   a_tensor_coefficient.exchange();

   int stencil_size = a_stencil_values.box().numPts();
   int block_entry_offset = a_block_column * stencil_size;
   int* entries = new int[stencil_size];
   for (int n=0; n<stencil_size; n++) entries[n] = n + block_entry_offset;

   const DisjointBoxLayout& grids = m_geometry.grids();
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = m_coord_sys_ptr->boundaries();

   LevelData<FArrayBox> diagonal_plot(grids, 1, IntVect::Zero);

   Vector< Vector<CoDim1Stencil> > codim1_stencils;
   Vector< Vector<CoDim2Stencil> > codim2_stencils;

   constructBoundaryStencils(a_fourth_order, a_bc, codim1_stencils, codim2_stencils );

#if 0
   // This is a newer implementation in which the loops over directions and
   // sides is outside the loop over cells
   
   // Set the intra-block matrix entries

   for (DataIterator dit(grids); dit.ok(); ++dit) {
     const Box& box = grids[dit];
     IntVect lower(box.loVect());
     IntVect upper(box.hiVect());

     diagonal_plot[dit].setVal(0.);

     int block_number = m_coord_sys_ptr->whichBlock(box);
     const RealVect& dx = (m_coord_sys_ptr->getCoordSys(block_number))->dx();
     const Box& domain_box = (m_coord_sys_ptr->mappingBlocks())[block_number];

     FluxBox& this_coef = a_tensor_coefficient[dit];

     const FArrayBox* alpha = NULL;
     if ( a_alpha_coefficient.isDefined() ) {
        alpha = &a_alpha_coefficient[dit];
     }

     const FArrayBox* beta = NULL;
     bool add_beta = a_beta_coefficient.isDefined();
     if ( add_beta ) {
        beta = &a_beta_coefficient[dit];
     }

     int nvalues = stencil_size * box.numPts();
     double* values = new double[nvalues];
     for (int i=0; i<nvalues; ++i) values[i] = 0.;

     const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];

     for (int dir=0; dir<SpaceDim; ++dir) {
        bool lo_block_interface;
        bool hi_block_interface;
        if (m_mblex_potential_Ptr) {
           lo_block_interface = this_block_boundaries[dir].isInterface();
           hi_block_interface = this_block_boundaries[dir + SpaceDim].isInterface();
        }
        else {
           lo_block_interface = false;
           hi_block_interface = false;
        }
        bool force_codim2_condense = lo_block_interface || hi_block_interface;

        for (SideIterator sit; sit.ok(); ++sit) {
           Side::LoHiSide side = sit();

           for (int dir2=0; dir2<SpaceDim; ++dir2) {

              int ii = 0;
              for (BoxIterator bit(box); bit.ok(); ++bit) {
                 IntVect iv = bit();

                 a_stencil_values.setVal(0.);

                 // Compute the stencil couplings to cell iv without regard to boundary conditions, skipping
                 // the contributions due to cell faces on block interfaces, which are handled later by
                 // addUnstructuredMatrixEntries().

                 bool at_lo_interface = (lo_block_interface && iv[dir] == domain_box.smallEnd(dir));
                 bool at_hi_interface = (hi_block_interface && iv[dir] == domain_box.bigEnd(dir));
                 bool update_rhs_from_bc_only = (at_lo_interface && side == Side::LoHiSide::Lo)
                                             || (at_hi_interface && side == Side::LoHiSide::Hi);

                 accumStencilMatrixEntries(iv, dir, side, dir2, this_coef, dx,
                                           a_fourth_order, a_stencil_values);

                 modifyStencilForBCs( codim1_stencils[block_number], codim2_stencils[block_number],
                                      iv, a_stencil_values, a_rhs_from_bc[dit],
                                      update_rhs_from_bc_only, force_codim2_condense );

                 if ( alpha ) {
                    double alpha_val = alpha->operator()(iv);
                    a_stencil_values *= alpha_val;
                    a_rhs_from_bc[dit](iv,0) *= alpha_val;
                 }

                 double* this_stencil_values = a_stencil_values.dataPtr(0);

                 if ( add_beta ) {
                    this_stencil_values[a_diagonal_offset] += beta->operator()(iv);
                 }
                 diagonal_plot[dit](iv,0) += this_stencil_values[a_diagonal_offset];

                 for (int jj=0; jj<stencil_size; jj++) {
                    values[ii+jj] += this_stencil_values[jj];
                 }

                 ii += stencil_size;

              } // box loop
              add_beta = false;
           }
        }
     }

     HYPRE_SStructMatrixAddToBoxValues(a_matrix, block_number, lower.dataPtr(), upper.dataPtr(),
                                       a_block_row, stencil_size, entries, values);

     delete [] values;
   }
#else
   // This is the original implementation in which the loops over directions and
   // sides is inside the loop over cells

   LevelData<FArrayBox> structured_values(grids, stencil_size, IntVect::Zero);
   FArrayBox tmp_stencil_values(a_stencil_values.box(), a_stencil_values.nComp());
   double * stencil_values = a_stencil_values.dataPtr(0);

   setZero(a_rhs_from_bc);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
     const Box & box = grids[dit];

     int block_number = m_coord_sys_ptr->whichBlock(box);
     const NewCoordSys& block_coord_sys = *(m_coord_sys_ptr->getCoordSys(block_number));
     RealVect dx = block_coord_sys.dx();

     const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];
     bool lo_block_interface[SpaceDim];
     bool hi_block_interface[SpaceDim];
     if (m_mblex_potential_Ptr) {
        for (int dir=0; dir<SpaceDim; ++dir) {
           lo_block_interface[dir] = this_block_boundaries[dir].isInterface();
           hi_block_interface[dir] = this_block_boundaries[dir + SpaceDim].isInterface();
        }
     }
     else {
        for (int dir=0; dir<SpaceDim; ++dir) {
           lo_block_interface[dir] = false;
           hi_block_interface[dir] = false;
        }
     }

     const Box& domainBox = (m_coord_sys_ptr->mappingBlocks())[block_number];

     FluxBox& this_coef = a_tensor_coefficient[dit];

     const FArrayBox* alpha = NULL;
     if ( a_alpha_coefficient.isDefined() ) {
        alpha = &a_alpha_coefficient[dit];
     }

     const FArrayBox* beta = NULL;
     if ( a_beta_coefficient.isDefined() ) {
        beta = &a_beta_coefficient[dit];
     }

     // Set the intra-block matrix entries

     BoxIterator bit(box);
     for (bit.begin();bit.ok();++bit) {
       IntVect iv = bit();

       a_stencil_values.setVal(0.);

       // Compute the stencil couplings to cell iv without regard to boundary conditions, skipping
       // the contributions due to cell faces on block interfaces, which are handled later by
       // addUnstructuredMatrixEntries().

       for (int dir=0; dir<SpaceDim; ++dir) {
          bool at_lo_interface = (lo_block_interface[dir] && iv[dir] == domainBox.smallEnd(dir));
          bool at_hi_interface = (hi_block_interface[dir] && iv[dir] == domainBox.bigEnd(dir));
          bool force_codim2_condense = lo_block_interface[dir] || hi_block_interface[dir];

          for (SideIterator sit; sit.ok(); ++sit) {
             Side::LoHiSide side = sit();

             bool update_rhs_from_bc_only = (at_lo_interface && side == Side::LoHiSide::Lo)
                || (at_hi_interface && side == Side::LoHiSide::Hi);

             for (int dir2=0; dir2<SpaceDim; ++dir2) {

                tmp_stencil_values.setVal(0.);

                accumStencilMatrixEntries(iv, dir, side, dir2, this_coef, dx,
                                          a_fourth_order, tmp_stencil_values);

                modifyStencilForBCs( codim1_stencils[block_number], codim2_stencils[block_number],
                                     iv, tmp_stencil_values, a_rhs_from_bc[dit],
                                     update_rhs_from_bc_only, force_codim2_condense );

                a_stencil_values += tmp_stencil_values;
             }
          }
       }

       if ( alpha ) {
          a_stencil_values *= alpha->operator()(iv);
          a_rhs_from_bc[dit](iv,0) *= alpha->operator()(iv);
       }

       if ( beta ) {
          stencil_values[a_diagonal_offset] += beta->operator()(iv);
       }

       for (int jj=0; jj<stencil_size; jj++) {
          structured_values[dit](iv,jj) = stencil_values[jj];
       }
     }
   }

   for (DataIterator dit(grids); dit.ok(); ++dit) {
     const Box & box = grids[dit];
     IntVect lower(box.loVect());
     IntVect upper(box.hiVect());

     int block_number = m_coord_sys_ptr->whichBlock(box);

     // Set the intra-block matrix entries

     int nvalues  = stencil_size * box.numPts();
     double * values = new double[nvalues];

     int ii = 0;
     int diagonal_index = 0;

     BoxIterator bit(box);
     for (bit.begin();bit.ok();++bit) {
       IntVect iv = bit();

       for (int jj=0; jj<stencil_size; jj++) {
          stencil_values[jj] = structured_values[dit](iv,jj);
       }

       for (int jj=0; jj<stencil_size; jj++) {
          values[ii+jj] = stencil_values[jj];
       }

       ii += stencil_size;
       diagonal_index++;

     }

     HYPRE_SStructMatrixSetBoxValues(a_matrix, block_number, lower.dataPtr(), upper.dataPtr(),
                                     a_block_row, stencil_size, entries, values);

     delete [] values;
   }

#endif

   delete [] entries;

   addUnstructuredMatrixEntries(a_block_row, a_block_column, a_alpha_coefficient, a_tensor_coefficient,
                                a_bc, a_stencil_values, a_fourth_order, a_unstructured_coupling,
                                codim1_stencils, codim2_stencils, a_matrix);
}


void
MBHypreSolver::extrapGhosts( const EllipticOpBC&   a_bc,
                             const bool            a_fourth_order,
                             LevelData<FArrayBox>& a_data ) const
{
   Vector< Vector<CoDim1Stencil> > codim1_stencils;
   Vector< Vector<CoDim2Stencil> > codim2_stencils;

   constructBoundaryStencils(a_fourth_order, a_bc, codim1_stencils, codim2_stencils );

   bool extrapolate_from_interior = true;
   bool include_bvs = true;
   accumPhysicalGhosts( codim1_stencils, codim2_stencils, extrapolate_from_interior, include_bvs, a_data );
}


void
MBHypreSolver::AMGSetup( const HYPRE_SStructMatrix& a_matrix )
{
   CH_TIME("MBHypreSolver::AMGSetup");

   if ( m_AMG_solver_allocated ) {
      HYPRE_BoomerAMGDestroy(m_par_AMG_solver);
   }
   
   HYPRE_BoomerAMGCreate(&m_par_AMG_solver);

   // Algorithm options
   if( m_AMG_old_setup ) {  // Restores defaults prior to the changes made in Hypre v2.11.2
      HYPRE_BoomerAMGSetOldDefault(m_par_AMG_solver);
      HYPRE_BoomerAMGSetRelaxOrder(m_par_AMG_solver, 1);
   } 
   else {
      HYPRE_BoomerAMGSetCoarsenType(m_par_AMG_solver, m_AMG_coarsen_type);
   }
   HYPRE_BoomerAMGSetStrongThreshold(m_par_AMG_solver, m_AMG_strong_threshold);

   HYPRE_ParCSRMatrix par_A;
   HYPRE_SStructMatrixGetObject(a_matrix, (void **) &par_A);

   HYPRE_ParVector par_b;  // not used, but needs to be passed to HYPRE_BoomerAMGSetup
   HYPRE_ParVector par_x;  // not used, but needs to be passed to HYPRE_BoomerAMGSetup

   HYPRE_BoomerAMGSetup(m_par_AMG_solver, par_A, par_b, par_x);

   m_AMG_solver_allocated = true;
}


void
MBHypreSolver::ILUSetup( const HYPRE_SStructMatrix& a_matrix )
{
   CH_TIME("MBHypreSolver::ILUSetup");
   if ( !m_matrix_finalized ) {
      MayDay::Error("MBHypreSolver::ILUSetup(): Hypre matrix has not been finalized");
   }

   if ( m_ILU_solver_allocated ) {
      HYPRE_ILUDestroy(m_par_ILU_solver);
   }
   
   HYPRE_ILUCreate(&m_par_ILU_solver);

   HYPRE_ILUSetType(m_par_ILU_solver, m_ILU_type);
   HYPRE_ILUSetLevelOfFill(m_par_ILU_solver, m_ILU_fill_level);
   HYPRE_ILUSetDropThreshold(m_par_ILU_solver, m_ILU_drop_threshold);
   HYPRE_ILUSetDropThresholdArray(m_par_ILU_solver, m_ILU_drop_threshold_array);
   HYPRE_ILUSetMaxNnzPerRow(m_par_ILU_solver, m_ILU_max_nnz_per_row);
   HYPRE_ILUSetMaxIter(m_par_ILU_solver, m_ILU_max_iter);
   HYPRE_ILUSetSchurMaxIter(m_par_ILU_solver, m_ILU_max_schur_iter);
   HYPRE_ILUSetTol(m_par_ILU_solver, m_ILU_tol);
   HYPRE_ILUSetPrintLevel(m_par_ILU_solver, m_ILU_print_level);

   HYPRE_ParCSRMatrix par_A;
   HYPRE_SStructMatrixGetObject(a_matrix, (void **) &par_A);

   HYPRE_ParVector par_b;  // not used, but needs to be passed to HYPRE_ILUSetup
   HYPRE_ParVector par_x;  // not used, but needs to be passed to HYPRE_ILUSetup

   HYPRE_ILUSetup(m_par_ILU_solver, par_A, par_b, par_x);

   m_ILU_solver_allocated = true;
}


void
MBHypreSolver::MGRSetup( const HYPRE_SStructMatrix& a_matrix )
{
   CH_TIME("MBHypreSolver::MGRSetup");
   if ( !m_matrix_finalized ) {
      MayDay::Error("MBHypreSolver::MGRSetup(): Hypre matrix has not been finalized");
   }

   if ( m_MGR_solver_allocated ) {
      if(m_mgr_num_cindexes)
         hypre_TFree(m_mgr_num_cindexes, HYPRE_MEMORY_HOST);
      m_mgr_num_cindexes = NULL;

      if(m_mgr_cindexes) {
         for (int i=0; i<m_mgr_nlevels; i++) {
            if(m_mgr_cindexes[i])
               hypre_TFree(m_mgr_cindexes[i], HYPRE_MEMORY_HOST);
         }
         hypre_TFree(m_mgr_cindexes, HYPRE_MEMORY_HOST);
         m_mgr_cindexes = NULL;
      }

      if (m_mgr_amg_solver) HYPRE_BoomerAMGDestroy(m_mgr_amg_solver);
      HYPRE_MGRDestroy(m_par_MGR_solver);
   }
   
   // Expert-only(!) settings
   m_mgr_nlevels = 1;
   m_mgr_num_cindexes = NULL;
   m_mgr_cindexes = NULL;
   HYPRE_Int mgr_bsize = 2;
   HYPRE_Int mgr_non_c_to_f = 1;
   HYPRE_Int mgr_frelax_method = 1;
   HYPRE_Int mgr_relax_type = 0;
   HYPRE_Int mgr_num_relax_sweeps = 2;
   HYPRE_Int mgr_interp_type = 3;
   HYPRE_Int mgr_num_interp_sweeps = 2;
   HYPRE_Int mgr_gsmooth_type = 0;
   HYPRE_Int mgr_num_gsmooth_sweeps = 0;
   HYPRE_Int mgr_restrict_type = 0;
   HYPRE_Int mgr_num_restrict_sweeps = 0;  

   HYPRE_MGRCreate(&m_par_MGR_solver);

   m_mgr_num_cindexes = hypre_CTAlloc(HYPRE_Int,  m_mgr_nlevels, HYPRE_MEMORY_HOST);
   for (int i=0; i<m_mgr_nlevels; i++) {
      /* assume 1 coarse index per level */
      m_mgr_num_cindexes[i] = 1;
   }
   m_mgr_cindexes = hypre_CTAlloc(HYPRE_Int*,  m_mgr_nlevels, HYPRE_MEMORY_HOST);
   for (int i=0; i<m_mgr_nlevels; i++) {
      m_mgr_cindexes[i] = hypre_CTAlloc(HYPRE_Int,  m_mgr_num_cindexes[i], HYPRE_MEMORY_HOST);
   }
   for (int i=0; i<m_mgr_nlevels; i++) {
      /* assume coarse point is at index 0 */
      m_mgr_cindexes[i][0] = m_MGR_cpoint;
   }

   HYPRE_ParCSRMatrix par_A;
   HYPRE_SStructMatrixGetObject(a_matrix, (void **) &par_A);

   /* set MGR data by block */
   HYPRE_BigInt rowstart = hypre_ParCSRMatrixFirstRowIndex(par_A);
   HYPRE_BigInt rowend = hypre_ParCSRMatrixLastRowIndex(par_A);
   HYPRE_Int fsize = (rowend - rowstart + 1)/2 ;
   HYPRE_BigInt next_block = rowstart + fsize;

   HYPRE_BigInt idx_array[2] = {rowstart,next_block};
   HYPRE_MGRSetCpointsByContiguousBlock( m_par_MGR_solver, mgr_bsize, m_mgr_nlevels,
                                         idx_array, m_mgr_num_cindexes, m_mgr_cindexes);

   /* set intermediate coarse grid strategy */
   HYPRE_MGRSetNonCpointsToFpoints(m_par_MGR_solver, mgr_non_c_to_f);
   /* set F relaxation strategy */
   HYPRE_MGRSetFRelaxMethod(m_par_MGR_solver, mgr_frelax_method);
   /* set relax type for single level F-relaxation and post-relaxation */
   HYPRE_MGRSetRelaxType(m_par_MGR_solver, mgr_relax_type);
   HYPRE_MGRSetNumRelaxSweeps(m_par_MGR_solver, mgr_num_relax_sweeps);
   /* set interpolation type */
   HYPRE_MGRSetRestrictType(m_par_MGR_solver, mgr_restrict_type);
   HYPRE_MGRSetNumRestrictSweeps(m_par_MGR_solver, mgr_num_restrict_sweeps);
   HYPRE_MGRSetInterpType(m_par_MGR_solver, mgr_interp_type);
   HYPRE_MGRSetNumInterpSweeps(m_par_MGR_solver, mgr_num_interp_sweeps);
   /* set print level */
   HYPRE_MGRSetPrintLevel(m_par_MGR_solver, m_MGR_print_level);
   /* set max iterations */
   //   HYPRE_MGRSetMaxIter(m_par_MGR_solver, max_iter);
   //   HYPRE_MGRSetTol(m_par_MGR_solver, tol);

   HYPRE_MGRSetGlobalsmoothType(m_par_MGR_solver, mgr_gsmooth_type);
   HYPRE_MGRSetMaxGlobalsmoothIters( m_par_MGR_solver, mgr_num_gsmooth_sweeps );
   //HYPRE_Int num_functions = 1;
   //HYPRE_MGRSetLevelFRelaxNumFunctions(m_par_MGR_solver, &num_functions);
   /* create AMG coarse grid solver */
   HYPRE_BoomerAMGCreate(&m_mgr_amg_solver);
   HYPRE_BoomerAMGSetCGCIts(m_mgr_amg_solver, 1);
   HYPRE_BoomerAMGSetInterpType(m_mgr_amg_solver, 0);
   HYPRE_BoomerAMGSetPostInterpType(m_mgr_amg_solver, 0);
   HYPRE_BoomerAMGSetCoarsenType(m_mgr_amg_solver, 6);
   HYPRE_BoomerAMGSetPMaxElmts(m_mgr_amg_solver, 0);
   HYPRE_BoomerAMGSetCycleType(m_mgr_amg_solver, 1);
   HYPRE_BoomerAMGSetFCycle(m_mgr_amg_solver, 0);
   HYPRE_BoomerAMGSetNumSweeps(m_mgr_amg_solver, 1);
   HYPRE_BoomerAMGSetRelaxType(m_mgr_amg_solver, 3);
   /*
     if (relax_down > -1)
     HYPRE_BoomerAMGSetCycleRelaxType(amg_solver, relax_down, 1);
     if (relax_up > -1)
     HYPRE_BoomerAMGSetCycleRelaxType(amg_solver, relax_up, 2);
     if (relax_coarse > -1)
     HYPRE_BoomerAMGSetCycleRelaxType(amg_solver, relax_coarse, 3);
   */
   HYPRE_BoomerAMGSetRelaxOrder(m_mgr_amg_solver, 1);
   HYPRE_BoomerAMGSetMaxLevels(m_mgr_amg_solver, 25);
   HYPRE_BoomerAMGSetSmoothType(m_mgr_amg_solver, 6);
   HYPRE_BoomerAMGSetSmoothNumSweeps(m_mgr_amg_solver, 1);
   HYPRE_BoomerAMGSetMaxIter(m_mgr_amg_solver, m_MGR_amg_max_iter);
   HYPRE_BoomerAMGSetTol(m_mgr_amg_solver, 0.0);
   HYPRE_BoomerAMGSetPrintLevel(m_mgr_amg_solver, m_MGR_amg_print_level);

   /* set the MGR coarse solver. Comment out to use default Coarse Grid solver in MGR */
   HYPRE_MGRSetCoarseSolver( m_par_MGR_solver, HYPRE_BoomerAMGSolve,
                             HYPRE_BoomerAMGSetup, m_mgr_amg_solver);
   // AFF solver
   /*
     HYPRE_Solver aff_solver;
     HYPRE_BoomerAMGCreate(&aff_solver);
     HYPRE_BoomerAMGSetCGCIts(aff_solver, 1);
     HYPRE_BoomerAMGSetInterpType(aff_solver, 0);
     HYPRE_BoomerAMGSetPostInterpType(aff_solver, 0);
     HYPRE_BoomerAMGSetCoarsenType(aff_solver, 6);
     HYPRE_BoomerAMGSetPMaxElmts(aff_solver, 0);
     HYPRE_BoomerAMGSetCycleType(aff_solver, 1);
     HYPRE_BoomerAMGSetFCycle(aff_solver, 0);
     HYPRE_BoomerAMGSetNumSweeps(aff_solver, 1);
     HYPRE_BoomerAMGSetRelaxType(aff_solver, 3);

     HYPRE_BoomerAMGSetRelaxOrder(aff_solver, 1);
     HYPRE_BoomerAMGSetMaxLevels(aff_solver, 25);
     HYPRE_BoomerAMGSetSmoothType(aff_solver, 6);
     HYPRE_BoomerAMGSetSmoothNumSweeps(aff_solver, 1);
     HYPRE_BoomerAMGSetMaxIter(aff_solver, 1);
     HYPRE_BoomerAMGSetTol(aff_solver, 0.0);
     HYPRE_BoomerAMGSetPrintLevel(aff_solver, 1);
     // set the MGR aff solver. Comment out to use internal solvers in MGR 
     HYPRE_MGRSetFSolver(m_par_MGR_solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, aff_solver);
   */

   HYPRE_ParVector par_b;  // apparently not used, but needs to be passed to HYPRE_MGRSetup
   HYPRE_ParVector par_x;  // apparently not used, but needs to be passed to HYPRE_MGRSetup

   HYPRE_MGRSetup(m_par_MGR_solver, par_A, par_b, par_x);

   m_MGR_solver_allocated = true;
}


void
MBHypreSolver::AMG( const HYPRE_SStructMatrix&  a_matrix,
                    const HYPRE_SStructVector&  a_b,
                    const double                a_tol,
                    const int                   a_max_iter,
                    const bool                  a_verbose,
                    const HYPRE_SStructVector&  a_x )
{
   CH_TIMERS("MBHypreSolver::AMG");

   HYPRE_ParCSRMatrix    par_A;
   HYPRE_SStructMatrixGetObject(a_matrix, (void **) &par_A);
   HYPRE_ParVector       par_b;
   HYPRE_SStructVectorGetObject(a_b, (void **) &par_b);
   HYPRE_ParVector       par_x;
   HYPRE_SStructVectorGetObject(a_x, (void **) &par_x);

   HYPRE_BoomerAMGSetTol(m_par_AMG_solver, a_tol);
   HYPRE_BoomerAMGSetMaxIter(m_par_AMG_solver, a_max_iter);

   HYPRE_BoomerAMGSolve(m_par_AMG_solver, par_A, par_b, par_x);

   int num_iterations;
   HYPRE_BoomerAMGGetNumIterations(m_par_AMG_solver, &num_iterations);

   double final_res_norm;
   HYPRE_BoomerAMGGetFinalRelativeResidualNorm(m_par_AMG_solver, &final_res_norm);

   if (a_verbose && procID()==0) {
      cout << "        --> AMG solver residual = " << final_res_norm << " after " 
           << num_iterations << " iterations" << endl;
   }

   if (m_hypre_object_type == HYPRE_PARCSR) {
      HYPRE_SStructVectorGather(a_x);
   }
}


void
MBHypreSolver::ILU( const HYPRE_SStructMatrix&  a_matrix,
                    const HYPRE_SStructVector&  a_b,
                    const double                a_tol,
                    const int                   a_max_iter,
                    const bool                  a_verbose,
                    const HYPRE_SStructVector&  a_x )
{
   CH_TIMERS("MBHypreSolver::ILU");
   CH_TIMER("ILU_solve",t_ILU_solve);

   HYPRE_ParCSRMatrix    par_A;
   HYPRE_ParVector       par_b;
   HYPRE_ParVector       par_x;

   HYPRE_SStructMatrixGetObject(a_matrix, (void **) &par_A);
   HYPRE_SStructVectorGetObject(a_b, (void **) &par_b);
   HYPRE_SStructVectorGetObject(a_x, (void **) &par_x);

   CH_START(t_ILU_solve);
   HYPRE_ILUSolve(m_par_ILU_solver, par_A, par_b, par_x);
   CH_STOP(t_ILU_solve);

   int num_iterations;
   HYPRE_ILUGetNumIterations(m_par_ILU_solver, &num_iterations);
   double final_res_norm;
   HYPRE_ILUGetFinalRelativeResidualNorm(m_par_ILU_solver, &final_res_norm);
   if (a_verbose && procID()==0) {
      cout << "        --> ILU solver residual = " << final_res_norm << " after " 
           << num_iterations << " iterations" << endl;
   }

   if (m_hypre_object_type == HYPRE_PARCSR) {
      HYPRE_SStructVectorGather(a_x);
   }
}


void
MBHypreSolver::MGR( const HYPRE_SStructMatrix&  a_matrix,
                    const HYPRE_SStructVector&  a_b,
                    const double                a_tol,
                    const int                   a_max_iter,
                    const bool                  a_verbose,
                    const HYPRE_SStructVector&  a_x )
{
   CH_TIMERS("MBHypreSolver::MGR");

   // Set the convergence criteria
   HYPRE_MGRSetMaxIter(m_par_MGR_solver, a_max_iter);
   HYPRE_MGRSetTol(m_par_MGR_solver, a_tol);
   
   HYPRE_ParCSRMatrix    par_A;
   HYPRE_SStructMatrixGetObject(a_matrix, (void **) &par_A);
   HYPRE_ParVector       par_b;
   HYPRE_SStructVectorGetObject(a_b, (void **) &par_b);
   HYPRE_ParVector       par_x;
   HYPRE_SStructVectorGetObject(a_x, (void **) &par_x);

   // Do the solve
   HYPRE_MGRSolve(m_par_MGR_solver, par_A, par_b, par_x);
      
   if ( a_verbose ) {
      // Get the convergence data
      int num_iterations;
      HYPRE_MGRGetNumIterations(m_par_MGR_solver, &num_iterations);

      HYPRE_Real final_res_norm;
      HYPRE_MGRGetFinalRelativeResidualNorm(m_par_MGR_solver, &final_res_norm);

      if ( procID()==0 ) {
         cout << "        --> MGR solver residual = " << final_res_norm << " after " 
              << num_iterations << " iterations" << endl;
      }
   }

   if (m_hypre_object_type == HYPRE_PARCSR) {
      HYPRE_SStructVectorGather(a_x);
   }
}


void
MBHypreSolver::AMG_preconditioned_GMRES( const HYPRE_SStructMatrix&  a_matrix,
                                         const HYPRE_SStructMatrix&  a_precond,
                                         const HYPRE_SStructVector&  a_b,
                                         const double                a_tol,
                                         const int                   a_max_iter,
                                         const double                a_amg_tol,
                                         const int                   a_amg_max_iter,
                                         const bool                  a_verbose,
                                         const HYPRE_SStructVector&  a_x ) const
{
   HYPRE_Solver          par_solver;
   HYPRE_Solver          par_precond;
   HYPRE_ParCSRMatrix    par_A;
   HYPRE_ParCSRMatrix    par_P;
   HYPRE_ParVector       par_b;
   HYPRE_ParVector       par_x;

   HYPRE_SStructMatrixGetObject(a_matrix, (void **) &par_A);
   HYPRE_SStructMatrixGetObject(a_precond, (void **) &par_P);
   HYPRE_SStructVectorGetObject(a_b, (void **) &par_b);
   HYPRE_SStructVectorGetObject(a_x, (void **) &par_x);

   HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &par_solver);

   /* set the GMRES paramaters */
   HYPRE_GMRESSetMaxIter(par_solver, a_max_iter);
   HYPRE_GMRESSetTol(par_solver, a_tol);
   HYPRE_GMRESSetPrintLevel(par_solver, 2);
   HYPRE_GMRESSetLogging(par_solver, 1);

   /* use BoomerAMG as preconditioner */
   HYPRE_BoomerAMGCreate(&par_precond);
   HYPRE_BoomerAMGSetCoarsenType(par_precond, 6);
   HYPRE_BoomerAMGSetOldDefault(par_precond);
   //   HYPRE_BoomerAMGSetStrongThreshold(par_precond, 0.25);
   HYPRE_BoomerAMGSetStrongThreshold(par_precond, 0.8);
   //   HYPRE_BoomerAMGSetTol(par_precond, a_amg_tol);
   HYPRE_BoomerAMGSetTol(par_precond, 1.e-6);
   HYPRE_BoomerAMGSetPrintLevel(par_precond, 3);
   //   HYPRE_BoomerAMGSetPrintFileName(par_precond, "ex9.out.log");
   //   HYPRE_BoomerAMGSetMaxIter(par_precond, a_amg_max_iter);
   HYPRE_BoomerAMGSetMaxIter(par_precond, 50);

   /* set the preconditioner */
   HYPRE_ParCSRGMRESSetPrecond(par_solver,
                               HYPRE_BoomerAMGSolve,
                               HYPRE_BoomerAMGSetup,
                               par_precond);

   HYPRE_GMRESSetup(par_solver, (HYPRE_Matrix)par_P, (HYPRE_Vector)par_b, (HYPRE_Vector)par_x);
   HYPRE_GMRESSolve(par_solver, (HYPRE_Matrix)par_A, (HYPRE_Vector)par_b, (HYPRE_Vector)par_x);

   // Get number of preconditioner iterations used
   int num_final_iterations;
   HYPRE_GMRESGetNumIterations(par_solver, &num_final_iterations);

   // Get final preconditioner residual norm
   double final_norm2;
   HYPRE_GMRESGetFinalRelativeResidualNorm(par_solver, &final_norm2);

   if (a_verbose && procID()==0) {
      cout << "      --> GMRES residual = " << final_norm2 << " after " << num_final_iterations << " iterations" << endl;
   }

   if (m_hypre_object_type == HYPRE_PARCSR) {
      HYPRE_SStructVectorGather(a_x);
   }
}


void
MBHypreSolver::copyToHypreVector( const LevelData<FArrayBox>&  a_in,
                                  HYPRE_SStructVector&         a_out ) const
{
   CH_TIME("MBHypreSolver::copyToHypreVector");

   const DisjointBoxLayout& grids = m_geometry.grids();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box& box = grids[dit];
      IntVect lower(box.loVect());
      IntVect upper(box.hiVect());

      int block_number = m_coord_sys_ptr->whichBlock(box);

      FArrayBox tmp(box,1);

      for (int var=0; var<m_nvar; ++var) {
         tmp.copy(a_in[dit], var, 0);
         HYPRE_SStructVectorSetBoxValues(a_out, block_number, lower.dataPtr(), upper.dataPtr(), var, tmp.dataPtr());
      }
   }
} 



void
MBHypreSolver::copyFromHypreVector( const HYPRE_SStructVector&  a_in,
                                    LevelData<FArrayBox>&       a_out ) const
{
   CH_TIME("MBHypreSolver::copyFromHypreVector");

   const DisjointBoxLayout& grids = m_geometry.grids();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box& box = grids[dit];
      IntVect lower(box.loVect());
      IntVect upper(box.hiVect());

      int block_number = m_coord_sys_ptr->whichBlock(box);

      FArrayBox tmp(box,1);

      for (int var=0; var<m_nvar; ++var) {
         HYPRE_SStructVectorGetBoxValues(a_in, block_number, lower.dataPtr(), upper.dataPtr(), var, tmp.dataPtr());
         a_out[dit].copy(tmp, 0, var);
      }
   }
} 



#include "NamespaceFooter.H"
