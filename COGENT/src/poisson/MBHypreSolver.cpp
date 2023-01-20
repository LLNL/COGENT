#include "MBHypreSolver.H"
#include "MBStencilIterator.H"
#include "BlockRegister.H"
#include "SparseCoupling.H"
#include "MBSolverF_F.H"
#include "MagGeom.H"

#include "_hypre_sstruct_mv.h"
//#include "_hypre_parcsr_ls.h"

#include "AMRIO.H"

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
     m_mgr_CF_indexes(NULL),
     m_matrix_initialized(false),
     m_matrix_finalized(false)
{
   m_rhs_from_bc.define(a_geom.grids(), a_nvar, IntVect::Zero);
   m_advection_rhs_from_bc.define(a_geom.grids(), a_nvar, IntVect::Zero);
   createHypreData();
}
      

MBHypreSolver::~MBHypreSolver()
{
   for (int n=0; n<m_advection_stencil_size; ++n) delete [] m_advection_offsets[n];
   delete [] m_advection_offsets;

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
      if ( a_pp.query("mgr_schur", m_MGR_schur) == 0 ) m_MGR_schur = false;
      if ( a_pp.query("mgr_schur_1", m_MGR_schur_1) == 0 ) m_MGR_schur_1 = false;
      if ( a_pp.query("print_level", m_MGR_print_level) == 0 ) m_MGR_print_level = 0;
      if ( a_pp.query("amg_print_level", m_MGR_amg_print_level) == 0 ) m_MGR_amg_print_level = 0;
      if ( a_pp.query("amg_tol", m_MGR_amg_tol) == 0 ) m_MGR_amg_tol = 0.;
      if ( a_pp.query("amg_max_iter", m_MGR_amg_max_iter) == 0 ) m_MGR_amg_max_iter = 1;
      if ( a_pp.query("amg_threshold", m_MGR_amg_threshold) == 0 ) m_MGR_amg_threshold = 0.25;
      if ( a_pp.query("amg_verbose", m_MGR_amg_verbose) == 0 ) m_MGR_amg_verbose = false;
      if ( a_pp.query("interp_type", m_MGR_interp_type) == 0 ) m_MGR_interp_type = 3;
      if ( a_pp.query("restrict_type", m_MGR_restrict_type) == 0 ) m_MGR_restrict_type = 0;

      if ( a_pp.query("relax_type", m_MGR_relax_type) != 0 ) {
         // For backward compatibility with old inputs, check to see if frelax_method has been specified.
         // If so, and relax_type has also been specified, make sure they are the same
         int frelax_method;
         if ( (a_pp.query("frelax_method", frelax_method) != 0) && frelax_method != m_MGR_relax_type ) {
            MayDay::Error("MBHypreSolver::setMethodParams: Different values of frelax_method and relax_type specified");
         }
      }
      else if ( a_pp.query("frelax_method", m_MGR_relax_type) != 0 ) {
         // If frelax_method was specified as 1, change it to 0
         if ( m_MGR_relax_type == 1 ) m_MGR_relax_type == 0; 
      }
      else {
         m_MGR_relax_type = 0;
      }

      // Override certain inputs that need to be set in a particular way when MGR is being used to solve
      // a block system by constructing a Schur complement, solving it, and then backsolving
      if ( m_MGR_schur && m_MGR_schur_1 ) {
         MayDay::Error("MBHypreSolver::setMethodParams: m_MGR_schur and m_MGR_schur_1 cannot both be true");
      }
      if ( m_MGR_schur ) {
         m_MGR_cpoint = 0;
      }
      else if ( m_MGR_schur_1 ) {
         m_MGR_cpoint = 1;
      }
      if ( m_MGR_schur || m_MGR_schur_1 ) {
         m_MGR_interp_type = 2;
         m_MGR_restrict_type = 0;
         m_MGR_relax_type = 0;
      }
      
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
   else if ( m_method == "AIR" ) {
      if ( a_pp.query("air_restrict_type", m_AIR_restrict_type) == 0 ) m_AIR_restrict_type = 1;
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
   else if ( m_method == "AMG" || m_method == "AIR" ) {
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
MBHypreSolver::dumpVector(const string& a_file_name) const
{
  HYPRE_SStructVectorPrint(a_file_name.c_str(), m_b, 0);
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
      AMGSetup(m_A, m_x, m_b);
   }
   else if ( m_method == "ILU" ) {
      ILUSetup(m_A);
   }
   else if ( m_method == "MGR" ) {
      MGRSetup(m_A);
   }
   else if ( m_method == "AIR" ) {
      AIRSetup(m_A);
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
MBHypreSolver::constructAdvectionOffsets()
{
   m_advection_stencil_size = 2*(SpaceDim)*m_advection_stencil_radius + 1;
   m_advection_offsets = new int*[m_advection_stencil_size];

   // Stencil center offset
   int n = 0;
   m_advection_offsets[n] = new int[CH_SPACEDIM];
   for (int dir=0; dir<CH_SPACEDIM; ++dir) {
      m_advection_offsets[n][dir] = 0;
   }

   for (int dir=0; dir<SpaceDim; ++dir) {

     for (int m=1; m<=m_advection_stencil_radius; ++m) {
       // Stencil offsets along positive leg
       n++;
       m_advection_offsets[n] = new int[CH_SPACEDIM];
       for (int tdir=0; tdir<CH_SPACEDIM; ++tdir) {
	 if (tdir == dir) {
           m_advection_offsets[n][tdir] = m;
	 }
         else {
           m_advection_offsets[n][tdir] = 0;
         }
       }
       // Stencil offsets along negative leg
       n++;
       m_advection_offsets[n] = new int[CH_SPACEDIM];
       for (int tdir=0; tdir<CH_SPACEDIM; ++tdir) {
         if (tdir == dir) {
           m_advection_offsets[n][tdir] = -m;
         }
         else {
           m_advection_offsets[n][tdir] = 0;
         }
       }
     }
   }
}


void
MBHypreSolver::createGraph( const Box& a_stencil_box,
                            const int  a_radius )
{
   // First construct the offsets for advection terms in the operator, if any.  N.B.: It is currently
   // assumed such offsets are a subset of the offsets computed below for the second-order
   // elliptic term.  It is also assumed that the advection operator is discretized using either
   // second-order, centered differences ("C2") or first-order upwind ("UW1").  If any other
   // scheme options are added, it is essential to ensure that the corresponding offsets remain
   // within those of the second-order term.  Otherwise, the code below will need to be generalized
   // compute the graph object from the union of the offsets of the advection and second-order terms.
   
   m_advection_stencil_radius = 1;
   constructAdvectionOffsets();

   // Now construct the offsets for the second-order term and the corresponding Hypre graph object

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

         if(m_mgr_CF_indexes)
            hypre_TFree(m_mgr_CF_indexes, HYPRE_MEMORY_HOST);
         m_mgr_CF_indexes = NULL;

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

int MBHypreSolver::findStructuredEntry( const int      a_block_column,
					const Box&     a_stencil_box,
					const IntVect& a_iv ) const
{
   /*
     The Hypre matrix coefficients corresponding to the equation centered at cell a_iv are
     identified by an entry number, which is returned by this function.  Regular "stencil"
     couplings to cells are enumerated first, followed by the unstructured couplings to
     cells in other blocks contained in a_unstructured_ivs.
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

   // bool found_entry = false;

   // int entry;
   // for (entry=0; entry<m_advection_stencil_size; ++entry) {
   //    IntVect offset(m_advection_offsets[entry]);
   //    offset += a_iv_center;
   //    if (offset == a_iv) {
   // 	 found_entry = true;
   // 	 break;
   //    }
   // }

   CH_assert(found_entry);


   return entry;
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

         const Box& domain_box = (m_coord_sys_ptr->mappingBlocks())[block_number];
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
               Box stencil_box_valid = stencil_box & domain_box;
               IntVectSet extra_block_ghosts = stencil_box & this_ghosts_fab;

               double alpha = alphaPtr? alphaPtr->operator()(iv): 1.;

               // Get all of the indices that are sparsely coupled to this cell
               const IntVectSet& unstructured_couplings = a_unstructured_coupling[dit](iv);

               for (int dir=0; dir<SpaceDim; ++dir) {
                  bool at_lo_interface = lo_block_interface[dir] && (iv[dir] == domain_box.smallEnd(dir));
                  bool at_hi_interface = hi_block_interface[dir] && (iv[dir] == domain_box.bigEnd(dir));
                  bool force_codim2_condense = lo_block_interface[dir] || hi_block_interface[dir];
               
                  bool at_lo_boundary = (!at_lo_interface && iv[dir] == domain_box.smallEnd(dir));
                  bool at_hi_boundary = (!at_hi_interface && iv[dir] == domain_box.bigEnd(dir));

                  int bc_type = EllipticOpBC::UNDEFINED;
                  if ( at_lo_boundary ) bc_type = a_bc.getBCType(block_number, dir, Side::LoHiSide::Lo);
                  if ( at_hi_boundary ) bc_type = a_bc.getBCType(block_number, dir, Side::LoHiSide::Hi);

                  for (SideIterator sit; sit.ok(); ++sit) {
                     Side::LoHiSide side = sit();

                     IntVect iv_face;
                     iv_face = iv;
                     iv_face.shift(dir,side);

                     // Halve the contributions on face contained in block interfaces to effect averaging
                     double fac;
                     if ( (at_lo_interface && side == Side::LoHiSide::Lo) ||
                          (at_hi_interface && side == Side::LoHiSide::Hi) ) {
                        fac = m_flux_average[dir]? 0.5: 1.0;
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
                                             iv, a_stencil_values, dummy, 0, false, force_codim2_condense, false );

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
                                 if ( (at_lo_interface && side == Side::LoHiSide::Lo) ||
                                      (at_hi_interface && side == Side::LoHiSide::Hi) ) {
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
                                    if ( (at_lo_interface && side == Side::LoHiSide::Lo) ||
                                         (at_hi_interface && side == Side::LoHiSide::Hi) ) {
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
MBHypreSolver::addAdvectionUnstructuredMatrixEntries( const int                                     a_block_row,
						      const int                                     a_block_column,
                                                      const LevelData<FArrayBox>&                   a_alpha_coefficient,
						      const LevelData<BaseFab<Vector<IntVect> > >&  a_structured_couplings,
						      const LevelData<BaseFab<Vector<Real> > >&     a_structured_weights,
						      const LayoutData< BaseFab<IntVectSet> >&      a_unstructured_couplings,
						      HYPRE_SStructMatrix&                          a_matrix ) const
{
   CH_TIMERS("MBHypreSolver::addAdvectionUnstructuredMatrixEntries()");
   CH_TIMER("add_stencil_contrib",t_add_stencil_contrib);
   CH_TIMER("add_extrablock_contrib",t_add_extrablock_contrib);

   if ( m_mblex_potential_Ptr ) {

      const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = m_coord_sys_ptr->boundaries();
      const LayoutData< RefCountedPtr< IVSFAB<MBStencil> > >& stencil = m_mblex_potential_Ptr->stencils();

      const DisjointBoxLayout & grids = m_geometry.gridsFull();

      BlockBaseRegister<BaseFab<SparseCoupling> > block_register(m_coord_sys_ptr, grids, 0);

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box & box = grids[dit];

         int block_number = m_coord_sys_ptr->whichBlock(box);

         const NewCoordSys* block_coord_sys = m_coord_sys_ptr->getCoordSys(block_number);
         RealVect dx = block_coord_sys->dx();

         const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];
         bool lo_block_interface[CH_SPACEDIM];
         bool hi_block_interface[CH_SPACEDIM];
         for (int dir=0; dir<SpaceDim; ++dir) {
            lo_block_interface[dir] = this_block_boundaries[dir].isInterface();
            hi_block_interface[dir] = this_block_boundaries[dir + SpaceDim].isInterface();
         }

         const FArrayBox* alphaPtr = NULL;
         if ( a_alpha_coefficient.isDefined() ) {
            alphaPtr = &a_alpha_coefficient[dit];
         }

         const BaseFab<Vector<IntVect> >& this_structured_couplings = a_structured_couplings[dit];
         const BaseFab<Vector<Real> >& this_structured_weights = a_structured_weights[dit];

         const Box& domainBox = (m_coord_sys_ptr->mappingBlocks())[block_number];
         if ( stencil[dit] != NULL ) {
            const IVSFAB<MBStencil>& this_stencil = *stencil[dit];

            BaseFab<SparseCoupling> data[CFG_DIM];
            for(int dir=0; dir<CFG_DIM; ++dir ) {
               if (m_flux_average[dir]) {
                  data[dir].define(surroundingNodes(box, dir), 1);
               }
            }

            // Here, we loop over all cells whose stencil is affected by a block boundary
            IntVectSet interblock_ivs = getInterBlockCoupledCells(block_number, m_advection_stencil_radius, box);
            for (IVSIterator ivsit(interblock_ivs); ivsit.ok(); ++ivsit) {
               IntVect iv = ivsit();

               double alpha = alphaPtr? alphaPtr->operator()(iv): 1.;

               // Get all of the indices that are sparsely coupled to this cell
               const IntVectSet& unstructured_couplings = a_unstructured_couplings[dit](iv);

               for (int dir=0; dir<CFG_DIM; ++dir) {
                  bool at_lo_interface = lo_block_interface[dir] && (iv[dir] == domainBox.smallEnd(dir));
                  bool at_hi_interface = hi_block_interface[dir] && (iv[dir] == domainBox.bigEnd(dir));

                  Box stencil_box(iv,iv);
                  stencil_box.grow(m_advection_stencil_radius);
                  Box stencil_box_valid = stencil_box & domainBox;

                  for (SideIterator sit; sit.ok(); ++sit) {
                     Side::LoHiSide side = sit();

                     if ( (side == Side::LoHiSide::Lo && at_lo_interface) ||
                          (side == Side::LoHiSide::Hi && at_hi_interface) ) {

                        IntVect iv_face;

                        // Halve the contributions on face contained in block interfaces to effect averaging
                        double fac;
                        if ( (at_lo_interface && side == Side::LoHiSide::Lo) ||
                             (at_hi_interface && side == Side::LoHiSide::Hi) ) {
                           fac = m_flux_average[dir]? 0.5: 1.0;
                           iv_face = iv;
                           iv_face.shift(dir,side);
                        }
                        else {
                           fac = 1.;
                        }

                        int coupling_comp = dir + ((1-sign(side))/2)*CFG_DIM;
                        for (int m=0; m<this_structured_couplings(iv,coupling_comp).size(); ++m) {
                           IntVect stencil_box_iv = this_structured_couplings(iv,coupling_comp)[m];
                           
                           double s = fac * this_structured_weights(iv,coupling_comp)[m];

                           if ( s != 0. ) {

                              if ( stencil_box_valid.contains(stencil_box_iv) ) {

                                 CH_START(t_add_stencil_contrib);
                                 // If the current side of cell iv is on a block interface, add the stencil
                                 // contribution to the BlockRegister data.
                                 // The minus sign accounts for the normal component, which is negated on
                                 // the adacent cell.  Also, since we have not
                                 // yet accounted for the regular stencil contributions, do that now.
                                 if ( (at_lo_interface && side == Side::LoHiSide::Lo) || 
                                      (at_hi_interface && side == Side::LoHiSide::Hi) ) {
                                    if (m_flux_average[dir]) {
                                       // Need to defer alpha multiplication until after flux averaging
                                       data[dir](iv_face).add(stencil_box_iv, -s);
                                    }
                                    
                                    int entry = findHypreEntry(a_block_column, stencil_box, unstructured_couplings, stencil_box_iv);

                                    double value = s * alpha;
                                    
                                    int num_entries = 1;
                                    HYPRE_SStructMatrixAddToValues(a_matrix, block_number, iv.dataPtr(),
                                                                   a_block_row, num_entries, &entry, &value);
                                 }
                                 CH_STOP(t_add_stencil_contrib);
                              }
                              else if ( (stencil_box_iv[dir] < domainBox.smallEnd(dir) 
                                         && lo_block_interface[dir]) ||
                                        (stencil_box_iv[dir] > domainBox.bigEnd(dir)
                                         && hi_block_interface[dir]) ) {

                                 CH_START(t_add_extrablock_contrib);

                                 const MBStencil& elements = this_stencil(stencil_box_iv, 0);
                                 int num_unstructured_entries = elements.size();
                                 double * values = new double[num_unstructured_entries];
                                 int * entries = new int[num_unstructured_entries];

                                 MBStencilIterator stit(elements);
                                 int k = 0;
                                 for (stit.begin(); stit.ok(); ++stit) {
                                    const MBStencilElement& stencilElement = stit();
                                    IntVect interp_cell = stencilElement.cell();

                                    entries[k] = findHypreEntry(a_block_column, stencil_box, unstructured_couplings, interp_cell);
                                    values[k] = s * stencilElement.weight();
                                    
                                    // If we are on a block interface, add the stencil contribution to the
                                    // BlockRegister data.  The minus sign accounts for the normal component,
                                    // which is negated on the adacent cell.  Need to defer multiplication by
                                    // alpha until after the flux averaging.
                                    if ( m_flux_average[dir] &&
                                         ((at_lo_interface && side == Side::LoHiSide::Lo) ||
                                          (at_hi_interface && side == Side::LoHiSide::Hi)) ) {
                                       data[dir](iv_face).add(interp_cell, -values[k]);
                                    }

                                    values[k] *= alpha;
                                    k++;
                                 }

                                 HYPRE_SStructMatrixAddToValues(a_matrix, block_number, iv.dataPtr(),
                                                                a_block_row, num_unstructured_entries, entries, values);

                                 delete [] entries;
                                 delete [] values;
                                 
                                 CH_STOP(t_add_extrablock_contrib);
                              }
                           }
                        }
                     }
                  }
               }
            }

            // Store the contributions at block interfaces
            for (int dir=0; dir<CFG_DIM; dir++) {
               if ( m_flux_average[dir] ) {
                  for (SideIterator sit; sit.ok(); ++sit) {
                     Side::LoHiSide side = sit();
                     if (block_register.hasInterface(dit(), dir, side)) {
                         block_register.store(data[dir], dit(), dir, side);
                     }
                  }
               }
            }
         }
      }

      // Exchange the block interface contributions and accumulate them to the matrix
      block_register.exchange();

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         int block_number = m_coord_sys_ptr->whichBlock(grids[dit]);

         const FArrayBox* alphaPtr = NULL;
         if ( a_alpha_coefficient.isDefined() ) {
            alphaPtr = &a_alpha_coefficient[dit];
         }

         for (SideIterator sit; sit.ok(); ++sit) {
            Side::LoHiSide side = sit();
            for (int dir=0; dir<CFG_DIM; dir++) {
               if ( m_flux_average[dir] && block_register.hasInterface(dit(), dir, side) ) {
                  Box faceBox = adjCellBox(grids[dit], dir, side, 1);
                  // if Lo, then shift +1; if Hi, then shift -1
                  faceBox.shiftHalf(dir, -sign(side));
                  BaseFab<SparseCoupling> otherFab(faceBox, 1);
                  block_register.fill(otherFab, dit(), dir, side, side);
                  BoxIterator bit(faceBox);
                  for (bit.begin(); bit.ok(); ++bit) {
                     IntVect iv_face = bit();
                     IntVect iv(iv_face);
                     iv.shift(dir,-side);

                     Box stencil_box(iv,iv);
                     stencil_box.grow(dir, m_advection_stencil_radius);
                     
                     const IntVectSet& unstructured_couplings = a_unstructured_couplings[dit](iv);

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

                  }
               }
            }
         }
      }
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

   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries 
      = m_coord_sys_ptr->boundaries();

   Vector< Vector<CoDim1Stencil> > codim1_stencils;
   Vector< Vector<CoDim2Stencil> > codim2_stencils;

   constructBoundaryStencils(a_fourth_order, a_bc, codim1_stencils, codim2_stencils);

   testMatrixConstruct( a_alpha_coefficient, 
                        a_tensor_coefficient,
                        a_beta_coefficient,
                        block_boundaries,
                        codim1_stencils,
                        codim2_stencils,
                        a_bc,
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

   Vector< Vector<CoDim1Stencil> > codim1_stencils;
   Vector< Vector<CoDim2Stencil> > codim2_stencils;

   constructBoundaryStencils(a_fourth_order, a_bc, codim1_stencils, codim2_stencils );

   LevelData<FluxBox> bc_factor;
   getNeumannNaturalFactor(a_bc, bc_factor);

#if 0
   // This is a newer implementation in which the loops over directions and
   // sides are outside the loop over cells
   
   // Set the intra-block matrix entries

   for (DataIterator dit(grids); dit.ok(); ++dit) {
     const Box& box = grids[dit];
     IntVect lower(box.loVect());
     IntVect upper(box.hiVect());

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

           int bc_type = a_bc.getBCType(block_number, dir, side);

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

                 bool at_lo_boundary = (!at_lo_interface && iv[dir] == domain_box.smallEnd(dir));
                 bool at_hi_boundary = (!at_hi_interface && iv[dir] == domain_box.bigEnd(dir));
                 bool at_boundary = at_lo_boundary || at_hi_boundary;

                 if (  !(bc_type == EllipticOpBC::NATURAL && at_boundary) ) {

                    accumStencilMatrixEntries(iv, dir, side, dir2, this_coef,
                                              dx, a_fourth_order, a_stencil_values);
                 }

                 modifyStencilForBCs( codim1_stencils[block_number], codim2_stencils[block_number],
                                      iv, a_stencil_values, a_rhs_from_bc[dit], a_block_row,
                                      update_rhs_from_bc_only, force_codim2_condense );

                 if ( alpha ) {
                    double alpha_val = alpha->operator()(iv);
                    a_stencil_values *= alpha_val;
                    a_rhs_from_bc[dit](iv, a_block_row) *= alpha_val;
                 }

                 double* this_stencil_values = a_stencil_values.dataPtr(0);

                 if ( add_beta ) {
                    this_stencil_values[a_diagonal_offset] += beta->operator()(iv);
                 }

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

   // We set only the current component to zero, however this will need to be re-thought
   // once we have inhomogenous BCs on the off-diagonal spots of the block matrix.
   for(DataIterator dit(grids); dit.ok(); ++dit)
     {
       a_rhs_from_bc[dit].setVal(0., a_block_row);
     }

   for (DataIterator dit(grids); dit.ok(); ++dit) {
     const Box & box = grids[dit];

     int block_number = m_coord_sys_ptr->whichBlock(box);
     const MagBlockCoordSys& block_coord_sys = *(MagBlockCoordSys*)(m_coord_sys_ptr->getCoordSys(block_number));
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

     const Box& domain_box = (m_coord_sys_ptr->mappingBlocks())[block_number];

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
          bool at_lo_interface = (lo_block_interface[dir] && iv[dir] == domain_box.smallEnd(dir));
          bool at_hi_interface = (hi_block_interface[dir] && iv[dir] == domain_box.bigEnd(dir));
          bool force_codim2_condense = lo_block_interface[dir] || hi_block_interface[dir];

          bool at_lo_boundary = (!at_lo_interface && iv[dir] == domain_box.smallEnd(dir));
          bool at_hi_boundary = (!at_hi_interface && iv[dir] == domain_box.bigEnd(dir));

          int bc_type = EllipticOpBC::UNDEFINED;
          if ( at_lo_boundary ) bc_type = a_bc.getBCType(block_number, dir, Side::LoHiSide::Lo);
          if ( at_hi_boundary ) bc_type = a_bc.getBCType(block_number, dir, Side::LoHiSide::Hi);

          for (SideIterator sit; sit.ok(); ++sit) {
             Side::LoHiSide side = sit();

             bool update_rhs_from_bc_only = (at_lo_interface && side == Side::LoHiSide::Lo)
                                         || (at_hi_interface && side == Side::LoHiSide::Hi);

             bool update_rhs = !((at_lo_boundary && 
                                  (bc_type == EllipticOpBC::NATURAL || bc_type == EllipticOpBC::NEUMANN) ) ||
                                 (at_hi_boundary &&
                                  (bc_type == EllipticOpBC::NATURAL || bc_type == EllipticOpBC::NEUMANN) ));

             for (int dir2=0; dir2<SpaceDim; ++dir2) {

                tmp_stencil_values.setVal(0.);

                accumStencilMatrixEntries(iv, dir, side, dir2, this_coef,
                                          dx, a_fourth_order, tmp_stencil_values);

                modifyStencilForBCs( codim1_stencils[block_number], codim2_stencils[block_number],
                                     iv, tmp_stencil_values, a_rhs_from_bc[dit], a_block_row,
                                     update_rhs_from_bc_only, force_codim2_condense, update_rhs );

                a_stencil_values += tmp_stencil_values;
             }

             // Add the right-hand side contributions from Neumann and natural boundary
             // conditions.  Contributions from other types of boundary conditions have
             // already been taken care of by the modifyStencilForBCs immediately above.
             if ( ((at_lo_boundary && side == Side::LoHiSide::Lo) || 
                   (at_hi_boundary && side == Side::LoHiSide::Hi)) &&
                  (bc_type == EllipticOpBC::NEUMANN || bc_type == EllipticOpBC::NATURAL) ) {
                IntVect iv_face = iv;
                if (side == Side::LoHiSide::Hi) iv_face[dir]++;
                
                double bv = getBV(a_bc, block_coord_sys, block_number, dir, side, iv_face);
                a_rhs_from_bc[dit](iv,0) += bv * bc_factor[dit][dir](iv_face,0);
             }
          }
       }

       if ( alpha ) {
          a_stencil_values *= alpha->operator()(iv);
          a_rhs_from_bc[dit](iv, a_block_row) *= alpha->operator()(iv);
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

     HYPRE_SStructMatrixAddToBoxValues(a_matrix, block_number, lower.dataPtr(), upper.dataPtr(),
                                       a_block_row, stencil_size, entries, values);

     delete [] values;
   }

#endif

   delete [] entries;

   addUnstructuredMatrixEntries(a_block_row, a_block_column, a_alpha_coefficient, a_tensor_coefficient,
                                a_bc, a_stencil_values, a_fourth_order, a_unstructured_coupling,
                                codim1_stencils, codim2_stencils, a_matrix);
}


void MBHypreSolver::addAdvectionMatrixBlockGeneral( const int                                     a_block_row,
                                                    const int                                     a_block_column,
                                                    const LevelData<FArrayBox>&                   a_alpha_coefficient,
                                                    const LevelData<FluxBox>&                     a_vector_coefficient,
                                                    const EllipticOpBC&                           a_bc)
{
  CH_TIMERS("MBHypreSolver::addAdvectionMatrixBlockGeneral()");
  const int coupling_comps = 2*SpaceDim;
  const IntVect coupling_ghosts = IntVect::Unit;
  const DisjointBoxLayout& grids = m_geometry.grids();
  LevelData<BaseFab<Vector<IntVect> > > couplings(grids, coupling_comps, coupling_ghosts);
  LevelData<BaseFab<Vector<Real> > > weights(grids, coupling_comps, coupling_ghosts);

  // We set m_advection_rhs_from_bc to zero here. Then constructAdvectionStencils puts
  // a 1 in the location where we want to keep the boundary values. This is for
  // zero-order extrapolation of the boundary values. For higher-order, something
  // different needs to be done.
  
  LevelData<FArrayBox> multFactor;
  multFactor.define(grids, m_nvar, IntVect::Zero);
  for(DataIterator dit(m_advection_rhs_from_bc.dataIterator()); dit.ok(); ++dit) {
    m_advection_rhs_from_bc[dit].setVal(0., a_block_column);
    multFactor[dit].setVal(0., a_block_column);
  }

  // Construct the stencils corresponding to second-order centered differencing.  First-order
  // upwind ("UW1") is also a valid option here.  
  constructAdvectionStencils("C2", a_vector_coefficient, couplings, weights, a_bc);

  addAdvectionMatrixBlockCoupled(a_block_row,
                                 a_block_column,
                                 a_alpha_coefficient,
                                 couplings,
                                 weights,
                                 a_bc,
                                 m_A_stencil_values);
}

// Alpha coefficient needs to contain a_shift * physical_volume / mapped_volume; where a_shift
// is the time-step size coefficient (1/dt)?
void MBHypreSolver::addAdvectionMatrixBlockCoupled( const int                                     a_block_row,
                                                    const int                                     a_block_column,
                                                    const LevelData<FArrayBox>&                   a_alpha_coefficient,
                                                    const LevelData<BaseFab<Vector<IntVect> > >&  a_structured_couplings,
                                                    const LevelData<BaseFab<Vector<Real> > >&     a_structured_weights,
                                                    const EllipticOpBC&                           a_bc,
                                                    const FArrayBox&                              a_stencil_values)
{
  CH_TIMERS("MBHypreSolver::addAdvectionMatrixBlockCoupled()");
  CH_TIMER("add_struct_matrix",t_add_struct_matrix);
  CH_TIMER("add_unstruct_matrix",t_add_unstruct_matrix);
  CH_TIMER("amg_setup",t_amg_setup);
  CH_TIMER("hypre_assemble",t_hypre_assemble);
  const DisjointBoxLayout& grids = m_geometry.gridsFull();

  bool lo_block_interface[SpaceDim];
  bool hi_block_interface[SpaceDim];
  bool lo_block_boundary[SpaceDim];
  bool hi_block_boundary[SpaceDim];

  Box stencil_offsets = a_stencil_values.box();
  int radius = (stencil_offsets.size(0)-1)/2;
  CH_assert(radius >= m_advection_stencil_radius);
  stencil_offsets -= radius*IntVect::Unit;

  CH_START(t_add_struct_matrix);
  const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = m_coord_sys_ptr->boundaries();
  for (DataIterator dit(grids); dit.ok(); ++dit) {
    const Box& box = grids[dit];
    int block_number = m_coord_sys_ptr->whichBlock(box);
    const Box& block_box = (m_coord_sys_ptr->mappingBlocks())[block_number];

    const FArrayBox* alpha = NULL;
    if ( a_alpha_coefficient.isDefined() ) {
      alpha = &a_alpha_coefficient[dit];
    }

    // const FArrayBox& this_volume = m_volume[dit];
    // const NewCoordSys& block_coord_sys = *(m_coord_sys_ptr->getCoordSys(block_number));

    const FArrayBox& alpha_coeff = a_alpha_coefficient[dit];
    const BaseFab<Vector<IntVect> >& this_couplings = a_structured_couplings[dit];
    const BaseFab<Vector<Real> >& this_weights = a_structured_weights[dit];

    const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];
    for (int dir=0; dir<SpaceDim; ++dir) {
      lo_block_interface[dir] = this_block_boundaries[dir].isInterface();
      hi_block_interface[dir] = this_block_boundaries[dir + SpaceDim].isInterface();
      lo_block_boundary[dir] = this_block_boundaries[dir].isDomainBoundary();
      hi_block_boundary[dir] = this_block_boundaries[dir + SpaceDim].isDomainBoundary();
    }

    // Set the intra-block matrix entries

    for (BoxIterator bit(box); bit.ok(); ++bit) {
      IntVect iv = bit();
      Box stencil_box = stencil_offsets + iv;

      // Initialize the matrix with a_shift times the identity.  The J (= physical cell volume
      // divided by mapped cell volume) factor is needed because the matrix being constructed
      // here is assumed to act on vectors in the physical frame, but the time-integration
      // shift is done in the computational frame
      // int center_entry = findStructuredEntry(a_block_column, stencil_box, iv);
      // double center_value = alpha_coeff(iv,0);
      // HYPRE_SStructMatrixSetValues(m_A, block_number, iv.dataPtr(), 0, 1, &center_entry,
      //                              &center_value);

      for (int dir=0; dir<SpaceDim; ++dir) {
        bool at_lo_interface = lo_block_interface[dir] && (iv[dir] == block_box.smallEnd(dir));
        bool at_hi_interface = hi_block_interface[dir] && (iv[dir] == block_box.bigEnd(dir));

        for (SideIterator sit; sit.ok(); ++sit) {
          Side::LoHiSide side = sit();

          // Couplings corresponding to faces at block interfaces are handled later
          if ( m_mblex_potential_Ptr == NULL || 
               (!(side == Side::LoHiSide::Lo && at_lo_interface) &&
                !(side == Side::LoHiSide::Hi && at_hi_interface) ) ) {

            int coupling_comp = dir + ((1-sign(side))/2)*(SpaceDim);
            for (int n=0; n<this_couplings(iv,coupling_comp).size(); ++n) {
              IntVect coupling = this_couplings(iv,coupling_comp)[n];
              double value = this_weights(iv,coupling_comp)[n]*alpha_coeff(iv,0);

              // Zero the couplings across physical boundaries
              if ( (coupling[dir] < block_box.smallEnd(dir) && lo_block_boundary[dir]) ||
                   (coupling[dir] > block_box.bigEnd(dir)   && hi_block_boundary[dir]) ) {
                value = 0.;
              }

              int entry = findStructuredEntry(a_block_column, stencil_box, coupling);
              HYPRE_SStructMatrixAddToValues(m_A, block_number, iv.dataPtr(), a_block_row,
                                             1, &entry, &value);
            }
          }

        }
      }
    }
  }
  CH_STOP(t_add_struct_matrix);

  CH_START(t_add_unstruct_matrix);
  addAdvectionUnstructuredMatrixEntries(a_block_row, a_block_column, a_alpha_coefficient, a_structured_couplings,
                                        a_structured_weights, m_A_unstructured_coupling, m_A);
  CH_STOP(t_add_unstruct_matrix);
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
MBHypreSolver::AMGSetup( const HYPRE_SStructMatrix&  a_matrix,
                         const HYPRE_SStructVector&  a_x,
                         const HYPRE_SStructVector&  a_b )
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

   HYPRE_ParVector par_b;
   HYPRE_SStructVectorGetObject(a_b, (void **) &par_b);

   HYPRE_ParVector par_x;
   HYPRE_SStructVectorGetObject(a_x, (void **) &par_x);

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

      if(m_mgr_CF_indexes)
         hypre_TFree(m_mgr_CF_indexes, HYPRE_MEMORY_HOST);
      m_mgr_CF_indexes = NULL;

      if (m_mgr_amg_solver) HYPRE_BoomerAMGDestroy(m_mgr_amg_solver);
      HYPRE_MGRDestroy(m_par_MGR_solver);
   }
   
   // Expert-only(!) settings
   m_mgr_nlevels = 1;
   m_mgr_num_cindexes = NULL;
   m_mgr_cindexes = NULL;
   HYPRE_Int mgr_bsize = 2;
   HYPRE_Int mgr_non_c_to_f = 1;
   HYPRE_Int mgr_num_relax_sweeps = 3;
   HYPRE_Int mgr_num_interp_sweeps = 2;
   HYPRE_Int mgr_gsmooth_type = 0;
   HYPRE_Int mgr_num_gsmooth_sweeps = 0;
   HYPRE_Int mgr_num_restrict_sweeps = 0;  

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
   HYPRE_ParVector par_b;  // apparently not used, but needs to be passed to HYPRE_*Setup
   HYPRE_ParVector par_x;  // apparently not used, but needs to be passed to HYPRE_*Setup

   HYPRE_MGRCreate(&m_par_MGR_solver);

   // Create and set the marker array defining the C and F points
   setMGRCFIndexes(&m_mgr_CF_indexes);
   HYPRE_MGRSetCpointsByPointMarkerArray(m_par_MGR_solver, mgr_bsize, m_mgr_nlevels,
                                         m_mgr_num_cindexes, m_mgr_cindexes, m_mgr_CF_indexes);

   // Set intermediate coarse grid strategy
   HYPRE_MGRSetNonCpointsToFpoints(m_par_MGR_solver, mgr_non_c_to_f);
   // Set relax type for single level F-relaxation and post-relaxation */
   HYPRE_MGRSetRelaxType(m_par_MGR_solver, m_MGR_relax_type);
   HYPRE_MGRSetNumRelaxSweeps(m_par_MGR_solver, mgr_num_relax_sweeps);
   // Set interpolation type
   HYPRE_MGRSetRestrictType(m_par_MGR_solver, m_MGR_restrict_type);
   HYPRE_MGRSetNumRestrictSweeps(m_par_MGR_solver, mgr_num_restrict_sweeps);
   HYPRE_MGRSetInterpType(m_par_MGR_solver, m_MGR_interp_type);
   HYPRE_MGRSetNumInterpSweeps(m_par_MGR_solver, mgr_num_interp_sweeps);
   // Set print level
   HYPRE_MGRSetPrintLevel(m_par_MGR_solver, m_MGR_print_level);

   HYPRE_MGRSetGlobalSmoothType(m_par_MGR_solver, mgr_gsmooth_type);
   HYPRE_MGRSetMaxGlobalSmoothIters( m_par_MGR_solver, mgr_num_gsmooth_sweeps );

   //========== Create AMG coarse grid solver ====================

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
   HYPRE_BoomerAMGSetTol(m_mgr_amg_solver, m_MGR_amg_tol);
   HYPRE_BoomerAMGSetPrintLevel(m_mgr_amg_solver, m_MGR_amg_print_level);
   HYPRE_BoomerAMGSetStrongThreshold(m_mgr_amg_solver, m_MGR_amg_threshold);

   /* set the MGR coarse solver. Comment out to use default Coarse Grid solver in MGR */
   HYPRE_MGRSetCoarseSolver( m_par_MGR_solver, HYPRE_BoomerAMGSolve,
                             HYPRE_BoomerAMGSetup, m_mgr_amg_solver);

   //========== MGR setup ====================

   HYPRE_MGRSetup(m_par_MGR_solver, par_A, par_b, par_x);

   m_MGR_solver_allocated = true;
}


void
MBHypreSolver::setMGRCFIndexes( HYPRE_Int** a_CF_indexes )
{
   const DisjointBoxLayout& grids = m_geometry.grids();

   int num_local_dof = 0;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      num_local_dof += grids[dit].numPts();
   }
   num_local_dof *= m_nvar;

   *a_CF_indexes = hypre_CTAlloc(HYPRE_Int, num_local_dof, HYPRE_MEMORY_HOST);

   // The index ordering is determined by the Hypre convention
   HYPRE_Int* ptr = *a_CF_indexes;
   for (int part=0; part<m_coord_sys_ptr->numBlocks(); ++part) {
      for (int var=0; var<m_nvar; ++var) {
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            const Box& box = grids[dit];
            int block_number = m_coord_sys_ptr->whichBlock(box);

            if ( block_number == part ) {
               for (BoxIterator bit(box); bit.ok(); ++bit) {
                  *ptr++ = var;
               }
            }
         }
      }
   }
} 


void MBHypreSolver::AIRSetup( const HYPRE_SStructMatrix& a_matrix )
{
   CH_TIME("MBHypreSolver::AIRSetup");

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

   HYPRE_BoomerAMGSetRestriction(m_par_AMG_solver, m_AIR_restrict_type);

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
   
   HYPRE_BoomerAMGSetAggNumLevels(m_par_AMG_solver, 0);

   HYPRE_BoomerAMGSetMaxCoarseSize(m_par_AMG_solver, 0);

   //   HYPRE_BoomerAMGSetStrongThreshold(m_par_AMG_solver, 0.25);
   //   HYPRE_BoomerAMGSetStrongThresholdR(m_par_AMG_solver, 0.25);
   //   HYPRE_BoomerAMGSetCoarsenType(m_par_AMG_solver, 6);  // Falgout coarsening
   //   HYPRE_BoomerAMGSetCycleRelaxType(m_par_AMG_solver, 3, 3); // hybrid Gauss-Seidel on coarsest level
   //   HYPRE_BoomerAMGSetCycleRelaxType(m_par_AMG_solver, 9, 3); // Gaussian elimination on coarsest level
   //   HYPRE_BoomerAMGSetCycleType(m_par_AMG_solver, 1);  // V = 1, W = 2
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
   if ( m_MGR_schur || m_MGR_schur_1 ) {
      HYPRE_MGRSetMaxIter(m_par_MGR_solver, 1);
   }
   else {
      HYPRE_MGRSetMaxIter(m_par_MGR_solver, a_max_iter);
   }
   HYPRE_MGRSetTol(m_par_MGR_solver, a_tol);
   
   HYPRE_ParCSRMatrix    par_A;
   HYPRE_SStructMatrixGetObject(a_matrix, (void **) &par_A);
   HYPRE_ParVector       par_b;
   HYPRE_SStructVectorGetObject(a_b, (void **) &par_b);
   HYPRE_ParVector       par_x;
   HYPRE_SStructVectorGetObject(a_x, (void **) &par_x);

   // Do the solve
   HYPRE_MGRSolve(m_par_MGR_solver, par_A, par_b, par_x);
      
   if ( m_MGR_amg_verbose ) {

      int num_iterations;
      HYPRE_BoomerAMGGetNumIterations(m_mgr_amg_solver, &num_iterations);

      HYPRE_Real final_res_norm;
      HYPRE_BoomerAMGGetFinalRelativeResidualNorm(m_mgr_amg_solver, &final_res_norm);

      if ( procID()==0 ) {
         cout << "        --> MGR AMG coarse solver residual = " << final_res_norm << " after "
              << num_iterations << " iterations" << endl;
      }

   }

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
