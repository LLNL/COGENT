#ifndef _MBHYPRESOLVER_H_
#define _MBHYPRESOLVER_H_

#include "MBSolver.H"
#include "EllipticOpBC.H"

// Hypre headers
#include "krylov.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_sstruct_mv.h"

#include "NamespaceHeader.H"


class MBHypreSolver
   : public MBSolver
{
public:

   /// Constructor with initialization.
   /**
    *  @param[in] geom the geometry
    *  @param[in] nvar number of equations/unknowns
    *  @param[in] discretization_order discretization order
    *  @param[in] mblex_ptr level exchange object pointer (leave NULL
    *                       unless you know why you don't want to)
    */
   MBHypreSolver( const MultiBlockLevelGeom&      geom,
                  const int                       nvar,
                  const int                       discretization_order,
                  MultiBlockLevelExchangeCenter*  mblex_ptr = NULL );

   /// Destructor.
   /**
    */
   virtual ~MBHypreSolver();

   virtual void multiplyMatrix( const LevelData<FArrayBox>& in,
                                LevelData<FArrayBox>&       out ) const;

   virtual void solve( const LevelData<FArrayBox>&  rhs,
                       LevelData<FArrayBox>&        solution,
                       bool                         a_homogeneous_bcs );

   virtual void dumpMatrix(const string& file_name) const;

   virtual void setMethodParams( const ParmParse&  pp );

   static const std::string pp_name;

   void extrapGhosts( const EllipticOpBC&   a_bc,
                      const bool            a_fourth_order,
                      LevelData<FArrayBox>& a_data ) const;

   virtual void initializeMatrix();

   virtual void finalizeMatrix();

   protected:

   void parseParameters( ParmParse& a_pp );

   virtual void constructMatrixGeneral( LevelData<FArrayBox>&  alpha_coefficient,
                                        LevelData<FluxBox>&    tensor_coefficient,
                                        LevelData<FArrayBox>&  beta_coefficient,
                                        const EllipticOpBC&    bc );

   virtual void constructMatrixBlockGeneral( const int              block_row,
                                             const int              block_column,
                                             LevelData<FArrayBox>&  alpha_coefficient,
                                             LevelData<FluxBox>&    tensor_coefficient,
                                             LevelData<FArrayBox>&  beta_coefficient,
                                             const EllipticOpBC&    bc );

   void createHypreData();

   virtual void createGraph( const Box& stencil_box,
                             const int  radius );

   void destroyHypreData();

   int findHypreEntry( const int          block_column,
                       const Box&         stencil_box,
                       const IntVectSet&  unstructured_ivs,
                       const IntVect&     iv ) const;

   void addUnstructuredGraphEntries( const int                                radius,
                                     const LayoutData< BaseFab<IntVectSet> >& unstructured_coupling,
                                     const int                                var,
                                     const int                                to_var,
                                     HYPRE_SStructGraph&                      graph ) const;

   void addUnstructuredMatrixEntries( const int                                 block_row,
                                      const int                                 block_column,
                                      const LevelData<FArrayBox>&               alpha_coefficient,
                                      const LevelData<FluxBox>&                 tensor_coeffcient,
                                      const EllipticOpBC&                       bc,
                                      FArrayBox&                                stencil_values,
                                      const bool                                fourthOrder,
                                      const LayoutData< BaseFab<IntVectSet> >&  unstructured_coupling,
                                      Vector< Vector<CoDim1Stencil> >&          codim1_stencils,
                                      Vector< Vector<CoDim2Stencil> >&          codim2_stencils,
                                      HYPRE_SStructMatrix&                      matrix ) const;

   void constructHypreMatrix( LevelData<FArrayBox>&               alpha_coefficient, 
                              LevelData<FluxBox>&                 tensor_coefficient,
                              LevelData<FArrayBox>&               beta_coefficient,
                              const EllipticOpBC&                 bc,
                              FArrayBox&                          stencil_values,
                              const int                           diagonal_offset,
                              LayoutData< BaseFab<IntVectSet> >&  unstructured_coupling,
                              const bool                          fourthOrder,
                              LevelData<FArrayBox>&               rhs_from_bc );

   void constructHypreMatrixBlock( const int                           block_row,
                                   const int                           block_column,
                                   LevelData<FArrayBox>&               alpha_coefficient, 
                                   LevelData<FluxBox>&                 tensor_coefficient,
                                   LevelData<FArrayBox>&               beta_coefficient,
                                   const EllipticOpBC&                 bc,
                                   HYPRE_SStructGraph&                 graph,
                                   FArrayBox&                          stencil_values,
                                   const int                           diagonal_offset,
                                   LayoutData< BaseFab<IntVectSet> >&  unstructured_coupling,
                                   const bool                          fourthOrder,
                                   HYPRE_SStructMatrix&                matrix,
                                   LevelData<FArrayBox>&               rhs_from_bc );

   int blockIndex(const int block_row, const int block_column) const
                 {return block_row*m_nvar + block_column;} 

   void AMGSetup( const HYPRE_SStructMatrix& matrix );

   void ILUSetup( const HYPRE_SStructMatrix& matrix );

   void MGRSetup( const HYPRE_SStructMatrix& matrix );

   void setMGRCFIndexes( HYPRE_Int** CF_indexes );

   void AMG( const HYPRE_SStructMatrix&  matrix,
             const HYPRE_SStructVector&  b,
             const double                tol,
             const int                   max_iter,
             const bool                  verbose,
             const HYPRE_SStructVector&  x );

   void ILU( const HYPRE_SStructMatrix&  matrix,
             const HYPRE_SStructVector&  b,
             const double                tol,
             const int                   max_iter,
             const bool                  verbose,
             const HYPRE_SStructVector&  x );

   void MGR( const HYPRE_SStructMatrix&  matrix,
             const HYPRE_SStructVector&  b,
             const double                tol,
             const int                   max_iter,
             const bool                  verbose,
             const HYPRE_SStructVector&  x );

   void AMG_preconditioned_GMRES( const HYPRE_SStructMatrix&  matrix,
                                  const HYPRE_SStructMatrix&  precond,
                                  const HYPRE_SStructVector&  b,
                                  const double                tol,
                                  const int                   max_iter,
                                  const double                amg_tol,
                                  const int                   amg_max_iter,
                                  const bool                  verbose,
                                  const HYPRE_SStructVector&  x ) const;

   void copyToHypreVector( const LevelData<FArrayBox>&  in,
                           HYPRE_SStructVector&         out ) const;

   void copyFromHypreVector( const HYPRE_SStructVector&  in,
                             LevelData<FArrayBox>&       out ) const;

   LevelData<FArrayBox> m_rhs_from_bc;

   // Hypre data structures

   bool m_hypre_allocated;
   HYPRE_SStructGrid m_grid;
   HYPRE_SStructStencil* m_A_stencil;
   HYPRE_SStructGraph m_A_graph;
   HYPRE_SStructMatrix m_A;
   mutable HYPRE_SStructVector m_b;
   mutable HYPRE_SStructVector m_x;
   int m_hypre_object_type;

   int m_nvar;
   FArrayBox m_A_stencil_values;
   int m_A_diagonal_offset;

   LayoutData< BaseFab<IntVectSet> > m_A_unstructured_coupling;

   // AMG parameters
   HYPRE_Solver m_par_AMG_solver;
   bool m_AMG_solver_allocated;
   double m_AMG_strong_threshold;
   int m_AMG_coarsen_type;
   bool m_AMG_old_setup;

   // ILU parameters
   HYPRE_Solver m_par_ILU_solver;
   bool m_ILU_solver_allocated;
   int m_ILU_type;
   int m_ILU_fill_level;
   double m_ILU_drop_threshold;
   double m_ILU_drop_threshold_array[CFG_DIM];
   int m_ILU_max_nnz_per_row;
   int m_ILU_max_iter;
   int m_ILU_max_schur_iter;
   double m_ILU_tol;
   int m_ILU_print_level;

   // MGR parameters
   HYPRE_Solver m_par_MGR_solver;
   bool m_MGR_solver_allocated;
   HYPRE_Solver m_MGR_aff_solver;
   bool m_MGR_aff_solver_allocated;
   int m_MGR_print_level;
   int m_MGR_amg_print_level;
   int m_MGR_cpoint;
   int m_MGR_amg_max_iter;
   bool m_MGR_schur;
   double m_MGR_amg_threshold;
   double m_MGR_frelax_amg_threshold;
   int m_MGR_frelax_method;
   int m_MGR_interp_type;
   int m_MGR_restrict_type;
   int m_MGR_relax_type;
   int m_MGR_aff_amg_max_iter;
   double m_MGR_aff_amg_tol;

   HYPRE_Int* m_mgr_num_cindexes;
   HYPRE_Int** m_mgr_cindexes;
   HYPRE_Int* m_mgr_CF_indexes;
   HYPRE_Int m_mgr_nlevels;
   HYPRE_Solver m_mgr_amg_solver;
   
   bool m_matrix_initialized;
   bool m_matrix_finalized;
};



#include "NamespaceFooter.H"

#endif
