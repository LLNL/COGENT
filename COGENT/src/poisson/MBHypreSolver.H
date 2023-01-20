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

   virtual void dumpVector(const string& file_name) const;

   virtual void setMethodParams( const ParmParse&  pp );

   void setBasicAMGMethod()
      {m_method = "AMG";
       m_AMG_coarsen_type = 6;
       m_AMG_strong_threshold = 0.25;
       m_AMG_old_setup = false;}

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
                                             const EllipticOpBC&    bc) override;

  virtual void addAdvectionMatrixBlockGeneral( const int                   a_block_row,
                                               const int                   a_block_column,
                                               const LevelData<FArrayBox>& a_alpha_coefficient,
                                               const LevelData<FluxBox>&   a_vector_coefficient,
                                               const EllipticOpBC&         a_bc) override;

   virtual void addAdvectionMatrixBlockCoupled( const int                                    a_block_row,
                                                const int                                    a_block_column,
                                                const LevelData<FArrayBox>&                  a_alpha_coefficient,
                                                const LevelData<BaseFab<Vector<IntVect> > >& a_structured_couplings,
                                                const LevelData<BaseFab<Vector<Real> > >&    a_structured_weights,
                                                const EllipticOpBC&                          a_bc,
                                                const FArrayBox&                             a_stencil_values);

   void createHypreData();

   void constructAdvectionOffsets();

   virtual void createGraph( const Box& stencil_box,
                             const int  radius );

   void destroyHypreData();

   int findStructuredEntry( const int      a_block_column,
                            const Box&     a_stencil_box,
                            const IntVect& a_iv) const;

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

   void addAdvectionUnstructuredMatrixEntries( const int                                     block_row,
					       const int                                     block_column,
                                               const LevelData<FArrayBox>&                   alpha_coefficient,
					       const LevelData<BaseFab<Vector<IntVect> > >&  structured_couplings,
					       const LevelData<BaseFab<Vector<Real> > >&     structured_weights,
					       const LayoutData< BaseFab<IntVectSet> >&      unstructured_couplings,
					       HYPRE_SStructMatrix&                          matrix ) const;

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

   void AMGSetup( const HYPRE_SStructMatrix& matrix,
                  const HYPRE_SStructVector& x,
                  const HYPRE_SStructVector& b);

   void ILUSetup( const HYPRE_SStructMatrix& matrix );

   void MGRSetup( const HYPRE_SStructMatrix& matrix );

   void setMGRCFIndexes( HYPRE_Int** CF_indexes );

   void AIRSetup( const HYPRE_SStructMatrix& matrix );

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
   LevelData<FArrayBox> m_advection_rhs_from_bc;

   // Hypre data structures

   int** m_advection_offsets;
   int m_advection_stencil_radius;
   int m_advection_stencil_size;
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
   int m_MGR_print_level;
   int m_MGR_amg_print_level;
   int m_MGR_cpoint;
   int m_MGR_amg_max_iter;
   double m_MGR_amg_tol;
   bool m_MGR_schur;
   bool m_MGR_schur_1;
   double m_MGR_amg_threshold;
   bool m_MGR_amg_verbose;
   int m_MGR_interp_type;
   int m_MGR_restrict_type;
   int m_MGR_relax_type;

   HYPRE_Int* m_mgr_num_cindexes;
   HYPRE_Int** m_mgr_cindexes;
   HYPRE_Int* m_mgr_CF_indexes;
   HYPRE_Int m_mgr_nlevels;
   HYPRE_Solver m_mgr_amg_solver;
   
   // AIR parameters
   HYPRE_Int m_AIR_restrict_type;  /* 0: P^T, 1: AIR, 2: AIR-2, >=3 nAIR */

   bool m_matrix_initialized;
   bool m_matrix_finalized;
};



#include "NamespaceFooter.H"

#endif
