#ifndef _GKVLASOVAMG_
#define _GKVLASOVAMG_

#include "ParmParse.H"
#include "PhaseGeom.H"
#include "SparseCoupling.H"
#include "IntVectSet.H"
#include "GKVlasovMBCoupling.H"

// Hypre headers
#include "krylov.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_sstruct_mv.h"

#include "NamespaceHeader.H"

class GKVlasovAMG
{
  public:

    GKVlasovAMG( const ParmParse&  pp,
                 const PhaseGeom&  geometry );

    ~GKVlasovAMG();

    void constructMatrix( const string&                                 face_avg_type,
                          const double                                  shift,
                          const LevelData<BaseFab<Vector<IntVect> > >&  structured_couplings,
                          const LevelData<BaseFab<Vector<Real> > >&     weights );

    void addUnstructuredMatrixEntries( const LevelData<BaseFab<Vector<IntVect> > >&  structured_couplings,
                                       const LevelData<BaseFab<Vector<Real> > >&     structured_weights,
                                       const LayoutData< BaseFab<IntVectSet> >&      unstructured_couplings,
                                       HYPRE_SStructMatrix&                          matrix ) const;

    int findStructuredEntry( const IntVect&  iv_center,
                             const IntVect&  iv ) const;
   
    int findEntry( const IntVect&     iv_center,
                   const IntVectSet&  unstructured_ivs,
                   const IntVect&     iv ) const;
   
    void multiplyMatrix( const LevelData<FArrayBox>&  in,
                         LevelData<FArrayBox>&        out ) const;

    void setConvergenceParams( const double  method_tol,
                               const int     method_max_iter,
                               const bool    method_verbose );

    void solve( const LevelData<FArrayBox>&  rhs,
                LevelData<FArrayBox>&        solution );
   
  protected:

    IntVectSet getInterBlockCoupledCells( const int   block_number,
                                          const int   radius,
                                          const Box&  box ) const;
   
    void getUnstructuredCouplings( int                                 radius,
                                   LayoutData< BaseFab<IntVectSet> >&  unstructured_couplings ) const;

    IntVectSet getUnstructuredCouplingsForCell( const Vector< BaseFab<IntVectSet>* >&  a_coupling,
                                                const Box&                             a_stencil_box,
                                                const int                              a_dir,
                                                const IntVect&                         a_iv ) const;
   
    void addUnstructuredGraphEntries( const int                                 radius,
                                      const LayoutData< BaseFab<IntVectSet> >&  unstructured_coupling,
                                      HYPRE_SStructGraph&                       graph ) const;
   
    void createHypreData();

    void AMGSetup( const HYPRE_SStructMatrix&  matrix );

    void AMG( const HYPRE_SStructMatrix&  matrix,
              const HYPRE_SStructVector&  b,
              const double                tol,
              const int                   max_iter,
              const bool                  verbose,
              const HYPRE_SStructVector&  x );

    void copyToHypreVector( const LevelData<FArrayBox>&  in,
                            HYPRE_SStructVector&         out ) const;

    void copyFromHypreVector( const HYPRE_SStructVector&  in,
                              LevelData<FArrayBox>&       out ) const;

  private:

   const PhaseGeom& m_geometry;
   LevelData<FArrayBox> m_volume;
   BlockBaseRegister<BaseFab<SparseCoupling> >* m_block_register;

   GKVlasovMBCoupling* m_mb_coupling; 

   // Hypre data structures

   bool m_hypre_allocated;
   int m_hypre_object_type;
   HYPRE_SStructGrid m_grid;
   HYPRE_SStructStencil m_A_stencil;
   HYPRE_SStructGraph m_A_graph;
   HYPRE_SStructMatrix m_A;
   mutable HYPRE_SStructVector m_b;
   mutable HYPRE_SStructVector m_x;

   HYPRE_Int m_restrict_type;  /* 0: P^T, 1: AIR, 2: AIR-2, >=3 nAIR */
   
   int m_stencil_radius;
   int m_stencil_size;
   int** m_offsets;

   LayoutData< BaseFab<IntVectSet> > m_A_unstructured_coupling;

   HYPRE_Solver m_par_AMG_solver;
   bool m_AMG_solver_allocated;

   // Convergence and verbosity parameters
   double m_tol;
   int m_max_iter;
   bool m_verbose;
   bool m_convergence_params_set;

   bool m_flux_average[CFG_DIM];
};


#include "NamespaceFooter.H"
#endif
