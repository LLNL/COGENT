#ifndef _GKVLASOVAMG_
#define _GKVLASOVAMG_

#include "ParmParse.H"
#include "MultiBlockCoordSys.H"

// Hypre headers
#include "krylov.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_sstruct_mv.h"

#include "NamespaceHeader.H"

class GKVlasovAMG
{
  public:

    GKVlasovAMG( const ParmParse& pp );

    ~GKVlasovAMG() {};

    void constructMatrix( const MultiBlockCoordSys*               coord_sys,
                          const string&                           face_avg_type,
                          LevelData<BaseFab<Vector<IntVect> > >&  couplings,
                          LevelData<BaseFab<Vector<Real> > >&     weights );

    int findHypreEntry( const Box&      stencil_box,
                        const IntVect&  iv ) const;
   
    void multiplyMatrix( const LevelData<FArrayBox>&  in,
                         LevelData<FArrayBox>&        out ) const;

    void setConvergenceParams( const double  method_tol,
                               const int     method_max_iter,
                               const bool    method_verbose );

    void solve( const LevelData<FArrayBox>&  rhs,
                LevelData<FArrayBox>&        solution,
                bool                         homogeneous_bcs );
   
  protected:

    void createHypreData( const DisjointBoxLayout&   grids,
                          const MultiBlockCoordSys*  coord_sys,
                          const string&              face_avg_type );

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

   const MultiBlockCoordSys* m_coord_sys_ptr;

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
   
   int m_A_diagonal_offset;
   int m_stencil_radius;

   HYPRE_Solver m_par_AMG_solver;
   bool m_AMG_solver_allocated;

   int m_discretization_order;

   // Convergence and verbosity parameters
   double m_tol;
   int m_max_iter;
   bool m_verbose;
   bool m_convergence_params_set;
};


#include "NamespaceFooter.H"
#endif
