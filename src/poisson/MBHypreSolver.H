#ifndef _MBHYPRESOLVER_H_
#define _MBHYPRESOLVER_H_

#include "MBSolver.H"
#include "PotentialBC.H"

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
    *  @param[in] pp the input database.
    *  @param[in] geom the geometry
    */
   MBHypreSolver( const MagGeom&   geom,
                  const int        discretization_order);
      
   /// Destructor.
   /**
    */
   virtual ~MBHypreSolver();

   virtual void multiplyMatrix( const LevelData<FArrayBox>& in,
                                LevelData<FArrayBox>&       out );

   virtual void updateMatrixDiagonals( const LevelData<FArrayBox>& diagonal_increment );

   void constructMatrix( LevelData<FluxBox>&   diffusion_coefficient,
                         const PotentialBC&    bc,
                         const bool            poloidally_average = false );

   void constructMatrix( LevelData<FArrayBox>& alpha_coefficient,
                         LevelData<FluxBox>&   diffusion_coefficient,
                         const PotentialBC&    bc,
                         const bool            poloidally_average = false );

   void constructMatrix( LevelData<FluxBox>&   diffusion_coefficient,
                         LevelData<FArrayBox>& beta_coefficient,
                         const PotentialBC&    bc,
                         const bool            poloidally_average = false );

   void constructMatrix( LevelData<FArrayBox>& alpha_coefficient,
                         LevelData<FluxBox>&   diffusion_coefficient,
                         LevelData<FArrayBox>& beta_coefficient,
                         const PotentialBC&    bc,
                         const bool            poloidally_average = false );

   void dumpMatrix(const string& file_name) const;

   void solveWithMultigrid( const LevelData<FArrayBox>& rhs,
                            const double                tol,
                            const int                   max_iter,
                            const bool                  verbose,
                            LevelData<FArrayBox>&       solution );

   /// Perform multigrid preconditioned GMRES linear system solve
   /**
    * Perform multigrid preconditioned GMRES lienar system solve
    *
    * @param[in]  rhs       Right-hand side vector
    * @param[out] solution  Solution vector
    */
   int solveWithGMRES(const LevelData<FArrayBox>& rhs,
                      const double                tol,
                      const int                   max_iter,
                      const double                amg_tol,
                      const int                   amg_max_iter,
                      const bool                  verbose,
                      LevelData<FArrayBox>&       solution);

   /// ParmParse object name
   /**
    * ParmParse object name.  Used to locate input file parameters of interest to
    * this class.
    */
   static char* pp_name;

protected:

   /// Create the Hypre data structures
   /**
    * Create the Hypre data structures
    */
   void createHypreData();

   /// Destroy the Hypre data structures
   /**
    * Destroy the Hypre data structures
    */
   void destroyHypreData();

   void addUnstructuredGraphEntries( const int                                radius,
                                     const LayoutData< BaseFab<IntVectSet> >& unstructured_coupling,
                                     HYPRE_SStructGraph&                      graph ) const;

   int findHypreEntry(const Box&        stencil_box,
                      const IntVectSet& unstructured_ivs,
                      const IntVect&    iv) const;

   void addUnstructuredMatrixEntries( const LevelData<FArrayBox>&              alpha_coefficient,
                                      const LevelData<FluxBox>&                diffusion_coeffcient,
                                      const PotentialBC&                       bc,
                                      FArrayBox&                               stencil_values,
                                      const bool                               fourthOrder,
                                      const bool                               dropOrder,
                                      const LayoutData< BaseFab<IntVectSet> >& unstructured_coupling,
                                      HYPRE_SStructMatrix&                     matrix ) const;

   void constructMatrix( LevelData<FArrayBox>&              alpha_coefficient, 
                         LevelData<FluxBox>&                diffusion_coefficient,
                         LevelData<FArrayBox>&              beta_coefficient,
                         const PotentialBC&                 bc,
                         HYPRE_SStructGraph&                graph,
                         const int                          stencil_size,
                         FArrayBox&                         stencil_values,
                         LevelData<FArrayBox>&              radial,
                         LevelData<FArrayBox>&              diagonal,
                         const int                          diagonal_offset,
                         LayoutData< BaseFab<IntVectSet> >& unstructured_coupling,
                         const bool                         average_poloidally,
                         const bool                         fourthOrder,
                         const bool                         dropOrder,
                         HYPRE_SStructMatrix&               matrix ) const;

   void MultigridSolve_Hypre( const HYPRE_SStructMatrix& matrix,
                              const HYPRE_SStructVector& b,
                              const double               tol,
                              const int                  max_iter,
                              const bool                 verbose,
                              const HYPRE_SStructVector& x ) const;

   int GMRESSolve_Hypre( const HYPRE_SStructMatrix&  matrix,
                         const HYPRE_SStructMatrix&  precond,
                         const HYPRE_SStructVector&  b,
                         const double                tol,
                         const int                   max_iter,
                         const double                amg_tol,
                         const int                   amg_max_iter,
                         const bool                  verbose,
                         const HYPRE_SStructVector&  x ) const;

   void copyToHypreVector( const LevelData<FArrayBox>& in,
                           HYPRE_SStructVector&        out );

   void copyFromHypreVector( const HYPRE_SStructVector& in,
                             LevelData<FArrayBox>&      out ) const;

   // Hypre data structures

   HYPRE_SStructGrid m_grid;
   HYPRE_SStructStencil m_A_stencil;
   HYPRE_SStructGraph m_A_graph;
   HYPRE_SStructMatrix m_A;
   HYPRE_SStructVector m_b;
   HYPRE_SStructVector m_x;
   bool m_hypre_allocated;
   int m_hypre_object_type;

   int m_A_stencil_size;
   int m_A_radius;
   FArrayBox m_A_stencil_values;
   int m_A_diagonal_offset;

   LayoutData< BaseFab<IntVectSet> > m_A_unstructured_coupling;
};



#include "NamespaceFooter.H"

#endif
