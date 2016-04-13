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
    *  @param[in] geom the geometry
    *  @param[in] discretization_ordr discretization order
    */
   MBHypreSolver( const MultiBlockLevelGeom&   geom,
                  const int                    discretization_order);
      
   /// Destructor.
   /**
    */
   virtual ~MBHypreSolver();

   virtual void constructMatrixGeneral( LevelData<FArrayBox>& alpha_coefficient,
                                        LevelData<FluxBox>&   tensor_coefficient,
                                        LevelData<FArrayBox>& beta_coefficient,
                                        const PotentialBC&    bc );

   virtual void multiplyMatrix( const LevelData<FArrayBox>& in,
                                LevelData<FArrayBox>&       out ) const;

   virtual void solve( const LevelData<FArrayBox>&  rhs,
                       LevelData<FArrayBox>&        solution,
                       bool                         a_homogeneous_bcs );

   virtual void dumpMatrix(const string& file_name) const;

   void setParams( const string& method,
                   const double  method_tol,
                   const int     method_max_iter,
                   const bool    method_verbose,
                   const string& precond_method,
                   const double  precond_tol,
                   const int     precond_max_iter,
                   const bool    precond_verbose );

   static const std::string pp_name;

protected:

   void createHypreData();

   void destroyHypreData();

   int findHypreEntry(const Box&        stencil_box,
                      const IntVectSet& unstructured_ivs,
                      const IntVect&    iv) const;

   void addUnstructuredGraphEntries( const int                                radius,
                                     const LayoutData< BaseFab<IntVectSet> >& unstructured_coupling,
                                     HYPRE_SStructGraph&                      graph ) const;

   void addUnstructuredMatrixEntries( const LevelData<FArrayBox>&               alpha_coefficient,
                                      const LevelData<FluxBox>&                 tensor_coeffcient,
                                      const PotentialBC&                        bc,
                                      FArrayBox&                                stencil_values,
                                      const bool                                fourthOrder,
                                      const bool                                dropOrder,
                                      const LayoutData< BaseFab<IntVectSet> >&  unstructured_coupling,
                                      Vector< Vector<CoDim1Stencil> >&          codim1_stencils,
                                      Vector< Vector<CoDim2Stencil> >&          codim2_stencils,
                                      HYPRE_SStructMatrix&                      matrix ) const;

   void constructHypreMatrix( LevelData<FArrayBox>&               alpha_coefficient, 
                              LevelData<FluxBox>&                 tensor_coefficient,
                              LevelData<FArrayBox>&               beta_coefficient,
                              const PotentialBC&                  bc,
                              HYPRE_SStructGraph&                 graph,
                              FArrayBox&                          stencil_values,
                              const int                           diagonal_offset,
                              LayoutData< BaseFab<IntVectSet> >&  unstructured_coupling,
                              const bool                          fourthOrder,
                              const bool                          dropOrder,
                              HYPRE_SStructMatrix&                matrix,
                              LevelData<FArrayBox>&               rhs_from_bc ) const;

   void AMG( const HYPRE_SStructMatrix&  matrix,
             const HYPRE_SStructVector&  b,
             const double                tol,
             const int                   max_iter,
             const bool                  verbose,
             const HYPRE_SStructVector&  x ) const;

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

   // Method, preconditioner and parameters

   string m_method;
   double m_method_tol;
   int m_method_max_iter;
   bool m_method_verbose;

   string m_precond_method;
   double m_precond_tol;
   int m_precond_max_iter;
   bool m_precond_verbose;

   bool m_params_set;

   // Hypre data structures

   bool m_hypre_allocated;
   HYPRE_SStructGrid m_grid;
   HYPRE_SStructStencil m_A_stencil;
   HYPRE_SStructGraph m_A_graph;
   HYPRE_SStructMatrix m_A;
   mutable HYPRE_SStructVector m_b;
   mutable HYPRE_SStructVector m_x;
   int m_hypre_object_type;

   FArrayBox m_A_stencil_values;
   int m_A_diagonal_offset;

   LayoutData< BaseFab<IntVectSet> > m_A_unstructured_coupling;
};



#include "NamespaceFooter.H"

#endif
