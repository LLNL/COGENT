#ifndef _MB_PETSC_SOLVER_H_
#define _MB_PETSC_SOLVER_H_

#include "MBSolver.H"
#include "PotentialBC.H"

// PETSc headers
#include <petscksp.h>

#include "NamespaceHeader.H"



class MBPETScSolver
   : public MBSolver
{
public:

   /// Constructor with initialization.
   /**
    *  @param[in] pp the input database.
    *  @param[in] volume cell volumes.
    *  @param[in] geom the geometry
    */
   MBPETScSolver( const MultiBlockLevelGeom&   geom,
                  const LevelData<FArrayBox>&  volume,
                  const int                    discretization_order);
      
   /// Destructor.
   /**
    */
   virtual ~MBPETScSolver();

   virtual void constructMatrixGeneral( LevelData<FArrayBox>& alpha_coefficient,
                                        LevelData<FluxBox>&   diffusion_coefficient,
                                        LevelData<FArrayBox>& beta_coefficient,
                                        const PotentialBC&    bc );

   virtual void multiplyMatrix( const LevelData<FArrayBox>& in,
                                LevelData<FArrayBox>&       out ) const;

   virtual void solve( const LevelData<FArrayBox>&  rhs,
                       LevelData<FArrayBox>&        solution,
                       bool                         homogeneous_bcs );

   virtual void dumpMatrix(const string& file_name) const;

   void setParams( const string& method,
                   const double  method_tol,
                   const int     method_max_iter,
                   const bool    method_verbose,
                   const string& precond_method,
                   const double  precond_tol,
                   const int     precond_max_iter,
                   const bool    precond_verbose );

protected:

   /// Create the PETSc data structures
   /**
    * Create the PETSc data structures
    */
   void createPETScData();

   /// Destroy the PETSc data structures
   /**
    * Destroy the PETSc data structures
    */
   void destroyPETScData();

   void createAMGSolver();

   void createGMRESSolver();

   void destroySolver();

   int numLocalGIDs() const {return m_ending_gid - m_starting_gid + 1;}

   int getGID( const IntVect& a_iv ) const;

   int getGID( const int block_number, const IntVect& iv ) const;

   void getNumNonzeros( const int                                a_radius,
                        const LayoutData< BaseFab<IntVectSet> >& a_unstructured_couplings,
                        int**                                    a_d_nnz,
                        int**                                    a_o_nnz ) const;

   void addUnstructuredMatrixEntries( const LevelData<FArrayBox>&              alpha_coefficient,
                                      const LevelData<FluxBox>&                diffusion_coeffcient,
                                      const PotentialBC&                       bc,
                                      FArrayBox&                               stencil_values,
                                      const bool                               fourthOrder,
                                      const bool                               dropOrder,
                                      const LayoutData< BaseFab<IntVectSet> >& unstructured_coupling,
                                      Vector< Vector<CoDim1Stencil> >&         codim1_stencils,
                                      Vector< Vector<CoDim2Stencil> >&         codim2_stencils,
                                      Mat&                                     matrix ) const;

   void constructPETScMatrix( LevelData<FArrayBox>&              alpha_coefficient, 
                              LevelData<FluxBox>&                diffusion_coefficient,
                              LevelData<FArrayBox>&              beta_coefficient,
                              const PotentialBC&                 bc,
                              FArrayBox&                         stencil_values,
                              const int                          diagonal_offset,
                              LayoutData< BaseFab<IntVectSet> >& unstructured_coupling,
                              const bool                         fourthOrder,
                              const bool                         dropOrder,
                              Mat&                               matrix,
                              LevelData<FArrayBox>&              rhs_from_bc ) const;

   void AMG( const Mat&   matrix,
             const Vec&   b,
             const double tol,
             const int    max_iter,
             const bool   verbose,
             const Vec&   x ) const;

   void AMG_Preconditioned_GMRES( const Mat&   matrix,
                                  const Vec&   b,
                                  const double tol,
                                  const int    max_iter,
                                  const double amg_tol,
                                  const int    amg_max_iter,
                                  const bool   verbose,
                                  const Vec&   x ) const;

   IntVect checkPeriodic() const;

   void copyToPETScVector( const LevelData<FArrayBox>& in,
                           Vec&                        out ) const;

   void copyFromPETScVector( const Vec&            in,
                             LevelData<FArrayBox>& out ) const;

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

   // PETSc data structures
   Mat m_A;
   Vec m_b;
   mutable Vec m_x;
   Vec m_in_vector;
   Vec m_out_vector;
   bool m_petsc_allocated;
   bool m_amg_solver_allocated;
   bool m_gmres_solver_allocated;
   int m_starting_gid;
   int m_ending_gid;
   PC m_pc;
   KSP m_amg_ksp;
   KSP m_gmres_ksp;

   FArrayBox m_A_stencil_values;
   int m_A_diagonal_offset;

   LayoutData< BaseFab<IntVectSet> > m_A_unstructured_coupling;


   Vector< Vector<Box> > m_block_boxes;
   Vector< Vector<int> > m_block_boxes_starting_gid;
   LevelData< BaseFab<int> > m_gids;
};



#include "NamespaceFooter.H"

#endif
