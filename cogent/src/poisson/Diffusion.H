#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include "FieldSolver.H"
#include "MultiBlockLevelExchangeAverage.H"
#ifdef with_petsc
#include "MBPETScSolver.H"
#else
#include "MBHypreSolver.H"
#endif

#include "MagFluxAlignedMBLEXCenter.H"

#include "GridFunction.H"
#include "GridFunctionLibrary.H"

#include "NamespaceHeader.H"

/// Diffusion operator class.
/**
 * This operator class solves the mapped gyrokinetic Poisson equation:
 *
 * \f[
 *  {\nabla} \cdot \left ( N^T D N/J  {\nabla} \Phi \right ) = \rho,
 * \f]
*/
class Diffusion
   : public FieldSolver
{
public:

   /// Constructor with initialization.
   /**
    *  @param[in] pp the input database.
    *  @param[in] geom the geometry data.
    */
   Diffusion( const ParmParse&  pp,
              const MagGeom&    geom );
      
   /// Destructor.
   /**
    */
   virtual ~Diffusion();

   virtual void updateBoundaries( const PotentialBC&  bc );

   /// Construct the operator coefficients
   /**
    * Computes the operator coefficients.
    *
    */

   virtual void setOperatorCoefficients( const double                       scalar_factor,
                                         const RefCountedPtr<GridFunction>  coefficient_function,
                                         const PotentialBC&                 bc,
                                         const bool                         update_preconditioner );

   void computeCoefficients( const double                       scalar_factor,
                             const RefCountedPtr<GridFunction>  coefficient_function,
                             LevelData<FluxBox>&                unmapped_coefficients,
                             LevelData<FluxBox>&                mapped_coefficients );

   virtual void setPreconditionerConvergenceParams( const double tol,
                                                    const int    max_iter,
                                                    const double precond_tol,
                                                    const int    precond_max_iter );

   virtual void solvePreconditioner( const LevelData<FArrayBox>& in,
                                     LevelData<FArrayBox>&       out );

   virtual void multiplyCoefficients( LevelData<FluxBox>& data,
                                      const bool mapped_coeff ) const;

   /// ParmParse object name
   /**
    * ParmParse object name.  Used to locate input file parameters of interest to
    * this class.
    */
   static const char* pp_name;
   
protected:

#ifdef with_petsc
   MBPETScSolver* m_preconditioner;
#else
   MBHypreSolver* m_preconditioner;
#endif

   LevelData<FluxBox> m_mapped_coefficients;
   LevelData<FluxBox> m_unmapped_coefficients;
   LevelData<FArrayBox> m_mapped_volume_reciprocal;
   
   MagFluxAlignedMBLEXCenter* m_mblx_ptr;
};



#include "NamespaceFooter.H"

#endif
