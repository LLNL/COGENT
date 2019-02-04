#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include "EllipticOp.H"
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
   : public EllipticOp
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

   virtual void updateBoundaries( const EllipticOpBC&  bc );

   /// Construct the operator coefficients
   /**
    * Computes the operator coefficients.
    *  @param[in] pointwise face-centered unmapped_coefficients.
    *  @param[in] pointwise face-centered mapped_coefficients.
    *  @param[in] boundary conditions
    *
    */

   virtual void setOperatorCoefficients( const LevelData<FluxBox>&           unmapped_coefficients_fc,
                                         const LevelData<FluxBox>&           mapped_coefficients_fc,
                                         const EllipticOpBC&                 bc );

   virtual void updateImExPreconditioner( const LevelData<FArrayBox>&      mshift,
                                          const EllipticOpBC&              bc );

   virtual void setPreconditionerConvergenceParams( const double  tol,
                                                    const int     max_iter,
                                                    const double  precond_tol,
                                                    const int     precond_max_iter );

   virtual void setImExPreconditionerConvergenceParams( const double  tol,
                                                        const int     max_iter,
                                                        const double  precond_tol,
                                                        const int     precond_max_iter );

   virtual void solvePreconditioner( const LevelData<FArrayBox>&  in,
                                     LevelData<FArrayBox>&        out );

   virtual void solveImExPreconditioner( const LevelData<FArrayBox>&  in,
                                         LevelData<FArrayBox>&        out );

   virtual void multiplyCoefficients( LevelData<FluxBox>&  data,
                                      const bool           mapped_coeff ) const;

   /// ParmParse object name
   /**
    * ParmParse object name.  Used to locate input file parameters of interest to
    * this class.
    */
   static const char* pp_name;
   
protected:

   bool m_coefficients_defined;

#ifdef with_petsc
   MBPETScSolver* m_preconditioner;
   MBPETScSolver* m_imex_preconditioner;
#else
   MBHypreSolver* m_preconditioner;
   MBHypreSolver* m_imex_preconditioner;
#endif

   LevelData<FluxBox> m_mapped_coefficients; //face-averaged
   LevelData<FluxBox> m_unmapped_coefficients; //face-centered
   LevelData<FArrayBox> m_volume_reciprocal;
   
   MagFluxAlignedMBLEXCenter* m_mblx_ptr;
};



#include "NamespaceFooter.H"

#endif
