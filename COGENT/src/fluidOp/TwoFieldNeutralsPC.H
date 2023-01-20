#ifndef _TWOFIELDNEUTRALSPC_H_
#define _TWOFIELDNEUTRALSPC_H_

#include "EllipticOp.H"
#include "GKPoisson.H"
#include "FluxSurface.H"
#include "MultiBlockLevelExchangeAverage.H"
#include "MBHypreSolver.H"

#include "MagFluxAlignedMBLEXCenter.H"

#include "GridFunction.H"
#include "GridFunctionLibrary.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "KineticSpecies.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include "NamespaceHeader.H"
namespace PS = PS_NAMESPACE;
namespace VEL = VEL_NAMESPACE;

/// Advection-Diffusion Two-Field Neutrals Operations
/**
 * This operator class solves the analytic Jacobian for the
 * Two Field Neutrals model. Primarily intended to be used
 * as a preconditioner.
 * This is a 2x2 matrix where the follwoing notation is used in the
 * code and comments:
 * 
 * P = | P_11 | P_12 |
 *     | P_21 | P_22 |
 *
 * and in written form:
 * P_11: topLeft
 * P_12: topRight
 * P_21: bottomLeft
 * P_22: bottomRight
 *
 * P_12/topRight is always zero.
*/
class TwoFieldNeutralsPC
   : public EllipticOp
{
public:

   /// Constructor with initialization.
   /**
    *  @param[in] pp the input database.
    *  @param[in] geom the geometry data.
    */
   TwoFieldNeutralsPC(const ParmParse&   pp,
				const MagGeom&     geom );
      
   /// Destructor.
   /**
    */
   virtual ~TwoFieldNeutralsPC();


   /// Updates the source term components of the analytic Jacobian preconditioner
   /**
    *  In the 2x2 block matrix, sets the three the source term components.
    *  The top-right component (corresponding to P_12) is always zero in the 
    *  analytic Jacobian and so is neglected.
    *
    *  @param[in] a_topLeft_beta The top-left coefficient (corresponding to P_11)
    *  @param[in] a_bottomLeft_beta The bottom-left coefficient (corresponding to P_21)
    *  @param[in] a_bottomRight_beta The bottom-right coefficient (corresponding to P_22)
    */
   virtual void updateSourceTermComponents( const LevelData<FArrayBox>& a_topLeft_beta,
                                            const LevelData<FArrayBox>& a_bottomLeft_beta,
                                            const LevelData<FArrayBox>& a_bottomRight_beta );

   /// Updates the top-left coefficients of the analytic Jacobian preconditioner, corresponding to P_11
   /**
    *  @param[in] a_unmapped_coeff The un-mapped version of the coefficient
    *  @param[in] a_mapped_coeff The mapped version of the coefficient
    */
   virtual void updateTopLeftCoefficients( const LevelData<FluxBox>& a_unmapped_coeff,
                                           const LevelData<FluxBox>& a_mapped_coeff);

   /// Updates the bottom-left coefficients of the analytic Jacobian preconditioner, corresponding to P_21
   /**
    *  @param[in] a_unmapped_coeff The un-mapped version of the coefficient
    *  @param[in] a_mapped_coeff The mapped version of the coefficient
    */
   virtual void updateBottomLeftCoefficients( const LevelData<FluxBox>& a_unmapped_coeff,
                                              const LevelData<FluxBox>& a_mapped_coeff);


   /// Updates the bottom-right coefficients of the analytic Jacobian preconditioner, corresponding to P_22
   /**
    *  @param[in] a_unmapped_coeff The un-mapped version of the coefficient
    *  @param[in] a_mapped_coeff The mapped version of the coefficient
    */
   virtual void updateBottomRightCoefficients( const LevelData<FluxBox>& a_unmapped_coeff,
                                               const LevelData<FluxBox>& a_mapped_coeff);

   /// Updates the preconditioner matrix once the coefficients have been individually updated
   /**
    *  @param[in] a_diffusion_bc Boundary conditions for diffusion coefficients
    *  @param[in] a_viscosity_bc Boundary conditions for viscosity coefficients
    */
   virtual void updatePreconditioner( const EllipticOpBC& a_diffusion_bc,
                                      const EllipticOpBC& a_viscosity_bc );

   /// Solves the preconditioner. A null-op, use the other solvePreconditioner function
   /**
    */
   virtual void solvePreconditioner( const LevelData<FArrayBox>& in, LevelData<FArrayBox>& out );

   /// Solves the preconditioner
   /**
    *  Solves the preconditioner matrix for P*z=RHS.
    *  P is established by setting the coefficients and source term components,
    *  then calling updatePreconditioner.
    *  Since this is a 2x2 matrix, two values of RHS are required, and two values
    *  of output are returned:
    *  The density RHS and result correspond to the top row of the matrix
    *  The parallel-momentum RHS and result correspond to the top row of the matrix
    *
    *  @param[in] a_rhs_density The RHS for the density equation
    *  @param[in] a_rhs_par_mom The RHS for the parallel-momentum equation
    *  @param[in] a_z_density The result from the density equation
    *  @param[in] a_z_par_mom The result from the parallel-momentum equation
    */
   virtual void solvePreconditioner( const LevelData<FArrayBox>& a_rhs_density,
				     const LevelData<FArrayBox>& a_rhs_par_mom,
				     LevelData<FArrayBox>&       a_z_density,
				     LevelData<FArrayBox>&       a_z_par_mom );


   /// Multiplies stored coeffiicents by the input data. A null-op for this preconditioner.
   /**
    */
   virtual void multiplyCoefficients( LevelData<FluxBox>& a_data,
                                      const bool          a_mapped_coeff ) const;

   /// Sets convergence criteria for the preconditioner
   /**
    */
   virtual void setImExPreconditionerConvergenceParams( const double  tol,
                                                        const int     max_iter,
                                                        const double  precond_tol,
                                                        const int     precond_max_iter );

private:

   bool m_diffusion_coeffs_defined;
   bool m_viscosity_coeffs_defined;
   bool m_bottomLeft_coeffs_defined;

   MBHypreSolver* m_imex_preconditioner;

   LevelData<FArrayBox> m_volume_reciprocal;

   bool m_precond_verbose;
   bool m_precond_precond_verbose;

  // These are the data structures for the preconditioner block-matrix.
  // Each block is a 2x2 matrix, there is no 'topRight' here because
  // those coefficients are all zero (i.e. P_{01} = 0).

  // The 'beta' are the source term components, cell values
   LevelData<FArrayBox> m_topLeft_beta;
   LevelData<FArrayBox> m_bottomLeft_beta;
   LevelData<FArrayBox> m_bottomRight_beta;

  // The 'tensor' are the diffusion and advection components, face values
   LevelData<FluxBox> m_topLeft_mapped_tensor;
   LevelData<FluxBox> m_bottomLeft_mapped_tensor;
   LevelData<FluxBox> m_bottomRight_mapped_tensor;

   LevelData<FluxBox> m_topLeft_unmapped_tensor;
   LevelData<FluxBox> m_bottomLeft_unmapped_tensor;
   LevelData<FluxBox> m_bottomRight_unmapped_tensor;

   MagFluxAlignedMBLEXCenter* m_mblx_ptr;
   

   int m_verbosity;
};



#include "NamespaceFooter.H"

#endif // _TWOFIELDNEUTRALSPC_H_
