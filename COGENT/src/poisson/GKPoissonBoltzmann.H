#ifndef _GKPOISSONBOLTZMANN_H_
#define _GKPOISSONBOLTZMANN_H_

#include "GKPoisson.H"
#include "FluxSurface.H"
#include "BoltzmannElectron.H"

#include "GKBoltzmannRelation.H"

#include "MBTridiagonalSolver.H"

#include "NamespaceHeader.H"

/*
  The Boltzmann prefactor strategy selects from among the following options for
  setting the prefactor n0 in the Boltzmann relation

      n_e = n0 * exp(phi/T_e)

   GLOBAL_NEUTRALITY:  At every time step, the averages of ne and ni over the
                       entire problem domain are made equal.
   GLOBAL_NEUTRALITY_INITIAL:  The averages of ne and ni over the entire problem
                       are made equal only at the initial time t=0.  The same
                       prefactor is used for all subsequent times.
   FS_NEUTRALITY:      At every time step, the flux surface averages of ne and ni
                       are made equal.
   FS_NEUTRALITY_GLOBAL_NI: At every time step, the numerator of the
                       prefactor is the average of ni over the entire domain
                       at the current time, and the denomininator is the
                       current flux surface average of exp(phi/Te).
   FS_NEUTRALITY_INITIAL_GLOBAL_NI: At every time step, the numerator of the
                       prefactor is the average of ni over the entire domain
                       at the intial time t=0, and the denomininator is the
                       current flux surface average of exp(phi/Te).
   FS_NEUTRALITY_INITIAL_FS_NI: At every time step, the numerator of the
                       prefactor is the flux surface average of ni at the intial
                       time t=0, and the denomininator is the current flux surface
                       average of exp(phi/Te).
*/

enum prefactor_strategy {GLOBAL_NEUTRALITY,
                         GLOBAL_NEUTRALITY_INITIAL,
                         FS_NEUTRALITY,
                         FS_NEUTRALITY_GLOBAL_NI,
                         FS_NEUTRALITY_INITIAL_GLOBAL_NI,
                         FS_NEUTRALITY_INITIAL_FS_NI};

class GKPoissonBoltzmann
   : public GKPoisson
{
public:


   GKPoissonBoltzmann( ParmParse&                  pp,
                       const MagGeom&              geom,
                       const Real                  larmor_number,
                       const Real                  debye_number,
                       const LevelData<FArrayBox>& initial_ion_charge_density );

   ~GKPoissonBoltzmann( );

   void computePotentialAndElectronDensity(
      LevelData<FArrayBox>&        phi,
      BoltzmannElectron&           ne,
      const LevelData<FArrayBox>&  ni,
      const EllipticOpBC&          bc,
      const bool                   first_step = false);

   virtual void setOperatorCoefficients( const LevelData<FArrayBox>&  ni,
                                         const EllipticOpBC&          bc ); 

   virtual void setOperatorCoefficients( const LevelData<FArrayBox>&  ion_mass_density,
                                         const EllipticOpBC&          bc,
                                         double&                      lo_value,
                                         double&                      hi_value,
                                         LevelData<FArrayBox>&        radial_gkp_divergence_average);

   static const char* pp_name;

   ///
   /**
      Given the current state of the residual the correction, apply your preconditioner to a_cor.
   */
   virtual void preCond( LevelData<FArrayBox>&       cor,
                         const LevelData<FArrayBox>& residual );

   ///
   /**
      In the context of solving L(phi) = rhs, set a_lhs = L(a_phi).  If a_homogeneous is true,
      evaluate the operator using homogeneous boundary conditions.
   */
   virtual void applyOp( LevelData<FArrayBox>& a_out,
                         const LevelData<FArrayBox>& a_in,
                         bool a_homogeneous = false);

   /// Compute the electron density ne = prefactor exp(phi/T)
   /** Compute the electron density ne = prefactor exp(phi/T).  The prefactor depends
    *  upon whether global or flux surface neutrality has been requested
    *
    * @param[in]   phi       Potential
    * @param[in]   ni        Ion density
    * @param[out]  ne        Electron density
    */
   void computeElectronDensity( const LevelData<FArrayBox>& phi,
                                const LevelData<FArrayBox>& ni,
                                BoltzmannElectron& ne );

   void computeBoundaryData( FArrayBox& a_inner_divertor_bvs,
                             FArrayBox& a_outer_divertor_bvs,
                             const BoltzmannElectron&      a_ne,
                             const LevelData<FArrayBox>& a_Zni,
                             const LevelData<FArrayBox>& a_Jpar);

   void computeSheathEllipticOpBC( const BoltzmannElectron&  ne,
                                   const double *            jpar_outer_plate, 
                                   const double *            jpar_inner_plate,
                                   const double *            bmag_outer_plate, 
                                   const double *            bmag_inner_plate,
                                   const int                 norm_dir );

   void getPhiTilde( const LevelData<FArrayBox>& Zni,
                     const BoltzmannElectron&    ne,
                     LevelData<FArrayBox>&       phi_tilde );

protected:

   /// Update the Boltzmann derivative contribution to the Jacobian system
   /** Update the contribution to the Jacobian system due to the derivative
    *  of the Boltzmann relation
    *
    * @param[in]  ne  Boltzmann electron species
    */
   void updateLinearSystem( const BoltzmannElectron& ne, const EllipticOpBC& bc );

   /// Solve the GK Poisson system
   /** Solve the GK Poisson system
    *
    * @param[in]      ni  Ion density
    * @param[in,out]  ne  Boltzmann electron species
    * @param[out]     phi Potential
    */
   void solve( const LevelData<FArrayBox>&  ni,
               const EllipticOpBC&          bc,
               BoltzmannElectron&           ne,
               LevelData<FArrayBox>&        phi );

   void solveRadial( const LevelData<FArrayBox>&  ni,
                     const EllipticOpBC&          bc,
                     BoltzmannElectron&           ne,
                     LevelData<FArrayBox>&        phi );

   void solveLinear( const LevelData<FArrayBox>&  ni,
                     const EllipticOpBC&          bc,
                     BoltzmannElectron&           ne,
                     LevelData<FArrayBox>&        phi );

   void solveSubspaceIteration( const LevelData<FArrayBox>&  ni,
                                const EllipticOpBC&          bc,
                                BoltzmannElectron&           ne,
                                LevelData<FArrayBox>&        phi );

   void solveSimpleBoltzmann( const LevelData<FArrayBox>&   ni,
                                BoltzmannElectron&          ne,
                                LevelData<FArrayBox>&       phi );

   /// Solve the Jacobian system Jz = r
   /** Solve the Jacobian system Jz = r
    *
    * @param[in]   rvec  Right-hand side vector
    * @param[in]   use_relative_tolerance  Solve to an absolute tolerance
    * @param[out]  zvec  Solution vector
    */
   void solveJacobian( const LevelData<FArrayBox>& rvec,
                       LevelData<FArrayBox>& zvec,
                       const bool use_absolute_tolerance = false );

   /// Compute residual = GKP(phi) - ni + ne
   /** Compute residual = GKP(phi) - ni + ne
    *
    * @param[in]   phi       Potential
    * @param[in]   ni        Ion density
    * @param[in]   ne        Electron density
    * @param[out]  residual  Residual
    */
   double computeResidual( const LevelData<FArrayBox>& ni,
                           const LevelData<FArrayBox>& phi,
                           BoltzmannElectron&          ne,
                           LevelData<FArrayBox>&       residual );

   void computeChargeDensity( const LevelData<FArrayBox>& Zni,
                              const LevelData<FArrayBox>& phi,
                              BoltzmannElectron&          ne,
                              LevelData<FArrayBox>&       charge_density );

   /// Compute the derivative of the electron density with respect to the potential phi
   /** Compute the derivative of the electron density with respect to the potential phi.
    *  Here, we ignore the possible dependence of the Boltzmann relation prefactor on
    *  phi (the contribution to the Jacobian made by those derivative is accounted for
    *  elsewhere), so the derivative returned by this function is just ne/Te
    *
    * @param[in]  ne           Boltzmann electron species
    * @param[out] boltz_deriv  Derivative of ne with respect to the potential phi
    */
   void computeBoltzmannDerivative( const BoltzmannElectron& ne,
                                    LevelData<FArrayBox>& boltz_deriv ) const;

   /// Compute the maximum neutrality error relative to the ion density on all flux surfaces.
   /** Compute the maximum neutrality error relative to the ion density on all flux surfaces.
    *  On each flux surface, we compute the total charge divided by the ion charge.  We then return
    *  the maximum absolute value over all of the flux surfaces.
    *
    *
    * @param[in]  ni           Ion density
    * @param[in]  ne           Boltzmann electron species
    * @param[out] error        L2Norm(ni-ne) / L2Norm(ni)
    */
   double fluxSurfaceNeutralityRelativeError( const LevelData<FArrayBox>& ni,
                                              const LevelData<FArrayBox>& ne ) const;

   /// Compute the global charge neutrality error relative to the total ion charge.
   /** Compute the global charge neutrality error relative to the total ion charge.  Specifically,
    *  compute the total charge in the domain divided by the total ion charge in the domain
    *
    * @param[in]  ni           Ion density
    * @param[in]  ne           Boltzmann electron species
    * @param[out] retval       Neutrality error
    */
   double globalNeutralityRelativeError( const LevelData<FArrayBox>& ni,
                                         const LevelData<FArrayBox>& ne ) const;

   /// Integrate over the domain
   /** Integrate over the domain
    *
    * @param[in]  data         Data to be integrated
    * @param[out] retval       Integral
    */
   double integrateMapped( const LevelData<FArrayBox>& data ) const;

   /// Average over the domain
   /** Average over the domain
    *
    * @param[in]  data         Data to be averaged
    * @param[out] retval       Average
    */
   double averageMapped( const LevelData<FArrayBox>& data ) const;

   void computePrefactorNumerator( const LevelData<FArrayBox>& ion_density );

   void applyP( LevelData<FArrayBox>& a_data ) const;

   void applyQ( LevelData<FArrayBox>& a_data ) const;

   void multiplyM( LevelData<FArrayBox>& a_data ) const;

   void divideM( LevelData<FArrayBox>& a_data ) const;

   void multiplyD( LevelData<FArrayBox>& a_data ) const;

   int getConsequtiveRadialIndex( const int a_mapped_index,
                                  const int a_block_number ) const;
   
   void initializeElectronTemperature(LevelData<FArrayBox>&    Te,
                                      const BoltzmannElectron& ne);

   LevelData<FArrayBox> m_M;
   LevelData<FArrayBox> m_boltzmann_prefactor;
   LevelData<FArrayBox> m_boltzmann_prefactor_saved_numerator;
   LevelData<FArrayBox> m_initial_ion_charge_density;
   LevelData<FArrayBox> m_D;
   LevelData<FArrayBox> m_Te;

   double m_nonlinear_relative_tolerance;
   double m_nonlinear_change_tolerance;
   double m_jacobian_solve_tolerance;

   int m_nonlinear_max_iterations;

   bool m_recompute_prefactor;

   FluxSurface m_flux_surface;

   GKBoltzmannRelation m_boltzmann_relation;

   prefactor_strategy m_prefactor_strategy;

   double m_initial_ion_charge_density_average;

   LevelData<FArrayBox> m_ni_fs_av_saved;

   bool m_radial_solve_only;
   bool m_subspace_iteration_solve;
   bool m_simple_boltzmann;
   bool m_linear_response;
   bool m_linear_solve;
   bool m_zero_initial_solution;

   int m_max_subspace_iterations;
   double m_subspace_iteration_tol;

   bool m_compute_core_prefactor;
   LevelData<FArrayBox> m_core_prefactor_numerator;

   MBTridiagonalSolver* m_precond_Psolver;
#ifdef with_petsc
   MBPETScSolver* m_precond_Qsolver;
#else
   MBHypreSolver* m_precond_Qsolver;
#endif

   double*  m_Zni_outer_plate;
   double*  m_Zni_inner_plate;
   double*  m_phi_outer_plate;
   double*  m_phi_inner_plate;

   // If true, preserves the initial ion density average
   // in the Newton residual calculation to enable solvability
   // for problems with a (constant) null space
   bool m_preserve_initial_ni_average;

   bool m_gkp_verbose;
};

#include "NamespaceFooter.H"

#endif   /* _GKPOISSONBOLTZMANN_H_ */
