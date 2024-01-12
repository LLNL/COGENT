#ifndef  _IDEALMHDOP_H_
#define  _IDEALMHDOP_H_

#include "FluidOpInterface.H"
#include "FluidVarBC.H"
#include "Diffusion.H"
#include "Poisson.H"
#include "SpaceUtils.H.multidim"

#include "EdgeDataBox.H"

#include "NamespaceHeader.H"

/**
 * multi-dim ideal MHD equation operator class.
 *
 * d(rho)/dt + div(mom) = 0
 * d(mom)/dt + div(mom*u + IP + pi) = JxB
 * d(ene)/dt + div((ene+P)*u + u*pi) = JdotE
 * dB/dt     + cur(E) = 0
 *
 * rho = M*n, M is mass n is density
 * mom_x = rho*ux
 * mom_y = rho*uy
 * mom_z = rho*uz
 * ene = 0.5*rho*u^2 + P/(g-1)
 * g = 1+2/degFree is adiabatic coefficient
 *
 * mu0*J = curl(B)
 * E = uxB
 *
 * pi is viscosity tensor
 * pi = -eta*W (W is rate of strain tensor)
 * W = grad(u) + grad(u)^T - 2/3I*div(u)
 *
 * using Braginskii convention for tensors
 * divpi_j = d(pi_jk)/dx_k
 * u*pi_k = u_j*pi_jk
 *
*/
class IdealMhdOp
   : public FluidOpInterface
{
public:

  /// Null constructor (used for state variable creation) 
  IdealMhdOp( const MagGeom&  geometry )
     : m_geometry(geometry) {};

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  IdealMhdOp( const string&   pp_str,
             const string&   species_name,
             const MagGeom&  geometry,
             const int       verbosity );

  /// Destructor.
  /**
   */
  virtual ~IdealMhdOp();

   /// Evaluates the RHS.
   /**
    *  Pure virtual function that will be resposible for evaluating
    *  field RHS on the given field.
    *
    *  @param[in,out] rhs           -  contribution to be added to rhs.
    *  @param[in] fields            -  current solution for fields.
    *  @param[in] fluids            -  current solution for fluids.
    *  @param[in] kinetic_specties  -  current solution for kinetic species.
    *  @param[in] component         -  component of the field vector being computed
    *  @param[in] time              -  the time at which the field RHS is to be evaluated
    */
   virtual void accumulateExplicitRHS(FluidSpeciesPtrVect&               rhs,
                                      const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                                      const FluidSpeciesPtrVect&         fluid_species,
                                      const PS::ScalarPtrVect&           scalars,
                                      const EMFields&                    EM_fields,
                                      const int                          component,
                                      const Real                         time);
   
   virtual void accumulateImplicitRHS(FluidSpeciesPtrVect&               rhs,
                                      const PS::KineticSpeciesPtrVect&   kinetic_species,
                                      const FluidSpeciesPtrVect&         fluid_species,
                                      const PS::ScalarPtrVect&           scalars,
                                      const EMFields&                    EM_fields,
                                      const int                          component,
                                      const Real                         time);
   
   virtual void defineBlockPC( std::vector<PS::Preconditioner<PS::ODEVector,PS::AppCtxt>*>&,
                               std::vector<PS::DOFList>&,
                               const PS::ODEVector&,
                               void*,
                               const std::string&,
                               const std::string&,
                               bool,
                               const FluidSpecies&,
                               const PS::GlobalDOFFluidSpecies&,
                               const int,
                               const int );

   virtual void updateBlockPC(  std::vector<PS::Preconditioner<PS::ODEVector,PS::AppCtxt>*>&,
                                const PS::KineticSpeciesPtrVect&,
                                const FluidSpeciesPtrVect&,
                                const Real,
                                const int,
                                const int,
                                const Real,
                                const bool,
                                const int );

   virtual void updatePCImEx( const FluidSpeciesPtrVect&        fluid_species,
                              const PS::KineticSpeciesPtrVect&  kinetic_species,
                              const double                      time,
                              const int                         step,
                              const int                         stage,
                              const double                      shift,
                              const int                         component,
                              const std::string& );

   virtual void solvePCImEx( FluidSpeciesPtrVect&              fluid_species_solution,
                             const PS::KineticSpeciesPtrVect&  kinetic_species_rhs,
                             const FluidSpeciesPtrVect&        fluid_species_rhs,
                             const std::string&,
                             const int                         component );

   virtual bool trivialSolutionOp() const {return true;}
   
   virtual void getMemberVar( LevelData<FArrayBox>&  a_var,
                              const CFGVars&         a_fluid_vars,
                              const string&          a_name ) const;

   virtual CFGVars* createStateVariable( const string&       a_pp_prefix,
                                         const std::string&  a_name,
                                         const MagGeom&      a_geometry,
                                         const IntVect&      a_ghost_vect ) const
   {
      return new FluidSpecies(a_pp_prefix, a_name, a_geometry, a_ghost_vect);
   }
   
   virtual void initialize( CFGVars&      a_species,
                            const double  a_time );

   virtual bool isInitializationConstrained(const FluidSpecies& fluid_phys);
   
   virtual void applyInitializationConstraints(FluidSpeciesPtrVect&               fluid_comp,
                                               FluidSpeciesPtrVect&               fluid_phys,
                                               const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                                               const EMFields&                    EM_fields,
                                               const int                          component,
                                               const double                       time );

   virtual void fillGhostCells( FluidSpecies&  species_phys,
                                const double   time );

   /// Compute a stable time step.
   /**
    * Computes and returns an estimate of the maximum stable time step.
    */
   Real computeDtExplicitTI( const FluidSpeciesPtrVect&  a_fluid_comp );
   Real computeDtImExTI( const FluidSpeciesPtrVect&  a_fluid_comp );
   
   Real computeDtMhd( const FluidSpeciesPtrVect&  a_fluid_comp );
   
   
   static const char* pp_name;

private:

   int m_verbosity;
   double m_courant_time_step;
   const MagGeom& m_geometry;
   int m_implicitMhd, m_implicitViscosity;
   std::string m_advScheme, m_advScheme_B;
   std::string m_species_name;
   int m_CGL;
   int m_magP;
   bool m_blankN, m_blankB;
   int m_thermPforE0, m_magPforE0;
   int m_initializeBfromVectorPotential;
   int m_initializeProfilesFromPSI;
   int m_initializeBackgroundPlusPert; 
   int m_addEnergyPerturbation;
   int m_perturbDensity;
   int m_perturbBtor;
   Real m_rbtor_fact;
   Real m_bpol_fact;

   std::string m_opt_string;
   int m_my_pc_idx;
      
   // allocate contains for metric stuff to transform from phys to mapped
   // (eventually these will be removed as function are moved to magGeom.cp)
   LevelData<FArrayBox> m_Jacobian;
   LevelData<FArrayBox> m_dXdxi;
   LevelData<FArrayBox> m_dxidX;
   LevelData<FArrayBox> m_Nmatrix;
   LevelData<FArrayBox> m_Xphys;
   LevelData<FArrayBox> m_g_y; // used only when solving virtual B
   bool m_twoDaxisymm;
   bool m_energy_is_thermal;
   bool m_energy_is_thermal_2;
   bool m_isothermal;
   bool m_enforceFluxBCs;
   bool m_use_perturbative_form;
   bool m_free_stream_preserving;
   
   // allocate containers for species fluid variables
   //
   LevelData<FArrayBox> m_Jrho_PC;
   LevelData<FArrayBox> m_shift_PC;
   double m_alpha_PC;
   LevelData<FluxBox> m_ellip_coeff, m_ellip_coeff_mapped;
   
   LevelData<FArrayBox> m_rhoDen_cc;
   LevelData<FArrayBox> m_momDen_cc;
   LevelData<FArrayBox> m_eneDen_cc;
   LevelData<FArrayBox> m_momDen_virtual_cc;
   LevelData<FArrayBox> m_By_phys_cc;
   //
   LevelData<FArrayBox>  m_By_contr_cc;
   LevelData<FArrayBox>  m_By_covar_cc;
   LevelData<FArrayBox>  m_B_covar_cc;
   LevelData<FArrayBox>  m_B0_phys, m_B1_phys;

   // allocate containers for derived variables
   //
   LevelData<FArrayBox> m_pressure, m_P0, m_N0, m_deltaN;
   LevelData<FArrayBox> m_temperature, m_T0, m_deltaT;
   LevelData<FArrayBox> m_Bpressure;
   LevelData<FArrayBox> m_velocity;
   LevelData<FluxBox> m_velocity_cf;
   LevelData<FArrayBox> m_velocity_virtual;
   
   LevelData<FArrayBox> m_velocity_contra, m_momDen_contra;
   
   // allocate containers for viscosity tensor
   //
   double m_relative_density_floor;
   double m_tau;
   double m_etaMin;
   LevelData<FArrayBox> m_etaVisc;
   LevelData<FArrayBox> m_momVisc_source;
   LevelData<FluxBox> m_etaVisc_cf;
   LevelData<FluxBox> m_JaW_cf;
   LevelData<FluxBox> m_m0JaFluxVisc_cf;
   LevelData<FluxBox> m_m1JaFluxVisc_cf;
   LevelData<FluxBox> m_m2JaFluxVisc_cf;
   LevelData<FluxBox> m_enJaFluxVisc_cf;
   LevelData<FluxBox> m_m0JaFluxVisc_norm;
   LevelData<FluxBox> m_m1JaFluxVisc_norm;
   LevelData<FluxBox> m_m2JaFluxVisc_norm;
   LevelData<FluxBox> m_momJaFluxVisc_norm;
   LevelData<FluxBox> m_enJaFluxVisc_norm;

   // allocate containers for fluid variable fluxes at cell center
   //
   LevelData<FArrayBox> m_rhoFlux_cc;
   LevelData<FArrayBox> m_mxFlux_cc;
   LevelData<FArrayBox> m_myFlux_cc;
   LevelData<FArrayBox> m_mzFlux_cc;
   LevelData<FArrayBox> m_enFlux_cc;
   LevelData<FArrayBox> m_enParFlux_cc;
   LevelData<FArrayBox> m_mvFlux_cc;
   
   // allocate containers for Cspeed
   //
   LevelData<FArrayBox>  m_Cspeed_cc;
   LevelData<FluxBox> m_CspeedL_norm;
   LevelData<FluxBox> m_CspeedR_norm;
   
   // allocate containers for normal fluxes on faces
   //
   LevelData<FluxBox> m_rhoFlux_norm;
   LevelData<FluxBox> m_mxFlux_norm;
   LevelData<FluxBox> m_myFlux_norm;
   LevelData<FluxBox> m_mzFlux_norm;
   LevelData<FluxBox> m_enFlux_norm;
   LevelData<FluxBox> m_enParFlux_norm;
   LevelData<FluxBox> m_momFlux_norm;
   LevelData<FluxBox> m_mvFlux_norm;
   LevelData<FArrayBox> m_div_IdentFlux_R;
   LevelData<FArrayBox> m_div_IdentFlux_Z;
   
   // allocate containers for magnetic field stuff
   LevelData<EdgeDataBox>  m_E0_ce;
   LevelData<FArrayBox>  m_JaJ0_cc, m_JaJ00, m_JaJ01;
   LevelData<FArrayBox>  m_Apar_cc;
   LevelData<FArrayBox>  m_E0_cc;
   LevelData<FArrayBox>  m_curlE0_virtual;
   LevelData<FluxBox>    m_curlE0_cf;
   LevelData<FluxBox>    m_B_contr_cf;
   LevelData<FArrayBox>  m_JaUdotGradP, m_JaPdivU;
   LevelData<FArrayBox>  m_JaJdotE;   // Jacobian*J\cdotE
   LevelData<FArrayBox>  m_JaJcrossB; // Jacobian*J\timesB
   LevelData<FArrayBox>  m_JaJcrossB_10; // Jacobian*J\timesB
   //
   LevelData<FArrayBox> m_JaJ0_virtual_cc;
   LevelData<FArrayBox> m_JaJcrossB_virtual;
   LevelData<FArrayBox> m_E0_virtual_cc;
   LevelData<NodeFArrayBox> m_E0_virtual_nc;
   
   // allocate containers for stuff used to enforce flux BCs
   LevelData<NodeFArrayBox> m_NodeBC_zeros;
   LevelData<EdgeDataBox> m_EdgeBC_zeros;
   LevelData<FluxBox> m_FluxBC_zeros;
   LevelData<FluxBox> m_rhoFluxBC_norm;
   LevelData<FluxBox> m_mxFluxBC_norm;
   LevelData<FluxBox> m_myFluxBC_norm;
   LevelData<FluxBox> m_mzFluxBC_norm;
   LevelData<FluxBox> m_enFluxBC_norm;
   LevelData<FluxBox> m_enParFluxBC_norm;
   LevelData<FluxBox> m_mvFluxBC_norm;
   LevelData<FluxBox> m_momFluxBC_norm;
   //
   
   mutable LevelData<FArrayBox> m_dummyDiv;
   mutable LevelData<FArrayBox> m_dummyFArray_oneComp;
   mutable LevelData<FArrayBox> m_dummyFArray_spaceDim;
   LevelData<FArrayBox> m_dummyDiv_mom, m_forceDensity, m_rhs_rho;
   
   mutable LevelData<FluxBox> m_dummyFlux_oneComp;
   mutable LevelData<FluxBox> m_dummyFlux_spaceDim;
   mutable LevelData<EdgeDataBox> m_dummyEdge_oneComp;
   
   mutable LevelData<NodeFArrayBox> m_dummyNode_oneComp;

   LevelData<FArrayBox> m_PSI;

   /// Parse parameters.
   /**
    * Private method to obtain control parameters from "field_model.field_name" section
    * of ParmParse database, where "field_name" is replaced by the field name (e.g., field_model.VectPotential)
    *
    * \param[in] ppcls field_model.field input database.
    */
   void parseParameters( ParmParse& pp );
   

   /// Print parameters.
   /**
    * Private method to print parameters.
    */
   void printParameters();
   
   Diffusion*  m_diffusionOp_mom0;
   Diffusion*  m_diffusionOp_mom1; 
   Diffusion*  m_diffusionOp_mom2; 
   Poisson*  m_diffusionOp_Apar; 
   EllipticOpBC*  m_bc_mom0;
   EllipticOpBC*  m_bc_mom1;
   EllipticOpBC*  m_bc_mom2;
   EllipticOpBC*  m_bc_Apar;
   std::vector<FluidVarBC*>  m_fluid_bc;
   std::vector<FluidVarBC*>  m_fluid_face_var_bc;
   
   int m_density_cbc;
   int m_momentumDensity_cbc;
   int m_energyDensity_cbc;
   int m_momentumDensity_virtual_cbc;
   int m_magneticField_virtual_cbc;
   int m_magneticField_fbc;
  
   void accumulateMhdRHS( FluidSpecies&  a_rhs_fluid,
                    const FluidSpecies&  a_soln_fluid,
                    const Real           a_time );
   
   void accumulateViscRHS( FluidSpecies&  a_rhs_fluid,
                     const FluidSpecies&  a_soln_fluid,
                     const int            a_accumulateMom,
                     const int            a_accumulateEne );
   
   void computeDiffusionCoefficients( LevelData<FluxBox>&  a_ellip_coeff,
                                      LevelData<FluxBox>&  a_ellip_coeff_mapped,
                                const LevelData<FluxBox>&  a_D_fc,
                                const Real                 a_time) const;
   
   void inspectEdgeDataBox(const LevelData<EdgeDataBox>& a_Edge,
                           const int                 a_dir) const;
   
   void convertPhysFluxToMappedFlux( FArrayBox& a_mapped_flux,
                               const FArrayBox& a_physical_flux,
                               const FArrayBox& a_N) const;

   void computeNTFfaceArea( LevelData<FluxBox>&    a_Flux_norm,
                      const LevelData<FArrayBox>&  a_Flux_mapped_cc,
                      const LevelData<FArrayBox>&  a_Jf_cc,
                      const LevelData<FArrayBox>&  a_Cspeed_cc,
                      const LevelData<FluxBox>&    a_CspeedR_norm,
                      const LevelData<FluxBox>&    a_CspeedL_norm ) const;
   
   void upwindFluxVector( LevelData<FluxBox>&    a_Flux_phys_fc,
                    const LevelData<FArrayBox>&  a_Flux_phys_cc,
                    const LevelData<FArrayBox>&  a_fun_cc,
                    const LevelData<FArrayBox>&  a_Cspeed_cc,
                    const LevelData<FluxBox>&    a_CspeedR_norm,
                    const LevelData<FluxBox>&    a_CspeedL_norm ) const;
   
   void computeIdealEatEdges( LevelData<EdgeDataBox>&  a_E_covar,
                        const LevelData<FArrayBox>&    a_V_phys_cc,
                        const LevelData<FArrayBox>&    a_B_contr_cc,
                        const LevelData<FArrayBox>&    a_Cspeed_cc ) const;
   
   void computeIdealEatEdges( LevelData<EdgeDataBox>&  a_E_covar,
                        const LevelData<FArrayBox>&    a_V_phys_cc,
                        const LevelData<FluxBox>&      a_B_contr,
                        const LevelData<FArrayBox>&    a_Cspeed_cc ) const;
   
   void computeIdealEatNodes( LevelData<NodeFArrayBox>&  a_E_covar,
                        const LevelData<FArrayBox>&      a_V_phys_cc,
                        const LevelData<FluxBox>&        a_B_contr,
                        const LevelData<FArrayBox>&      a_Cspeed_cc ) const;
   
   void setCourantTimeStep(const LevelData<FArrayBox>& a_Cspeed);

   void defineLevelDatas( const DisjointBoxLayout&  a_grids,
                          const IntVect&            a_ghostVect );
   
   void setCellCenterValues( const FluidSpecies&  a_soln_fluid );
   
   void setCellCenterFluxes( const FluidSpecies&  a_soln_fluid );
   
   void setCellCenterFluxesCGL( const FluidSpecies&  a_soln_fluid );

   void setMagneticFieldTerms( const FluidSpecies&  a_soln_fluid );
   
   void setMagneticFieldTerms_2Dvirtual( const FluidSpecies&  a_soln_fluid );
   
   void setMagneticFieldTerms_2DinPlane( const FluidSpecies&  a_soln_fluid );
   
   void setMagneticFieldTerms_2Dfull( const FluidSpecies&  a_soln_fluid );
   
   void setMagneticFieldTerms_3D( const FluidSpecies&  a_soln_fluid );

   void setFaceCenteredFluxes( const FluidSpecies&  a_soln_fluid );
   
   void setFaceCenteredFluxes_freestreamPreserving( const FluidSpecies&  a_soln_fluid );
   
   void setFaceCenteredFluxes_ZIP( const FluidSpecies&  a_soln_fluid );

   void setFaceCenteredFlux( LevelData<FluxBox>&    a_Flux_norm,
                       const LevelData<FArrayBox>&  a_Flux_phys_cc );
   
   void enforceFluxBCs( const FluidSpecies&  a_soln_fluid,
                        const Real           a_time );
  
   void setCurlE( const FluidSpecies&  a_soln_fluid,
                  const Real           a_time );
   
   void setViscosityCoefficient( LevelData<FArrayBox>&  a_etaVisc,
                                 LevelData<FluxBox>&    a_etaVisc_cf,
                           const FluidSpecies&          a_soln_fluid ); 
   
   void computeViscosityEnergyFlux( LevelData<FluxBox>&  a_UJaPi_Flux_phys,
                              const LevelData<FluxBox>&  a_velocity_phys );
   
   void computeViscSourceAxisymm( LevelData<FArrayBox>&  a_Wthth,
                            const LevelData<FArrayBox>&  a_velocity_phys,
                            const LevelData<FluxBox>&    a_velocity_phys_cf );
   
   void updateRHSs_visc( FluidSpecies&  a_rhs_fluid,
                   const FluidSpecies&  a_soln_fluid,
                   const int            a_accumulateMom,
                   const int            a_accumulateEne );

   void updateRHSs( FluidSpecies&  a_rhs_fluid,
              const FluidSpecies&  a_soln_fluid );
   
   void setMappedCspeed( const double  a_gamma,
                         const int     a_Ptherm,
                         const int     a_Pmag );

   void computePhysicalDivergence( LevelData<FArrayBox>&  a_divF_phys,
                             const LevelData<FluxBox>&    a_F_mapped,
                             const int                    a_mapped ) const;
   
   void getE0_2Dvirtual( LevelData<FArrayBox>&  a_E0, 
                   const FluidSpecies&          a_fluid_species ) const;
   
   void getE0_2DinPlane( LevelData<FArrayBox>&  a_E0, 
                   const FluidSpecies&          a_fluid_species ) const;
   
   void getE0_2Dfull( LevelData<FArrayBox>&  a_E0, 
                const FluidSpecies&          a_fluid_species ) const;
   
   void getE0_3D( LevelData<FArrayBox>&  a_E0, 
            const FluidSpecies&          a_fluid_species ) const;
   
   void getJ0_2Dvirtual( LevelData<FArrayBox>&  a_J0, 
                   const FluidSpecies&          a_fluid_species ) const;
   
   void getJ0_2DinPlane( LevelData<FArrayBox>&  a_J0, 
                   const FluidSpecies&          a_fluid_species ) const;
   
   void getJ0_2Dfull( LevelData<FArrayBox>&  a_J0, 
                const FluidSpecies&          a_fluid_species ) const;
   
   void getJ0_3D( LevelData<FArrayBox>&  a_J0, 
            const FluidSpecies&          a_fluid_species ) const;
   
   void getJ0_3D2( LevelData<FArrayBox>&  a_J0, 
            const FluidSpecies&          a_fluid_species ) const;
   
   void computeCrossProduct( LevelData<FArrayBox>&  a_A,
                       const LevelData<FArrayBox>&  a_B, 
                       const LevelData<FArrayBox>&  a_C,
                       const int                    a_virtual_case ) const;
   
   void computeCrossProduct( LevelData<FArrayBox>&  a_A_inPlane,
                       const LevelData<FArrayBox>&  a_B_inPlane, 
                       const LevelData<FArrayBox>&  a_C_inPlane,
                       const LevelData<FArrayBox>&  a_B_virtual, 
                       const LevelData<FArrayBox>&  a_C_virtual ) const;
   
};


#include "NamespaceFooter.H"

#endif
