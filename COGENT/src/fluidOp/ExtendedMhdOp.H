#ifndef  _EXTENDEDMHDOP_H_
#define  _EXTENDEDMHDOP_H_

#include "FluidOpInterface.H"
#include "FluidVarBC.H"
#include "Diffusion.H"
#include "SpaceUtils.H"
#include "ConstFact.H"

#include "EdgeDataBox.H"
#include "NodeFArrayBox.H"

#include "NamespaceHeader.H"

/**
 * multi-dim ideal MHD equation operator class.
 *
 * d(rho)/dt + div(mom) = 0
 * d(mom)/dt + div(mom*u + IP + pii) = JxB
 * d(enei)/dt + div((enei+Pi)*u  + u*pii ) = nUdotE + Qie 
 * d(enee)/dt + div((enee+Pe)*ue + ue*pie) = JdotE - nUdotE - Qie
 * dB/dt     + cur(E) = 0
 * 
 * epsilon*dJ/dt = n*[E-E0-lambdai/n*(J0xB-nabla(Pe))] - n*eta*J
 * delta*dE/dt   = J0 - J
 *
 * rho = M*n, M is mass n is density
 * mom_x = rho*ux
 * mom_y = rho*uy
 * mom_z = rho*uz
 * ene = 0.5*rho*u^2 + P/(g-1)
 * g = 1+2/degFree is adiabatic coefficient
 *
 * mu0*J0 = curl(B)
 * E0 = uxB
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
class ExtendedMhdOp
   : public FluidOpInterface
{
public:

  /// Null constructor (used for state variable creation) 
  ExtendedMhdOp( const MagGeom&  geometry )
     : m_geometry(geometry) {};

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  ExtendedMhdOp( const string&   pp_str,
                 const string&   species_name,
                 const MagGeom&  geometry,
                 const double    Nscale,
                 const double    Tscale,
                 const double    Xscale,
                 const double    Mscale,
                 const int       verbosity );

  /// Destructor.
  /**
   */
  virtual ~ExtendedMhdOp();

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
                                      const EField&                      E_field,
                                      const int                          component,
                                      const Real                         time);
   
   virtual void accumulateImplicitRHS(FluidSpeciesPtrVect&               rhs,
                                      const PS::KineticSpeciesPtrVect&   kinetic_species,
                                      const FluidSpeciesPtrVect&         fluid_species,
                                      const PS::ScalarPtrVect&           scalars,
                                      const EField&                      E_field,
                                      const int                          component,
                                      const Real                         time);
   
   virtual void defineBlockPC( std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>&,
                               std::vector<PS::DOFList>&,
                               const PS::GKVector&,
                               PS::GKOps&,
                               const std::string&,
                               const std::string&,
                               bool,
                               const FluidSpecies&,
                               const PS::GlobalDOFFluidSpecies&,
                               int );

   virtual void updateBlockPC(  std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>&,
                                const PS::KineticSpeciesPtrVect&,
                                const FluidSpeciesPtrVect&,
                                const Real,
                                const Real,
                                const bool,
                                const int );

   virtual void updatePCImEx( const FluidSpeciesPtrVect&        fluid_species,
                              const PS::KineticSpeciesPtrVect&  kinetic_species,
                              const double                      time,
                              const double                      shift,
                              const int                         component );

   virtual void solvePCImEx( FluidSpeciesPtrVect&              fluid_species_solution,
                             const PS::KineticSpeciesPtrVect&  kinetic_species_rhs,
                             const FluidSpeciesPtrVect&        fluid_species_rhs,
                             const int                         component );

   virtual void postTimeEval( FluidSpecies&  a_species_comp,
                              FluidSpecies&  a_species_phys,
                              const Real     a_dt,
                              const Real     a_time,
                              const int      a_stage );

   virtual void postStageEval( FluidSpecies&  a_species_comp,
                               FluidSpecies&  a_species_phys,
                               const Real     a_dt,
                               const Real     a_time,
                               const int      a_stage );
      
   //virtual LevelData<FArrayBox>* getVarPtr( const string&  a_var_name )
   //{ return m_deltaN_Ptr; }
   
   virtual void getMemberVar( LevelData<FArrayBox>&  a_Var,
                        const FluidSpecies&          a_fluid_species,
                        const string&                a_name ) const;

   virtual bool trivialSolutionOp() const {return true;}

   virtual CFGVars* createStateVariable( const string&       a_pp_prefix,
                                         const std::string&  a_name,
                                         const MagGeom&      a_geometry,
                                         const IntVect&      a_ghost_vect ) const
   {
      return new FluidSpecies(a_pp_prefix, a_name, a_geometry, a_ghost_vect);
   }
   
   virtual void initialize( CFGVars&      a_species,
                            const double  a_time );

   virtual void initializeWithBC( FluidSpecies&  species_comp,
                                  FluidSpecies&  species_phys,
                            const double         time );

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
   const MagGeom& m_geometry;
   
   bool m_useRelaxedHeatFlux;
   bool m_useRelaxedViscosity;
   bool m_useRelaxedOhmsLaw;
   double m_tauiMin;
   double m_taueMin;
   double m_etaResMin;
   double m_etaMin;
   double m_etaVis_ele;
   double m_viscLimiter;
   double m_gyroViscLimiter;
   double m_gyroHeatViscCoeff;
   double m_kappaMin;
   double m_heatLimiter;
   double m_gyroHeatLimiter;
   double m_xionFactor;
   double m_xeleFactor;
   
   std::string m_advScheme;
   double m_meRel;
   double m_epRel;
   bool m_useJ0forJcc;
   bool m_fullyImplicitJ;
   int m_iterHallMax;
   double m_HallResTol;
   std::string m_advSchemeHall;
   //bool m_useFiniteIonCorrections;
   bool m_modelDriftTerms;
   bool m_initializeJfromAmperesLaw;
   bool m_initializeEfromOhmsLaw;
   std::string m_species_name;
   double m_courant_time_step;
   std::string m_opt_string;
   int m_my_pc_idx;

   Copier m_exchangeCopier;
   
   double m_mM;
   double m_wpescale_Hz;
   double m_wpiscale_Hz;
   double m_wcescale_Hz;
   double m_wciscale_Hz;
   double m_Tscale_eV;
   double m_tscale_s;
   double m_tauescale1eV_s;
   double m_tauiscale1eV_s;

   double m_Ve0sq, m_Vi0sq;
   double m_lambdai0;
   double m_Le0sq;
   double m_epsilon0;
   double m_delta0;
   double m_taui0, m_taue0;

   // allocate contains for metric stuff to transform from phys to mapped
   // (eventually these will be removed as function are moved to magGeom.cp)
   LevelData<FArrayBox> m_Jacobian;
   LevelData<FluxBox> m_Jacobian_cf;
   //LevelData<FArrayBox> m_dXdxi;
   //LevelData<FArrayBox> m_dxidX;
   LevelData<FArrayBox> m_Nmatrix;
   LevelData<FArrayBox> m_Xphys;
   LevelData<FArrayBox> m_g_y; // used only when solving virtual B
   LevelData<FArrayBox> m_cellVol;
   bool m_twoDaxisymm;
   
   // allocate containers for viscosity tensor
   //
   LevelData<FArrayBox> m_etaVisC, m_etaVisG, m_JaPion;
   LevelData<FArrayBox> m_momVisc_source;
   LevelData<FluxBox> m_etaVisC_cf, m_etaVisG_cf;
   LevelData<FluxBox> m_DiffVx_cf, m_DiffVz_cf;
   LevelData<FArrayBox> m_DiffVth;
   LevelData<FluxBox> m_JaW_cf, m_JaPion_cf;
   LevelData<FluxBox> m_m0JaFluxVisc_cf;
   LevelData<FluxBox> m_m1JaFluxVisc_cf;
   LevelData<FluxBox> m_enJaFluxVisc_cf;
   LevelData<FluxBox> m_m0JaFluxVisc_norm;
   LevelData<FluxBox> m_m1JaFluxVisc_norm;
   LevelData<FluxBox> m_momJaFluxVisc_norm;
   LevelData<FluxBox> m_enJaFluxVisc_norm;
   //
   LevelData<FArrayBox> m_qi_wedge_cc;
   LevelData<FluxBox> m_JaWq_cf;   
   
   // containers for viscosity in Ohm's law
   //
   LevelData<NodeFArrayBox> m_gradVepara_cn;   
   LevelData<FArrayBox> m_gradVeperp_cc;   
   LevelData<FArrayBox> m_gradVepara_cc;   
   LevelData<NodeFArrayBox> m_divVe_cn;   
   LevelData<EdgeDataBox> m_divVe_ce;   
   LevelData<NodeFArrayBox> m_Pie_diag_cn, m_Pie_diag_cn0, m_Pie_diag_cn_old, m_Pie_diagD_cn;   
   LevelData<FArrayBox> m_Pie_xz_cc, m_Pie_xz_cc0, m_Pie_xz_cc_old, m_Pie_xzD_cc;   
   LevelData<EdgeDataBox> m_Pie_yy_ce, m_Pie_yy_ce0, m_Pie_yy_ce_old, m_Pie_yyD_ce;   
   LevelData<EdgeDataBox> m_divPie_ce;   
   LevelData<FluxBox> m_enEleFluxVisc_cf;
   LevelData<FluxBox> m_enEleFluxVisc_norm;

   // below are for relaxation scheme
   //
   LevelData<FluxBox> m_rhoDen_cf;
   LevelData<FluxBox> m_m0JaFluxVisc_cf0;
   LevelData<FluxBox> m_m1JaFluxVisc_cf0;
   LevelData<FArrayBox> m_momVisc_source0;
   LevelData<FluxBox> m_m0JaFluxVisc_cf_old;
   LevelData<FluxBox> m_m1JaFluxVisc_cf_old;
   LevelData<FArrayBox> m_momVisc_source_old;
   
   // containers for heat flux
   //
   LevelData<FArrayBox> m_Tele_cc, m_Tion_cc;
   LevelData<FluxBox> m_Tele_cf, m_Tion_cf;
   LevelData<FluxBox> m_Dqe, m_Dqi;
   //
   LevelData<FluxBox> m_qe_old, m_qe_cf, m_qe_norm, m_qe0, m_qe_wedge;
   LevelData<FluxBox> m_qi_old, m_qi_cf, m_qi_norm, m_qi0, m_qi_wedge;
   LevelData<FArrayBox> m_kappae, m_kappae_wedge;
   LevelData<FArrayBox> m_kappai, m_kappai_wedge;
   LevelData<FluxBox> m_kappae_cf, m_kappae_wedge_cf;
   LevelData<FluxBox> m_kappai_cf, m_kappai_wedge_cf;
  
   ///////////////////////
   
   // cell center vars
   //
   LevelData<FArrayBox> m_rhoDen_cc;
   LevelData<FArrayBox> m_momDen_cc;
   LevelData<FArrayBox> m_eneDeni_cc, m_eneDene_cc;
   LevelData<FArrayBox> m_momDen_virtual_cc;
   LevelData<FArrayBox> m_magField_virtual_cc;
   
   // allocate containers for derived variables
   //
   LevelData<FArrayBox> m_pressure;
   LevelData<FArrayBox> m_pressure_ele;
   LevelData<FArrayBox> m_pressure_ion;
   LevelData<FArrayBox> m_Bpressure;
   LevelData<FArrayBox> m_velocity;
   LevelData<FluxBox> m_velocity_cf;
   LevelData<FArrayBox> m_velocity_virtual;
   
   // allocate containers for fluid variable fluxes at cell center
   //
   LevelData<FArrayBox> m_rhoFlux_cc;
   LevelData<FArrayBox> m_mxFlux_cc;
   LevelData<FArrayBox> m_myFlux_cc;
   LevelData<FArrayBox> m_mzFlux_cc;
   LevelData<FArrayBox> m_enFluxi_cc;
   LevelData<FArrayBox> m_enFluxe_cc;
   LevelData<FArrayBox> m_mvFlux_cc;
   
   // allocate containers for Cspeed
   //
   LevelData<FArrayBox>  m_Cspeed_cc;
   LevelData<FluxBox> m_CspeedL_norm;
   LevelData<FluxBox> m_CspeedR_norm;
   LevelData<FArrayBox> m_CLight_cc;
   
   // allocate containers for normal fluxes on faces
   //
   LevelData<FluxBox> m_rhoFlux_norm;
   LevelData<FluxBox> m_mxFlux_norm;
   LevelData<FluxBox> m_myFlux_norm;
   LevelData<FluxBox> m_mzFlux_norm;
   LevelData<FluxBox> m_enFluxi_norm;
   LevelData<FluxBox> m_enFluxe_norm;
   LevelData<FluxBox> m_enFluxeHall_norm;
   LevelData<FluxBox> m_momFlux_norm;
   LevelData<FluxBox> m_mvFlux_norm;
   LevelData<FArrayBox> m_div_IdentFlux_R;
   LevelData<FArrayBox> m_div_IdentFlux_Z;
 
   ////////////////////////
  
   // allocate containers for magnetic field stuff
   //
   //LevelData<NodeFArrayBox> virtualEleField;
   LevelData<NodeFArrayBox> m_divE;
   LevelData<EdgeDataBox> m_deltaN_ce;
   LevelData<FArrayBox> m_deltaN_cc;
   //LevelData<FArrayBox>* m_deltaN_Ptr = &m_deltaN_cc;
   //
   LevelData<EdgeDataBox>  m_E0_ce;
   LevelData<EdgeDataBox>  m_J0_ce;
   LevelData<FArrayBox>  m_J0_cc;
   LevelData<FArrayBox>  m_E0_cc;
   LevelData<FArrayBox>  m_curlE;
   LevelData<FArrayBox>  m_By_contr_cc;
   LevelData<FArrayBox>  m_By_covar_cc;
   LevelData<FArrayBox>  m_JaSeneDeni;   // ion energy den source
   LevelData<FArrayBox>  m_JaSeneDene;   // ele energy den source
   LevelData<FArrayBox>  m_JaNUdotE;     // Jacobian*N*U\cdotE
   LevelData<FArrayBox>  m_JaJdotE;      // Jacobian*J\cdotE
   LevelData<FArrayBox>  m_J0crossBv_cc; // Jacobian*J\timesB
   
   // allocate containers for stuff used to enforce flux BCs
   //
   LevelData<FluxBox> m_rhoFlux_ce;
   LevelData<FluxBox> m_mxFlux_ce;
   LevelData<FluxBox> m_myFlux_ce;
   LevelData<FluxBox> m_mzFlux_ce;
   LevelData<FluxBox> m_enFluxi_ce;
   LevelData<FluxBox> m_enFluxe_ce;
   LevelData<FluxBox> m_mvFlux_ce;
   
   ///////////////////
   
   LevelData<EdgeDataBox> m_EdgeBC_zeros;
   //LevelData<EdgeDataBox> m_EdgeBC_res;
   LevelData<FluxBox> m_FluxBC_zeros;
   LevelData<FluxBox> m_rhoFluxBC_norm;
   LevelData<FluxBox> m_mxFluxBC_norm;
   LevelData<FluxBox> m_myFluxBC_norm;
   LevelData<FluxBox> m_mzFluxBC_norm;
   LevelData<FluxBox> m_enFluxiBC_norm;
   LevelData<FluxBox> m_enFluxeBC_norm;
   LevelData<FluxBox> m_mvFluxBC_norm;
   LevelData<FluxBox> m_momFluxBC_norm;
   //
   LevelData<FArrayBox> m_dummyDiv, m_dummyFArray_SpaceDim;
   LevelData<FArrayBox> m_dummyDiv_mom, m_dummyFArray_oneComp;
   LevelData<FluxBox> m_dummyFlux_oneComp, m_dummyFlux_SpaceDim;


   // allocate containers for things on edges
   //
   LevelData<EdgeDataBox> m_NonEdges;         // density on cell edges
   LevelData<EdgeDataBox> m_Eold;     // EonEdge at previous time-step
   LevelData<EdgeDataBox> m_Jold;     // JonEdge at previous time-step
   LevelData<FArrayBox> m_Ecc;   // electric field at cell center (phys)
   LevelData<FArrayBox> m_Jcc;   // current density at cell center
   LevelData<FArrayBox> m_Jcc_res;

   /////////////////////

   LevelData<FArrayBox>    m_VhallBy_cc; // Vhall*By covar at cell center
   LevelData<FluxBox>      m_VhallBy_cf; // Vhall*By covar at face center
   LevelData<EdgeDataBox>  m_Ehall_ce;   // Ehall = -VhallxB covar at cell edges
   LevelData<EdgeDataBox>  m_gfun_covar;      // covar g-function in Hall equation
   LevelData<EdgeDataBox>  m_gfun_covar_op;   // covar g-fun on opposite edges
   LevelData<FArrayBox>    m_Vhall_cc;   // Vhall = -lambdai*J0/N contravar at cell center
   LevelData<FArrayBox>    m_VhallDrift_contra; // VhallDrift = lambdai*Bcrossgrad(N)/N^2 contravar at cell center
   LevelData<FluxBox>      m_Vhall_cf;   // Vhall = -lambdai*J0/N contravar at face center
   LevelData<FArrayBox>    m_betaHall_cc;   // |Vhall|/cvac
   LevelData<EdgeDataBox>  m_betaHall_ce;   // |Vhall|/cvac
   LevelData<FArrayBox>    m_Ve_cc;      // Ve = U - lambdai*J0/N at cell center
   LevelData<EdgeDataBox>  m_Ve_ce;      // Ve = U - lambdai*J0/N contravar at cell edges
   LevelData<FluxBox>      m_Ve_cf;      // Ve = U - lambdai*J0/N contravar at face edges
  
   LevelData<EdgeDataBox>  m_etaHall_ce;    // Ehall = Li/N*JxB = J x etaHall
   LevelData<FluxBox>      m_etaHall_cf;    // Ehall = Li/N*JxB = J x etaHall

   LevelData<FArrayBox> m_gradPe_covar_cc;
   LevelData<FArrayBox> m_gradTeeff_covar_cc;
   LevelData<FArrayBox> m_gradT_cc;
   LevelData<EdgeDataBox> m_gradPe_covar_ce;  // covar grad(Pe) on cell edges

   ////////////////////
   
   LevelData<FArrayBox> m_taui, m_taue;
   LevelData<FArrayBox> m_etaRes, m_Qie;
   LevelData<EdgeDataBox> m_etaRes_ce;
   LevelData<FArrayBox> m_xion, m_xele;
   

   /// Parse parameters.
   /**
    * Private method to obtain control parameters from "field_model.field_name" section
    * of ParmParse database, where "field_name" is replaced by the field name 
    * (e.g., field_model.VectPotential)
    *
    * \param[in] ppcls field_model.field input database.
    */
   void parseParameters( ParmParse& pp, 
                         const double, const double,
                         const double, const double );
   

   /// Print parameters.
   /**
    * Private method to print parameters.
    */
   void printParameters();
   
   //FluidVarBC*       m_fluid_variable_bc;
   std::vector<FluidVarBC*>  m_fluid_bc;
   std::vector<FluidVarBC*>  m_fluid_face_var_bc;
   std::vector<FluidVarBC*>  m_fluid_edge_var_bc;
  
   void accumulateMhdRHS( FluidSpecies&  a_rhs_fluid,
                    const FluidSpecies&  a_soln_fluid,
                    const Real           a_time );
   
   void accumulateMaxwellRHS( FluidSpecies&  a_rhs_fluid,
                        const FluidSpecies&  a_soln_fluid,
                        const Real           a_time );
   
   void accumulateViscRHS( FluidSpecies&  a_rhs_fluid,
                     const FluidSpecies&  a_soln_fluid );
   
   void accumulateHeatRHS( FluidSpecies&  a_rhs_fluid,
                     const FluidSpecies&  a_soln_fluid );
   
   void convertPhysFluxToMappedFlux(LevelData<FArrayBox>& a_mapped_flux,
                              const LevelData<FArrayBox>& a_physical_flux,
                              const LevelData<FArrayBox>& a_N) const;
   
   void convertPhysFluxToMappedFlux( FArrayBox& a_mapped_flux,
                               const FArrayBox& a_physical_flux,
                               const FArrayBox& a_N) const;

   void computeNTFfaceArea( LevelData<FluxBox>&    a_Flux_norm,
                      const LevelData<FArrayBox>&  a_Flux_cc_mapped,
                      const LevelData<FArrayBox>&  a_Jf_cc,
                      const LevelData<FArrayBox>&  a_Cspeed_cc,
                      const LevelData<FluxBox>&    a_CspeedR_norm,
                      const LevelData<FluxBox>&    a_CspeedL_norm ) const;
   
   void computeIdealEatEdges( LevelData<EdgeDataBox>&  a_Edge_covar,
                        const LevelData<FArrayBox>&    a_V_phys_cc,
                        const LevelData<FArrayBox>&    a_B_contr_cc,
                        const LevelData<FArrayBox>&    a_Cspeed_cc ) const;
   
   void computeHallEatEdges( LevelData<EdgeDataBox>&  a_etaHall_ce,
                             LevelData<EdgeDataBox>&  a_Ehall_covar,
                       const LevelData<FArrayBox>&    a_J_contr_cc,
                       const LevelData<FArrayBox>&    a_B_contr_cc,
                       const int                      a_NoUpdateBonFaces );
   
   void setHallDriftSpeed( const LevelData<FArrayBox>&  a_N,
                           const LevelData<FArrayBox>&  a_B_phys );

   void computeElectronVelForVisc( LevelData<EdgeDataBox>&  a_Ve_en,
                             const int                      a_method,
                             const LevelData<EdgeDataBox>&  a_J_ce,
                             const LevelData<EdgeDataBox>&  a_N_ce,
                             const LevelData<FArrayBox>&    a_N_cc );

   void advanceElectronViscFluxes( LevelData<NodeFArrayBox>&  a_Pie_diag_cn,
                                   LevelData<FArrayBox>&      a_Pie_xz_cc,
                                   LevelData<EdgeDataBox>&    a_Pie_yy_ce,
                             const LevelData<EdgeDataBox>&    a_Ve_ce,
                             const int                        a_stage,
                             const double                     a_time,
                             const double                     a_dt );

   void computeStrainTensorStag( LevelData<NodeFArrayBox>&  a_Pie_diag_cn,
                                 LevelData<FArrayBox>&      a_Pie_xz_cc,
                                 LevelData<EdgeDataBox>&    a_Pie_yy_ce,
                           const LevelData<EdgeDataBox>&    a_J_ce );

   void computeDivOfElectronVisc( LevelData<EdgeDataBox>&    a_divPie_ce,
                            const LevelData<NodeFArrayBox>&  a_Pie_para_cn,
                            const LevelData<FArrayBox>&      a_Pie_xz_cc,
                            const LevelData<EdgeDataBox>&    a_Pie_yy_ce );

   void addDivPieToGradPele( LevelData<EdgeDataBox>&  a_gradPe_ce,
                       const LevelData<EdgeDataBox>&  a_divPie_ce );
   
   void setCourantTimeStep(const LevelData<FArrayBox>& a_Cspeed);

   void defineLevelDatas( const DisjointBoxLayout&  a_grids,
                          const IntVect&            a_ghostVect );
   
   void setCellCenterValues( const FluidSpecies&  a_soln_fluid );
   
   void setCellCenterFluxes( const FluidSpecies&  a_soln_fluid );
  
   void addHallEnergyFlux( const FluidSpecies&  a_soln_fluid );

   void setMagneticFieldTerms( const FluidSpecies&  a_soln_fluid );

   void setFaceCenteredFluxes( const FluidSpecies&  a_soln_fluid );
   
   void enforceFluxBCs( FluidSpecies&  a_rhs_fluid,
                  const FluidSpecies&  a_soln_fluid,
                  const Real           a_time );

   void updateCollisionTerms( const FluidSpecies&          a_soln_fluid,
                              const LevelData<FArrayBox>&  a_Pele, 
                              const LevelData<FArrayBox>&  a_Pion );
   
   void setHeatFluxCoefficients( const double                 a_gamma,
                                 const LevelData<FArrayBox>&  a_taue, 
                                 const LevelData<FArrayBox>&  a_taui, 
                                 const LevelData<FArrayBox>&  a_xele, 
                                 const LevelData<FArrayBox>&  a_xion, 
                                 const LevelData<FArrayBox>&  a_Pele, 
                                 const LevelData<FArrayBox>&  a_Pion );

   void computeCollHeatFlux( LevelData<FluxBox>&    a_q,
                             LevelData<FluxBox>&    a_Dq,
                       const LevelData<FArrayBox>&  a_Temp_cc,
                       const LevelData<FluxBox>&    a_Temp_cf,
                       const LevelData<FluxBox>&    a_Dens_cf,
                       const LevelData<FluxBox>&    a_kappa_cf,
                       const int                    a_charge,
                       const double                 a_tauMin,
                       const double                 a_limiter,
                       const double                 a_V0sq,
                       const double                 a_gamma );
   
   void setIonViscosityCoefficients( const LevelData<FArrayBox>&  a_taui, 
                                     const LevelData<FArrayBox>&  a_xion, 
                                     const LevelData<FArrayBox>&  a_Pion ); 
 
   void computeViscFluxes( LevelData<FluxBox>&    a_m0JaFluxVisc_cf,
                           LevelData<FluxBox>&    a_m1JaFluxVisc_cf,
                           LevelData<FArrayBox>&  a_momVisc_source,
                     const LevelData<FluxBox>&    a_velocity_phys_cf,
                     const LevelData<FArrayBox>&  a_velocity_phys,
                     const LevelData<FluxBox>&    a_etaVisC_cf,
                     const LevelData<FArrayBox>&  a_etaVisC );
   
   void computeViscEnergyFlux( LevelData<FluxBox>&  a_UJaPi_Flux_phys,
                         const LevelData<FluxBox>&  a_JaPix_Flux_phys,
                         const LevelData<FluxBox>&  a_JaPiy_Flux_phys,
                         const LevelData<FluxBox>&  a_velocity_phys );

   void computeEleViscEnergyFlux( LevelData<FluxBox>&        a_UePie_Flux_phys,
                            const LevelData<FArrayBox>&      a_Pie_xz_cc,
                            const LevelData<NodeFArrayBox>&  a_Pie_diag_cn,
                            const LevelData<FluxBox>&        a_velocity_ele_cf );
   
   void computeEleVelOnCells( LevelData<FArrayBox>&  a_Ue_cc,
                        const LevelData<FArrayBox>&  a_J_cc,
                        const LevelData<FArrayBox>&  a_N_cc,
                        const LevelData<FArrayBox>&  a_U_cc );
   
   void computeViscSourceAxisymm( LevelData<FArrayBox>&  a_Wthth,
                            const LevelData<FArrayBox>&  a_velocity_phys,
                            const LevelData<FluxBox>&    a_velocity_phys_cf );
  
   void addGyroViscToFluxes( LevelData<FluxBox>&    a_m0JaFluxVisc_cf,
                             LevelData<FluxBox>&    a_m1JaFluxVisc_cf,
                       const LevelData<FArrayBox>&  a_velocity_phys,
                       const LevelData<FluxBox>&    a_etaVisG_cf );

   void addGyroHeatToViscFluxes( LevelData<FluxBox>&    a_m0JaFluxVisc_cf,
                                 LevelData<FluxBox>&    a_m1JaFluxVisc_cf,
                           const double                 a_gamma,
                           const LevelData<FArrayBox>&  a_qwedge_phys,
                           const LevelData<FluxBox>&    a_etaVisG_cf );

   void computeGyroHeatFluxPhys( LevelData<FArrayBox>&  a_q_wedge,
                           const LevelData<FArrayBox>&  a_Temp,
                           const LevelData<FArrayBox>&  a_kappa_wedge,
                           const double                 a_charge );

   void addGyroHeatToFlux( LevelData<FluxBox>&    a_q,
                     const LevelData<FArrayBox>&  a_Temp_cc,
                     const LevelData<FluxBox>&    a_Temp_cf,
                     const LevelData<FluxBox>&    a_Dens_cf,
                     const LevelData<FluxBox>&    a_kappa_cf,
                     const double                 a_charge,
                     const double                 a_limiter,
                     const double                 a_V0sq );

   void addGyroHeatToSource( LevelData<FArrayBox>&  a_RHS,
                       const LevelData<FArrayBox>&  a_Temp_cc,
                       const LevelData<FArrayBox>&  a_Dens_cc,
                       const LevelData<FArrayBox>&  a_kappaG_cc,
                       const double                 a_charge,
                       const double                 a_limiter,
                       const double                 a_V0sq );
   
   void updateRHSs_heat( FluidSpecies&  a_rhs_fluid,
                   const FluidSpecies&  a_soln_fluid );
   
   void updateRHSs_visc( FluidSpecies&  a_rhs_fluid,
                   const FluidSpecies&  a_soln_fluid );

   void updateRHSs( FluidSpecies&  a_rhs_fluid,
              const FluidSpecies&  a_soln_fluid );
   
   void updateMaxwellRHSs( FluidSpecies&  a_rhs_fluid,
                     const FluidSpecies&  a_soln_fluid );
   
   void setMappedCspeed( const double  a_gamma,
                         const int     a_Ptherm,
                         const int     a_Pmag );
   
   void setMappedCLight();

   void computeGfunEdges( LevelData<EdgeDataBox>&  a_gfun_covar,
                    const double                   a_dt,
                    const LevelData<EdgeDataBox>&  a_E0_covar,
                    const LevelData<EdgeDataBox>&  a_gradPe_covar,
                    const LevelData<EdgeDataBox>&  a_betaHall,
                    const LevelData<EdgeDataBox>&  a_J0_covar,
                    const LevelData<EdgeDataBox>&  a_NonEdges,
                    const LevelData<EdgeDataBox>&  a_Eold,
                    const LevelData<EdgeDataBox>&  a_Jold ) const;

   void advanceEandJ( LevelData<EdgeDataBox>&  soln_Ecovar,
                      LevelData<EdgeDataBox>&  soln_Jcovar,
                const double                   a_dt,
                const LevelData<EdgeDataBox>&  a_E0_covar,
                const LevelData<EdgeDataBox>&  a_eta_ce,
                const LevelData<EdgeDataBox>&  a_Ehall_covar,
                const LevelData<EdgeDataBox>&  a_etaHall_ce,
                const LevelData<EdgeDataBox>&  a_gfun_covar,
                const LevelData<EdgeDataBox>&  a_gfun_covar_op,
                const LevelData<EdgeDataBox>&  a_gradPe_covar,
                const LevelData<EdgeDataBox>&  a_J0_covar,
                const LevelData<EdgeDataBox>&  a_NonEdges,
                const LevelData<EdgeDataBox>&  a_Eold,
                const LevelData<EdgeDataBox>&  a_Jold,
                const int                      a_dir,
                const int                      a_method ) const;

};


#include "NamespaceFooter.H"

#endif
