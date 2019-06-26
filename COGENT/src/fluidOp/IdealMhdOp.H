#ifndef  _IDEALMHDOP_H_
#define  _IDEALMHDOP_H_

#include "FluidOpInterface.H"
#include "FluidVarBC.H"
#include "Diffusion.H"
#include "SpaceUtils.H"

#include "EdgeDataBox.H"

#include "NamespaceHeader.H"

/**
 * multi-dim ideal MHD equation operator class.
 *
 * d(rho)/dt + div(m) = 0
 * d(mom)/dt + div(u*m) = JxB - grad(P)
 * d(ene)/dt + div(u*(ene+P)) = JdotE
 * dB/dt     + cur(E) = 0
 *
 * rho = M*n, M is mass n is density
 * mom_x = ux*n
 * mom_y = uy*n
 * mom_z = uz*n
 * ene = 0.5*n*u^2 + P/(g-1)
 *
 * mu0*J = cur(B)
 * E = uxB
 *
*/
class IdealMhdOp
   : public FluidOpInterface
{
public:

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

   virtual bool trivialSolutionOp() const {return true;}

   virtual CFGVars* createStateVariable( const string&       a_pp_prefix,
                                         const std::string&  a_name,
                                         const MagGeom&      a_geometry,
                                         const IntVect&      a_ghost_vect ) const
   {
      return new FluidSpecies(a_pp_prefix, a_name, a_geometry, a_ghost_vect);
   }

   virtual void fillGhostCells( FluidSpecies&  species_phys,
                                const double   time );

   /// Compute a stable time step.
   /**
    * Computes and returns an estimate of the maximum stable time step.
    */
   Real computeDtExplicitTI( const FluidSpeciesPtrVect&  fluid_species );
   //Real computeDtImExTI( const FluidSpeciesPtrVect&  fluid_species );
   
   
   static const char* pp_name;

private:

   int m_verbosity;
   double m_courant_time_step;
   const MagGeom& m_geometry;
   int m_is_time_implicit;
   std::string m_advScheme;
   std::string m_species_name;

   std::string m_opt_string;
   int m_my_pc_idx;
      
   // allocate contains for metric stuff to transform from phys to mapped
   // (eventually these will be removed as function are moved to magGeom.cp)
   LevelData<FArrayBox> m_Jacobian;
   LevelData<FArrayBox> m_Nmatrix;
   LevelData<FArrayBox> m_X_phys;
   LevelData<FArrayBox> m_g_y; // used only when solving virtual B
   LevelData<FArrayBox> m_cellVol;
   bool m_twoDaxisymm;
   
   // allocate containers for species fluid variables
   //
   LevelData<FArrayBox> m_rhoDen_cc;
   LevelData<FArrayBox> m_momDen_cc;
   LevelData<FArrayBox> m_eneDen_cc;
   LevelData<FArrayBox> m_momDen_virtual_cc;
   LevelData<FArrayBox> m_magField_virtual_cc;
   
   // allocate containers for derived variables
   //
   LevelData<FArrayBox> m_pressure;
   LevelData<FArrayBox> m_Bpressure;
   LevelData<FArrayBox> m_velocity;
   LevelData<FArrayBox> m_velocity_virtual;
   
   // allocate containers for fluid variable fluxes at cell center
   //
   LevelData<FArrayBox> m_rhoFlux_cc;
   LevelData<FArrayBox> m_mxFlux_cc;
   LevelData<FArrayBox> m_myFlux_cc;
   LevelData<FArrayBox> m_mzFlux_cc;
   LevelData<FArrayBox> m_enFlux_cc;
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
   LevelData<FluxBox> m_momFlux_norm;
   LevelData<FluxBox> m_mvFlux_norm;
   LevelData<FArrayBox> m_div_IdentFlux_R;
   LevelData<FArrayBox> m_div_IdentFlux_Z;
   
   // allocate containers for magnetic field stuff
   //
   //LevelData<NodeFArrayBox> virtualEleField;
   LevelData<EdgeDataBox>  m_E0_ce;
   LevelData<EdgeDataBox>  m_J0_ce;
   LevelData<FArrayBox>  m_J0_cc;
   LevelData<FArrayBox>  m_E0_cc;
   LevelData<FArrayBox>  m_curlE0;
   LevelData<FArrayBox>  m_By_contr_cc;
   LevelData<FArrayBox>  m_By_covar_cc;
   LevelData<FArrayBox>  m_JJdotE; // Jacobian*J\cdotE
   LevelData<FArrayBox>  m_JJcrossB; // Jacobian*J\timesB
   
   // allocate containers for stuff used to enforce flux BCs
   //
   LevelData<FluxBox> m_rhoFlux_ce;
   LevelData<FluxBox> m_mxFlux_ce;
   LevelData<FluxBox> m_myFlux_ce;
   LevelData<FluxBox> m_mzFlux_ce;
   LevelData<FluxBox> m_enFlux_ce;
   LevelData<FluxBox> m_mvFlux_ce;
   //
   LevelData<EdgeDataBox> m_EdgeBC_zeros;
   LevelData<FluxBox> m_FluxBC_zeros;
   LevelData<FluxBox> m_rhoFluxBC_norm;
   LevelData<FluxBox> m_mxFluxBC_norm;
   LevelData<FluxBox> m_myFluxBC_norm;
   LevelData<FluxBox> m_mzFluxBC_norm;
   LevelData<FluxBox> m_enFluxBC_norm;
   LevelData<FluxBox> m_mvFluxBC_norm;
   LevelData<FluxBox> m_momFluxBC_norm;
   //

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
   
   Diffusion*        m_diffusion_op;
   EllipticOpBC*     m_bc;
   //FluidVarBC*       m_fluid_variable_bc;
   std::vector<FluidVarBC*>  m_fluid_bc;
  
   RefCountedPtr<GridFunction> m_D_rad;
   RefCountedPtr<GridFunction> m_D_par;
   RefCountedPtr<GridFunction> m_D_perp;
   RefCountedPtr<GridFunction> m_electron_dens;
   RefCountedPtr<GridFunction> m_electron_temp;

   virtual void accumulateRHS(FluidSpeciesPtrVect&               rhs,
                              const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                              const FluidSpeciesPtrVect&         fluid_species,
                              const PS::ScalarPtrVect&           scalars,
                              const EField&                      E_field,
                              const int                          component,
                              const Real                         time);
   
   void computeDiffusionCoefficients(LevelData<FluxBox>& D_tensor,
                                     LevelData<FluxBox>& D_tensor_mapped,
                                     const Real          time) const;
   
   void inspectFArrayBox(const LevelData<FArrayBox>& a_F0,
                         const int                   a_comp) const;
   
   void inspectFluxBox(const LevelData<FluxBox>& a_Flux,
                       const int                 a_dir) const;
   
   void inspectEdgeDataBox(const LevelData<EdgeDataBox>& a_Edge,
                           const int                 a_dir) const;
   
   
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
   
   void setCourantTimeStep(const LevelData<FArrayBox>& a_Cspeed);

   void defineLevelDatas( const DisjointBoxLayout&  a_grids,
                          const IntVect&            a_ghostVect );
   
   void setCellCenterValues( const FluidSpecies&  a_soln_fluid );
   
   void setCellCenterFluxes( const FluidSpecies&  a_soln_fluid );

   void setMagneticFieldTerms( const FluidSpecies&  a_soln_fluid );

   void setFaceCenteredFluxes( const FluidSpecies&  a_soln_fluid );
   
   void enforceFluxBCs( FluidSpecies&  a_rhs_fluid,
                  const FluidSpecies&  a_soln_fluid,
                  const Real           a_time );
   
   void updateRHSs( FluidSpecies&  a_rhs_fluid,
              const FluidSpecies&  a_soln_fluid );
   
   void setCspeed( const LevelData<FArrayBox>&  a_density,
                   const double                 a_gamma );
   

};


#include "NamespaceFooter.H"

#endif
