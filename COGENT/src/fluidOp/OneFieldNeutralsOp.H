#ifndef  _ONEFIELDNEUTRALSOP_H_
#define  _ONEFIELDNEUTRALSOP_H_

#include "FluidOpInterface.H"
#include "FluidVarBC.H"
#include "Diffusion.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "GKVlasov.H"
#include "FluidOpPreconditioner.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include <string>

#include "NamespaceHeader.H"

/**
 * This class computes the RHS for a single-field neutral (gas) diffusion model:
 *
 *    drho_g/dt = \nabla( D \nabla(rho_g T_g)) - nu_iz * rho_g + Source
 *
 * rho_g is the mass density (rho_g = n_g *m_g)
 * D = 1/(m_g nu_cx) is the diffusion tensor coefficient
 * nu_cx(n_i, T_i) is the charge-exchange coefficient
 * nu_iz(n_e,T_e) is the ionization coefficinet
 * Source is a prescribed source function.
 *
 * T_g=T_i is assumed (and the default value is unity).
 * T_i and n_e and can be either a prescribed function
 * or computed self-consistently from ion kinetic species.
 *
 * Recycling BC can be used with this model to set the
 * onput neutral particle flux as the recycling coefficient
 * times the ion outflow flux.
 *
*/

class OneFieldNeutralsOp
   : public FluidOpInterface
{
public:
   
  /// Null constructor (used for state variable creation) 
  OneFieldNeutralsOp( const MagGeom&  geometry )
     : m_geometry(geometry) {};
   
  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  OneFieldNeutralsOp(const string&   pp_str,
                     const string&   species_name,
                     const MagGeom&  geometry,
                     const Real      larmor,
                     const int       verbosity );

  /// Destructor.
  /**
   */
  virtual ~OneFieldNeutralsOp();

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
                              const int                         component);

   virtual void solvePCImEx( FluidSpeciesPtrVect&              fluid_species_solution,
                             const PS::KineticSpeciesPtrVect&  kinetic_species_rhs,
                             const FluidSpeciesPtrVect&        fluid_species_rhs,
                             const int                         component );

   virtual void preOpEval( const PS::KineticSpeciesPtrVect&   kinetic_species,
                           const FluidSpeciesPtrVect&         fluid_species,
                           const PS::ScalarPtrVect&           scalars,
                           const EField&                      E_field,
                           const double                       time );
   
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
   Real computeDtImExTI( const FluidSpeciesPtrVect&  fluid_species );
   
   static const char* pp_name;
   
protected:
   
   virtual void getPlasmaParameters(const Real time);
   
   void computeDiffusionWithRecyclingBC(LevelData<FArrayBox>&        out,
                                        FluidVarBC&                  recycl_var_bc,
                                        const LevelData<FArrayBox>&  in,
                                        const Diffusion&             diffusion_op,
                                        const LevelData<FluxBox>&    ion_normal_flux,
                                        const Real                   time) const;
   
   virtual void computeDiffusionCoefficients(LevelData<FluxBox>& D_tensor,
                                             LevelData<FluxBox>& D_tensor_mapped,
                                             const Real          mass,
                                             const Real          time) const;
   
   /// Computes ion mass density
   void computeIonChargeDensity(LevelData<FArrayBox>&             mass_density,
                                const PS::KineticSpeciesPtrVect&  species ) const;

   /// Computes normal ion particle flux integrated over a face
   void computeIntegratedNormalIonParticleFlux(LevelData<FluxBox>&                ion_particle_flux,
                                               const PS::KineticSpeciesPtrVect&   species,
                                               const EField&                      E_field,
                                               const Real&                        time) const;

   
   void computeIonizationRate(LevelData<FArrayBox>& ionization_rate,
                              const Real            time) const;
   
   void computeChargeExchangeRate(LevelData<FArrayBox>&   chx_rate,
                                  const Real              time) const;
   
   void convertCellToFace(LevelData<FluxBox>&         a_face_data,
                          const LevelData<FArrayBox>& a_cell_data) const;
   
   void setZero( LevelData<FArrayBox>& data ) const;
   void setZero( LevelData<FluxBox>& data ) const;
   
   IntVect m_ghostVect;

   LevelData<FluxBox> m_ion_normal_flux;
   LevelData<FArrayBox> m_ne, m_Te, m_Tg;
   RefCountedPtr<GridFunction> m_Tg_func;

   Real m_larmor;

   bool m_consistent_diffusion;
   bool m_consistent_density;
   bool m_consistent_temperature;
   
private:
   
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
   
   virtual void accumulateRHS(FluidSpeciesPtrVect&               rhs,
                              const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                              const FluidSpeciesPtrVect&         fluid_species,
                              const PS::ScalarPtrVect&           scalars,
                              const EField&                      E_field,
                              const int                          component,
                              const Real                         time);
   
   
   PS::GKVlasov*  m_vlasov;
   
   Diffusion*        m_diffusion_op;
   EllipticOpBC*     m_bc;
   FluidVarBC*       m_fluid_variable_bc;
  
   RefCountedPtr<GridFunction> m_D_rad_func;
   RefCountedPtr<GridFunction> m_D_par_func;
   RefCountedPtr<GridFunction> m_D_perp_func;
   RefCountedPtr<GridFunction> m_fixed_source_func;
   RefCountedPtr<GridFunction> m_ne_func;
   RefCountedPtr<GridFunction> m_Te_func;

   LevelData<FArrayBox> m_volume;
   LevelData<FArrayBox> m_D_rad, m_D_par, m_D_perp;
   
   const MagGeom& m_geometry;
   
   int m_verbosity;

   bool m_uedge_chx;
   bool m_dens_recycling_bc;
   bool m_is_time_implicit;
   bool m_first_call;
   bool m_include_ionization;
   
   std::string m_opt_string;
   int m_my_pc_idx;
   
};


#include "NamespaceFooter.H"

#endif
