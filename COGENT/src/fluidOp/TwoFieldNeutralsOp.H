#ifndef  _TWOFIELDNEUTRALSOP_H_
#define  _TWOFIELDNEUTRALSOP_H_

#include "FluidOpInterface.H"
#include "FluidVarBC.H"
#include "Diffusion.H"
#include "OneFieldNeutralsOp.H"
#include "TwoFieldNeutralsPC.H"
#include "SpaceUtils.H.multidim"

#include "NamespaceHeader.H"

/**
 * Single fluid operator class.
*/
class TwoFieldNeutralsOp
   : public OneFieldNeutralsOp
{
public:

  /// Null constructor (used for state variable creation) 
  TwoFieldNeutralsOp( const MagGeom&  geometry )
     :OneFieldNeutralsOp(geometry),
      m_geometry(geometry) {};

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  TwoFieldNeutralsOp(const string&   pp_str,
                     const string&   species_name,
                     const MagGeom&  geometry,
                     const Real      a_larmor,
                     const int       verbosity );

  /// Destructor.
  /**
   */
  virtual ~TwoFieldNeutralsOp();

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
   Real computeDt( const FluidSpeciesPtrVect&  fluid_species );
   
   /// Set LevelData<FArrayBox> object for plotting
   virtual void getMemberVar(LevelData<FArrayBox>&   a_Var,
                             const FluidSpecies&     a_fluid_species,
                             const string&           a_name ) const;
   
   static const char* pp_name;

protected:
   
   virtual void getPlasmaParameters(const Real time);

   virtual void computeDiffusionCoefficients(LevelData<FluxBox>& D_tensor,
                                             LevelData<FluxBox>& D_tensor_mapped,
                                             const Real          mass,
                                             const Real          time) const;
   
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

   void accumulateDiffusiveTerms(FluidSpecies&        fluid_rhs,
                                 const FluidSpecies&  fluid_soln,
                                 const Real           time) const;
   
   void accumulateParAdvectionTerms(FluidSpecies&        fluid_rhs,
                                    const FluidSpecies&  fluid_soln,
                                    const Real           time) const;

   
   void accumulatePerpAdvectionTerms(FluidSpecies&        fluid_rhs,
                                     const FluidSpecies&  fluid_soln,
                                     const Real           time) const;
   
   void accumulateSourceTerms(FluidSpecies&        fluid_rhs,
                              const FluidSpecies&  fluid_soln,
                              const Real           time) const;
   
   void computeAdvectionRHS(LevelData<FArrayBox>&        div,
                            const LevelData<FArrayBox>&  dist,
                            const LevelData<FluxBox>&    velocity,
                            const bool                   recycling_bc,
                            const bool                   homogeneous_flux_bc,
                            const Real                   time) const;
   
   void computeFlux(LevelData<FluxBox>&       flux,
                    const LevelData<FluxBox>& faceDist,
                    const LevelData<FluxBox>& velocity,
                    const bool                fourth_order) const;
   
   void computeParallelPressureGradient(LevelData<FArrayBox>&         par_grad,
                                        const LevelData<FArrayBox>&   soln_density,
                                        const Real                    mass) const;
   
   void computePerpPressureGradient(LevelData<FluxBox>&           perp_grad,
                                    const LevelData<FArrayBox>&   soln_density,
                                    const Real                    mass) const;
   
   
   void computeViscosityCoefficients(LevelData<FluxBox>&         ellip_coeff,
                                     LevelData<FluxBox>&         ellip_coeff_mapped,
                                     const LevelData<FArrayBox>& density,
                                     const Real                  mass,
                                     const Real                  time) const;
   
   void computeAtomicRates(LevelData<FArrayBox>&   ionization_rate,
                           LevelData<FArrayBox>&   recombination_rate,
                           LevelData<FArrayBox>&   chx_rate,
                           const Real              time) const;

   void computeIntegratedNormalPrescribedFlux(LevelData<FluxBox>&         flux,
                                              const LevelData<FArrayBox>& ne,
                                              const LevelData<FArrayBox>& Vi) const;
   

   const MagGeom& m_geometry;
   
   int m_verbosity;
   int m_my_pc_idx;

   Real m_bckgr_density;

   bool m_is_time_implicit;
   bool m_include_advection_bc;
   bool m_dens_recycling_bc;
   bool m_mom_recycling_bc;
   bool m_consistent_velocity;
   bool m_fixed_recycling;
   bool m_first_call;
   bool m_extrapolate_density;

   std::string m_advScheme;
   std::string m_opt_string;
   
   bool m_analytic_pc; // Setting for enabling the analytic preconditioner
   
   TwoFieldNeutralsPC* m_twofieldneutrals_pc; // The analytic preconditioner

   Diffusion*        m_diffusion_op;
   EllipticOpBC*     m_diffusion_bc;

   Diffusion*        m_viscosity_op;
   mutable EllipticOpBC*     m_viscosity_bc;

   std::vector<FluidVarBC*>  m_fluid_bc;
   
   LevelData<FArrayBox> m_ionization_rate;
   LevelData<FArrayBox> m_recombination_rate;
   LevelData<FArrayBox> m_chx_rate;
   LevelData<FArrayBox> m_ion_parallel_vel;
   LevelData<FluxBox>   m_ion_parallel_vel_face;
   LevelData<FArrayBox> m_pc_density;
   
   RefCountedPtr<GridFunction> m_ion_parallel_vel_func;
   
};


#include "NamespaceFooter.H"

#endif
