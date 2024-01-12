#ifndef  _VORTICITYOP_H_
#define  _VORTICITYOP_H_

#undef ALFVEN_WAVE_TEST

#include "FluidOpInterface.H"
#include "GKFluidOp.H"
#include "FluidOpUtils.H"
#include "GKVorticity.H"
#include "Potentials.H"
#include "EMFields.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "GKVlasov.H"
#include "FluidOpPreconditioner.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include <string>
#include "NamespaceHeader.H"

/**
 * VorticityOp operator class.
*/
class VorticityOp
   : public FluidOpInterface
{
public:

  /// Null constructor (used for state variable creation) 
  VorticityOp( const MagGeom&  geometry )
     : m_geometry(geometry),
       m_vlasov(NULL),
       m_gyropoisson_op(NULL),
       m_imex_pc_op(NULL),
       m_par_cond_op(NULL),
       m_parallel_current_divergence_op(NULL),
       m_vorticity_diffusion_op(NULL),
       m_diamagnetic_correction_op(NULL),
       m_potential_bcs(NULL),
       m_fluid_op_utils(NULL) {};
   
  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  VorticityOp( const ParmParse&   pp,
               const MagGeom&     geometry,
               const double       larmor,
               const int          verbosity );

  /// Destructor.
  /**
   */
  virtual ~VorticityOp();

   /// Evaluates the RHS.
   /**
    *  Pure virtual function that will be resposible for evaluating
    *  field RHS on the given field.
    *
    *  @param[in,out] rhs           -  contribution to be added to rhs.
    *  @param[in] fields            -  current solution for fields.
    *  @param[in] fluids            -  current solution for fluids.
    *  @param[in] kinetic_specties  -  current solution for kinetic species.
    *  @param[in] fieldVecComp      -  component of the field vector to which operator is applied.
    *  @param[in] time              -  the time at which the field RHS is to be evaluated
    */
   virtual void accumulateExplicitRHS( FluidSpeciesPtrVect&               rhs,
                                       const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                                       const FluidSpeciesPtrVect&         fluid_species,
                                       const PS::ScalarPtrVect&           scalars,
                                       const EMFields&                    EM_fields,
                                       const int                          fluid_vec_comp,
                                       const Real                         time );

   virtual void accumulateImplicitRHS( FluidSpeciesPtrVect&               rhs,
                                       const PS::KineticSpeciesPtrVect&   kinetic_species,
                                       const FluidSpeciesPtrVect&         fluid_species,
                                       const PS::ScalarPtrVect&           scalars,
                                       const EMFields&                    EM_fields,
                                       const int                          fluid_vec_comp,
                                       const Real                         time );

   virtual void preSolutionOpEval( const PS::KineticSpeciesPtrVect&   kinetic_species,
                                   const FluidSpeciesPtrVect&         fluid_species,
                                   const PS::ScalarPtrVect&           scalars,
                                   const EMFields&                    EM_fields,
                                   const double                       time );

   virtual void preOpEval( const PS::KineticSpeciesPtrVect&   kinetic_species,
                           const FluidSpeciesPtrVect&         fluid_species,
                           const PS::ScalarPtrVect&           scalars,
                           const EMFields&                    EM_fields,
                           const double                       time );

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

   virtual void evalSolutionOp( FluidSpeciesPtrVect&               rhs,
                                const PS::KineticSpeciesPtrVect&   kinetic_species_comp,
                                const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                                const FluidSpeciesPtrVect&         fluid_species_comp,
                                const FluidSpeciesPtrVect&         fluid_species_phys,
                                const PS::ScalarPtrVect&           scalars,
                                const int                          component,
                                const Real                         time );

   virtual void solveSolutionPC( FluidSpeciesPtrVect&              a_fluid_species_solution,
                                 const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                                 const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                                 const int                         a_component );

   virtual void updatePCImEx( const FluidSpeciesPtrVect&       fluid_species,
                              const PS::KineticSpeciesPtrVect& kinetic_species,
                              const double                     time,
                              const int                        step,
                              const int                        stage,
                              const double                     shift,
                              const int                        component,
                              const std::string&               op_name );

   virtual void solvePCImEx( FluidSpeciesPtrVect&              fluid_species_solution,
                             const PS::KineticSpeciesPtrVect&  kinetic_species_rhs,
                             const FluidSpeciesPtrVect&        fluid_species_rhs,
                             const std::string&                op_name,
                             const int                         component );

   virtual void postTimeEval( FluidSpecies&  species_comp,
                              FluidSpecies&  species_phys,
                              const Real     dt,
                              const Real     time,
                              const int      stage );
   
   virtual bool trivialSolutionOp() const {return false;}

   virtual void setStepConstKinCoeff(const bool a_step_const_kin_coeff)
   {
     if (a_step_const_kin_coeff) {
       if (!procID()) {
         std::cout << "VorticityOp: skipping preconditioner updates ";
         std::cout << "at time integration stages because kinetic coeffs ";
         std::cout << "are assumed constant in a step.\n";
       }
       m_update_pc_skip_stage_i = true;
       m_update_pc_skip_stage_e = true;
     }
   }

   virtual CFGVars* createStateVariable( const string&       a_pp_prefix,
                                         const std::string&  a_name,
                                         const MagGeom&      a_geometry,
                                         const IntVect&      a_ghost_vect ) const
   {
      return new Potentials(a_pp_prefix, a_name, a_geometry, a_ghost_vect, false);
   }

   virtual void fillGhostCells( FluidSpecies&  species_phys,
                                const double   time );

   // Compute the divergence of the perpedicular ion current due solely to magnetic drifts
   void computeDivPerpIonMagCurrentDensity( LevelData<FArrayBox>&             div_Jperp_mag,
                                            const EMFields&                   EM_fields,
                                            const PS::KineticSpeciesPtrVect&  species,
                                            const Real&                       time );

   void computeDivPerpIonExBCurrentDensity( LevelData<FArrayBox>&             div_Jperp,
                                            const EMFields&                   EM_fields,
                                            const PS::KineticSpeciesPtrVect&  species,
                                            const LevelData<FArrayBox>&       ion_charge_density,
                                            const Real&                       time);

   void computeDivPerpElectronMagCurrentDensity(LevelData<FArrayBox>&             div_Jperp,
                                                const LevelData<FArrayBox>&       ion_charge_density,
                                                const LevelData<FArrayBox>&       negative_vorticity,
                                                const Real&                       time);

  
   void  addDiamagneticCorrection(LevelData<FArrayBox>&              rhs,
                                  const LevelData<FArrayBox>&        soln,
                                  const EMFields&                    EM_fields);
   
   void addDivJperpTerm(LevelData<FArrayBox>&              rhs,
                        const LevelData<FArrayBox>&        soln,
                        const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                        const EMFields&                    EM_fields,
                        const Real                         time);
   
   void addReynoldsStressTerm(LevelData<FArrayBox>&              rhs,
                              const LevelData<FArrayBox>&        soln,
                              const EMFields&                    EM_fields);

   void addStabilizationTerms(LevelData<FArrayBox>&              rhs,
                              const LevelData<FArrayBox>&        soln);



   static const char* pp_name;

protected:

   /// Computes ion mass density
   void computeIonMassDensity(LevelData<FArrayBox>&             ion_mass_density,
                              const PS::KineticSpeciesPtrVect&  kinetic_species_phys,
                              const FluidSpeciesPtrVect&        fluid_species_comp) const;

   /// Computes ion charge density
   void computeIonChargeDensity(LevelData<FArrayBox>&            ion_charge_density,
                                const PS::KineticSpeciesPtrVect& kinetic_species_phys,
                                const FluidSpeciesPtrVect&       fluid_species_comp ) const;

   /// Computes ion mass over charge density
    void computeIonMassOverChargeDensity(LevelData<FArrayBox>&            ion_charge_density,
                                         const PS::KineticSpeciesPtrVect& kinetic_species_phys,
                                         const FluidSpeciesPtrVect&       fluid_species_comp ) const;
   
   void setCoreBC( const double   core_inner_bv,
                   const double   core_outer_bv,
                   EllipticOpBC&  bc ) const ;

   
   void setZero( LevelData<FArrayBox>& data ) const;
   void setZero( LevelData<FluxBox>& data ) const;
   
   GKPoisson*    m_gyropoisson_op;
   EllipticOpBC* m_potential_bcs;
   
   GKVorticity*    m_vorticity_diffusion_op;
   RefCountedPtr<EllipticOpBC>  m_vorticity_bcs;
   
   GKVorticity*   m_parallel_current_divergence_op;
   RefCountedPtr<EllipticOpBC>  m_parallel_current_divergence_op_bcs;
   
   Diffusion* m_diamagnetic_correction_op;
   RefCountedPtr<EllipticOpBC>  m_diamagnetic_correction_op_bcs;
   
   FluidOpUtils*      m_fluid_op_utils;
   RefCountedPtr<EllipticOpBC>  m_hyperviscosity_op_bcs;
   
   const MagGeom& m_geometry;

   LevelData<FArrayBox> m_ion_charge_density;
   LevelData<FArrayBox> m_ion_mass_density;
   LevelData<FArrayBox> m_electron_temperature;
   LevelData<FArrayBox> m_negativeDivJpar;
   LevelData<FluxBox> m_pol_diffusion_face;
   LevelData<FluxBox> m_perp_diffusion_face;
   LevelData<FArrayBox> m_divJperp_mag_e;
   LevelData<FArrayBox> m_divJperp_mag_i;
   
   Vector<int> m_harmonic_filtering;
   
   int m_verbosity;
   int m_it_counter;
   
   bool m_reynolds_stress;
   bool m_ExB_current_model;
   bool m_remove_axisymmetric_phi;
   bool m_second_order;
   bool m_low_pollution;
   bool m_include_pol_den_correction;
   bool m_include_diam_correction;
   bool m_use_vorticity_bcs;
   bool m_include_diffusion;
   bool m_include_hyperviscosity;
   bool m_use_hyperviscosity_bcs;
   bool m_suppress_non_zonal_comp_at_boundaries;
   bool m_include_boundary_relaxation;
   bool m_include_stabilization;
   bool m_minimal_implicit_model;
   
   std::string m_advection_scheme;
   
private:
   
   void initializeDiffusionProfiles(LevelData<FluxBox>& pol_diffusion_face,
                                    LevelData<FluxBox>& perp_diffusion_face);
   
   void initializeDiamagneticCoeff(LevelData<FluxBox>& D_tensor,
                                   LevelData<FluxBox>& D_tensor_mapped);
   

   /// Parse parameters.
   void parseParameters( const ParmParse& pp );
   
   /// Print parameters.
   void printParameters();
  
   GKVorticity*    m_par_cond_op;
   GKVorticity*    m_imex_pc_op;

   PS::GKVlasov*  m_vlasov;

   std::string m_opt_string;
   std::vector<int> m_my_pc_idx_i;
   std::vector<std::string> m_pc_names;

   RefCountedPtr<EllipticOpBC> m_imex_pc_op_bcs;
   RefCountedPtr<EllipticOpBC> m_gyropoisson_op_bcs;
   RefCountedPtr<EllipticOpBC> m_par_cond_op_bcs;

   LevelData<FArrayBox> m_volume;
   LevelData<FluxBox> m_face_bXE_integrals;
   
   bool m_sigma_div_e_coefs_set;
   bool m_multiphysics_precond;
   bool m_include_pol_den_correction_precond;
   bool m_include_diffusion_precond;
   bool m_include_pol_den_correction_to_pe;
   bool m_consistent_lower_bc_only;
   bool m_consistent_upper_bc_only;
   bool m_par_cond_op_coefs_defined;
   bool m_update_pc_skip_stage_e;
   bool m_update_pc_skip_stage_i;
   
   int m_my_pc_idx_e;
   int m_num_pc_rhs;
   int m_update_pc_freq_e;
   int m_update_pc_freq_i;

   double m_larmor;

   RefCountedPtr<GridFunction> m_electron_temperature_func;
   RefCountedPtr<GridFunction> m_pol_diffusion_func;
   RefCountedPtr<GridFunction> m_perp_diffusion_func;
};


#include "NamespaceFooter.H"

#endif
