#ifndef _GKVLASOV_H_
#define _GKVLASOV_H_

#include "KineticFunction.H"
#include "TiDefinitions.H"
#include "GlobalDOF.H"
#include "AppCtxt.H"
#include "MultiPhysicsPreconditioner.H"
#include "GKVlasovAMG.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "GridFunction.H"
#include "GridFunctionLibrary.H"
#include "EField.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

class ODEVector;

/**
 * Gyrokinetic Vlasov operator class.
*/
class GKVlasov
{
public:

   /// Constructor with initialization.
   /**
    *  @param[in] pp the input database.
    *  @param[in] larmor_number the Larmor number for the operator.
    */
   GKVlasov( ParmParse&                       pp,
             const Real                       larmor_number );

   GKVlasov( ParmParse&                       pp,
             const Real                       larmor_number,
             const std::vector<std::string>&  species_name_list );

   /// Destructor.
   /**
    */
   virtual ~GKVlasov();

   virtual void accumulateRHS( GKRHSData&                             rhs,
                               const KineticSpeciesPtrVect&           kinetic_phys,
                               const CFG::LevelData<CFG::FArrayBox>&  phi,
                               const CFG::EField&                     E_field,
                               const bool                             implicit,
                               const Real&                            time );

   /// Compute a stable time step.
   /**
    * Computes and returns an estimate of the maximum stable time step.
    */

   Real computeDtExplicitTI(const CFG::EField&           Efield_cfg,
                            const KineticSpeciesPtrVect& soln,
                            const Real&                  time );

   Real computeDtImExTI( const CFG::EField&           Efield_cfg,
                         const KineticSpeciesPtrVect& soln,
                         const Real&                  time );

   /// Compute the time scale.
   /**
    * Computes and returns an estimate of the Vlasov time scale.
    */
   Real computeTimeScale( const CFG::EField&           EField_cfg,
                          const KineticSpeciesPtrVect& soln,
                          const Real&                  time );


   static const char* pp_name;
   
   
   /// Evaluates the RHS of the semi-discrete Vlasov operator.
   /**
    * @param[out] rhs data holder for rhs (without ghost cells).
    * @param[in] soln current solution with which to evaluate RHS.
    * @param[in] time current time.
    */

   // Operator evaluation using a phase space velocity computed directly from
   // a phase-space injected electric field (Efield), as opposed to the
   // divergence-free method
   void evalRHSExplicit( KineticSpecies&                        rhs,
                         const KineticSpecies&                  soln,
                         const CFG::LevelData<CFG::FArrayBox>&  phi,
                         const CFG::EField&                     Efield,
                         const int                              velocity_option,
                         const Real                             time );

   void evalRHSImplicit( KineticSpecies&                        rhs,
                         const KineticSpecies&                  soln,
                         const CFG::LevelData<CFG::FArrayBox>&  phi,
                         const CFG::EField&                     Efield,
                         const int                              velocity_option,
                         const Real                             time );

   void evalRHS(         KineticSpecies&                        rhs,
                         const KineticSpecies&                  soln,
                         const CFG::LevelData<CFG::FArrayBox>&  phi,
                         const CFG::EField&                     Efield,
                         const int                              velocity_option,
                         const Real                             time );


   void computePhysicalFlux(  LevelData<FluxBox>&                    flux_phys,
                              const KineticSpecies&                  soln_species,
                              const CFG::LevelData<CFG::FArrayBox>&  phi,
                              const CFG::EField&                     E_field,
                              const int                              velocity_option,
                              const Real                             time,
                              const bool                             apply_axisym_corr);

   void computeIntegratedFluxNormals(LevelData<FluxBox>&                    flux_normal,
                                     const KineticSpecies&                  soln_species,
                                     const CFG::LevelData<CFG::FArrayBox>&  phi,
                                     const CFG::LevelData<CFG::FArrayBox>&  Efield_cell,
                                     const CFG::LevelData<CFG::FArrayBox>&  phi_node,
                                     const bool                             fourth_order_Efield,
                                     const int                              velocity_option,
                                     const Real                             time );
   
   
   void computeIntegratedMomentFluxNormals(CFG::LevelData<CFG::FluxBox>&          flux_normal,
                                           const KineticSpecies&                  soln_species,
                                           const CFG::LevelData<CFG::FArrayBox>&  phi,
                                           const CFG::EField&                     E_field,
                                           const int                              velocity_option,
                                           const string&                          moment_type,
                                           const Real                             time);
   
   void defineMultiPhysicsPC(std::vector<Preconditioner<ODEVector,AppCtxt>*>&,
                             std::vector<DOFList>&,
                             const KineticSpeciesPtrVect&,
                             const GlobalDOFKineticSpeciesPtrVect&,
                             const ODEVector&,
                             void*,    
                             const std::string&,
                             const std::string&,
                             bool,
                             const int);

   void updateMultiPhysicsPC(std::vector<Preconditioner<ODEVector,AppCtxt>*>&,
                             const KineticSpeciesPtrVect&,
                             const CFG::LevelData<CFG::FArrayBox>&,
                             const CFG::EField&,
                             const GlobalDOFKineticSpeciesPtrVect&,
                             const Real, 
                             const int,
                             const int,
                             const Real, 
                             const bool);

   void defineBlockPC( std::vector<Preconditioner<ODEVector,AppCtxt>*>&,
                       std::vector<DOFList>&,
                       const ODEVector&,
                       void*,    
                       const std::string&,
                       const std::string&,
                       bool,
                       const KineticSpecies&,
                       const GlobalDOFKineticSpecies&,
                       const int,
                       const int);

   void updateBlockPC( std::vector<Preconditioner<ODEVector,AppCtxt>*>&,
                       const KineticSpecies&,
                       const CFG::LevelData<CFG::FArrayBox>&  phi,
                       const CFG::EField&                     E_field,
                       const GlobalDOFKineticSpecies&,
                       const Real,
                       const int,
                       const int,
                       const Real,
                       const bool,
                       const int );
   
   bool computeCfl() const {return m_compute_cfl;}

   inline void implicitSpecies(std::vector<std::string>& a_names) const
   {
     a_names.clear();
     std::map<std::string, tiOpType>::const_iterator it;
     for (it = m_ti_op_type.begin(); it != m_ti_op_type.end(); ++it) {
       if (it->second == _implicit_op_) {
         a_names.push_back(it->first);
       }
     }
   }

   void testPC( const GKVlasovAMG*         pc,
                const PhaseGeom&           geometry,
                const LevelData<FluxBox>&  velocity,
                const Real                 shift );

   void solvePCImEx(Preconditioner<ODEVector,AppCtxt>* const,
                    KineticSpeciesPtrVect&,
                    const KineticSpeciesPtrVect&,
                    const CFG::FluidSpeciesPtrVect&,
                    int ) const;
   
   void constructStencils( const LevelData<FluxBox>&               face_vel,
                           const PhaseGeom&                        geometry,
                           LevelData<BaseFab<Vector<IntVect> > >&  couplings,
                           LevelData<BaseFab<Vector<Real> > >&     weights ) const;

  protected:

   void computeAdvectionFluxPhys( const LevelData<FArrayBox>& dist_fn,
                                  const LevelData<FluxBox>&   velocity,
                                  LevelData<FluxBox>&         flux,
                                  const PhaseGeom&            phase_geom);

   void computeAdvectionFluxMapped( const LevelData<FArrayBox>& dist_fn,
                                    const LevelData<FluxBox>&   velocity_normal,
                                    LevelData<FluxBox>&         flux_normal,
                                    const PhaseGeom&            phase_geom);
   
   void applyMappedLimiter( LevelData<FluxBox>&         facePhi,
                            const LevelData<FArrayBox>& cellPhi,
                            const LevelData<FluxBox>&   faceVel,
                            const PhaseGeom&            geom,
                            const bool                  mapped_velocity );

   Real computeMappedDtSpecies(const LevelData<FluxBox>& faceVel,
                               const PhaseGeom&          geom,
                               Real                      cfl);

   Real computeMappedDtSpeciesFromNormals(const LevelData<FluxBox>& velocity_normal,
                                          const PhaseGeom&          geom,
                                          Real                      cfl);

   Real computeMappedTimeScaleSpecies(const LevelData<FluxBox>& faceVel,
                                      const PhaseGeom&          geom);

   Real computeMappedTimeScaleSpeciesFromNormals(const LevelData<FluxBox>& velocity_normal,
                                                 const PhaseGeom&          geom);

   //Computes delta_dfn (does not contain either B_star or J)
   void computeDeltaF(const KineticSpecies&         a_soln_species,
                      const LevelData<FArrayBox>&   a_dfn_no_bstar,
                      LevelData<FArrayBox>&         a_delta_dfn,
                      LevelData<FArrayBox>&         a_maxwellian_dfn);
   
   double globalMax(const double data) const;
   
   void testZeroDivergence(const LevelData<FluxBox>&  velocity,
                           const PhaseGeom&           geometry) const;

   void testEnergyConservation(const LevelData<FArrayBox>&            rhs,
                               const KineticSpecies&                  soln_species,
                               const CFG::LevelData<CFG::FArrayBox>&  phi,
                               const CFG::EField&                     E_field);
  

   /// Initializes the kinetic species data.
   /**
    * Working through the vector, initializes each KineticSpecies with
    * the associated initial data.
    *
    * @param[out] soln solution to initialize.
    * @param[in] time time at which to initialize.
    */
   void initialize( KineticSpeciesPtrVect& soln,
                    const Real      time );

   inline void computeFaceAverages( LevelData<FluxBox>&         face_data,
                                    const LevelData<FArrayBox>& cell_data,
                                    const bool                  second_order = false ) const;

   
   void getEFieldwFLR(  LevelData<FluxBox>&,
                        const KineticSpecies&,
                        const CFG::LevelData<CFG::FArrayBox>&,
                        const CFG::EField& );

   Real m_larmor_number;

   typedef enum {INVALID=-1, PPM, C2, UW1, UW3, UW5, WENO5, BWENO, NUM_FLUX } FluxType;
   FluxType m_face_avg_type;
   std::string m_method;

   static Real s_stability_bound[NUM_FLUX];

   // where we save the advective dt computed during the RHS computation.
   Real m_dt_dim_factor;
   
   bool m_verbose;
   bool m_time_step_diagnostics;
   bool m_subtract_maxwellian;
   bool m_update_maxwellian;
   bool m_compute_cfl;
   bool m_first_call;
   
   LevelData<FArrayBox> m_F0;
   LevelData<FArrayBox> m_maxwellian_dfn;
   LevelData<FArrayBox> m_delta_dfn;
   LevelData<FArrayBox> m_dfn_no_bstar;

   LevelData<FluxBox> m_flux;
   LevelData<FluxBox> m_flux_normal;
   LevelData<FluxBox> m_velocity;
   LevelData<FluxBox> m_velocity_normal;
   LevelData<FluxBox> m_delta_flux;
   LevelData<FluxBox> m_delta_flux_normal;
   
   LevelData<FluxBox> m_face_dist;
   RefCountedPtr<KineticFunction> m_ref_func;

   std::map<std::string, tiOpType> m_ti_op_type;
   std::map<std::string, int> m_my_pc_idx;

   bool m_precond_build_test;
   string m_precond_face_avg_type;

   bool m_include_FLR_effects;
   RefCountedPtr<CFG::GridFunction> m_FLR_switch_func;
   LevelData<FluxBox> m_FLR_bc_factor;
   int m_FLR_bc_buffer;
   bool m_plot_FLR_switch;

   int m_update_pc_freq;
   bool m_update_pc_skip_stage;

	bool m_useSG; // Whether to use sparse grid-friendly stencils when computing fluxes, etc.
};

#include "NamespaceFooter.H"

#endif
