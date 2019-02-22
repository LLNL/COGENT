#ifndef _GKVLASOV_H_
#define _GKVLASOV_H_

#include "KineticFunction.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "EField.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

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
   GKVlasov( ParmParse&  pp,
             const Real  larmor_number);

   /// Destructor.
   /**
    */
   virtual ~GKVlasov() {};

   virtual void accumulateRHS( GKRHSData&                             rhs,
                               const KineticSpeciesPtrVect&           kinetic_phys,
                               const CFG::LevelData<CFG::FArrayBox>&  phi,
                               const CFG::EField&                     E_field,
                               const Real&                            time );

   /// Compute a stable time step.
   /**
    * Computes and returns an estimate of the maximum stable time step.
    */

   Real computeDt( const CFG::EField&           Efield_cfg,
                   const KineticSpeciesPtrVect& soln );

   /// Compute the time scale.
   /**
    * Computes and returns an estimate of the Vlasov time scale.
    */
   Real computeTimeScale( const CFG::EField&           EField_cfg,
                          const KineticSpeciesPtrVect& soln );

   void computeFluxNormal( const LevelData<FArrayBox>& dist_fn,
                           const LevelData<FluxBox>&   velocity_normal,
                           LevelData<FluxBox>&         flux_normal,
                           const PhaseGeom&            phase_geom);

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
   void evalRHS( KineticSpecies&                        rhs,
                 const KineticSpecies&                  soln,
                 const CFG::LevelData<CFG::FArrayBox>&  phi,
                 const CFG::EField&                     Efield,
                 const Real                             time );

   // Operator evaluation using a divergence-free phase space velocity computed from
   // the cell-centered field (Efield_cell) and nodal potential (phi_node)
   void evalRHS( KineticSpecies&                        rhs,
                 const KineticSpecies&                  soln,
                 const CFG::LevelData<CFG::FArrayBox>&  phi,
                 const CFG::LevelData<CFG::FArrayBox>&  Efield_cell,
                 const CFG::LevelData<CFG::FArrayBox>&  phi_node,
                 const bool                             fourth_order_Efield,
                 const int                              velocity_option,
                 const Real                             time );

   void evalRHS( KineticSpecies&                        rhs,
                 const KineticSpecies&                  soln,
                 const CFG::LevelData<CFG::FArrayBox>&  Efield_cell,
                 const CFG::LevelData<CFG::FArrayBox>&  phi_node,
                 const bool                             fourth_order_Efield,
                 const int                              velocity_option,
                 const Real                             time );

protected:

   void computeFlux( const LevelData<FArrayBox>& dist_fn,
                     const LevelData<FluxBox>&   velocity,
                     LevelData<FluxBox>&         flux,
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

   Real m_larmor_number;

   typedef enum {INVALID=-1, PPM, UW1, UW3, UW5, WENO5, BWENO, NUM_FLUX } FluxType;
   FluxType m_face_avg_type;

   static Real s_stability_bound[NUM_FLUX];

   // where we save the advective dt computed during the RHS computation.
   Real m_dt_dim_factor;
   
   bool m_verbose;
   bool m_time_step_diagnostics;
   bool m_subtract_maxwellian;
   bool m_update_maxwellian;
   
   LevelData<FArrayBox> m_F0;
   RefCountedPtr<KineticFunction> m_ref_func;
};

#include "NamespaceFooter.H"

#endif
