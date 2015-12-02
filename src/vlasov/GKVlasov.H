#ifndef _GKVLASOV_H_
#define _GKVLASOV_H_

#include "KineticSpecies.H"

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
   virtual ~GKVlasov();

   /// Evaluates the RHS of the semi-discrete Vlasov operator.
   /**
    * @param[out] rhs data holder for rhs (without ghost cells).
    * @param[in] soln current solution with which to evaluate RHS.
    * @param[in] time current time.
    */
   virtual void evalRHS( KineticSpecies&           rhs,
                         const KineticSpecies&     soln,
                         const LevelData<FluxBox>& Efield,
                         const Real                time );

   virtual void evalRHS( KineticSpecies&                 rhs_species,
                         double&                         lo_value,
                         double&                         hi_value,
                         CFG::LevelData<CFG::FArrayBox>& radial_flux_divergence_average,
                         const KineticSpecies&           soln_species,
                         const LevelData<FluxBox>&       Efield,
                         const Real                      time );

   /// Compute a stable time step.
   /**
    * Computes and returns an estimate of the maximum stable time step.
    */
   Real computeDt( const LevelData<FluxBox>& Efield,
                   const KineticSpeciesPtrVect& soln );

   Real computeMappedDtSpecies(const LevelData<FluxBox>& faceVel,
                               const PhaseGeom&          geom,
                               Real                      cfl);

   void applyMappedLimiter( LevelData<FluxBox>&         facePhi,
                            const LevelData<FArrayBox>& cellPhi,
                            const LevelData<FluxBox>&   faceVel,
                            const PhaseGeom&            geom);

  static const char* pp_name;

  void computeFlux( const LevelData<FArrayBox>& dist_fn,
                    const LevelData<FluxBox>&   velocity,
                    LevelData<FluxBox>&         flux,
                    const PhaseGeom&            phase_geom);


  /// Computes FS_averaged particle flux normalized by the shell_volume (i.e., <Flux*grad(Psi)>/Shell_volume).
  /**
   * Computes <grad(psi) * Flux>/(one-cell wide shell volume) 
   * radial_flux_divergence_average does not have ghost cells, and stores the data in the radial direction  as follows
   * radial_flux_divergence_average[i] = <grad(psi) *Flux[i]>/shell_volume[i]
   * lo_value = radial_flux_divergence_average[0] = <grad(psi) * Flux[0]>/shell_volume[0]
   * hi_value = <grad(psi) * Flux[nrad+1]>/shell_volume[nrad]
   * n_rad number of the core cells in the radial direction 
   * <...> denotes the shell-average
   */

  void computeRadialFluxDivergence(const PhaseGeom&                geometry,
                                   LevelData<FluxBox>&             flux,
                                   double                          mass,
                                   double                          charge,
                                   double&                         lo_value,
                                   double&                         hi_value,
                                   CFG::LevelData<CFG::FArrayBox>& radial_flux_divergence_average) const;

private:

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
   Real m_saved_dt;
   Real m_dt_dim_factor;
   
   bool m_verbose;
   bool m_time_step_diagnostics;
};

#include "NamespaceFooter.H"

#endif
