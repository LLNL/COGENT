#ifndef _GKVLASOVAMPERE_H_
#define _GKVLASOVAMPERE_H_

#include "GKVlasov.H"

#include "NamespaceHeader.H"

/**
 * Gyrokinetic Vlasov operator class for the Ampere option.
*/
class GKVlasovAmpere
   : public GKVlasov
{
public:

   /// Constructor with initialization.
   /**
    *  @param[in] pp the input database.
    *  @param[in] larmor_number the Larmor number for the operator.
    */
   GKVlasovAmpere( ParmParse&                       pp,
                   const Real                       larmor_number,
                   const bool                       self_consistent_bcs_only = false );

   GKVlasovAmpere( ParmParse&                       pp,
                   const Real                       larmor_number,
                   const std::vector<std::string>&  species_name_list,
                   const bool                       self_consistent_bcs_only = false );

   /// Destructor.
   /**
    */
   virtual ~GKVlasovAmpere() {};

   virtual void accumulateRHS( GKRHSData&                             rhs,
                               const KineticSpeciesPtrVect&           kinetic_phys,
                               const CFG::LevelData<CFG::FArrayBox>&  phi,
                               const CFG::EField&                     E_field,
                               const bool                             implicit,
                               const Real&                            time );

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

   void evalRHSExplicit( KineticSpecies&                        rhs_species,
                         double&                                lo_value,
                         double&                                hi_value,
                         CFG::LevelData<CFG::FArrayBox>&        radial_flux_divergence_average,
                         const KineticSpecies&                  soln_species,
                         const CFG::LevelData<CFG::FArrayBox>&  phi,
                         const CFG::EField&                     E_field,
                         const int                              velocity_option,
                         const Real                             time );

   void evalRHSImplicit( KineticSpecies&                        rhs_species,
                         double&                                lo_value,
                         double&                                hi_value,
                         CFG::LevelData<CFG::FArrayBox>&        radial_flux_divergence_average,
                         const KineticSpecies&                  soln_species,
                         const CFG::LevelData<CFG::FArrayBox>&  phi,
                         const CFG::EField&                     E_field,
                         const int                              velocity_option,
                         const Real                             time );

   void evalRHS(         KineticSpecies&                        rhs_species,
                         double&                                lo_value,
                         double&                                hi_value,
                         CFG::LevelData<CFG::FArrayBox>&        radial_flux_divergence_average,
                         const KineticSpecies&                  soln_species,
                         const CFG::LevelData<CFG::FArrayBox>&  phi,
                         const CFG::EField&                     E_field,
                         const int                              velocity_option,
                         const Real                             time );


protected:
   
   /// Computes FS_averaged particle flux normalized by the shell_volume (i.e., <Flux*grad(Psi)>/Shell_volume).
   /**
    * Computes <grad(psi) * Flux> * S/(one-cell wide shell volume) 
    * radial_flux_divergence_average does not have ghost cells, and stores the data in the radial direction  as follows
    * radial_flux_divergence_average[i] = <grad(psi) *Flux[i]> * S[i]/shell_volume[i]
    * lo_value = radial_flux_divergence_average[0] = <grad(psi) * Flux[0]> * S[0]/shell_volume[0]
    * hi_value = <grad(psi) * Flux[nrad+1]> * S[nrad]/shell_volume[nrad]
    * n_rad number of the core cells in the radial direction 
    * <...> denotes the shell-average, S is the surface area of a flux surface
    * <grad(psi) * Flux> * S is the total flux through the flux-surface
    */

   void computeRadialFluxDivergence(const RefCountedPtr<PhaseGeom>&  geometry,
                                    LevelData<FluxBox>&              flux,
                                    double                           mass,
                                    double                           charge,
                                    double&                          lo_value,
                                    double&                          hi_value,
                                    CFG::LevelData<CFG::FArrayBox>&  radial_flux_divergence_average) const;

   void computeRadialFluxDivergence(const RefCountedPtr<PhaseGeom>&  geometry,
                                    const LevelData<FluxBox>&        flux,
                                    const KineticSpecies&            species,
                                    double&                          lo_value,
                                    double&                          hi_value) const;

   
   void computeRadialNormalFluxDivergence(const RefCountedPtr<PhaseGeom>&  geometry,
                                          LevelData<FluxBox>&              flux,
                                          double                           mass,
                                          double                           charge,
                                          double&                          lo_value,
                                          double&                          hi_value,
                                          CFG::LevelData<CFG::FArrayBox>&  radial_flux_divergence_average) const;

   void computeRadialNormalFluxDivergence(const RefCountedPtr<PhaseGeom>&  geometry,
                                          const LevelData<FluxBox>&        flux,
                                          const KineticSpecies&            species,
                                          double&                          lo_value,
                                          double&                          hi_value) const;
   bool m_self_consistent_bcs_only;

   CFG::LevelData<CFG::FArrayBox> m_total_radial_flux_divergence_average;
   double m_lo_radial_flux_divergence_average;
   double m_hi_radial_flux_divergence_average;
   
   mutable LevelData<FArrayBox> m_volume;
   mutable LevelData<FluxBox> m_bdry_flux;
   mutable LevelData<FluxBox> m_bdry_flux_normal;
   mutable LevelData<FArrayBox> m_phase_divergence_bdry;
};

#include "NamespaceFooter.H"

#endif
