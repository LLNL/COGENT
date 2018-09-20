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
   GKVlasovAmpere( ParmParse&  pp,
                   const Real  larmor_number,
                   const bool  self_consistent_bcs_only = false );

   /// Destructor.
   /**
    */
   virtual ~GKVlasovAmpere() {};

   virtual void accumulateRHS( GKRHSData&                    rhs,
                               const KineticSpeciesPtrVect&  kinetic_phys,
                               const CFG::EField&            E_field,
                               const Real&                   time );

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
   void evalRHS( KineticSpecies&                  rhs_species,
                 double&                          lo_value,
                 double&                          hi_value,
                 CFG::LevelData<CFG::FArrayBox>&  radial_flux_divergence_average,
                 const KineticSpecies&            soln_species,
                 const CFG::EField&               E_field,
                 const Real                       time );

   // Operator evaluation using a divergence-free phase space velocity computed from
   // the cell-centered field (Efield_cell) and nodal potential (phi_node)
   void evalRHS( KineticSpecies&                        rhs_species,
                 double&                                lo_value,
                 double&                                hi_value,
                 CFG::LevelData<CFG::FArrayBox>&        radial_flux_divergence_average,
                 const KineticSpecies&                  soln_species,
                 const CFG::LevelData<CFG::FArrayBox>&  Efield_cell,
                 const CFG::LevelData<CFG::FArrayBox>&  phi_node,
                 const bool                             fourth_order_Efield,
                 const int                              velocity_option,
                 const Real                             time );

protected:

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

   void computeRadialNormalFluxDivergence(const PhaseGeom&                geometry,
                                          LevelData<FluxBox>&             flux,
                                          double                          mass,
                                          double                          charge,
                                          double&                         lo_value,
                                          double&                         hi_value,
                                          CFG::LevelData<CFG::FArrayBox>& radial_flux_divergence_average) const;

   bool m_self_consistent_bcs_only;

   CFG::LevelData<CFG::FArrayBox> m_total_radial_flux_divergence_average;
   double m_lo_radial_flux_divergence_average;
   double m_hi_radial_flux_divergence_average;
};

#include "NamespaceFooter.H"

#endif
