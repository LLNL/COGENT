#ifndef _GKSYSPARAM_H_
#define _GKSYSPARAM_H_

#include <string>
#include <vector>

#include "Dimensions.H"
#include "Directions.H"
#include "REAL.H"

#include "NamespaceHeader.H"

class GKSystemParameters
{
  public:

    GKSystemParameters();

    void readParams(const std::string& a_parsekey = "gksystem");

    inline const std::string& parseKey() const { return m_parsekey; }

    inline int kineticGhosts() const { return m_kinetic_ghosts; }
    inline int fluidGhosts() const { return m_fluid_ghosts; }
    inline const std::string& magGeomType() const { return m_mag_geom_type; }
    inline int verbosity() const { return m_verbosity; }

    inline const std::vector<int>&  numCells() const { return m_num_cells; }
    inline const std::vector<int>&  numVelocityCells() const { return m_num_velocity_cells; }
    inline const std::vector<bool>& isPeriodic() const { return m_is_periodic; }
    inline const std::vector<int>&  cfgDecomp() const { return m_configuration_decomposition; }
    inline const std::vector<int>&  velDecomp() const { return m_velocity_decomposition; }
    inline const std::vector<int>&  phaseDecomp() const { return m_phase_decomposition; }

    inline const std::string& tiClass() const { return m_ti_class; }
    inline const std::string& tiMethod() const { return m_ti_method; }

    inline const std::vector<std::string>& diagnosticsDfnPlots() const
    {
      return m_diagnostics_selected_dfn_plots;
    }
    inline const std::vector<std::string>& diagnosticsCfgFieldVars() const
    {
      return m_diagnostics_selected_cfg_field_vars;
    }
    inline const std::vector<std::string>& diagnosticsCfgVarsKinSpec() const
    {
      return m_diagnostics_selected_cfg_vars_kin_spec;
    }
    inline const std::vector<std::string>& diagnosticsCfgVarsTotalKinSpec() const
    {
      return m_diagnostics_selected_cfg_vars_total_kin_spec;
    }
    inline bool diagnosticsFluids() const { return m_diagnostics_fluids; }

    inline const std::vector<int> fixedPlotIndices() const { return m_fixed_plotindices; }
    inline void fixedPlotIndices(const std::vector<int>& a_indices) { m_fixed_plotindices = a_indices; }

    inline bool enforceStagePositivity() const { return m_enforce_stage_positivity; }
    inline bool enforceStepPositivity() const { return m_enforce_step_positivity; }

    inline int positivityNIter() const { return m_positivity_n_iter; }
    inline bool positivityVerbose() const { return m_positivity_verbose; }
 
    inline bool enforceStepFloor() const { return m_enforce_step_floor; }
    inline Real floorValue() const { return m_floor_value; }
    inline bool absoluteFloor() const { return m_absolute_floor; }
 
    inline bool oldVorticityModel() const { return m_old_vorticity_model; }
    inline bool fixedEField() const { return m_fixed_efield; }
    inline bool zeroEField() const { return m_zero_efield; }

    inline bool imexTransportImplicit() const { return m_transport_imex_implicit; }
    inline bool imexNeutralsImplicit() const { return m_neutrals_imex_implicit; }

    inline bool skipEFieldStageUpdate() const { return m_skip_efield_stage_update; }

    inline bool consistentPotentialBCs() const { return m_consistent_potential_bcs; }

    inline bool useVlasovAmpereBCs() const { return m_use_vlasov_ampere_bcs; }
    inline bool ESolExtrapolation() const { return m_Esol_extrapolation; }
    inline bool dealignmentCorrections() const { return m_dealignment_corrections; }

    inline bool ampereLaw() const { return m_ampere_law; }
    inline bool ampereColdElectrons() const { return m_ampere_cold_electrons; }
    inline bool amperePostStepUpdate() const { return m_ampere_post_step_update; }

    inline bool transportModelOn() const { return m_transport_model_on; }
    inline bool neutralsModelOn() const { return m_neutrals_model_on; }

    inline bool enforceQuasineutrality() const { return m_enforce_quasineutrality; }

    inline bool stepConstKineticCoeffFluidOp() const { return m_step_const_kin_coeff_fluidop; }

    inline void numCells(const std::vector<int>& a_num_cells)
    {
      CH_assert(a_num_cells.size() == PDIM);
      m_num_cells = a_num_cells;
    }

    inline void cfgDecomp( const std::vector<int>& a_decomp )
    {
      CH_assert(a_decomp.size() == CFG_DIM);
      m_configuration_decomposition = a_decomp;
    }

    inline void velDecomp( const std::vector<int>& a_decomp )
    {
      CH_assert(a_decomp.size() == VEL_DIM);
      m_velocity_decomposition = a_decomp;
    }

    inline void phaseDecomp( const std::vector<int>& a_decomp )
    {
      CH_assert(a_decomp.size() == PDIM);
      m_phase_decomposition = a_decomp;
    }

  protected:

    std::string m_parsekey;

    int m_kinetic_ghosts;
    int m_fluid_ghosts;
    std::string m_mag_geom_type;
    int m_verbosity;

    std::vector<int>  m_num_cells;
    std::vector<int>  m_num_velocity_cells;
    std::vector<bool> m_is_periodic;
    std::vector<int>  m_configuration_decomposition;
    std::vector<int>  m_velocity_decomposition;
    std::vector<int>  m_phase_decomposition;

    std::string m_ti_class;
    std::string m_ti_method;

    std::vector<std::string> m_diagnostics_all_dfn_plots;
    std::vector<std::string> m_diagnostics_selected_dfn_plots;

    std::vector<std::string> m_diagnostics_all_cfg_field_vars;
    std::vector<std::string> m_diagnostics_selected_cfg_field_vars;

    std::vector<std::string> m_diagnostics_all_cfg_vars_kin_spec;
    std::vector<std::string> m_diagnostics_selected_cfg_vars_kin_spec;

    std::vector<std::string> m_diagnostics_all_cfg_vars_total_kin_spec;
    std::vector<std::string> m_diagnostics_selected_cfg_vars_total_kin_spec;

    bool m_diagnostics_fluids;

    std::vector<int> m_fixed_plotindices;

    bool m_enforce_stage_positivity;
    bool m_enforce_step_positivity;

    int m_positivity_n_iter;
    bool m_positivity_verbose;

    bool m_enforce_step_floor;
    Real m_floor_value;
    bool m_absolute_floor;

    bool m_old_vorticity_model;
    bool m_fixed_efield;
    bool m_zero_efield;

    bool m_transport_imex_implicit;
    bool m_neutrals_imex_implicit;

    bool m_skip_efield_stage_update;
    bool m_consistent_potential_bcs;
    bool m_use_vlasov_ampere_bcs;
    bool m_Esol_extrapolation;
    bool m_dealignment_corrections;

    bool m_ampere_law;
    bool m_ampere_cold_electrons;
    bool m_ampere_post_step_update;

    bool m_transport_model_on;
    bool m_neutrals_model_on;

    bool m_enforce_quasineutrality;
    bool m_step_const_kin_coeff_fluidop;
};

#include "NamespaceFooter.H"
#endif
