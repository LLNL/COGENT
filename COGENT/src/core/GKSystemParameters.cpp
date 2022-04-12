#include "ParmParse.H"
#include "GKSystemParameters.H"

#include "NamespaceHeader.H"

GKSystemParameters::GKSystemParameters()
{
  m_parsekey = "invalid";

  m_kinetic_ghosts = 4;
  m_fluid_ghosts = 4;
  m_mag_geom_type = "none";
  m_verbosity = 0;

  m_num_cells.clear();
  m_num_velocity_cells.clear();
  m_is_periodic.clear();
  m_configuration_decomposition.clear();
  m_velocity_decomposition.clear();
  m_phase_decomposition.clear();

  m_ti_class = "rk";
  m_ti_method = "4";

  m_diagnostics_all_dfn_plots.clear();
  m_diagnostics_all_dfn_plots.push_back("dfn");
  m_diagnostics_all_dfn_plots.push_back("deltaF");
  m_diagnostics_all_dfn_plots.push_back("dfn_at_mu");
  m_diagnostics_all_dfn_plots.push_back("vpartheta");
  m_diagnostics_all_dfn_plots.push_back("bfvpartheta");
  m_diagnostics_all_dfn_plots.push_back("rtheta");
  m_diagnostics_all_dfn_plots.push_back("vparmu");

  m_diagnostics_all_cfg_field_vars.clear();
  m_diagnostics_all_cfg_field_vars.push_back("potential");
  m_diagnostics_all_cfg_field_vars.push_back("potential_non_zonal");
  m_diagnostics_all_cfg_field_vars.push_back("efield");
  m_diagnostics_all_cfg_field_vars.push_back("ExBdata");

  m_diagnostics_all_cfg_vars_kin_spec.clear();
  m_diagnostics_all_cfg_vars_kin_spec.push_back("density");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("momentum");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("ParallelMomentum");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("PoloidalMomentum");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("ParallelVelocity");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("energyDensity");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("kineticEnergyDensity");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("parallelEnergyDensity");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("perpEnergyDensity");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("pressure");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("parallelPressure");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("perpPressure");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("gradPoverN");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("parallelHeatFlux");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("totalParallelHeatFlux");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("temperature");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("parallelTemperature");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("perpTemperature");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("fourthMoment");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("ParticleFlux");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("HeatFlux");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("ExBiHeatFlux");
  m_diagnostics_all_cfg_vars_kin_spec.push_back("ExBKineticEnergyFlux");

  m_diagnostics_all_cfg_vars_total_kin_spec.clear();
  m_diagnostics_all_cfg_vars_total_kin_spec.push_back("total_density");

  m_diagnostics_selected_dfn_plots.clear();
  m_diagnostics_selected_cfg_field_vars.clear();
  m_diagnostics_selected_cfg_vars_kin_spec.clear();
  m_diagnostics_selected_cfg_vars_total_kin_spec.clear();
  m_diagnostics_fluids = false;

  m_fixed_plotindices.clear();

  m_enforce_stage_positivity = false;
  m_enforce_step_positivity = false;

  m_positivity_n_iter = 5;
  m_positivity_verbose = false;

  m_enforce_step_floor = false;
  m_floor_value = 0.0;
  m_absolute_floor = false;

  m_old_vorticity_model = false;
  m_fixed_efield = false;
  m_no_efield = false;

  m_transport_imex_implicit = false;
  m_neutrals_imex_implicit = false;

  m_skip_efield_stage_update = false;
  m_consistent_potential_bcs = false;
  m_use_vlasov_ampere_bcs = true;
  m_Esol_extrapolation = true;
  m_dealignment_corrections = false;

  m_ampere_law = false;
  m_ampere_cold_electrons = false;
  m_ampere_post_step_update = false;

  m_transport_model_on = false;
  m_neutrals_model_on = false;

  m_enforce_quasineutrality = false;
  m_step_const_kin_coeff_fluidop = false;
}

void GKSystemParameters::readParams(const std::string&  a_parsekey)
{
  ParmParse ppgksys(a_parsekey.c_str());

  ppgksys.query( "kinetic_ghost_width", m_kinetic_ghosts );
  ppgksys.query( "fluid_ghost_width", m_fluid_ghosts );

  ppgksys.query( "verbosity", m_verbosity );
  CH_assert( m_verbosity >= 0 );

  if (ppgksys.contains("magnetic_geometry_mapping")) {
    m_parsekey = a_parsekey;
    ppgksys.get( "magnetic_geometry_mapping", m_mag_geom_type );
  }

  if (m_mag_geom_type != "none") {

    if ( m_mag_geom_type != "SingleNull" && m_mag_geom_type != "SNCore") {

      if (ppgksys.contains("num_cells")) {
        m_num_cells.resize( PDIM );
        for (int i=0; i<PDIM; ++i) m_num_cells[i] = 0;
        ppgksys.getarr( "num_cells", m_num_cells, 0, PDIM );
        for (int i=0; i<PDIM; ++i) CH_assert( m_num_cells[i]>0 );
      }
  
      if (ppgksys.contains("is_periodic")) {
        m_is_periodic.resize(PDIM);
        std::vector<int> isPeriodic( PDIM ); 
        ppgksys.getarr( "is_periodic", isPeriodic, 0, PDIM );
        for (int dim=0; dim<SpaceDim; dim++)  {
          m_is_periodic[dim] = (isPeriodic[dim] == 1);
        }
      }
      
      if (ppgksys.contains("configuration_decomp")) {
        m_configuration_decomposition.resize( CFG_DIM );
        for (int i=0; i<CFG_DIM; ++i) m_configuration_decomposition[i] = 0;
        ppgksys.getarr( "configuration_decomp", 
                          m_configuration_decomposition, 
                          0, 
                          CFG_DIM );
        for (int i=0; i<CFG_DIM; ++i) {
          CH_assert( m_configuration_decomposition[i]>0 );
        }
      }
      
      if (ppgksys.contains("phase_decomp")) {
        m_phase_decomposition.resize( PDIM );
        for (int i=0; i<PDIM; ++i) m_phase_decomposition[i] = 0;
        ppgksys.getarr( "phase_decomp", m_phase_decomposition, 0, PDIM );
        for (int i=0; i<PDIM; ++i) CH_assert( m_phase_decomposition[i]>0 );
      }
    
    } else {
      
      if (ppgksys.contains("num_velocity_cells")) {
        m_num_velocity_cells.resize( VEL_DIM );
        for (int i=0; i<VEL_DIM; ++i) m_num_velocity_cells[i] = 0;
        ppgksys.getarr( "num_velocity_cells", 
                        m_num_velocity_cells, 
                        0, 
                        VEL_DIM );
        for (int i=0; i<VEL_DIM; ++i) CH_assert( m_num_velocity_cells[i]>0 );
      }
      
    }
  }
  
  if (ppgksys.contains("velocity_decomp")) {
    m_velocity_decomposition.resize( VEL_DIM );
    for (int i=0; i<VEL_DIM; ++i) m_velocity_decomposition[i] = 0;
    ppgksys.getarr( "velocity_decomp", 
                      m_velocity_decomposition, 
                      0, 
                      VEL_DIM );
    for (int i=0; i<VEL_DIM; ++i) CH_assert( m_velocity_decomposition[i]>0 );
  }


  ppgksys.query("ti_class",m_ti_class);
  ppgksys.query("ti_method",m_ti_method);

  for (int i=0; i<m_diagnostics_all_dfn_plots.size(); i++) {
    std::string var = m_diagnostics_all_dfn_plots[i];
    std::string key = "hdf_" + var;
    bool flag = false;
    ppgksys.query(key.c_str(), flag);
    if (flag) m_diagnostics_selected_dfn_plots.push_back(var);
  }
  for (int i=0; i<m_diagnostics_all_cfg_field_vars.size(); i++) {
    std::string var = m_diagnostics_all_cfg_field_vars[i];
    std::string key = "hdf_" + var;
    bool flag = false;
    ppgksys.query(key.c_str(), flag);
    if (flag) m_diagnostics_selected_cfg_field_vars.push_back(var);
  }
  for (int i=0; i<m_diagnostics_all_cfg_vars_kin_spec.size(); i++) {
    std::string var = m_diagnostics_all_cfg_vars_kin_spec[i];
    std::string key = "hdf_" + var;
    bool flag = false;
    ppgksys.query(key.c_str(), flag);
    if (flag) m_diagnostics_selected_cfg_vars_kin_spec.push_back(var);
  }
  for (int i=0; i<m_diagnostics_all_cfg_vars_total_kin_spec.size(); i++) {
    std::string var = m_diagnostics_all_cfg_vars_total_kin_spec[i];
    std::string key = "hdf_" + var;
    bool flag = false;
    ppgksys.query(key.c_str(), flag);
    if (flag) m_diagnostics_selected_cfg_vars_total_kin_spec.push_back(var);
  }
  ppgksys.query("hdf_fluids",m_diagnostics_fluids);
  
  m_fixed_plotindices.resize( 5 );
  for (int i=0; i<5; ++i) m_fixed_plotindices[i]=0;
  ppgksys.queryarr("fixed_plot_indices",m_fixed_plotindices,0,5);
  
  bool enforce_positivity(false);
  if (ppgksys.contains("enforce_positivity")) {
    ppgksys.get("enforce_positivity", enforce_positivity);
  }
  
  if (enforce_positivity) {

    std::string ptype("stage");
    m_enforce_stage_positivity = true;
    if (ppgksys.contains("enforce_positivity_type")) {
      ppgksys.get("enforce_positivity_type", ptype);
      if (ptype=="step") {
        m_enforce_step_positivity = true;
        m_enforce_stage_positivity = false;
      } else if (ptype!="stage") {
        MayDay::Error("Invalid positivity enforcement type");
      }
    }
    
    m_positivity_n_iter = 5;
    if (ppgksys.contains("max_positivity_iter")) {
      ppgksys.get("max_positivity_iter", m_positivity_n_iter);
    }
    
    m_positivity_verbose = false;
    if (ppgksys.contains("positivity_verbose_output")) {
      ppgksys.get("positivity_verbose_output", m_positivity_verbose);
    }

  }
  
  bool enforce_floor(false);
  if (ppgksys.contains("enforce_floor")) {
    ppgksys.get("enforce_floor", enforce_floor);
  }
  
  if (enforce_floor) {
    
    m_enforce_step_floor = true;
    
    m_floor_value = 0.;
    ppgksys.query( "floor_value", m_floor_value );
    
    m_absolute_floor = true;
    ppgksys.query( "absolute_floor", m_absolute_floor );
    
  }
  
  if (ppgksys.contains("old_vorticity_model")) {
    ppgksys.get("old_vorticity_model", m_old_vorticity_model);
  } else {
    m_old_vorticity_model = false;
  }

  ppgksys.query( "zero_efield", m_no_efield );  // zero_efield is deprecated
  ppgksys.query( "no_efield", m_no_efield );
  if ( m_no_efield ) {
     m_fixed_efield = true;
  }
  else {
     ppgksys.query( "fixed_efield", m_fixed_efield ); 
  }

  ppgksys.query( "imex_transport_implicit", m_transport_imex_implicit );
  ppgksys.query( "imex_neutrals_implicit", m_neutrals_imex_implicit );

  ppgksys.query("skip_efield_stage_update", m_skip_efield_stage_update);
  ppgksys.query( "consistent_potential_bcs", m_consistent_potential_bcs );
  ppgksys.query( "use_vlasov_ampere_bcs", m_use_vlasov_ampere_bcs );
  ppgksys.query( "extrapolated_sol_efield", m_Esol_extrapolation );
  ppgksys.query( "efield_dealignment_corrections", m_dealignment_corrections );

  ppgksys.query( "ampere_law", m_ampere_law );
  if (m_ampere_law) {
    m_consistent_potential_bcs = true;
    ppgksys.query( "ampere_cold_electrons", m_ampere_cold_electrons );
  }

  ppgksys.query( "ampere_post_step_update", m_ampere_post_step_update );

  ppgksys.query("transport_model_on", m_transport_model_on);
  ppgksys.query("neutrals_model_on", m_neutrals_model_on);

  ppgksys.query( "gksystem.enforce_quasineutrality", m_enforce_quasineutrality );
  ppgksys.query("step_const_coef", m_step_const_kin_coeff_fluidop);

  return;
}

#include "NamespaceFooter.H"

