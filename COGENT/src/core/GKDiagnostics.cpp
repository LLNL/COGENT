#include "GKDiagnostics.H"

#define VLASOV_FLUXES

#include "NamespaceHeader.H"

void GKDiagnostics::getCfgVar(  CFG::LevelData<CFG::FArrayBox>& a_var,
                                const std::string&              a_var_name,
                                const bool&                     a_opt_arg ) const
{
  if (a_var_name == "potential") {
    getPotential( a_var, a_opt_arg );
  } else if (a_var_name == "efield") {
    getEField( a_var );
  } else if (a_var_name == "ExBdata") {
    getExBData( a_var );
  } else {
    MayDay::Error("Invalid var_name");
  }
  return;
}

void GKDiagnostics::getCfgVar( CFG::LevelData<CFG::FArrayBox>&  a_var,
                               const std::string&               a_var_name,
                               const KineticSpecies&            a_species) const
{
  if (a_var_name == "density") {
    getChargeDensity( a_var, a_species );
  } else if (a_var_name == "fluidVelocity") {
    getFluidVelocity( a_var, a_species );
  } else if (a_var_name == "parallelParticleFlux") {
    getParallelParticleFlux( a_var, a_species );
  } else if (a_var_name == "poloidalParticleFlux") {
    getPoloidalParticleFlux( a_var, a_species );
  } else if (a_var_name == "parallelVelocity") {
    getParallelVelocity( a_var, a_species );
  } else if (a_var_name == "energyDensity") {
    getEnergyDensity( a_var, a_species );
  } else if (a_var_name == "kineticEnergyDensity") {
    getKineticEnergyDensity( a_var, a_species );
  } else if (a_var_name == "parallelEnergyDensity") {
    getParallelEnergyDensity( a_var, a_species );
  } else if (a_var_name == "perpEnergyDensity") {
    getPerpEnergyDensity( a_var, a_species );
  } else if (a_var_name == "pressure") {
    getPressure( a_var, a_species );
  } else if (a_var_name == "parallelPressure") {
    getParallelPressure( a_var, a_species );
  } else if (a_var_name == "perpPressure") {
    getPerpPressure( a_var, a_species );
  } else if (a_var_name == "radialForceBalance") {
    getRadialForceBalance( a_var, a_species );
  } else if (a_var_name == "temperature") {
    getTemperature( a_var, a_species );
  } else if (a_var_name == "parallelTemperature") {
    getParallelTemperature( a_var, a_species );
  } else if (a_var_name == "perpTemperature") {
    getPerpTemperature( a_var, a_species );
  } else if (a_var_name == "parallelHeatFlux") {
    getParallelHeatFlux( a_var, a_species );
  } else if (a_var_name == "totalParallelHeatFlux") {
    getTotalParallelHeatFlux( a_var, a_species );
  } else if (a_var_name == "fourthMoment") {
    getFourthMoment( a_var, a_species );
  } else if (a_var_name == "radialCurrent") {
    getRadialCurrent( a_var, a_species );
  } else if (a_var_name == "radialParticleFlux") {
    getRadialParticleFlux( a_var, a_species );
  } else if (a_var_name == "radialHeatFlux") {
    getRadialHeatFlux( a_var, a_species );
  } else if (a_var_name == "ExBHeatFlux") {
    getExBHeatFlux( a_var, a_species );
  } else if (a_var_name == "ExBKineticEnergyFlux") {
    getExBKineticEnergyFlux( a_var, a_species );
  } else {
    MayDay::Error("Invalid var_name");
  }
  return;
}

void GKDiagnostics::getCfgVar( CFG::LevelData<CFG::FArrayBox>&  a_var,
                               const std::string&               a_var_name,
                               const KineticSpeciesPtrVect&     a_species) const
{
  if (a_var_name == "total_density") {
    getChargeDensity( a_var, a_species );
  } else {
    MayDay::Error("Invalid var_name");
  }
  return;
}

void GKDiagnostics::getPhaseVar(  LevelData<FArrayBox>& a_var,
                                  const std::string&    a_var_name,
                                  const KineticSpecies& a_species) const
{
  if (a_var_name == "dfn") {
    getDistributionFunction( a_var, a_species );
  } else if (a_var_name == "bstar_par") {
    getBStarParallel( a_var, a_species );
  } else if (a_var_name == "deltaF") {
    getDeltaF( a_var, a_species );
  } else if (a_var_name == "dfn_at_mu") {
    getDistributionFunction( a_var, a_species );
  } else if (a_var_name == "vpartheta") {
    getDistributionFunction( a_var, a_species );
  } else if (a_var_name == "bfvpartheta") {
    getBStarParallelDfn( a_var, a_species );
  } else if (a_var_name == "rtheta") {
    getDistributionFunction( a_var, a_species );
  } else if (a_var_name == "vparmu") {
    getDistributionFunction( a_var, a_species );
  } else {
    MayDay::Error("Invalid var_name");
  }
  return;
}

void GKDiagnostics::getPotential( CFG::LevelData<CFG::FArrayBox>& a_phi,
                                  const bool& a_non_zonal_comp ) const
{
  const CFG::LevelData<CFG::FArrayBox>& phi(m_ops->getPhi());
  if (!a_non_zonal_comp){
    a_phi.define(phi);
    for (auto dit(phi.dataIterator()); dit.ok(); ++dit) {
      a_phi[dit].copy(phi[dit]);
    }
  } else {
    const PhaseGeom& phase_geometry( *m_phase_geometry );
    const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
    a_phi.define( mag_geom.gridsFull(), 2, CFG::IntVect::Zero );
    
    CFG::LevelData<CFG::FArrayBox> phi_axisym;
    phi_axisym.define(phi);
    
    CFG::FluxSurface flux_surface(mag_geom);
    flux_surface.averageAndSpread(phi, phi_axisym);
    for (auto dit(phi.dataIterator()); dit.ok(); ++dit) {
      a_phi[dit].copy(phi[dit],0,0,1);
      a_phi[dit].copy(phi[dit],0,1,1);
      a_phi[dit].minus(phi_axisym[dit],0,1,1);
    }
  }
  return;
}

void GKDiagnostics::getEField( CFG::LevelData<CFG::FArrayBox>& a_Efield) const
{
  
  /*
   Computes e_r, b_phi, b_pol (in 3D) or e_r, b_pol (in 2D)
   projections of Efield
   */
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
  const CFG::DisjointBoxLayout& grids = mag_geom.gridsFull();
  
#if CFG_DIM==2
  CFG::LevelData<CFG::FArrayBox> E_3D;
  m_ops->getEField(E_3D);
  a_Efield.define(grids, 2, CFG::IntVect::Zero);
  mag_geom.projectPoloidalVector(E_3D, a_Efield);
#else
  m_ops->getEField(a_Efield);
#endif
  
  // Does nothing for the slab case
  mag_geom.convertToFSProjections(a_Efield);

  return;
}


void GKDiagnostics::getExBData( CFG::LevelData<CFG::FArrayBox>& a_ExB_data ) const
{
  /*
   Compute ExB flow velocity and its shear

   First component: e_r * V_ExB  (radial ExB velocity)

   Second component: b_pol * V_ExB (poloidal ExB velocity)

   Third component: e_r * grad <b_pol * V_ExB>
                    (poloidal velocity shear)

   Fourth component: -<e_r * grad <b_pol * V_ExB> >
                    (FS-averaged poloidal velocity shear)
   
   */
  
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );

  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
  const CFG::DisjointBoxLayout& grids = mag_geom.gridsFull();

  const CFG::LevelData<CFG::FArrayBox>& phi(m_ops->getPhi());
  const CFG::PhiOps& phi_ops(m_ops->getPhiOps());
  
  CFG::LevelData<CFG::FArrayBox> Efield( grids, 3, CFG::IntVect::Unit);
  phi_ops.computeField( phi, Efield );
   
  CFG::LevelData<CFG::FArrayBox> ExB_drift( grids, 3, CFG::IntVect::Zero);
  mag_geom.computeEXBDrift(Efield, ExB_drift);
   
  const double larmor_number( m_units->larmorNumber() );
  for (CFG::DataIterator dit(grids); dit.ok(); ++dit) {
    ExB_drift[dit].mult(larmor_number);
  }
   
  a_ExB_data.define( grids, 4, CFG::IntVect::Zero);
  CFG::LevelData<CFG::FArrayBox> tmp( grids, 1, CFG::IntVect::Zero);
   
  // Compute e_r (outward radial) component of ExB velocity
  mag_geom.computeRadialProjection(tmp, ExB_drift, false);
  for (CFG::DataIterator dit(grids); dit.ok(); ++dit) {
    a_ExB_data[dit].copy(tmp[dit],0,0,1);
  }

  // Compute b_pol component of ExB velocity
  mag_geom.computePoloidalProjection(tmp, ExB_drift);

  // Compute FS shell-averaged b_pol component of ExB velocity, <(ExB)pol>
  CFG::FluxSurface flux_surface(mag_geom);
  CFG::LevelData<CFG::FArrayBox> tmp_axisym( grids, 1, CFG::IntVect::Zero);
  flux_surface.averageAndSpread(tmp, tmp_axisym);
   
  for (CFG::DataIterator dit(grids); dit.ok(); ++dit) {
    a_ExB_data[dit].copy(tmp_axisym[dit],0,1,1);
  }

  // Compute e_r (outward radial) component of grad(<(ExB)pol>)
  CFG::LevelData<CFG::FArrayBox> ExB_shear( grids, 3, CFG::IntVect::Unit);
  phi_ops.computeField(tmp_axisym, ExB_shear);
  
  mag_geom.computeRadialProjection(tmp, ExB_shear, false);
  for (CFG::DataIterator dit(grids); dit.ok(); ++dit) {
    a_ExB_data[dit].copy(tmp[dit],0,2,1);
    //Multiply by -1 to compensate for the minus sign
    // in computeField calculation that returns -grad(phi)
    a_ExB_data[dit].mult(-1.0,2,1);
  }

  // Compute flux surface (shell) average of the -<ExB> shear, -<e_r*grad(<(ExB)pol>)>
  CFG::LevelData<CFG::FArrayBox> negative_ExB_shear_axisym( grids, 1, CFG::IntVect::Zero);
  flux_surface.averageAndSpread(tmp, negative_ExB_shear_axisym);
  for (CFG::DataIterator dit(grids); dit.ok(); ++dit) {
    a_ExB_data[dit].copy(negative_ExB_shear_axisym[dit],0,3,1);
    a_ExB_data[dit].mult(-1.0,3,1);
  }

  return;
}

void GKDiagnostics::getDistributionFunction(LevelData<FArrayBox>& a_dfn,
                                            const KineticSpecies& a_soln_species ) const
{
  const PhaseGeom& species_geometry = a_soln_species.phaseSpaceGeometry();
  const LevelData<FArrayBox>& soln_dfn(a_soln_species.distributionFunction());
  a_dfn.define( species_geometry.gridsFull(), 1, IntVect::Zero );
  for (auto dit(a_dfn.dataIterator()); dit.ok(); ++dit) {
     a_dfn[dit].copy(soln_dfn[dit]);
  }

  species_geometry.divideBStarParallel( a_dfn );
  return;
}

void GKDiagnostics::getBStarParallelDfn(LevelData<FArrayBox>& a_dfn,
                                        const KineticSpecies& a_soln_species ) const
{
  const PhaseGeom& species_geometry = a_soln_species.phaseSpaceGeometry();
  const LevelData<FArrayBox>& soln_dfn(a_soln_species.distributionFunction());
  a_dfn.define( species_geometry.gridsFull(), 1, IntVect::Zero );
  for (auto dit(a_dfn.dataIterator()); dit.ok(); ++dit) {
     a_dfn[dit].copy(soln_dfn[dit]);
  }
  return;
}

void GKDiagnostics::getBStarParallel( LevelData<FArrayBox>& a_var,
                                      const KineticSpecies& a_soln_species ) const
{
  const PhaseGeom& species_geometry = a_soln_species.phaseSpaceGeometry();
  a_var.define( species_geometry.gridsFull(), 1, IntVect::Zero );
  species_geometry.getBStarParallel(a_var);
}

void GKDiagnostics::getDeltaF(LevelData<FArrayBox>& a_var,
                              const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  CFG::LevelData<CFG::FArrayBox> v_parallel( mag_geom.gridsFull(),
                                            1,
                                            CFG::IntVect::Zero );
  a_soln_species.parallelVelocity( v_parallel );
  
  CFG::LevelData<CFG::FArrayBox> density(  mag_geom.gridsFull(), 
                                           1, 
                                           CFG::IntVect::Zero );
  a_soln_species.numberDensity( density );

  CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.gridsFull(), 
                                           1, 
                                           CFG::IntVect::Zero );
  a_soln_species.pressure(pressure, v_parallel);

  for (auto dit(density.dataIterator()); dit.ok(); ++dit) {
    pressure[dit].divide(density[dit]);
  }

  const PhaseGeom& species_geometry = a_soln_species.phaseSpaceGeometry();
  a_var.define( species_geometry.gridsFull(), 1, IntVect::Zero );
  DeltaFKernel<FArrayBox> DeltaF_Kernel(density, pressure, v_parallel);
  DeltaF_Kernel.eval(a_var, a_soln_species);

  return;
}

void GKDiagnostics::getChargeDensity( CFG::LevelData<CFG::FArrayBox>& a_var,
                                      const KineticSpecies& a_soln_species ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   a_var.define( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.chargeDensity( a_var );

   return;
}

void GKDiagnostics::getFluidVelocity(CFG::LevelData<CFG::FArrayBox>& a_var,
                                     const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const double larmor_number( m_units->larmorNumber() );

  const CFG::EMFields& EM_fields = m_ops->getEMFields();
  const LevelData<FluxBox>& E_field_tmp = EM_fields.getEField();
  const LevelData<FluxBox>& Apar_derivatives_tmp = EM_fields.getAparDerivs();

  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 3, CFG::IntVect::Zero );
  a_soln_species.fluidVelocity( a_var, E_field_tmp, EM_fields, Apar_derivatives_tmp, larmor_number );

  mag_geom.convertToFSProjections(a_var);
  
  return;
}

void GKDiagnostics::getParallelParticleFlux(CFG::LevelData<CFG::FArrayBox>& a_var,
                                            const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.parallelParticleFlux(a_var);
}

void GKDiagnostics::getPoloidalParticleFlux(CFG::LevelData<CFG::FArrayBox>& a_var,
                                            const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const double larmor_number( m_units->larmorNumber() );

  const CFG::EMFields& EM_fields = m_ops->getEMFields();
  const LevelData<FluxBox>& E_field_tmp = EM_fields.getEField();
  const LevelData<FluxBox>& Apar_derivatives_tmp = EM_fields.getAparDerivs();

  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.poloidalParticleFlux( a_var, E_field_tmp, EM_fields, Apar_derivatives_tmp, larmor_number );
  
  return;
}

void GKDiagnostics::getParallelVelocity(  CFG::LevelData<CFG::FArrayBox>& a_var,
                                          const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.parallelVelocity( a_var );
}

void GKDiagnostics::getEnergyDensity( CFG::LevelData<CFG::FArrayBox>& a_var,
                                      const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 3, CFG::IntVect::Zero );
  
  CFG::LevelData<CFG::FArrayBox> energy_density_loc( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.energyDensity( energy_density_loc, m_ops->getPhi() );

  for (auto dit(a_var.dataIterator()); dit.ok(); ++dit) {
    a_var[dit].copy(energy_density_loc[dit], 0, 0, 1);
  }

  // Get flux-surface shell-averaged and shell-integrated components
  CFG::FluxSurface flux_surface(mag_geom);
  CFG::LevelData<CFG::FArrayBox> energy_density_tmp( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  
  flux_surface.averageAndSpread(energy_density_loc, energy_density_tmp);
  for (auto dit(a_var.dataIterator()); dit.ok(); ++dit) {
    a_var[dit].copy(energy_density_tmp[dit], 0, 1, 1);
  }
  
  flux_surface.integrateAndSpread(energy_density_loc, energy_density_tmp);
  for (auto dit(a_var.dataIterator()); dit.ok(); ++dit) {
    a_var[dit].copy(energy_density_tmp[dit], 0, 2, 1);
  }
}

void GKDiagnostics::getKineticEnergyDensity(CFG::LevelData<CFG::FArrayBox>& a_var,
                                            const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 3, CFG::IntVect::Zero );
  
  CFG::LevelData<CFG::FArrayBox> energy_density_loc( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.kineticEnergyDensity( energy_density_loc );

  for (auto dit(a_var.dataIterator()); dit.ok(); ++dit) {
    a_var[dit].copy(energy_density_loc[dit], 0, 0, 1);
  }

  // Get flux-surface shell-averaged and shell-integrated components
  CFG::FluxSurface flux_surface(mag_geom);
  CFG::LevelData<CFG::FArrayBox> energy_density_tmp( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  
  flux_surface.averageAndSpread(energy_density_loc, energy_density_tmp);
  for (auto dit(a_var.dataIterator()); dit.ok(); ++dit) {
    a_var[dit].copy(energy_density_tmp[dit], 0, 1, 1);
  }
  
  flux_surface.integrateAndSpread(energy_density_loc, energy_density_tmp);
  for (auto dit(a_var.dataIterator()); dit.ok(); ++dit) {
    a_var[dit].copy(energy_density_tmp[dit], 0, 2, 1);
  }
}

void GKDiagnostics::getParallelEnergyDensity( CFG::LevelData<CFG::FArrayBox>& a_var,
                                              const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 3, CFG::IntVect::Zero );
  
  CFG::LevelData<CFG::FArrayBox> energy_density_loc( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.parallelEnergyDensity( energy_density_loc );

  for (auto dit(a_var.dataIterator()); dit.ok(); ++dit) {
    a_var[dit].copy(energy_density_loc[dit], 0, 0, 1);
  }

  // Get flux-surface shell-averaged and shell-integrated components
  CFG::FluxSurface flux_surface(mag_geom);
  CFG::LevelData<CFG::FArrayBox> energy_density_tmp( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  
  flux_surface.averageAndSpread(energy_density_loc, energy_density_tmp);
  for (auto dit(a_var.dataIterator()); dit.ok(); ++dit) {
    a_var[dit].copy(energy_density_tmp[dit], 0, 1, 1);
  }
  
  flux_surface.integrateAndSpread(energy_density_loc, energy_density_tmp);
  for (auto dit(a_var.dataIterator()); dit.ok(); ++dit) {
    a_var[dit].copy(energy_density_tmp[dit], 0, 2, 1);
  }
}

void GKDiagnostics::getPerpEnergyDensity( CFG::LevelData<CFG::FArrayBox>& a_var,
                                          const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.perpEnergyDensity( a_var );
}

void GKDiagnostics::getPressure(  CFG::LevelData<CFG::FArrayBox>& a_var,
                                  const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
  a_var.define( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.pressure( a_var );
  return;
}

void GKDiagnostics::getParallelPressure(CFG::LevelData<CFG::FArrayBox>& a_var,
                                        const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
  a_var.define( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.parallelPressure( a_var );
  return;
}

void GKDiagnostics::getPerpPressure(  CFG::LevelData<CFG::FArrayBox>& a_var,
                                      const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.perpPressure( a_var );
  return;
}

void GKDiagnostics::getRadialForceBalance(CFG::LevelData<CFG::FArrayBox>& a_var,
                                          const KineticSpecies& a_soln_species ) const
{
  /*
   Returns terms of the radial force balance
   projected onto the outer normal (e_r)
   First component = - grad(P_perp)/qn * e_r
   Second component = 1/c * (Vphi x Bpol) * e_r
   Third component = 1/c * (Vpol x Bphi) * e_r
   Fourth component = - grad(P_perp)/qn * e_r + E_r
                      + 1/c * (Vphi x Bpol) * e_r
                      + 1/c * (Vpol x Bphi) * e_r
   */
  
  CH_assert( m_phase_geometry != NULL );
  const double larmor_number( m_units->larmorNumber() );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
  const CFG::DisjointBoxLayout& cfg_grids = mag_geom.gridsFull();
  
  // Diagnostics quantity
  a_var.define( cfg_grids, 4, CFG::IntVect::Zero );
  
  // Get number density
  CFG::LevelData<CFG::FArrayBox> density(cfg_grids, 1, CFG::IntVect::Zero );
  a_soln_species.numberDensity( density );

  // Get perpendicular pressure
  CFG::LevelData<CFG::FArrayBox> perp_pressure(cfg_grids, 1, CFG::IntVect::Zero );
                                              
  a_soln_species.perpPressure(perp_pressure);
  
  CFG::LevelData<CFG::FArrayBox> perp_pressure_wghosts(cfg_grids, 1, 2*CFG::IntVect::Unit );
  for (auto dit(cfg_grids.dataIterator()); dit.ok(); ++dit) {
    perp_pressure_wghosts[dit].copy(perp_pressure[dit]);
  }
  
  mag_geom.fillInternalGhosts(perp_pressure_wghosts);
  
  int order = 2;
  int nghosts = 2;
  mag_geom.extrapolateAtPhysicalBoundaries(perp_pressure_wghosts, order, nghosts);
  
  // Get mapped pressure gradient
  CFG::LevelData<CFG::FArrayBox> pressure_gradient_mapped( cfg_grids, CFG_DIM, CFG::IntVect::Unit );
  mag_geom.computeMappedGradient(perp_pressure_wghosts, pressure_gradient_mapped, order);
  
  // Get physical pressure gradient
  CFG::LevelData<CFG::FArrayBox> pressure_gradient(cfg_grids, CFG_DIM, CFG::IntVect::Unit);
  
#if CFG_DIM==2
  mag_geom.unmapPoloidalGradient(pressure_gradient_mapped, pressure_gradient);
#else
  mag_geom.unmapGradient(pressure_gradient_mapped, pressure_gradient);
#endif
  
  // Get physical pressure gradient
  CFG::LevelData<CFG::FArrayBox> pressure_gradient_r(cfg_grids, 1, CFG::IntVect::Zero);
  mag_geom.computeRadialProjection(pressure_gradient_r, pressure_gradient, false);
  
  Real charge = a_soln_species.charge();
  
  // Fill the first component: -grad(P_perp)/qn
  for (auto dit(cfg_grids.dataIterator()); dit.ok(); ++dit) {
    a_var[dit].copy(pressure_gradient_r[dit],0,0,1);
    a_var[dit].divide(density[dit], 0, 0, 1);
    a_var[dit].divide(-charge, 0, 1);
  }
  
  // Get fluid velocity (Vgc+c/qn*curlM)
  CFG::LevelData<CFG::FArrayBox> velocity(cfg_grids, 3, CFG::IntVect::Zero);
  
  const CFG::EMFields& EM_fields = m_ops->getEMFields();
  const LevelData<FluxBox>& E_field_inj = EM_fields.getEField();
  const LevelData<FluxBox>& Apar_derivatives_inj = EM_fields.getAparDerivs();

  a_soln_species.fluidVelocity( velocity, E_field_inj, EM_fields, Apar_derivatives_inj, larmor_number );

  // Compute e_r, b_pol, and b_phi projections
  mag_geom.convertToFSProjections(velocity);

  // Compute b_phi x b_pol for sign
  const CFG::LevelData<CFG::FArrayBox>& Bvec = mag_geom.getCCBField();
  
  CFG::LevelData<CFG::FArrayBox> B_cyl(cfg_grids, 3, CFG::IntVect::Zero);
  if (CFG_DIM == 2)
    for (auto dit(cfg_grids.dataIterator()); dit.ok(); ++dit) {
      B_cyl[dit].copy(Bvec[dit]);
    }
  else {
    mag_geom.convertCartesianToCylindrical(B_cyl, Bvec);
  }

  CFG::LevelData<CFG::FArrayBox> B_pol(cfg_grids, 3, CFG::IntVect::Zero);
  CFG::LevelData<CFG::FArrayBox> B_phi(cfg_grids, 3, CFG::IntVect::Zero);
  for (auto dit(cfg_grids.dataIterator()); dit.ok(); ++dit) {
    B_pol[dit].setVal(0.);
    B_pol[dit].copy(B_cyl[dit],0,0,1);
    B_pol[dit].copy(B_cyl[dit],2,2,1);
    
    B_phi[dit].setVal(0.);
    B_phi[dit].copy(B_cyl[dit],1,1,1);
  }
  
  CFG::LevelData<CFG::FArrayBox> BphixBpol(cfg_grids, 3, CFG::IntVect::Zero);
  mag_geom.crossProduct(BphixBpol, B_phi, B_pol);
  
  // Get the radial projection of BphixBpol (for subsequent sign)
  CFG::LevelData<CFG::FArrayBox> sign(cfg_grids, 3, CFG::IntVect::Zero);
  mag_geom.computeRadialProjection(sign, BphixBpol, false);
  
  // Get Bpol_mag for VphixBpol term
  CFG::LevelData<CFG::FArrayBox> Bpol_mag(cfg_grids, 1, CFG::IntVect::Zero);
  mag_geom.getBpoloidalMag(Bpol_mag);
  
  
  // Compute Er
  CFG::LevelData<CFG::FArrayBox> E_field;
  m_ops->getEField( E_field );
  
  CFG::LevelData<CFG::FArrayBox> E_r(cfg_grids, 1, CFG::IntVect::Zero);
  mag_geom.computeRadialProjection(E_r, E_field, false);
  
  // Fill the rest of raial force balance components
  for (auto dit(cfg_grids.dataIterator()); dit.ok(); ++dit) {
    CFG::BoxIterator bit(cfg_grids[dit]);
    for (bit.begin();bit.ok();++bit) {
      CFG::IntVect iv = bit();
      
      // Compute Vphi x Bpol
      a_var[dit](iv,1) = velocity[dit](iv,1) * Bpol_mag[dit](iv,0);
      a_var[dit](iv,1) *= sign[dit](iv,0)/abs(sign[dit](iv,0))/larmor_number;
      
      // Compute Vpol x Bphi
      a_var[dit](iv,2) = velocity[dit](iv,2) * abs(B_phi[dit](iv,1));
      a_var[dit](iv,2) *= -sign[dit](iv,0)/abs(sign[dit](iv,0))/larmor_number;

      // Compute radial force balance
      double force_balance = a_var[dit](iv,0)+a_var[dit](iv,1)+a_var[dit](iv,2)+E_r[dit](iv,0);
      a_var[dit](iv,3) = force_balance;
    }
  }

  return;
}

void GKDiagnostics::getTemperature( CFG::LevelData<CFG::FArrayBox>& a_var,
                                    const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.temperature( a_var );
  return;
}

void GKDiagnostics::getParallelTemperature( CFG::LevelData<CFG::FArrayBox>& a_var,
                                            const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.parallelTemperature( a_var );
  return;
}

void GKDiagnostics::getPerpTemperature( CFG::LevelData<CFG::FArrayBox>& a_var,
                                        const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.perpTemperature( a_var );
  return;
}

void GKDiagnostics::getParallelHeatFlux(  CFG::LevelData<CFG::FArrayBox>& a_var,
                                          const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
  
  CFG::LevelData<CFG::FArrayBox> v_parallel( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.parallelVelocity( v_parallel );
  
  a_var.define( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.parallelHeatFlux( a_var, v_parallel );
  return;
}

void GKDiagnostics::getTotalParallelHeatFlux( CFG::LevelData<CFG::FArrayBox>& a_var,
                                              const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
  
  CFG::LevelData<CFG::FArrayBox> zero_velocity(  mag_geom.gridsFull(), 
                                                 1, 
                                                 CFG::IntVect::Zero );
  for (auto dit(zero_velocity.dataIterator()); dit.ok(); ++dit){
     zero_velocity[dit].setVal(0.0);
  }

  a_var.define(  mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.parallelHeatFlux( a_var, zero_velocity );
  return;
}

void GKDiagnostics::getFourthMoment(  CFG::LevelData<CFG::FArrayBox>& a_var,
                              const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
  
  CFG::LevelData<CFG::FArrayBox> density( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.numberDensity( density );
  
  CFG::LevelData<CFG::FArrayBox> v_parallel( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.parallelVelocity( v_parallel );
  
  CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.pressure(pressure, v_parallel);

  a_var.define( mag_geom.gridsFull(), 1, CFG::IntVect::Zero);
  a_soln_species.fourthMoment( a_var );

  CFG::LevelData<CFG::FArrayBox> temp;
  temp.define(pressure);

  for (auto dit(a_var.dataIterator()); dit.ok(); ++dit) {
    a_var[dit].divide(temp[dit]);    // fourthMom/Pressure
    temp[dit].divide(density[dit]);  // Pressure/Density
    a_var[dit].divide(temp[dit]);    // fourthMom/(N*T^2)
    a_var[dit].mult(4.0/15.0);       // should be unity for Maxwellian!!!
  }

  return;
}

void GKDiagnostics::getRadialCurrent(CFG::LevelData<CFG::FArrayBox>& a_var,
                                     const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
  a_var.define( mag_geom.gridsFull(), 2, CFG::IntVect::Zero );
  
  computeRadialFlux(a_var,
                    a_soln_species,
                    m_ops->getPhi(),
                    m_ops->getEMFields(),
                    PhaseGeom::FULL_VELOCITY,
                    ChargeDensityKernel<FluxBox>());

  return;
}

void GKDiagnostics::getRadialParticleFlux(CFG::LevelData<CFG::FArrayBox>& a_var,
                                          const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
  
  a_var.define( mag_geom.gridsFull(), 2, CFG::IntVect::Zero );

#ifdef VLASOV_FLUXES
  computeRadialFlux(a_var,
                    a_soln_species,
                    m_ops->getPhi(),
                    m_ops->getEMFields(),
                    PhaseGeom::FULL_VELOCITY,
                    DensityKernel<FluxBox>());
#else
  LevelData<FluxBox> E_field_tmp;
  m_ops->getEField( E_field_tmp );
  a_soln_species.radialParticleFlux( a_var,E_field_tmp,PhaseGeom::FULL_VELOCITY );
#endif
  
  return;
}

void GKDiagnostics::getRadialHeatFlux( CFG::LevelData<CFG::FArrayBox>& a_var,
                                      const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 2, CFG::IntVect::Zero );
  
#ifdef VLASOV_FLUXES
    computeRadialFlux(a_var,
                      a_soln_species,
                      m_ops->getPhi(),
                      m_ops->getEMFields(),
                      PhaseGeom::FULL_VELOCITY,
                      EnergyKernel<FluxBox>(m_ops->getPhi()));
#else
  LevelData<FluxBox> E_field_tmp;
  m_ops->getEField( E_field_tmp );
  
  LevelData<FArrayBox> phi_injected_tmp;
  phase_geometry.injectConfigurationToPhase( m_ops->getPhi(), phi_injected_tmp );

  a_soln_species.radialHeatFlux( a_var,
                                E_field_tmp,
                                phi_injected_tmp,
                                PhaseGeom::FULL_VELOCITY );
#endif
  return;
}

void GKDiagnostics::getExBHeatFlux( CFG::LevelData<CFG::FArrayBox>& a_var,
                                    const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 2, CFG::IntVect::Zero );
  
#ifdef VLASOV_FLUXES
    computeRadialFlux(a_var,
                      a_soln_species,
                      m_ops->getPhi(),
                      m_ops->getEMFields(),
                      PhaseGeom::EXB_DRIFT_VELOCITY,
                      EnergyKernel<FluxBox>(m_ops->getPhi()));
#else
  LevelData<FluxBox> E_field_tmp;
  m_ops->getEField( E_field_tmp );
  
  LevelData<FArrayBox> phi_injected_tmp;
  phase_geometry.injectConfigurationToPhase( m_ops->getPhi(), phi_injected_tmp );
  
  a_soln_species.radialHeatFlux(a_var,
                                E_field_tmp,
                                phi_injected_tmp,
                                PhaseGeom::EXB_DRIFT_VELOCITY );
#endif
  
  return;
}

void GKDiagnostics::getExBKineticEnergyFlux(  CFG::LevelData<CFG::FArrayBox>& a_var,
                                              const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
  
  a_var.define( mag_geom.gridsFull(), 2, CFG::IntVect::Zero );
  
#ifdef VLASOV_FLUXES
    computeRadialFlux(a_var,
                      a_soln_species,
                      m_ops->getPhi(),
                      m_ops->getEMFields(),
                      PhaseGeom::EXB_DRIFT_VELOCITY,
                      KineticEnergyKernel<FluxBox>());
#else
  LevelData<FluxBox> E_field_tmp;
  m_ops->getEField( E_field_tmp );

  LevelData<FArrayBox> phi_injected_tmp;
  CFG::LevelData<CFG::FArrayBox> zero_phi(m_ops->getPhi().getBoxes(),
                                 1,
                                 m_ops->getPhi().ghostVect() );
  for (auto dit(zero_phi.dataIterator()); dit.ok(); ++dit) {
    zero_phi[dit].setVal(0.);
  }
  phase_geometry.injectConfigurationToPhase( zero_phi, phi_injected_tmp );

  a_soln_species.radialHeatFlux(a_var,
                                E_field_tmp,
                                phi_injected_tmp,
                                PhaseGeom::EXB_DRIFT_VELOCITY );
#endif
  
  return;
}

void GKDiagnostics::getChargeDensity( CFG::LevelData<CFG::FArrayBox>& a_var,
                                      const KineticSpeciesPtrVect& a_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );

  CFG::LevelData<CFG::FArrayBox> charge_density( mag_geom.gridsFull(), 
                                                 1, 
                                                 CFG::IntVect::Zero );

  for (auto dit(a_var.dataIterator()); dit.ok(); ++dit) {
    a_var[dit].setVal(0.);
  }

  for (int species(0); species<a_species.size(); species++) {
    KineticSpecies& soln_species( *(a_species[species]) );
    soln_species.chargeDensity(charge_density);

    for (auto dit(a_var.dataIterator()); dit.ok(); ++dit) {
      a_var[dit] += charge_density[dit];
    }
  }

  if (m_ops->usingBoltzmannElectrons()) {
    const CFG::LevelData<CFG::FArrayBox>& ne( m_ops->boltzmannElectronNumberDensity() );
    for (auto dit(a_var.dataIterator()); dit.ok(); ++dit) {
      a_var[dit] -= ne[dit];
    }
  }

  return;
}

void GKDiagnostics::writeDfnPlots(const std::vector<std::string>& a_plotnames,
                                  const std::vector<std::string>& a_filenames,
                                  const KineticSpecies&           a_species,
                                  const Real                      a_cur_time,
                                  const int                       a_radial_index,
                                  const int                       a_poloidal_index,
                                  const int                       a_toroidal_index,
                                  const int                       a_vpar_index,
                                  const int                       a_mu_index ) const
{
  const PhaseGeom& species_geometry = a_species.phaseSpaceGeometry();

  for (int i = 0; i < a_plotnames.size(); i++) {

    std::string plotname = a_plotnames[i];

    LevelData<FArrayBox> var;
    getPhaseVar(  var, plotname, a_species );

    if (plotname == "dfn") {
      plotPhaseVar( var, 
                    species_geometry, 
                    a_filenames[i], 
                    a_cur_time );
    } else if (plotname == "bstar_par") {
      plotPhaseVar( var, 
                    species_geometry, 
                    a_filenames[i], 
                    a_cur_time );
    } else if (plotname == "deltaF") {
      plotPhaseVar( var, 
                    species_geometry, 
                    a_filenames[i], 
                    a_cur_time );
    } else if (plotname == "dfn_at_mu") {
      plotPhaseVarAtMu( var, 
                        species_geometry, 
                        a_filenames[i], 
                        a_mu_index,
                        a_cur_time );
    } else if (plotname == "vpartheta") {
      plotPhaseVarVParTheta(  var, 
                              species_geometry, 
                              a_filenames[i], 
                              a_radial_index,
                              a_toroidal_index,
                              a_mu_index,
                              a_cur_time );
    } else if (plotname == "bfvpartheta") {
      plotPhaseVarVParTheta(  var, 
                              species_geometry, 
                              a_filenames[i], 
                              a_radial_index,
                              a_toroidal_index,
                              a_mu_index,
                              a_cur_time );
    } else if (plotname == "rtheta") {
      plotPhaseVarAtVelIndex( var, 
                              species_geometry, 
                              a_filenames[i], 
                              a_vpar_index,
                              a_mu_index,
                              a_cur_time );
    } else if (plotname == "vparmu") {
      plotPhaseVarAtCfgIndex( var, 
                              species_geometry, 
                              a_filenames[i], 
                              a_radial_index,
                              a_poloidal_index,
                              a_toroidal_index,
                              a_cur_time );
    }
  }
  return;
}

void GKDiagnostics::defineFluidCellVar( CFG::LevelData<CFG::FArrayBox>& a_var,
                                        const CFG::FluidSpecies& a_species,
                                        const string& a_name ) const
{
  a_var.define( a_species.cell_var(a_name) );
  return;
}

void GKDiagnostics::defineFluidFaceVarAtCell( CFG::LevelData<CFG::FArrayBox>& a_var,
                                              const CFG::FluidSpecies& a_fluid_species,
                                              const string& a_var_name ) const
{
  // convert contravarient vars on faces to physical vars at cells
  const CFG::IntVect& ghostVect 
     = a_fluid_species.face_var(a_var_name).ghostVect();
  a_var.define( a_fluid_species.face_var(a_var_name).getBoxes(), 
                CFG_DIM, ghostVect);
  return;   
}

void GKDiagnostics::defineFluidEdgeVarAtCell( CFG::LevelData<CFG::FArrayBox>& a_var,
                                              const CFG::FluidSpecies& a_fluid_species,
                                              const string& a_var_name ) const
{
  // convert covariant vars on edges to physical vars at cells
  const CFG::IntVect& ghostVect 
     = a_fluid_species.edge_var(a_var_name).ghostVect();
  a_var.define(  a_fluid_species.edge_var(a_var_name).getBoxes(), 
                 CFG_DIM, ghostVect );
  return;
}

void GKDiagnostics::defineFluidOpMember(CFG::LevelData<CFG::FArrayBox>& a_var,
                                        const CFG::FluidSpecies& a_fluid_species,
                                        const std::string& a_member_var_name ) const
{
  m_fluid_op->getMemberVarForPlotting( a_var, 
                                       a_fluid_species, 
                                       a_member_var_name);
 
}

void GKDiagnostics::getFluidCellVar(  CFG::LevelData<CFG::FArrayBox>&  a_var,
                                      const CFG::CFGVars&              a_fluid_vars,
                                      const string&                    a_name ) const
{
  const CFG::LevelData<CFG::FArrayBox>& fluid_var(a_fluid_vars.cell_var(a_name));
  a_var.define( a_fluid_vars.cell_var(a_name) );
  for (auto dit(a_var.dataIterator()); dit.ok(); ++dit) {
    a_var[dit].copy(fluid_var[dit]);
  }
  return;
}

void GKDiagnostics::getFluidFaceVarAtCell(  CFG::LevelData<CFG::FArrayBox>&  a_var,
                                            const CFG::CFGVars&              a_fluid_vars,
                                            const string&                    a_var_name ) const
{
  // convert contravarient vars on faces to physical vars at cells
  const CFG::IntVect& ghostVect 
     = a_fluid_vars.face_var(a_var_name).ghostVect();
  a_var.define( a_fluid_vars.face_var(a_var_name).getBoxes(), 
                CFG_DIM, ghostVect);
  a_fluid_vars.interpFaceVarToCell(a_var, a_var_name);
  
  if(!a_fluid_vars.plotMappedFaceVars()) {
     const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
     mag_geom.convertPhysToContravar(a_var,1);
  }

  return;   
}

void GKDiagnostics::getFluidEdgeVarAtCell(  CFG::LevelData<CFG::FArrayBox>&  a_var,
                                            const CFG::CFGVars&              a_fluid_vars,
                                            const string&                    a_var_name ) const
{
  // convert covariant vars on edges to physical vars at cells
  const CFG::IntVect& ghostVect 
     = a_fluid_vars.edge_var(a_var_name).ghostVect();
  a_var.define(  a_fluid_vars.edge_var(a_var_name).getBoxes(), 
                 CFG_DIM, ghostVect );
  a_fluid_vars.interpEdgeVarToCell(a_var, a_var_name);
  
  if(!a_fluid_vars.plotMappedEdgeVars()) {
     const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
     mag_geom.convertPhysToCovar(a_var,1);
  }

  return;
}

void GKDiagnostics::getFluidNodeVarAtCell(  CFG::LevelData<CFG::FArrayBox>&  a_var,
                                            const CFG::CFGVars&              a_fluid_vars,
                                            const string&                    a_var_name ) const
{
  // convert var on node to var at cells
  const CFG::IntVect& ghostVect
     = a_fluid_vars.node_var(a_var_name).ghostVect();
  a_var.define(  a_fluid_vars.node_var(a_var_name).getBoxes(),
                 1, ghostVect );
  a_fluid_vars.interpNodeVarToCell(a_var, a_var_name);

  return;
}

void GKDiagnostics::getFluidOpMember( CFG::LevelData<CFG::FArrayBox>&  a_var,
                                      const CFG::CFGVars&              a_fluid_vars,
                                      const std::string&               a_member_var_name ) const
{
   m_fluid_op->getMemberVarForPlotting(a_var, a_fluid_vars, a_member_var_name);
}


void GKDiagnostics::computeRadialFlux(CFG::LevelData<CFG::FArrayBox>&        a_RadialFlux,
                                      const KineticSpecies&                  a_soln_species,
                                      const CFG::LevelData<CFG::FArrayBox>&  a_phi,
                                      const CFG::EMFields&                   a_EM_fields,
                                      const int                              a_velocity_option,
                                      const Kernel<FluxBox>&                 a_kernel) const
{
  /*
    Computes radial FS-averagead and FS-integrated fluxes of
    the provided kernel function (e.g., for the density kernel
    we obtain particle flux). This calculation returns fluxes
    that are consistent with those obtained by the
    vlasov discretization. That is all low-pollution options
    are taken into account.
   */

  // Get geometry
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
  
  // Create vlasov class using input file initializations
  ParmParse pp_vlasov("gkvlasov");
  const double larmor_number( m_units->larmorNumber() );
  GKVlasov vlasov(pp_vlasov, larmor_number);
  
  // Create a copy of the solution species with ghosts, which are needed for flux calculations
  int order = (phase_geometry.secondOrder()) ? 2 : 4;
  int ghosts = (order == 2) ? 2 : 3;
  KineticSpeciesPtr species_with_ghosts( a_soln_species.clone( ghosts*IntVect::Unit, true ) );
  
  LevelData<FArrayBox>& dfn_with_ghosts( species_with_ghosts->distributionFunction() );
   
  // Fill interior ghosts; have to use species-specific geometry
  const PhaseGeom& species_geometry = a_soln_species.phaseSpaceGeometry();
  species_geometry.fillInternalGhosts(dfn_with_ghosts);
  
  // Extrapolate to physical ghosts
  int nghost = 2;
  species_geometry.extrapolateAtPhysicalBoundaries(dfn_with_ghosts, order, nghost);
  
  // Get normal component of a moment flux integrated over a cell face
  CFG::LevelData<CFG::FluxBox> flux_normal(mag_geom.grids(), 1, CFG::IntVect::Unit);
      
  // The last artument is time, which is not presently
  // used to compute velocities.
  vlasov.computeIntegratedMomentFluxNormals(flux_normal,
                                            *species_with_ghosts,
                                            a_EM_fields,
                                            a_velocity_option,
                                            a_kernel,
                                            0.);
   
  // FS averaging can only be done for an FArrayBox, so first need to interpolate to cell-centeres
  // Becasue this is used only for diagnostics purposes, use slow C++ code
  CFG::LevelData<CFG::FArrayBox> flux_normal_cc(mag_geom.grids(), 1, CFG::IntVect::Zero);
  for (CFG::DataIterator dit(flux_normal_cc.dataIterator()); dit.ok(); ++dit) {
    CFG::BoxIterator bit(flux_normal_cc[dit].box());
    for (bit.begin(); bit.ok(); ++bit) {
      CFG::IntVect iv = bit();
      CFG::IntVect iv_shift = iv;
      iv_shift[RADIAL_DIR] += 1;
      flux_normal_cc[dit](iv,0) = 0.5 * ( flux_normal[dit][RADIAL_DIR](iv,0)
                                        + flux_normal[dit][RADIAL_DIR](iv_shift,0) );
    }
  }
  
  // Becasue the flux_normal object already contains the face element area, need to divide it out
  CFG::LevelData<CFG::FArrayBox> face_areas(mag_geom.grids(), CFG_DIM, CFG::IntVect::Unit);
  mag_geom.getPointwiseFaceAreas(face_areas);
  for (CFG::DataIterator dit(face_areas.dataIterator()); dit.ok(); ++dit) {
    flux_normal_cc[dit].divide(face_areas[dit],RADIAL_DIR, 0, 1);
  }
 
  // Create FluxSurface class performing FS-averaging (not shell averaging)
  CFG::FluxSurface flux_surface(mag_geom, false);

  // Create a temporary
  CFG::LevelData<CFG::FArrayBox> flux_tmp(mag_geom.grids(), 1, CFG::IntVect::Zero);

  // Get flux-surfaced averaged component
  flux_surface.averageAndSpread(flux_normal_cc, flux_tmp);
  for (CFG::DataIterator dit(a_RadialFlux.dataIterator()); dit.ok(); ++dit) {
    a_RadialFlux[dit].copy(flux_tmp[dit], 0, 0, 1);
  }
  
  // Get flux-surface integrated component
  flux_surface.integrateAndSpread(flux_normal_cc, flux_tmp);
  for (CFG::DataIterator dit(a_RadialFlux.dataIterator()); dit.ok(); ++dit) {
    a_RadialFlux[dit].copy(flux_tmp[dit], 0, 1, 1);
  }
}

#include "NamespaceFooter.H"
