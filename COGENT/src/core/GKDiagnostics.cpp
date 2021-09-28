#include "GKDiagnostics.H"

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

void GKDiagnostics::getCfgVar(  CFG::LevelData<CFG::FluxBox>& a_var,
                                const std::string&            a_var_name ) const
{
  if (a_var_name == "efield") {
    getEField( a_var );
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
  } else if (a_var_name == "momentum") {
    getMomentum( a_var, a_species );
  } else if (a_var_name == "ParallelMomentum") {
    getParallelMomentum( a_var, a_species );
  } else if (a_var_name == "PoloidalMomentum") {
    getPoloidalMomentum( a_var, a_species );
  } else if (a_var_name == "ParallelVelocity") {
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
  } else if (a_var_name == "gradPoverN") {
    getGradPoverN( a_var, a_species );
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
  } else if (a_var_name == "ParticleFlux") {
    getParticleFlux( a_var, a_species );
  } else if (a_var_name == "HeatFlux") {
    getHeatFlux( a_var, a_species );
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

  CFG::LevelData<CFG::FArrayBox> density(  mag_geom.gridsFull(), 
                                           1, 
                                           CFG::IntVect::Zero );
  a_soln_species.numberDensity( density );

  CFG::LevelData<CFG::FArrayBox> ParallelMom(  mag_geom.gridsFull(), 
                                               1, 
                                               CFG::IntVect::Zero );
  a_soln_species.ParallelMomentum( ParallelMom );

  for (auto dit(density.dataIterator()); dit.ok(); ++dit) {
    ParallelMom[dit].divide(density[dit]);
  }
  CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.gridsFull(), 
                                           1, 
                                           CFG::IntVect::Zero );
  a_soln_species.pressureMoment(pressure, ParallelMom);

  for (auto dit(density.dataIterator()); dit.ok(); ++dit) {
    pressure[dit].divide(density[dit]);
  }

  const PhaseGeom& species_geometry = a_soln_species.phaseSpaceGeometry();
  a_var.define( species_geometry.gridsFull(), 1, IntVect::Zero );
  DeltaFKernel<FArrayBox> DeltaF_Kernel(density, pressure, ParallelMom);
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

void GKDiagnostics::getMomentum(  CFG::LevelData<CFG::FArrayBox>& a_var,
                                  const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const double larmor_number( m_units->larmorNumber() );

  LevelData<FluxBox> E_field_tmp;
  m_ops->getEField( E_field_tmp );

  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), CFG_DIM, CFG::IntVect::Zero );
  a_soln_species.momentumDensity( a_var, E_field_tmp, larmor_number );

  return;
}

void GKDiagnostics::getParallelMomentum( CFG::LevelData<CFG::FArrayBox>& a_var,
                                         const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.ParallelMomentum(a_var);
}

void GKDiagnostics::getPoloidalMomentum(  CFG::LevelData<CFG::FArrayBox>& a_var,
                                          const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const double larmor_number( m_units->larmorNumber() );

  LevelData<FluxBox> E_field_tmp;
  m_ops->getEField( E_field_tmp );

  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.PoloidalMomentum( a_var, E_field_tmp, larmor_number );
  
  return;
}

void GKDiagnostics::getParallelVelocity(  CFG::LevelData<CFG::FArrayBox>& a_var,
                                          const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.ParallelVelocity( a_var );
}

void GKDiagnostics::getEnergyDensity( CFG::LevelData<CFG::FArrayBox>& a_var,
                                      const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 3, CFG::IntVect::Zero );
  
  CFG::LevelData<CFG::FArrayBox> energy_density_loc( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.energyMoment( energy_density_loc, m_ops->getPhi() );

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
  a_soln_species.kineticEnergyMoment( energy_density_loc );

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
  a_soln_species.pressureMoment( a_var );
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

void GKDiagnostics::getGradPoverN(CFG::LevelData<CFG::FArrayBox>& a_var,
                                  const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  CFG::LevelData<CFG::FArrayBox> density(  mag_geom.gridsFull(), 
                                           1, 
                                           CFG::IntVect::Zero );
  CFG::LevelData<CFG::FArrayBox> ParallelMom(  mag_geom.gridsFull(), 
                                               1, 
                                               CFG::IntVect::Zero );
  a_soln_species.numberDensity( density );
  a_soln_species.ParallelMomentum( ParallelMom );

  for (auto dit(density.dataIterator()); dit.ok(); ++dit) {
    ParallelMom[dit].divide(density[dit]);
  }
  CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.gridsFull(), 
                                           1, 
                                           CFG::IntVect::Zero );
  a_soln_species.pressureMoment(pressure, ParallelMom);

  const CFG::DisjointBoxLayout& cfg_grids = mag_geom.gridsFull();
  CFG::LevelData<CFG::FArrayBox> pressure_grown( cfg_grids, 
                                                 1, 
                                                 2*CFG::IntVect::Unit );
  for (auto dit(pressure.dataIterator()); dit.ok(); ++dit) {
    pressure_grown[dit].copy(pressure[dit]);
    const CFG::MagBlockCoordSys& coord_sys 
        = mag_geom.getBlockCoordSys(cfg_grids[dit]);
    const CFG::ProblemDomain& domain = coord_sys.domain();
    secondOrderCellExtrapAtDomainBdry(  pressure_grown[dit], 
                                        cfg_grids[dit], 
                                        domain  );
  }
  mag_geom.fillInternalGhosts(pressure_grown);

#if CFG_DIM==2
  CFG::LevelData<CFG::FArrayBox> pressure_gradient_mapped( cfg_grids, 
                                                           CFG_DIM, 
                                                           CFG::IntVect::Unit );
  mag_geom.computeMappedPoloidalGradientWithGhosts(pressure_grown, 
                                                   pressure_gradient_mapped, 
                                                   2);
  
  CFG::LevelData<CFG::FArrayBox> pressure_gradient(cfg_grids, 
                                                   CFG_DIM, 
                                                   CFG::IntVect::Unit);
  mag_geom.unmapPoloidalGradient(pressure_gradient_mapped, 
                                 pressure_gradient);

  a_var.define( cfg_grids, CFG_DIM, CFG::IntVect::Zero );
  for (auto dit(cfg_grids.dataIterator()); dit.ok(); ++dit) {
     mag_geom.getBlockCoordSys(cfg_grids[dit]).computePsiThetaProjections(pressure_gradient[dit]);
     a_var[dit].copy(pressure_gradient[dit]);
     for (int n=0; n<a_var.nComp(); ++n) {
        a_var[dit].divide(density[dit], 0, n);
     }
  }
  
#else
  MayDay::Error("GKDiagnostics::gradPover this diagnostcs is not implemented in 3D");
#endif
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
  
  CFG::LevelData<CFG::FArrayBox> density( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  CFG::LevelData<CFG::FArrayBox> ParallelMom( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.numberDensity( density );
  a_soln_species.ParallelMomentum( ParallelMom );
  
  for (auto dit(density.dataIterator()); dit.ok(); ++dit) {
     ParallelMom[dit].divide(density[dit]);
  }
  a_var.define( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
  a_soln_species.parallelHeatFluxMoment( a_var, ParallelMom );
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
  a_soln_species.parallelHeatFluxMoment( a_var, zero_velocity );
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

  CFG::LevelData<CFG::FArrayBox> ParallelMom( mag_geom.gridsFull(), 
                                 1, 
                                 CFG::IntVect::Zero );
  a_soln_species.ParallelMomentum( ParallelMom );
  for (auto dit(density.dataIterator()); dit.ok(); ++dit) {
    ParallelMom[dit].divide(density[dit]);
  }

  CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.gridsFull(), 
                                           1, 
                                           CFG::IntVect::Zero );
  a_soln_species.pressureMoment(pressure, ParallelMom);

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

void GKDiagnostics::getParticleFlux( CFG::LevelData<CFG::FArrayBox>& a_var,
                                     const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );

  LevelData<FluxBox> E_field_tmp;
  m_ops->getEField( E_field_tmp );

  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
  
  a_var.define( mag_geom.gridsFull(), 2, CFG::IntVect::Zero );
  a_soln_species.ParticleFlux( a_var,E_field_tmp,PhaseGeom::FULL_VELOCITY );
  return;
}

void GKDiagnostics::getHeatFlux(  CFG::LevelData<CFG::FArrayBox>& a_var,
                                  const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  LevelData<FluxBox> E_field_tmp;
  LevelData<FArrayBox> phi_injected_tmp;

  m_ops->getEField( E_field_tmp );
  phase_geometry.injectConfigurationToPhase( m_ops->getPhi(), phi_injected_tmp );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 2, CFG::IntVect::Zero );
  a_soln_species.HeatFlux( a_var, 
                           E_field_tmp, 
                           phi_injected_tmp, 
                           PhaseGeom::FULL_VELOCITY );
  return;
}

void GKDiagnostics::getExBHeatFlux( CFG::LevelData<CFG::FArrayBox>& a_var,
                                    const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  LevelData<FluxBox> E_field_tmp;
  LevelData<FArrayBox> phi_injected_tmp;

  m_ops->getEField( E_field_tmp );
  phase_geometry.injectConfigurationToPhase( m_ops->getPhi(), phi_injected_tmp );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 2, CFG::IntVect::Zero );
  a_soln_species.HeatFlux( a_var, 
                           E_field_tmp, 
                           phi_injected_tmp, 
                           PhaseGeom::EXB_DRIFT_VELOCITY );
  return;
}

void GKDiagnostics::getExBKineticEnergyFlux(  CFG::LevelData<CFG::FArrayBox>& a_var,
                                              const KineticSpecies& a_soln_species ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  LevelData<FluxBox> E_field_tmp;
  LevelData<FArrayBox> phi_injected_tmp;

  m_ops->getEField( E_field_tmp );
  
  CFG::LevelData<CFG::FArrayBox> zero_phi(m_ops->getPhi().getBoxes(), 
                                 1, 
                                 m_ops->getPhi().ghostVect() );
   for (auto dit(zero_phi.dataIterator()); dit.ok(); ++dit) {
     zero_phi[dit].setVal(0.);
   }
  
  phase_geometry.injectConfigurationToPhase( zero_phi, phi_injected_tmp );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

  a_var.define( mag_geom.gridsFull(), 2, CFG::IntVect::Zero );
  a_soln_species.HeatFlux( a_var, 
                           E_field_tmp, 
                           phi_injected_tmp, 
                           PhaseGeom::EXB_DRIFT_VELOCITY );
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

void GKDiagnostics::getFluidCellVar(  CFG::LevelData<CFG::FArrayBox>& a_var,
                                      const CFG::FluidSpecies& a_species,
                                      const string& a_name ) const
{
  const CFG::LevelData<CFG::FArrayBox>& fluid_var(a_species.cell_var(a_name));
  a_var.define( a_species.cell_var(a_name) );
  for (auto dit(a_var.dataIterator()); dit.ok(); ++dit) {
    a_var[dit].copy(fluid_var[dit]);
  }
  return;
}

void GKDiagnostics::getFluidFaceVarAtCell(  CFG::LevelData<CFG::FArrayBox>& a_var,
                                            const CFG::FluidSpecies& a_fluid_species,
                                            const string& a_var_name ) const
{
  // convert contravarient vars on faces to physical vars at cells
  const CFG::IntVect& ghostVect 
     = a_fluid_species.face_var(a_var_name).ghostVect();
  a_var.define( a_fluid_species.face_var(a_var_name).getBoxes(), 
                CFG_DIM, ghostVect);
  a_fluid_species.interpFaceVarToCell(a_var, a_var_name);
  
  const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
  mag_geom.convertPhysToContravar(a_var,1); // 1 is for inverse

  return;   
}

void GKDiagnostics::getFluidEdgeVarAtCell(  CFG::LevelData<CFG::FArrayBox>& a_var,
                                            const CFG::FluidSpecies& a_fluid_species,
                                            const string& a_var_name ) const
{
  // convert covariant vars on edges to physical vars at cells
  const CFG::IntVect& ghostVect 
     = a_fluid_species.edge_var(a_var_name).ghostVect();
  a_var.define(  a_fluid_species.edge_var(a_var_name).getBoxes(), 
                 CFG_DIM, ghostVect );
  a_fluid_species.interpEdgeVarToCell(a_var, a_var_name);
  
  const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
  mag_geom.convertPhysToCovar(a_var,1); // 1 is for inverse

  return;
}

void GKDiagnostics::getFluidOpMember( CFG::LevelData<CFG::FArrayBox>& a_var,
                                      const CFG::FluidSpecies& a_fluid_species,
                                      const std::string& a_member_var_name ) const
{
  m_fluid_op->getMemberVarForPlotting( a_var, 
                                       a_fluid_species, 
                                       a_member_var_name);
 
}

#include "NamespaceFooter.H"
