#include "FourthOrderUtil.H"
#include "KineticSpecies.H"
#include "MomentOp.H"
#include "Misc.H"
#include "Directions.H"
#include "CFArrayBox.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FluxSurface.H"
#include "FourthOrderUtil.H"
#include "MagGeom.H"
#include "MagBlockCoordSys.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

KineticSpecies::KineticSpecies(
         const string&                   a_name,
         const Real                      a_mass,
         const Real                      a_charge,
         const RefCountedPtr<PhaseGeom>& a_geometry,
         const bool                      a_is_gyrokinetic
   )
  : m_geometry( a_geometry ),
    m_name( a_name ),
    m_mass( a_mass ),
    m_charge( a_charge ),
    m_moment_op( MomentOp::instance() ),
    m_velocity_option(PhaseGeom::FULL_VELOCITY),
    m_gyroavg_op(NULL),
    m_is_gyrokinetic(a_is_gyrokinetic)
{
    int ghost = (m_geometry->secondOrder()) ? 0 : 1;
    const DisjointBoxLayout& dbl = m_geometry->gridsFull();
    m_velocity.define(dbl, SpaceDim, ghost * IntVect::Unit);
}


KineticSpecies::KineticSpecies( const KineticSpecies& a_foo )
  : m_geometry( a_foo.m_geometry ),
    m_name( a_foo.m_name ),
    m_mass( a_foo.m_mass ),
    m_charge( a_foo.m_charge ),
    m_moment_op( MomentOp::instance() ),
    m_velocity_option(a_foo.m_velocity_option),
    m_is_gyrokinetic( a_foo.m_is_gyrokinetic )
{
   gyroaverageOp(a_foo.gyroaverageOp());
   m_dist_func.define( a_foo.m_dist_func );
   m_velocity.define( a_foo.m_velocity );
}

void KineticSpecies::numberDensity( CFG::LevelData<CFG::FArrayBox>& a_rho ) const
{
   m_moment_op.compute( a_rho, *this, DensityKernel<FArrayBox>() );
}


void KineticSpecies::massDensity( CFG::LevelData<CFG::FArrayBox>& a_rho ) const
{
   m_moment_op.compute( a_rho, *this, MassDensityKernel<FArrayBox>() );
}


void KineticSpecies::chargeDensity( CFG::LevelData<CFG::FArrayBox>& a_rho ) const
{
   m_moment_op.compute( a_rho, *this, ChargeDensityKernel<FArrayBox>() );
}

void KineticSpecies::gyroaveragedChargeDensity( CFG::LevelData<CFG::FArrayBox>& a_rho ) const
{
   CH_assert(isGyrokinetic());
   KineticSpecies tmp(*this);
   m_gyroavg_op->applyOp(tmp.distributionFunction(), (*this).distributionFunction());
   m_moment_op.compute( a_rho, tmp, ChargeDensityKernel<FArrayBox>() );
}


void KineticSpecies::fluidVelocity(CFG::LevelData<CFG::FArrayBox>& a_velocity,
                                   const LevelData<FluxBox>& a_field,
                                   const CFG::EMFields&      a_EM_fields,
                                   const LevelData<FluxBox>& a_Apar_derivatives,
                                   const double larmor  ) const
{
   /*
      Computes Vfluid = V_gc + c/qn * curl (M * b)
      M = int( -mu * dfn * d^3v) is the magnetization
      make use of curl (M * b) = M*curlb + gradM x b
    
      This function returns R Phi Z components in 2D
      and X, Y, Z components in 3D
    */
   
   CH_assert(a_velocity.nComp() == 3);
   
   const CFG::MagGeom& mag_geom = m_geometry->magGeom();
   const CFG::MagCoordSys& coords = *mag_geom.getCoordSys();

   CFG::LevelData<CFG::FArrayBox> gc_particle_flux;
   gc_particle_flux.define(a_velocity);

#if CFG_DIM == 2
   CFG::LevelData<CFG::FArrayBox> gc_flow_rz(a_velocity.disjointBoxLayout(), 2, a_velocity.ghostVect());
   m_moment_op.compute( gc_flow_rz,
                       *this,
                        GuidingCenterVelocityKernel<FArrayBox>(a_field, a_EM_fields, a_Apar_derivatives, PhaseGeom::DIAGNOSTICS) );

   CFG::LevelData<CFG::FArrayBox> gc_flow_phi(a_velocity.disjointBoxLayout(), 2, a_velocity.ghostVect());
   m_moment_op.compute( gc_flow_phi,
                       *this,
                        GuidingCenterVelocityKernel<FArrayBox>(a_field, a_EM_fields, a_Apar_derivatives, PhaseGeom::TOROIDAL_VELOCITY) );
   
   for (CFG::DataIterator dit(gc_particle_flux.dataIterator()); dit.ok(); ++dit) {
      gc_particle_flux[dit].copy(gc_flow_rz[dit],0,0,1);
      gc_particle_flux[dit].copy(gc_flow_phi[dit],0,1,1);
      gc_particle_flux[dit].copy(gc_flow_rz[dit],1,2,1);
   }
#else
   m_moment_op.compute( gc_particle_flux, *this,
                        GuidingCenterVelocityKernel<FArrayBox>(a_field, a_EM_fields, a_Apar_derivatives, PhaseGeom::DIAGNOSTICS) );
#endif
   
   // Compute magnetization M
   CFG::LevelData<CFG::FArrayBox> magnetization(a_velocity.disjointBoxLayout(), 1,
                                                a_velocity.ghostVect());
   m_moment_op.compute( magnetization, *this, MagnetizationKernel<FArrayBox>() );

   CFG::LevelData<CFG::FArrayBox> magnetization_grown(a_velocity.disjointBoxLayout(), 1,
                                                      a_velocity.ghostVect()+2*CFG::IntVect::Unit);

   // Extrapolate to physical ghosts
   for (CFG::DataIterator dit(magnetization.dataIterator()); dit.ok(); ++dit) {
     magnetization_grown[dit].copy(magnetization[dit]);
   }
   mag_geom.fillInternalGhosts(magnetization_grown);

   int order = 2;
   int nghosts = 2;
   mag_geom.extrapolateAtPhysicalBoundaries(magnetization_grown, order, nghosts);
   
   // Compute magnetization gradient, gradM
   CFG::LevelData<CFG::FArrayBox> gradM_mapped(a_velocity.disjointBoxLayout(), 3,
                                               a_velocity.ghostVect()+CFG::IntVect::Unit);
   
#if CFG_DIM == 2
   CFG::LevelData<CFG::FArrayBox> gradM_mapped_poloidal(a_velocity.disjointBoxLayout(), 2,
                                               a_velocity.ghostVect()+CFG::IntVect::Unit);
   
   mag_geom.computeMappedGradient(magnetization_grown, gradM_mapped_poloidal, order);
   mag_geom.injectPoloidalVector(gradM_mapped_poloidal, gradM_mapped);
   
#else
   mag_geom.computeMappedGradient(magnetization_grown, gradM_mapped, order);
#endif

   CFG::LevelData<CFG::FArrayBox> gradM_phys(a_velocity.disjointBoxLayout(), 3,
                                             a_velocity.ghostVect()+CFG::IntVect::Unit);
   
   mag_geom.unmapGradient(gradM_mapped, gradM_phys);
   
   // Get magnetic geometry data
   const CFG::LevelData<CFG::FArrayBox>& bunit = mag_geom.getCCBFieldDir();
   const CFG::LevelData<CFG::FArrayBox>& curlb = mag_geom.getCCCurlBFieldDir();

   // Create objects to store gradMxb and Mcurlb
   CFG::LevelData<CFG::FArrayBox> gradMxb(a_velocity.disjointBoxLayout(), 3, a_velocity.ghostVect());
   CFG::LevelData<CFG::FArrayBox> Mcurlb(a_velocity.disjointBoxLayout(), 3, a_velocity.ghostVect());
   
   // Compute gradM x b
   mag_geom.crossProduct(gradMxb, gradM_phys, bunit);
   
   // Compute n*Vflow
   for (CFG::DataIterator dit(a_velocity.dataIterator()); dit.ok(); ++dit) {

      // Compute M*curlb
      for (int nComp=0; nComp<3; ++nComp) {
         Mcurlb[dit].copy(magnetization[dit],0,nComp,1);
         Mcurlb[dit].mult(curlb[dit],nComp,nComp,1);
      }
      
      // Perform normalization
      double fac = larmor / 2.0 / m_charge;
      gradMxb[dit].mult(fac);
      Mcurlb[dit].mult(fac);
      
      // Collect all terms
      a_velocity[dit].copy(gc_particle_flux[dit],0,0,3);
      a_velocity[dit].plus(Mcurlb[dit],0,0,3);
      a_velocity[dit].plus(gradMxb[dit],0,0,3);
   }
   
   // Divide by density to get the flow velocity
   CFG::LevelData<CFG::FArrayBox> density(a_velocity.disjointBoxLayout(), 1,
                                          a_velocity.ghostVect());
   m_moment_op.compute( density, *this, DensityKernel<FArrayBox>() );

   for (CFG::DataIterator dit(density.dataIterator()); dit.ok(); ++dit) {
      for (int nComp=0; nComp<3; ++nComp) {
         a_velocity[dit].divide(density[dit], 0, nComp, 1);
      }
   }
}

void KineticSpecies::parallelVelocity( CFG::LevelData<CFG::FArrayBox>& a_parallel_vel ) const
{
   CFG::LevelData<CFG::FArrayBox> density;
   density.define(a_parallel_vel);
   m_moment_op.compute( density, *this, DensityKernel<FArrayBox>() );

   m_moment_op.compute( a_parallel_vel, *this, ParallelVelKernel<FArrayBox>() );

   for (CFG::DataIterator dit(density.dataIterator()); dit.ok(); ++dit) {
      a_parallel_vel[dit].divide(density[dit]);
   }
}


void KineticSpecies::parallelParticleFlux( CFG::LevelData<CFG::FArrayBox>& a_particle_flux ) const
{
   m_moment_op.compute( a_particle_flux, *this, ParallelVelKernel<FArrayBox>() );
}

void KineticSpecies::poloidalParticleFlux(CFG::LevelData<CFG::FArrayBox>& a_particle_flux,
                                          const LevelData<FluxBox>& a_field,
                                          const CFG::EMFields&      a_EM_fields,
                                          const LevelData<FluxBox>& a_Apar_derivatives,
                                          const double a_larmor  ) const
{


   const CFG::MagGeom& mag_geom = m_geometry->magGeom();
   const CFG::MagCoordSys& coords = *mag_geom.getCoordSys();

   CFG::LevelData<CFG::FArrayBox> gc_particle_flux;
   gc_particle_flux.define(a_particle_flux);
   m_moment_op.compute( gc_particle_flux, *this, GuidingCenterPoloidalVelKernel<FArrayBox>(a_field, a_EM_fields, a_Apar_derivatives) );
   
   CFG::LevelData<CFG::FArrayBox> magnetization;
   magnetization.define(a_particle_flux);
   m_moment_op.compute( magnetization, *this, MagnetizationKernel<FArrayBox>() );

   CFG::LevelData<CFG::FArrayBox> magnetization_grown(a_particle_flux.disjointBoxLayout(),1,
                                                      a_particle_flux.ghostVect()+2*CFG::IntVect::Unit);
   
   //Extrapolating magnetization moment into two layers of ghost cells
   for (CFG::DataIterator dit(magnetization.dataIterator()); dit.ok(); ++dit) {
      int block_number( coords.whichBlock( mag_geom.gridsFull()[dit] ) );
      const CFG::MagBlockCoordSys& block_coord_sys = (const CFG::MagBlockCoordSys&)(*(coords.getCoordSys(block_number)));
      const CFG::ProblemDomain& block_domain = block_coord_sys.domain();
      magnetization_grown[dit].copy(magnetization[dit]);
      fourthOrderCellExtrapAtDomainBdry(magnetization_grown[dit], block_domain, mag_geom.gridsFull()[dit]);
   }
   mag_geom.fillInternalGhosts(magnetization_grown);
   
   //Compute radial component of the magnetization gradient (part of curlM calculation)
   CFG::LevelData<CFG::FArrayBox> gradM_mapped(a_particle_flux.disjointBoxLayout(), 2,
                                               a_particle_flux.ghostVect()+CFG::IntVect::Unit);
  
   mag_geom.computeMappedPoloidalGradientWithGhosts(magnetization_grown,gradM_mapped,2);
   
   CFG::LevelData<CFG::FArrayBox> gradM_phys;
   gradM_phys.define(gradM_mapped);
   mag_geom.unmapPoloidalGradient(gradM_mapped, gradM_phys);
   
   CFG::LevelData<CFG::FArrayBox> gradM_r(gradM_phys.disjointBoxLayout(),1, gradM_phys.ghostVect());
   mag_geom.computeRadialProjection(gradM_r, gradM_phys);
   
   //Get magnetic geometry data
   const CFG::LevelData<CFG::FArrayBox>& bunit = mag_geom.getCCBFieldDir();
   const CFG::LevelData<CFG::FArrayBox>& curlb = mag_geom.getCCCurlBFieldDir();
   CFG::LevelData<CFG::FArrayBox> curlb_pol(curlb.disjointBoxLayout(),1,curlb.ghostVect());
   mag_geom.computePoloidalProjection(curlb_pol, curlb);

   
   double fac = a_larmor / 2.0 / m_charge;

   for (CFG::DataIterator dit(a_particle_flux.dataIterator()); dit.ok(); ++dit) {
      
      //Compute curl of Magnetization vector
      curlb_pol[dit].mult(magnetization[dit]);
      curlb_pol[dit].mult(fac);
      gradM_r[dit].mult(bunit[dit],1,0,1);
      gradM_r[dit].mult(fac);
      
      //Add the poloidal component of the guiding center velocity and the curl of magnetization
      a_particle_flux[dit].copy(gc_particle_flux[dit]);
      a_particle_flux[dit].plus(curlb_pol[dit]);
      a_particle_flux[dit].plus(gradM_r[dit]);
   }
}

void KineticSpecies::radialParticleFlux(CFG::LevelData<CFG::FArrayBox>& a_particle_flux,
                                        const LevelData<FluxBox>&       a_field,
                                        const CFG::EMFields&            a_EM_fields,
                                        const LevelData<FluxBox>&       a_Apar_derivatives,
                                        const int                       a_velocity_option) const

{
   m_moment_op.compute( a_particle_flux, *this, RadialParticleFluxKernel<FArrayBox>(a_field, a_EM_fields, a_Apar_derivatives, a_velocity_option) );

   const CFG::MagGeom& mag_geom = m_geometry->magGeom();
   CFG::FluxSurface flux_surface(mag_geom, false);
   CFG::LevelData<CFG::FArrayBox> particle_flux_tmp(mag_geom.grids(), 1, CFG::IntVect::Zero);
   CFG::LevelData<CFG::FArrayBox> particle_loc_tmp(mag_geom.grids(), 1, CFG::IntVect::Zero);

   // Get flux-surfaced averaged component
   for (CFG::DataIterator dit(a_particle_flux.dataIterator()); dit.ok(); ++dit) {
     particle_loc_tmp[dit].copy(a_particle_flux[dit], 0, 0, 1);
   }
   flux_surface.averageAndSpread(particle_loc_tmp, particle_flux_tmp);
   
   for (CFG::DataIterator dit(a_particle_flux.dataIterator()); dit.ok(); ++dit) {
     a_particle_flux[dit].copy(particle_flux_tmp[dit], 0, 0, 1);
   }
   
   // Get flux-surface integrated component
   for (CFG::DataIterator dit(a_particle_flux.dataIterator()); dit.ok(); ++dit) {
     particle_loc_tmp[dit].copy(a_particle_flux[dit], 1, 0, 1);
   }
   flux_surface.integrateAndSpread(particle_loc_tmp, particle_flux_tmp);
   
   for (CFG::DataIterator dit(a_particle_flux.dataIterator()); dit.ok(); ++dit) {
     a_particle_flux[dit].copy(particle_flux_tmp[dit], 0, 1, 1);
   }
}

void KineticSpecies::radialHeatFlux(CFG::LevelData<CFG::FArrayBox>& a_heat_flux,
                                    const LevelData<FluxBox>&       a_field,
                                    const CFG::EMFields&            a_EM_fields,
                                    const LevelData<FluxBox>&       a_Apar_derivatives,
                                    const LevelData<FArrayBox>&     a_phi,
                                    const int                       a_velocity_option) const

{
   m_moment_op.compute( a_heat_flux, *this, RadialHeatFluxKernel<FArrayBox>(a_field, a_EM_fields, a_Apar_derivatives, a_phi, a_velocity_option) );

   const CFG::MagGeom& mag_geom = m_geometry->magGeom();
   CFG::FluxSurface flux_surface(mag_geom, false);
   CFG::LevelData<CFG::FArrayBox> heat_flux_tmp(mag_geom.grids(), 1, CFG::IntVect::Zero);
   CFG::LevelData<CFG::FArrayBox> heat_loc_tmp(mag_geom.grids(), 1, CFG::IntVect::Zero);

   // Get flux-surfaced averaged component
   for (CFG::DataIterator dit(a_heat_flux.dataIterator()); dit.ok(); ++dit) {
     heat_loc_tmp[dit].copy(a_heat_flux[dit], 0, 0, 1);
   }
   flux_surface.averageAndSpread(heat_loc_tmp, heat_flux_tmp);
  
   for (CFG::DataIterator dit(a_heat_flux.dataIterator()); dit.ok(); ++dit) {
     a_heat_flux[dit].copy(heat_flux_tmp[dit], 0, 0, 1);
   }
  
   // Get flux-surface integrated component
   for (CFG::DataIterator dit(a_heat_flux.dataIterator()); dit.ok(); ++dit) {
     heat_loc_tmp[dit].copy(a_heat_flux[dit], 1, 0, 1);
   }
   flux_surface.integrateAndSpread(heat_loc_tmp, heat_flux_tmp);
  
   for (CFG::DataIterator dit(a_heat_flux.dataIterator()); dit.ok(); ++dit) {
     a_heat_flux[dit].copy(heat_flux_tmp[dit], 0, 1, 1);
   }
}

void KineticSpecies::parallelHeatFlux(CFG::LevelData<CFG::FArrayBox>& a_parallel_heat_flux,
                                      CFG::LevelData<CFG::FArrayBox>& a_vparshift ) const
{
   m_moment_op.compute( a_parallel_heat_flux, *this, ParallelHeatFluxKernel<FArrayBox>(a_vparshift) );
}

void KineticSpecies::pressure(CFG::LevelData<CFG::FArrayBox>& a_pressure,
                              CFG::LevelData<CFG::FArrayBox>& a_vparshift ) const
{
   m_moment_op.compute( a_pressure, *this, PressureKernel<FArrayBox>(a_vparshift) );
}

void KineticSpecies::pressure( CFG::LevelData<CFG::FArrayBox>& a_pressure ) const
{
   CFG::LevelData<CFG::FArrayBox> vpar;
   vpar.define( a_pressure );
   this->parallelVelocity(vpar);

   m_moment_op.compute( a_pressure, *this, PressureKernel<FArrayBox>(vpar) );
}

void KineticSpecies::parallelPressure( CFG::LevelData<CFG::FArrayBox>& a_pressure,
                                     CFG::LevelData<CFG::FArrayBox>& a_vparshift ) const
{
   m_moment_op.compute( a_pressure, *this, ParallelPressureKernel<FArrayBox>(a_vparshift) );
}

void KineticSpecies::parallelPressure( CFG::LevelData<CFG::FArrayBox>& a_pressure ) const
{
   CFG::LevelData<CFG::FArrayBox> vpar;
   vpar.define( a_pressure );
   this->parallelVelocity(vpar);

   m_moment_op.compute( a_pressure, *this, ParallelPressureKernel<FArrayBox>(vpar) );
}

void KineticSpecies::perpPressure( CFG::LevelData<CFG::FArrayBox>& a_pressure ) const
{
   m_moment_op.compute( a_pressure, *this, PerpPressureKernel<FArrayBox>() );
}

void KineticSpecies::temperature( CFG::LevelData<CFG::FArrayBox>& a_temperature ) const
{
   CFG::LevelData<CFG::FArrayBox> density;
   density.define( a_temperature );
   m_moment_op.compute( density, *this, DensityKernel<FArrayBox>() );

   CFG::LevelData<CFG::FArrayBox> vpar;
   vpar.define( a_temperature );
   m_moment_op.compute( vpar, *this, ParallelVelKernel<FArrayBox>() );

   for (CFG::DataIterator dit(vpar.dataIterator()); dit.ok(); ++dit) {
      vpar[dit].divide(density[dit]);
   }
   
   m_moment_op.compute( a_temperature, *this, PressureKernel<FArrayBox>(vpar) );
   for (CFG::DataIterator dit(a_temperature.dataIterator()); dit.ok(); ++dit) {
      a_temperature[dit].divide(density[dit]);
   }
}

void KineticSpecies::parallelTemperature( CFG::LevelData<CFG::FArrayBox>& a_temperature ) const
{
   CFG::LevelData<CFG::FArrayBox> density;
   density.define( a_temperature );
   m_moment_op.compute( density, *this, DensityKernel<FArrayBox>() );

   CFG::LevelData<CFG::FArrayBox> vpar;
   vpar.define( a_temperature );
   m_moment_op.compute( vpar, *this, ParallelVelKernel<FArrayBox>() );

   for (CFG::DataIterator dit(vpar.dataIterator()); dit.ok(); ++dit) {
      vpar[dit].divide(density[dit]);
   }
   
   m_moment_op.compute( a_temperature, *this, ParallelPressureKernel<FArrayBox>(vpar) );
   for (CFG::DataIterator dit(a_temperature.dataIterator()); dit.ok(); ++dit) {
      a_temperature[dit].divide(density[dit]);
   }
}

void KineticSpecies::perpTemperature( CFG::LevelData<CFG::FArrayBox>& a_temperature ) const
{
   CFG::LevelData<CFG::FArrayBox> density;
   density.define( a_temperature );
   m_moment_op.compute( density, *this, DensityKernel<FArrayBox>() );
   m_moment_op.compute( a_temperature, *this, PerpPressureKernel<FArrayBox>() );

   for (CFG::DataIterator dit(a_temperature.dataIterator()); dit.ok(); ++dit) {
      a_temperature[dit].divide(density[dit]);
   }
}

void KineticSpecies::energyDensity( CFG::LevelData<CFG::FArrayBox>& a_energy,
                                    const CFG::LevelData<CFG::FArrayBox>& a_phi) const
{
   m_moment_op.compute( a_energy, *this, EnergyKernel<FArrayBox>(a_phi) );
}

void KineticSpecies::kineticEnergyDensity( CFG::LevelData<CFG::FArrayBox>& a_energy) const
{
   m_moment_op.compute( a_energy, *this, KineticEnergyKernel<FArrayBox>() );
}

void KineticSpecies::kineticEnergyDensity( CFG::LevelData<CFG::FArrayBox>&  a_energy,
                                           const LevelData<FArrayBox>&      a_function) const
{
   m_moment_op.compute( a_energy, *this, a_function, KineticEnergyKernel<FArrayBox>() );
}

void KineticSpecies::perpEnergyDensity( CFG::LevelData<CFG::FArrayBox>& a_energy ) const
{
   m_moment_op.compute( a_energy, *this, PerpEnergyKernel<FArrayBox>() );
}

void KineticSpecies::parallelEnergyDensity( CFG::LevelData<CFG::FArrayBox>& a_energy ) const
{
   m_moment_op.compute( a_energy, *this, ParallelEnergyKernel<FArrayBox>() );
}

void KineticSpecies::fourthMoment( CFG::LevelData<CFG::FArrayBox>& a_fourth ) const
{
   m_moment_op.compute( a_fourth, *this, FourthMomentKernel<FArrayBox>() );
}


void
KineticSpecies::computeFSavgMaxwellian( LevelData<FArrayBox>&  a_F0 ) const
{
   
   /* Computes Vpar-unshifted FS averaged maxwellian fit. No B_star_par or J factors are included */
   const CFG::MagGeom& mag_geom( m_geometry->magGeom() );

   IntVect ghost_vect = a_F0.ghostVect();
   CFG::IntVect ghost_cfg;
   ghost_cfg[RADIAL_DIR] = ghost_vect[RADIAL_DIR];
   ghost_cfg[POLOIDAL_DIR] = ghost_vect[POLOIDAL_DIR];
   
   CFG::LevelData<CFG::FArrayBox> density( mag_geom.grids(), 1, CFG::IntVect::Zero );
   m_moment_op.compute( density, *this, DensityKernel<FArrayBox>() );

   CFG::LevelData<CFG::FArrayBox> parallelMom( mag_geom.grids(), 1, ghost_cfg );
   for (CFG::DataIterator dit(density.dataIterator()); dit.ok(); ++dit) {
      parallelMom[dit].setVal(0.0);
   }
   CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.grids(), 1, CFG::IntVect::Zero );
   m_moment_op.compute( pressure, *this, PressureKernel<FArrayBox>(parallelMom) );
   

   for (CFG::DataIterator dit(density.dataIterator()); dit.ok(); ++dit) {
      pressure[dit].divide(density[dit]);
   }
   
   CFG::FluxSurface fs(mag_geom);
   
   CFG::LevelData<CFG::FArrayBox> fs_average_density(mag_geom.grids(), 1, CFG::IntVect::Zero);
   fs.averageAndSpread(density, fs_average_density);
   
   CFG::LevelData<CFG::FArrayBox> fs_average_parMom(mag_geom.grids(), 1, CFG::IntVect::Zero);
   fs.averageAndSpread(parallelMom, fs_average_parMom);
   
   CFG::LevelData<CFG::FArrayBox> fs_average_temperature(mag_geom.grids(), 1, CFG::IntVect::Zero);
   fs.averageAndSpread(pressure, fs_average_temperature);
   
   
   MaxwellianKernel<FArrayBox> maxwellian(fs_average_density,
                                          fs_average_temperature,
                                          fs_average_parMom);

   const DisjointBoxLayout& dbl( a_F0.getBoxes() );
   LevelData<FArrayBox> dfn_no_ghost( dbl, 1, IntVect::Zero );

   maxwellian.eval(dfn_no_ghost, *this);

   for (DataIterator dit(a_F0.dataIterator()); dit.ok(); ++dit) {
     a_F0[dit].copy(dfn_no_ghost[dit]);
   }

   for (DataIterator dit(a_F0.dataIterator()); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& coord_sys = m_geometry->getBlockCoordSys(dbl[dit]);
      const ProblemDomain& domain = coord_sys.domain();
      fourthOrderCellExtrapAtDomainBdry(a_F0[dit], domain, dbl[dit]);
   }
   m_geometry->fillInternalGhosts(a_F0);

   m_geometry->divideBStarParallel( a_F0 );

}


bool KineticSpecies::isSpecies( const string& a_name ) const
{
   if (m_name == a_name) return true;
   return false;
}


const KineticSpecies& KineticSpecies::operator=( const KineticSpecies& a_rhs )
{
   if (&a_rhs != this)
   {
      m_name = a_rhs.m_name;
      m_mass = a_rhs.m_mass;
      m_charge = a_rhs.m_charge;
      m_is_gyrokinetic = a_rhs.m_is_gyrokinetic;
      m_gyroavg_op = a_rhs.m_gyroavg_op;
      m_dist_func.define( a_rhs.m_dist_func );
   }
   return *this;
}


void KineticSpecies::copy( const KineticSpecies& a_rhs )
{
   if (&a_rhs != this)
   {
      m_name = a_rhs.m_name;
      m_mass = a_rhs.m_mass;
      m_charge = a_rhs.m_charge;
      m_is_gyrokinetic = a_rhs.m_is_gyrokinetic;
      m_gyroavg_op = a_rhs.m_gyroavg_op;

      DataIterator dit( m_dist_func.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         m_dist_func[dit].copy( a_rhs.m_dist_func[dit] );
      }
   }
}


void KineticSpecies::zeroData()
{
   DataIterator dit( m_dist_func.dataIterator() );
   for (dit.begin(); dit.ok(); ++dit) {
      m_dist_func[dit].setVal( 0.0 );
   }
}


void KineticSpecies::addData( const KineticSpecies& a_rhs,
                              const Real a_factor )
{
   try {
      DataIterator dit( m_dist_func.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         m_dist_func[dit].plus( a_rhs.m_dist_func[dit], a_factor );
      }
   }
   catch (std::bad_cast) {
      MayDay::Error( "Invalid SpeciesModel passed to KineticSpecies::addData!" );
   }
}

bool KineticSpecies::conformsTo( const KineticSpecies& a_rhs,
                                 const bool a_include_ghost_cells ) const
{
   try {
      const LevelData<FArrayBox>& thisData = m_dist_func;
      const LevelData<FArrayBox>& rhsData = a_rhs.m_dist_func;

      const DisjointBoxLayout& thisBoxes = thisData.disjointBoxLayout();
      const DisjointBoxLayout& rhsBoxes = rhsData.disjointBoxLayout();

      bool status( true );
      status &= thisBoxes.compatible( rhsBoxes );
      status &= ( thisData.nComp() == rhsData.nComp() );

      if ( a_include_ghost_cells) {
         status &= ( thisData.ghostVect() == rhsData.ghostVect() );
      }

      return status;
   }
   catch (std::bad_cast) {
      MayDay::Error( "Invalid SpeciesModel passed to KineticSpecies::comformsTo!" );
   }
   return false;
}

RefCountedPtr<KineticSpecies>
KineticSpecies::clone( const IntVect ghostVect,
                       const bool copy_soln_data ) const
{
   RefCountedPtr<KineticSpecies> result
      = RefCountedPtr<KineticSpecies>(
              new KineticSpecies( m_name, m_mass, m_charge, m_geometry, m_is_gyrokinetic ) );

   result->m_dist_func.define( m_dist_func.disjointBoxLayout(),
                               m_dist_func.nComp(),
                               ghostVect );

   if (m_is_gyrokinetic) {
     result->gyroaverageOp(m_gyroavg_op);
   }

   if (copy_soln_data) {
      LevelData<FArrayBox>& result_dfn = result->m_dist_func;
      DataIterator dit( result_dfn.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         result_dfn[dit].copy( m_dist_func[dit] );
      }
   }

   return result;
}


const DisjointBoxLayout& KineticSpecies::getGhostDBL() const
{
   if ( m_ghost_dbl.size() == 0 ) {

      CH_assert(m_dist_func.isDefined());
      const DisjointBoxLayout& grids = m_dist_func.disjointBoxLayout();
      const IntVect& ghosts = m_dist_func.ghostVect();

      m_ghost_dbl.deepCopy(grids);

      const PhaseCoordSys& coord_sys = m_geometry->phaseCoordSys();

      for (int block=0; block<coord_sys.numBlocks(); ++block) {
         const PhaseBlockCoordSys& block_coord_sys = (const PhaseBlockCoordSys&)(*coord_sys.getCoordSys(block));
         const ProblemDomain& domain = block_coord_sys.domain();
         const Box& domain_box = domain.domainBox();

         for (int i=0; i<m_ghost_dbl.rawPtr()->size(); ++i) {
            Box& new_box = *const_cast<Box*>(&((*m_ghost_dbl.rawPtr())[i].box));
            Box original_box = new_box;

            if ( domain_box.contains(original_box) ) {
               for (int dir=0; dir<CFG_DIM; ++dir) {
                  if ( !domain.isPeriodic(dir) ) {
                     if ( original_box.smallEnd(dir) == domain_box.smallEnd(dir) ) {
                        new_box.growLo(dir,ghosts[dir]);
                     }
                     if ( original_box.bigEnd(dir) == domain_box.bigEnd(dir) ) {
                        new_box.growHi(dir,ghosts[dir]);
                     }
                  }
               }
            }
         }
      }

      m_ghost_dbl.closeNoSort();
   }

   return m_ghost_dbl;
}


Real KineticSpecies::maxValue() const
{
   const DisjointBoxLayout& grids( m_dist_func.getBoxes() );
   Real local_maximum( -CH_BADVAL );
   DataIterator dit( m_dist_func.dataIterator() );
   for (dit.begin(); dit.ok(); ++dit) {
      Box box( grids[dit] );
      Real box_max( m_dist_func[dit].max( box ) );
      local_maximum = Max( local_maximum, box_max );
   }

   Real maximum( local_maximum );
#ifdef CH_MPI
   MPI_Allreduce( &local_maximum, &maximum, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
#endif

   return maximum;
}

Real KineticSpecies::minValue() const
{
   const DisjointBoxLayout& grids( m_dist_func.getBoxes() );
   Real local_minimum( CH_BADVAL );
   IntVect local_index(IntVect::Zero);
   DataIterator dit( m_dist_func.dataIterator() );
   for (dit.begin(); dit.ok(); ++dit) {
      Box box( grids[dit] );
      Real box_min( m_dist_func[dit].min( box ) );
      IntVect box_min_loc(IntVect::Zero);
      if (box_min<0.0) {
         box_min_loc = m_dist_func[dit].minIndex( box );
         pout() << box_min_loc << ":  " << box_min << endl;
      }
      if (box_min < local_minimum) {
	local_minimum = box_min;
	local_index = box_min_loc;  
      }
   }

   Real minimum(local_minimum);
   Real minimum_global(local_minimum);
   IntVect index(local_index);
   struct {
     double val;
     int rank;
   } pair_in, pair_out;


#ifdef CH_MPI
   pair_in.val = local_minimum;
   pair_in.rank	= procID();

   MPI_Allreduce( &pair_in, &pair_out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD );
   minimum_global = pair_out.val;
   MPI_Bcast(index.dataPtr(), SpaceDim, MPI_INT, pair_out.rank, MPI_COMM_WORLD);
   
   MPI_Allreduce( &local_minimum, &minimum, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD );
   
#endif

   pout() << index << ":  " << minimum_global << endl;
   return minimum;
}


void KineticSpecies::computeVelocity(LevelData<FluxBox>&       a_velocity,
                                     const CFG::EMFields&      a_EMfields,
                                     const bool                a_include_FLR,
                                     const int                 a_velocity_option,
                                     const Real&               a_time,
                                     const bool                a_apply_axisymmetric_correction) const
{
   CH_TIMERS("KineticSpecies::computeVelocity");
   CH_TIMER("copy_velocity_fort", t_copy_velocity_fort);

   CH_assert(a_velocity.ghostVect() <= m_velocity.ghostVect());
   CH_assert(a_velocity.nComp() == m_velocity.nComp());
   
   m_geometry->updateVelocities(a_EMfields,
                                m_velocity,
                                a_velocity_option,
                                (isGyrokinetic() && a_include_FLR) );

   CH_START(t_copy_velocity_fort);
   for (DataIterator dit(m_velocity.dataIterator()); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         SpaceUtils::copy(a_velocity[dit][dir],m_velocity[dit][dir]);
      }
   }
   CH_STOP(t_copy_velocity_fort);

   // Correct for the extra R factor added to the 4th-order implementation
   // of metric coefficients to discretely preserve free-stream
   if ( a_apply_axisymmetric_correction) {
      m_geometry->applyAxisymmetricCorrection(a_velocity);
   }
   
   if ( !m_geometry->secondOrder() ) {
      // Convert face-centered values to face averages on valid cell faces.
      // The values on ghost cell faces remain second-order.
      fourthOrderAverage(a_velocity);
      a_velocity.exchange();
   }
}


void KineticSpecies::computeMappedVelocityNormals(LevelData<FluxBox>& a_velocity_normal,
                                                  const CFG::EMFields& a_EM_fields,
                                                  const bool a_fourth_order_Efield,
                                                  const int  a_velocity_option) const
{
   const DisjointBoxLayout& dbl( m_dist_func.getBoxes() );
   a_velocity_normal.define( dbl, 1, IntVect::Unit );
   m_geometry->updateVelocityNormals( a_EM_fields,
                                      a_fourth_order_Efield, 
                                      a_velocity_normal, 
                                      a_velocity_option );
}


void KineticSpecies::computeMappedVelocity(LevelData<FluxBox>& a_velocity,
                                           //                                           const LevelData<FluxBox>& a_E_field,
                                           const CFG::EMFields&      a_EMfields,
                                           const LevelData<FluxBox>& a_Apar_derivatives,
                                           const bool a_include_FLR,
                                           const Real&  a_time) const
{
   
   CH_TIME("KineticSpecies::computeMappedVelocity");

   computeVelocity( a_velocity, 
                    a_EMfields,
                    a_include_FLR, 
                    PhaseGeom::FULL_VELOCITY, 
                    a_time, 
                    false);
   m_geometry->multNTransposePointwise( a_velocity );

}


#include "NamespaceFooter.H"
