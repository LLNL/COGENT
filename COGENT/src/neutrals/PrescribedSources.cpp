#include <math.h>

#include "PrescribedSources.H"
#include "Directions.H"
#include "ConstFact.H"
#include "MomentOp.H"
#include "inspect.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FluxSurface.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H" 


PrescribedSources::PrescribedSources( ParmParse& a_ppntr, const int a_verbosity )
   : m_verbosity(a_verbosity),
     m_include_relaxation(false),
     m_include_thermalization(false),
     m_include_heat_src(false),
     m_it_counter(0),
     m_update_freq(1),
     m_diagnostics(false)

{
   parseParameters( a_ppntr );

   if (m_verbosity>0) {
      printParameters();
   }
   m_first_call = true;
}


PrescribedSources::~PrescribedSources()
{
}


void PrescribedSources::evalNtrRHS(KineticSpecies&                   a_rhs_species,
                                   const KineticSpeciesPtrVect&      a_soln,
                                   const CFG::FluidSpeciesPtrVect&   a_fluid_species_phys,
                                   const int                         a_species,
                                   const Real                        a_time )
{
  /*
    
   Add prescribed sources to the RHS
   
   (1) Heat source (from Sarazin et al Nucl. Fusion 50 (2010) 054004)
   
   (2) Relaxation source of the form nu(X)*(JBf-JBF_ref), where F_ref is
       a prescribed reference function

   (3) Thermalization source of the form nu(X)*(JBf-JBF_maxw), where F_maxw
       is a maxwellian with prescribed N and Vpar, and FS-averaged T value 
       that provides energy conservation
   */
   
  CH_TIME("PrescribedSources::evalNtrRHS");
  
  //Get geometry
  LevelData<FArrayBox>& rhs_dfn( a_rhs_species.distributionFunction() );
  const DisjointBoxLayout& grids = rhs_dfn.getBoxes();
  
  //Get normalization factors
  if (m_first_call) {
    getNormalizationFactors();
  }
  
  // Create source objects
  if (!m_heat_src.isDefined()) m_heat_src.define(grids, 1, IntVect::Zero);
  if (!m_relaxation_src.isDefined()) m_relaxation_src.define(grids, 1, IntVect::Zero);
  if (!m_thermalization_src.isDefined()) m_thermalization_src.define(grids, 1, IntVect::Zero);
  
  // Compute heat source
  Real mass = a_rhs_species.mass();
  if (m_include_heat_src) {
    computeHeatSrc(m_heat_src, a_rhs_species, mass, a_time);
  }
  
  // Compute thermalization source
  if (m_include_thermalization){
    const KineticSpecies& soln_species( *(a_soln[a_species]) );
    computeThermalizationSrc(m_thermalization_src, soln_species, a_time);
  }
  
  // Compute relaxation source
  if (m_include_relaxation){
    const KineticSpecies& soln_species( *(a_soln[a_species]) );
    computeRelaxationSrc(m_relaxation_src, soln_species, a_time);
  }
  
  //Add source terms to the GKSystem RHS
  for (DataIterator dit(rhs_dfn.dataIterator()); dit.ok(); ++dit) {
    if (m_include_heat_src)       rhs_dfn[dit].plus( m_heat_src[dit] );
    if (m_include_thermalization) rhs_dfn[dit].plus( m_thermalization_src[dit] );
    if (m_include_relaxation)     rhs_dfn[dit].plus( m_relaxation_src[dit] );
  }

  // Diagnostics
  if (m_first_call && m_diagnostics) {
    LevelData<FArrayBox> src(grids,1,IntVect::Zero);
    for (DataIterator dit(src.dataIterator()); dit.ok(); ++dit) {
      src[dit].setVal(0.);
      if (m_include_heat_src)       src[dit].plus(m_heat_src[dit]);
      if (m_include_thermalization) src[dit].plus( m_thermalization_src[dit] );
      if (m_include_relaxation)     src[dit].plus( m_relaxation_src[dit] );
    }
    diagnostics(src, a_rhs_species, a_time);
  }
  
  m_first_call = false;
}


void PrescribedSources::computeHeatSrc(LevelData<FArrayBox>&  a_src,
                                       const KineticSpecies&  a_kinetic_species,
                                       const Real&            a_mass,
                                       const Real&            a_time)
{
  /*
   Computes heat source (from Sarazin et al Nucl. Fusion 50 (2010) 054004)
  */
  
  if (m_first_call) {
    
    // Get spatial profile function
    const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
    const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
    CFG::LevelData<CFG::FArrayBox> spatial_profile_cfg(mag_geom.grids(), 1, CFG::IntVect::Zero);
    m_heat_src_profile_func->assign( spatial_profile_cfg, mag_geom, a_time);
    
    // Get spatial profile normalization
    Real spatial_norm = mag_geom.volumeIntegrate(spatial_profile_cfg);
  
    // Inject spatial profile into the phase space
    phase_geom.injectConfigurationToPhase(spatial_profile_cfg, m_heat_src_profile_inj);

    // Get injected B
    const LevelData<FArrayBox>& B_injected = phase_geom.getBFieldMagnitude();
        
    // Compute phase shape of the heat source (NB:: doing so in cpp is not
    // as efficient as in Fortran, but we only do it once)
    const DisjointBoxLayout& grids = a_src.disjointBoxLayout();
    
    // Create a temporary
    LevelData<FArrayBox> src_tmp(grids, 2, IntVect::Zero);

    // Iterate over patches
    for (DataIterator dit(a_src.dataIterator()); dit.ok(); ++dit) {
      
      FArrayBox& this_src = src_tmp[dit];
      const FArrayBox& this_B = B_injected[dit];
      const FArrayBox& this_spatial_profile = m_heat_src_profile_inj[dit];
      
      const Box& Bbox = this_B.box();
      int vp_index = Bbox.smallEnd(VPARALLEL_DIR);
      int mu_index = Bbox.smallEnd(MU_DIR);

      // Get the physical velocity coordinates for this part of phase space
      const FArrayBox& velocityRealCoords = phase_geom.getVelocityRealCoords()[dit];

      // Interate over box indicies
      BoxIterator bit(this_src.box());
      for (bit.begin(); bit.ok(); ++bit) {
        IntVect iv = bit();
        IntVect iv_inj = iv;
        iv_inj[VPARALLEL_DIR] = vp_index;
        iv_inj[MU_DIR] = mu_index;
         
        double v_parallel = velocityRealCoords(iv,0);
        double mu = velocityRealCoords(iv,1);
        double v_perp2 =  mu * this_B(iv_inj) ;
        double energy = 0.5 * (a_mass * v_parallel * v_parallel + v_perp2);
        double energy_norm = energy / m_T_src;

        double fac = this_spatial_profile(iv_inj) / spatial_norm;
        fac *= pow(a_mass, 3.0/2.0) / pow(m_T_src, 5.0/2.0);
        fac *= m_power * m_heat_src_norm;
        
        this_src(iv,0) = fac * energy_norm * exp(-energy_norm);
        this_src(iv,1) = fac * 3.0/2.0 * exp(-energy_norm);
      }
    }
    
    phase_geom.multBStarParallel(src_tmp);
     
    // Renormalize src_tmp component_1 to enforce particle conservation
    // (we can do this more elegantly later, by casting the source expression to conservative form)
    CFG::LevelData<CFG::FArrayBox> particle_moment_cfg( mag_geom.grids(), 2, CFG::IntVect::Zero );
    
    MomentOp& moment_op = MomentOp::instance();
    moment_op.compute(particle_moment_cfg, a_kinetic_species, src_tmp, DensityKernel<FArrayBox>());

    LevelData<FArrayBox> particle_moment;
    phase_geom.injectAndExpandConfigurationToPhase(particle_moment_cfg, particle_moment);

    for (DataIterator dit(a_src.dataIterator()); dit.ok(); ++dit) {
      BoxIterator bit(a_src[dit].box());
      for (bit.begin(); bit.ok(); ++bit) {
        IntVect iv = bit();
        double denom = particle_moment[dit](iv,1);
        if (denom > 1.0e-10) {
          src_tmp[dit](iv,1) *= particle_moment[dit](iv,0)/denom;
        }
        a_src[dit](iv) = src_tmp[dit](iv,0) - src_tmp[dit](iv,1);
      }
    }
    phase_geom.multJonValid(a_src);
  }
}

void PrescribedSources::computeThermalizationSrc(LevelData<FArrayBox>&  a_src,
                                                 const KineticSpecies&  a_kinetic_species,
                                                 const Real&            a_time)
{
  
  // Get geometry objects
  const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
  const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
  const DisjointBoxLayout& dbl = a_src.getBoxes();
  
  if (m_first_call) {
    
    // Get src spatial profiles
    CFG::LevelData<CFG::FArrayBox> spatial_profile_cfg(mag_geom.grids(), 1, CFG::IntVect::Zero);
    m_thermalization_profile_func->assign( spatial_profile_cfg, mag_geom, a_time);
    phase_geom.injectAndExpandConfigurationToPhase(spatial_profile_cfg, m_thermalization_profile_inj);
    
    m_maxw_density.define(mag_geom.grids(), 1, CFG::IntVect::Zero);
    m_maxw_density_func->assign( m_maxw_density, mag_geom, a_time);

    m_maxw_vparallel.define(mag_geom.grids(), 1, CFG::IntVect::Zero);
    m_maxw_vparallel_func->assign( m_maxw_vparallel, mag_geom, a_time);
    
    m_maxw_dfn.define(dbl, 1, IntVect::Zero);
  }

  if ((m_it_counter % m_update_freq == 0) || m_first_call) {
    CFG::LevelData<CFG::FArrayBox> kin_energ(mag_geom.grids(), 1, CFG::IntVect::Zero);
    a_kinetic_species.kineticEnergyMoment(kin_energ);

    CFG::FluxSurface flux_surface(mag_geom);
    CFG::LevelData<CFG::FArrayBox> fs_kin_energ(mag_geom.grids(), 1, CFG::IntVect::Zero);
    flux_surface.averageAndSpread(kin_energ, fs_kin_energ);
  
    CFG::LevelData<CFG::FArrayBox> fs_average_temperature(mag_geom.grids(), 1, CFG::IntVect::Zero);

    CFG::DataIterator cfg_dit = m_maxw_density.dataIterator();
    for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
      fs_average_temperature[cfg_dit].copy(fs_kin_energ[cfg_dit]);
      fs_average_temperature[cfg_dit].divide(m_maxw_density[cfg_dit]);
      fs_average_temperature[cfg_dit].mult(2.0/3.0);
    }

    MaxwellianKernel<FArrayBox> maxwellian(m_maxw_density, fs_average_temperature, m_maxw_vparallel);
    maxwellian.eval(m_maxw_dfn, a_kinetic_species);

    //Enforce energy conservation
    CFG::LevelData<CFG::FArrayBox> maxw_kin_energ(mag_geom.grids(), 1, CFG::IntVect::Zero);
    CFG::LevelData<CFG::FArrayBox> fs_maxw_kin_energ(mag_geom.grids(), 1, CFG::IntVect::Zero);
    MomentOp& moment_op = MomentOp::instance();
    moment_op.compute(maxw_kin_energ, a_kinetic_species, m_maxw_dfn, KineticEnergyKernel<FArrayBox>());
    flux_surface.averageAndSpread(maxw_kin_energ, fs_maxw_kin_energ);

    LevelData<FArrayBox> fs_maxw_kin_energ_inj;
    LevelData<FArrayBox> fs_kin_energ_inj;
    phase_geom.injectAndExpandConfigurationToPhase(fs_maxw_kin_energ,fs_maxw_kin_energ_inj);
    phase_geom.injectAndExpandConfigurationToPhase(fs_kin_energ,fs_kin_energ_inj);
    for (DataIterator dit(a_src.dataIterator()); dit.ok(); ++dit) {
      m_maxw_dfn[dit].mult(fs_kin_energ_inj[dit]);
      m_maxw_dfn[dit].divide(fs_maxw_kin_energ_inj[dit]);
    }
  
    phase_geom.multJonValid(m_maxw_dfn);
  }
  
  const LevelData<FArrayBox>& soln_dfn( a_kinetic_species.distributionFunction() );
  for (DataIterator dit(a_src.dataIterator()); dit.ok(); ++dit) {
    a_src[dit].copy(soln_dfn[dit]);
  }
  phase_geom.multJonValid(a_src);

  for (DataIterator dit(a_src.dataIterator()); dit.ok(); ++dit) {
    a_src[dit].minus(m_maxw_dfn[dit]);
    a_src[dit].mult(m_thermalization_profile_inj[dit]);
    a_src[dit].mult(-1.0);
  }
}

void PrescribedSources::computeRelaxationSrc(LevelData<FArrayBox>&  a_src,
                                             const KineticSpecies&  a_kinetic_species,
                                             const Real&            a_time)
{
  
  // Get geometry objects
  const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
  const DisjointBoxLayout& dbl = a_src.getBoxes();
  
  if (m_first_call) {
    
    //Define the reference function object
    m_ref_dfn.define(dbl, 1, IntVect::Zero);
    
    // Get src spatial profiles
    const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
    CFG::LevelData<CFG::FArrayBox> spatial_profile_cfg(mag_geom.grids(), 1, CFG::IntVect::Zero);
    m_relaxation_profile_func->assign( spatial_profile_cfg, mag_geom, a_time);
    phase_geom.injectAndExpandConfigurationToPhase(spatial_profile_cfg, m_relaxation_profile_inj);

    // Create reference (J*Bstar_par*ref_function) distribution if defined
    KineticSpeciesPtr ref_species( a_kinetic_species.clone( IntVect::Zero, false ) );
    m_ref_func->assign( *ref_species, a_time );
    const LevelData<FArrayBox>& ref_dfn( ref_species->distributionFunction() );
    m_ref_dfn.define(ref_dfn);
  }
  
  const LevelData<FArrayBox>& soln_dfn( a_kinetic_species.distributionFunction() );
  for (DataIterator dit(a_src.dataIterator()); dit.ok(); ++dit) {
    a_src[dit].copy(soln_dfn[dit]);
  }
  phase_geom.multJonValid(a_src);

  for (DataIterator dit(a_src.dataIterator()); dit.ok(); ++dit) {
    a_src[dit].minus(m_ref_dfn[dit]);
    a_src[dit].mult(m_relaxation_profile_inj[dit]);
    a_src[dit].mult(-1.0);
  }
}

void PrescribedSources::getNormalizationFactors()

{
  //Get normalization parameters (units)
  double N, T, L, M;
  ParmParse ppunits( "units" );
  ppunits.get("number_density",N);  //[m^{-3}]
  ppunits.get("temperature",T);     //[eV]
  ppunits.get("length",L);          //[m]
  ppunits.get("mass",M);            //[a.u.]
  
  double pi = Constants::PI;
  
  // Get heat source normalization for the power input in MW
  m_heat_src_norm = 6.38e20/(3.0*sqrt(2.0*pi));
  m_heat_src_norm *= sqrt(M)/(pow(L,2)*pow(T,3.0/2.0)*N);
}

inline
void PrescribedSources::parseParameters( ParmParse& a_pp )
{
  KineticFunctionLibrary* phase_library = KineticFunctionLibrary::getInstance();
  CFG::GridFunctionLibrary* grid_library = CFG::GridFunctionLibrary::getInstance();
  std::string function_name;
  
  a_pp.query( "include_relaxation", m_include_relaxation );
  a_pp.query( "include_heat_source", m_include_heat_src );
  a_pp.query( "include_thermalization", m_include_thermalization );

  a_pp.query("update_cls_freq", m_update_freq );
  
  if (m_include_relaxation) {
    
    if (a_pp.contains("ref_function")) {
      a_pp.get("ref_function", function_name );
      m_ref_func = phase_library->find( function_name );
    }
    else {
      MayDay::Error("PrescribedSources:: either ref_function must be specified for the relaxation model ");
    }
    if (a_pp.contains("relaxation_profile_function")) {
      a_pp.get("relaxation_profile_function", function_name );
      m_relaxation_profile_func = grid_library->find( function_name );
    }
    else {
      MayDay::Error("PrescribedSources:: relaxation profile must be specified ");
    }
  }
  
  if (m_include_thermalization) {
    
    if (a_pp.contains("maxw_density") && a_pp.contains("maxw_vparallel")) {
      a_pp.get("maxw_density", function_name );
      m_maxw_density_func = grid_library->find( function_name );

      a_pp.get("maxw_vparallel", function_name );
      m_maxw_vparallel_func = grid_library->find( function_name );
    }
    else {
      MayDay::Error("PrescribedSources:: ref_density and ref_vparallel must be specified for the thermalization model ");
    }

    if (a_pp.contains("thermalization_profile_function")) {
      a_pp.get("thermalization_profile_function", function_name );
      m_thermalization_profile_func = grid_library->find( function_name );
    }
    else {
      MayDay::Error("PrescribedSources:: thermalization profile must be specified ");
    }
  }
  
  if (m_include_heat_src) {
    
    if (a_pp.contains("heat_source_profile_function")) {
      a_pp.get("heat_source_profile_function", function_name );
      m_heat_src_profile_func = grid_library->find( function_name );
    }
    else {
      MayDay::Error("PrescribedSources:: heat source spatial profile must be specified ");
    }
    
    if (a_pp.contains("source_temperature")) {
      a_pp.get( "source_temperature", m_T_src );
    }
    else{
      MayDay::Error("PrescribedSources:: source temperature must be specified ");
    }

    if (a_pp.contains("source_power_MW")) {
      a_pp.get( "source_power_MW", m_power );
    }
    else{
      MayDay::Error("PrescribedSources:: source power must be specified ");
    }
  }
  
  a_pp.query( "plot_src_data", m_diagnostics );
}


inline
void PrescribedSources::printParameters()
{
   if (procID()==0) {
      std::cout << "PrescribedSources neutral parameters:" << std::endl;
      if (m_include_heat_src) {
	std::cout << "  source_power_MW  =  " << m_power;
	std::cout << "  source_temperature = " << m_T_src;
	std::cout << "  heat source spatial function:" << std::endl;
	m_heat_src_profile_func->printParameters();
      }
      if (m_include_relaxation) m_relaxation_profile_func->printParameters();
      if (m_include_thermalization) m_thermalization_profile_func->printParameters();
   }
}

Real PrescribedSources::computeDtExplicitTI(const KineticSpeciesPtrVect& a_soln, const int a_species)
{
  return TimeScale(a_soln, a_species);
}

Real PrescribedSources::computeDtImExTI(const KineticSpeciesPtrVect& a_soln, const int a_species)
{
  return TimeScale(a_soln, a_species);
}

Real PrescribedSources::TimeScale(const KineticSpeciesPtrVect& a_soln, const int a_species)
{
  return DBL_MAX;
}

void
PrescribedSources::preTimeStep(const KineticSpeciesPtrVect& a_soln_mapped,
                               const int a_species,
                               const Real a_time,
                               const KineticSpeciesPtrVect& a_soln_physical )
{
   m_it_counter+=1;
}

void PrescribedSources::diagnostics(const LevelData<FArrayBox>& a_rhs,
                                    const KineticSpecies&       a_rhs_species,
                                    const double                a_time) const
{
  //Get geometry
  const PhaseGeom& phase_geom = a_rhs_species.phaseSpaceGeometry();
  const CFG::MagGeom& mag_geom( phase_geom.magGeom() );

  //Get moment operator
  MomentOp& moment_op = MomentOp::instance();
  
  //Plot particle source
  CFG::LevelData<CFG::FArrayBox> particle_src( mag_geom.grids(), 1, CFG::IntVect::Zero );
  moment_op.compute(particle_src, a_rhs_species, a_rhs, DensityKernel<FArrayBox>());
  phase_geom.plotConfigurationData( "particle_src", particle_src, a_time );
  
  //Plot parallel momentum source
  CFG::LevelData<CFG::FArrayBox> parMom_src( mag_geom.grids(), 1, CFG::IntVect::Zero );
  moment_op.compute(parMom_src, a_rhs_species, a_rhs, ParallelMomKernel<FArrayBox>());
  phase_geom.plotConfigurationData( "parMom_src", parMom_src, a_time );
  
  //Plot energy source
  CFG::LevelData<CFG::FArrayBox> energy_src( mag_geom.grids(), 1, CFG::IntVect::Zero );
  moment_op.compute(energy_src, a_rhs_species, a_rhs, KineticEnergyKernel<FArrayBox>());

  CFG::LevelData<CFG::FArrayBox> energy_src_fs( mag_geom.grids(), 1, CFG::IntVect::Zero );
  CFG::FluxSurface flux_surface(mag_geom);
  mag_geom.divideJonValid(energy_src);
  flux_surface.averageAndSpread(energy_src, energy_src_fs);
  
  phase_geom.plotConfigurationData( "energy_src_fs", energy_src_fs, a_time );
}

#include "NamespaceFooter.H"
