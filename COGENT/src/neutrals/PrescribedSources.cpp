#include <math.h>

#include "PrescribedSources.H"
#include "Directions.H"
#include "ConstFact.H"
#include "MomentOp.H"
#include "inspect.H"

#include "NamespaceHeader.H" 


PrescribedSources::PrescribedSources( ParmParse& a_ppntr, const int a_verbosity )
   : m_verbosity(a_verbosity),
     m_thermal_bath(false),
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

   Heat source (from Sarazin et al Nucl. Fusion 50 (2010) 054004)
    
   */
   
  CH_TIME("PrescribedSources::evalNtrRHS");
  
  //Get geometry
  LevelData<FArrayBox>& rhs_dfn( a_rhs_species.distributionFunction() );
  const DisjointBoxLayout& grids = rhs_dfn.getBoxes();
  
  //Get normalization factors
  if (m_first_call) {
    getNormalizationFactors();
  }
  
  // Compute Heat source
  if (!m_heat_src.isDefined()) m_heat_src.define(grids, 1, IntVect::Zero);
  Real mass = a_rhs_species.mass();
  if (!m_thermal_bath) {
    computeHeatSrc(m_heat_src, a_rhs_species, mass, a_time);
  }
  else {
    const KineticSpecies& soln_species( *(a_soln[a_species]) );
    computeThermalBathSrc(m_heat_src, soln_species, a_time);
  }
  
  
  //Add energy source term the GKSystem RHS
  for (DataIterator dit(rhs_dfn.dataIterator()); dit.ok(); ++dit) {
      rhs_dfn[dit].plus( m_heat_src[dit] );
  }
  
  if (m_first_call && m_diagnostics) {
    diagnostics(m_heat_src, a_rhs_species, a_time);
  }
  
  m_first_call = false;
}


void PrescribedSources::computeHeatSrc(LevelData<FArrayBox>&  a_heat_src,
                                       const KineticSpecies&  a_kinetic_species,
                                       const Real&            a_mass,
                                       const Real&            a_time)
{
  
  if (m_first_call) {
    
    // Get spatial profile function
    const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
    const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
    CFG::LevelData<CFG::FArrayBox> spatial_profile_cfg(mag_geom.grids(), 1, CFG::IntVect::Zero);
    m_heat_src_spatial_func->assign( spatial_profile_cfg, mag_geom, a_time);
    
    // Get spatial profile normalization
    Real spatial_norm = mag_geom.volumeIntegrate(spatial_profile_cfg);
  
    // Inject spatial profile into the phase space
    phase_geom.injectConfigurationToPhase(spatial_profile_cfg, m_spatial_profile_inj);

    // Get injected B
    const LevelData<FArrayBox>& B_injected = phase_geom.getBFieldMagnitude();
        
    // Compute phase shape of the heat source (NB:: doing so in cpp is not
    // as efficient as in Fortran, but we only do it once)
    const DisjointBoxLayout& grids = a_heat_src.disjointBoxLayout();
    
    // Create a temporary
    LevelData<FArrayBox> src_tmp(grids, 2, IntVect::Zero);

    // Iterate over patches
    for (DataIterator dit(a_heat_src.dataIterator()); dit.ok(); ++dit) {
      
      FArrayBox& this_src = src_tmp[dit];
      const FArrayBox& this_B = B_injected[dit];
      const FArrayBox& this_spatial_profile = m_spatial_profile_inj[dit];
      
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

    for (DataIterator dit(a_heat_src.dataIterator()); dit.ok(); ++dit) {
      BoxIterator bit(a_heat_src[dit].box());
      for (bit.begin(); bit.ok(); ++bit) {
        IntVect iv = bit();
        double denom = particle_moment[dit](iv,1);
        if (denom > 1.0e-10) {
          src_tmp[dit](iv,1) *= particle_moment[dit](iv,0)/denom;
        }
        a_heat_src[dit](iv) = src_tmp[dit](iv,0) - src_tmp[dit](iv,1);
      }
    }
    phase_geom.multJonValid(a_heat_src);
  }
}

void PrescribedSources::computeThermalBathSrc(LevelData<FArrayBox>&  a_heat_src,
                                              const KineticSpecies&  a_kinetic_species,
                                              const Real&            a_time)
{
  
  const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
  
  if (m_first_call) {
    // Get spatial profile src function
    const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
    CFG::LevelData<CFG::FArrayBox> spatial_profile_cfg(mag_geom.grids(), 1, CFG::IntVect::Zero);
    m_heat_src_spatial_func->assign( spatial_profile_cfg, mag_geom, a_time);
    phase_geom.injectAndExpandConfigurationToPhase(spatial_profile_cfg, m_spatial_profile_inj);

    // Create reference (J*Bstar_par*ref_function) distribution
    KineticSpeciesPtr ref_species( a_kinetic_species.clone( IntVect::Zero, false ) );
    m_ref_func->assign( *ref_species, a_time );
    const LevelData<FArrayBox>& ref_dfn( ref_species->distributionFunction() );
    m_ref_dfn.define(ref_dfn);
  }
  
  const LevelData<FArrayBox>& soln_dfn( a_kinetic_species.distributionFunction() );
  for (DataIterator dit(a_heat_src.dataIterator()); dit.ok(); ++dit) {
    a_heat_src[dit].copy(soln_dfn[dit]);
  }
  phase_geom.multJonValid(a_heat_src);

  for (DataIterator dit(a_heat_src.dataIterator()); dit.ok(); ++dit) {
    a_heat_src[dit].minus(m_ref_dfn[dit]);
    a_heat_src[dit].mult(m_spatial_profile_inj[dit]);
    a_heat_src[dit].mult(-1.0);
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

  CFG::GridFunctionLibrary* grid_library = CFG::GridFunctionLibrary::getInstance();
  std::string grid_function_name;
  
  a_pp.query( "thermal_bath", m_thermal_bath );
   
   
  if (a_pp.contains("heat_source_spatial_function")) {
    a_pp.get("heat_source_spatial_function", grid_function_name );
    m_heat_src_spatial_func = grid_library->find( grid_function_name );
  }
  else {
    MayDay::Error("PrescribedSources:: heat source spatial profile must be specified ");
  }
  
  if (m_thermal_bath) {
    KineticFunctionLibrary* library = KineticFunctionLibrary::getInstance();
    std::string function_name;

    if (a_pp.contains("ref_function")) {
      a_pp.get("ref_function", function_name );
      m_ref_func = library->find( function_name );
    }
    else {
      MayDay::Error("PrescribedSources:: ref function must be specified for thermal bath model ");
    }
  }
  
  else {
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
      std::cout << "  source_power_MW  =  " << m_power;
      std::cout << "  source_temperature = " << m_T_src;
      std::cout << "  heat source spatial function:" << std::endl;
      m_heat_src_spatial_func->printParameters();
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
  phase_geom.plotConfigurationData( "energy_src", energy_src, a_time );
}

#include "NamespaceFooter.H"
