#include <math.h>

#include "ConsDragDiff.H"
#include "Directions.H"
#include "PhaseGeom.H"
#include "PhaseBlockCoordSys.H"
#include "ConstFact.H"
#include "CollisionsF_F.H"

#include "MomentOp.H"
#include "Kernels.H"

#include "NamespaceHeader.H" //Should be the last one


ConsDragDiff::ConsDragDiff( const string& a_species_name, const string& a_ppcls_str, const int a_verbosity )
    : m_cls_only(false),
      moment_op( MomentOp::instance() ),
      m_first_call(true),
      m_it_counter(0),
      m_update_freq(1),
      m_time_implicit(true),
      m_diagnostics(false),
      m_verbosity(a_verbosity),
      m_cls_freq_func(NULL)
{
   m_ppcls = ParmParse(a_ppcls_str.c_str());
   m_species_name = a_species_name;
  
   ParseParameters(m_ppcls);
}

ConsDragDiff::~ConsDragDiff()
{
}

void ConsDragDiff::evalClsRHS( KineticSpeciesPtrVect&        a_rhs,
                               const KineticSpeciesPtrVect&  a_soln,
                               const int                     a_species,
                               const int                     a_species_bkgr,
                               const Real                    a_time )
{
  /*
    Evaluates the like-like species fully conservative drag-diffusion
    velocity space collision operator (a.k.a. Lenard & Bernstain 1958 model):
    df/dt_coll=div(flux_coll), where flux_coll = coll_freq*[(v-U)*f+vth*dfdv].
    Note that this operator exactly conservs density, momentum,
    and thermal energy while also relaxing the species to a Maxwellian state.
  */

  CH_TIME("ConsDragDiff::evalClsRHS");
  
  // get vlasov RHS for the current species
  KineticSpecies& rhs_species( *(a_rhs[a_species]) );
  LevelData<FArrayBox>& rhs_dfn = rhs_species.distributionFunction();

  // get solution distribution function (f*J*Bstarpar) for the current species
  const KineticSpecies& soln_species( *(a_soln[a_species]) );
  const LevelData<FArrayBox>& soln_fBJ = soln_species.distributionFunction();
  double mass = soln_species.mass();
  const DisjointBoxLayout& dbl = soln_fBJ.getBoxes();

  // get coordinate system parameters
  const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();
  const CFG::MagGeom & mag_geom = phase_geom.magGeom();
  
  // copy soln_fBJ so can perform exchange to ghost cells
  // we need to do so because we pass here a computational dfn
  // that does not have ghost cells filled
  if (!m_fBJ_vel_ghost.isDefined()) {
    IntVect velGhost = 2*IntVect::Unit;
    for (int dir=0; dir<CFG_DIM; dir++) {
      velGhost[dir] = 0;
    }
    m_fBJ_vel_ghost.define(dbl, 1, velGhost);
  }

  for (DataIterator dit(soln_fBJ.dataIterator()); dit.ok(); ++dit) {
    m_fBJ_vel_ghost[dit].copy(soln_fBJ[dit], dbl[dit]);
  }
  //since we only need ghost information in velocity space only
  // use simple exchange, instead of fillInternalGhosts
  m_fBJ_vel_ghost.exchange();
  
  // get injected magnetic field
  const LevelData<FArrayBox>& inj_B = phase_geom.getBFieldMagnitude();
  
  // get injected cell-centered velocity spatial normalization
  const int use_spatial_vel_norm = phase_geom.spatialVelNorm();
  const LevelData<FArrayBox>& inj_vel_norm = phase_geom.getVelNorm();
      
  // get problem domain and number of vpar and mu cells
  const ProblemDomain& phase_domain = phase_geom.domain();
  const Box& domain_box = phase_domain.domainBox();
  int num_vp_cells = domain_box.size(VPARALLEL_DIR);
  int num_mu_cells = domain_box.size(MU_DIR);

  // create the temporary fluxes needed to compute the conservative
  // mean velocity and conservative temperature
  LevelData<FluxBox> fluxes(dbl, 3, IntVect::Zero);
  DataIterator dit = fluxes.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    // get phase space dx
    const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(dbl[dit]);
    const RealVect& phase_dx =  block_coord_sys.dx();
    
    const FArrayBox& fBJ_on_patch = m_fBJ_vel_ghost[dit];
    const FArrayBox& B_on_patch   = inj_B[dit];
        
    const FArrayBox* this_velnormptr;

    // Create a dummy to be used in Fortran if no velocity normalization
    FArrayBox dummy(B_on_patch.box(),2);
                
    if (use_spatial_vel_norm == 1) {
      this_velnormptr = &(inj_vel_norm[dit]);
    }
    else {
      this_velnormptr = &dummy;
    }
        
    for (int dir=0; dir<SpaceDim; dir++)
    {
      FORT_EVAL_CONSDRAGDIFF_FLUX(CHF_BOX(fluxes[dit][dir].box()),
                                  CHF_FRA(fluxes[dit][dir]),
                                  CHF_CONST_FRA1(fBJ_on_patch,0),
                                  CHF_CONST_FRA1(B_on_patch,0),
                                  CHF_CONST_FRA((*this_velnormptr)),
                                  CHF_CONST_REAL(mass),
                                  CHF_CONST_REALVECT(phase_dx),
                                  CHF_CONST_INT(dir),
                                  CHF_CONST_INT(num_vp_cells),
                                  CHF_CONST_INT(num_mu_cells),
                                  CHF_CONST_INT(use_spatial_vel_norm));

    }
  }

  // compute the divergence of the velocity space fluxes
  LevelData<FArrayBox> Jpsi(dbl, 3, IntVect::Zero);
  phase_geom.mappedGridDivergenceFromFluxNormals(Jpsi, fluxes);

  // take momentum needed to compute conservative Upar and Temperature
  CFG::LevelData<CFG::FArrayBox> vpar_moms(mag_geom.grids(), 3, CFG::IntVect::Zero);
  CFG::LevelData<CFG::FArrayBox> pres_moms(mag_geom.grids(), 3, CFG::IntVect::Zero);

  //need actual Upar for pressure moment
  if (!m_dens.isDefined()) m_dens.define(mag_geom.grids(), 1, CFG::IntVect::Zero);
  if (!m_ushift.isDefined()) m_ushift.define(mag_geom.grids(), 1, CFG::IntVect::Zero);

  // optimize performace by cashing m_ushift variable and only updating it
  // every so often. In the results we degrade energy conservation a bit,
  // becasue the mean velocity used in the calculation of the distribution temeprature
  // will be slightly out of sink with the m_ushift variable that is used to evaluate
  // the "conservative" value of the temeprature (see Justin's thesis for details).
  // This should not affect momentum conservation.
  if ((m_it_counter % m_update_freq == 0) || m_first_call) {
     soln_species.numberDensity(m_dens);
     soln_species.parallelParticleFlux(m_ushift);
     CFG::DataIterator cfg_dit = m_ushift.dataIterator();
     for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit)
     {
        m_ushift[cfg_dit].divide(m_dens[cfg_dit]);
     }
  }

  moment_op.compute(vpar_moms, soln_species, Jpsi, ParallelVelKernel<FArrayBox>());
  moment_op.compute(pres_moms, soln_species, Jpsi, PressureKernel<FArrayBox>(m_ushift));
  LevelData<FArrayBox> inj_vpar_moms;
  LevelData<FArrayBox> inj_pres_moms;
  phase_geom.injectConfigurationToPhase(vpar_moms, inj_vpar_moms);
  phase_geom.injectConfigurationToPhase(pres_moms, inj_pres_moms);

  // create conservative mean parallel velocity and Temperature
  LevelData<FArrayBox> Upar(inj_vpar_moms.getBoxes(), 1, IntVect::Zero);
  LevelData<FArrayBox> Temp(inj_vpar_moms.getBoxes(), 1, IntVect::Zero);
  DataIterator mdit = Upar.dataIterator();
  for (mdit.begin(); mdit.ok(); ++mdit)
  {
    const FArrayBox& this_vpar_moms = inj_vpar_moms[mdit];
    const FArrayBox& this_pres_moms = inj_pres_moms[mdit];

    FORT_EVAL_CONS_UPAR_TEMP( CHF_BOX(Upar[mdit].box()),
                              CHF_FRA1(Upar[mdit],0),
                              CHF_FRA1(Temp[mdit],0),
                              CHF_CONST_FRA(this_vpar_moms),
                              CHF_CONST_FRA(this_pres_moms) );
  }

  
  // get collision frequency
  if (m_first_call) m_cls_freq.define(dbl, 1, IntVect::Zero);
    
  if (m_cls_freq_func != NULL && m_first_call) {
    CFG::LevelData<CFG::FArrayBox> cls_freq_cfg(mag_geom.grids(), 1, CFG::IntVect::Zero);
    m_cls_freq_func->assign( cls_freq_cfg, mag_geom, a_time);
    phase_geom.injectAndExpandConfigurationToPhase(cls_freq_cfg, m_cls_freq);
  }
  
  if (m_cls_freq_func == NULL) {
    if ((m_it_counter % m_update_freq == 0) || m_first_call) {
      mag_geom.divideJonValid(m_dens);
      LevelData<FArrayBox> dens_inj;
      phase_geom.injectConfigurationToPhase(m_dens, dens_inj);

      CFG::LevelData<CFG::FArrayBox> temperature(mag_geom.grids(), 1, CFG::IntVect::Zero);
      soln_species.temperature(temperature);
      mag_geom.divideJonValid(temperature);
      
      LevelData<FArrayBox> temp_inj;
      phase_geom.injectConfigurationToPhase(temperature, temp_inj);
      
      double mass = soln_species.mass();
      double charge = soln_species.charge();
          
      computeSelfConsistFreq(m_cls_freq, dens_inj, temp_inj, mass, charge);
    }
  }
    
  // create conservative drag-diffusion collisional RHS
  LevelData<FArrayBox> rhs_cls(dbl, 1, IntVect::Zero);
  DataIterator cdit = rhs_cls.dataIterator();
  for (cdit.begin(); cdit.ok(); ++cdit)
  {
    const FArrayBox& this_Upar = Upar[cdit];
    const FArrayBox& this_Temp = Temp[cdit];
    const FArrayBox& this_Jpsi = Jpsi[cdit];

    FORT_EVAL_CONSDRAGDIFF_RHS( CHF_BOX(rhs_cls[cdit].box()),
                                CHF_FRA1(rhs_cls[cdit],0),
                                CHF_CONST_FRA1(m_cls_freq[cdit],0),
                                CHF_CONST_FRA1(this_Upar,0),
                                CHF_CONST_FRA1(this_Temp,0),
                                CHF_CONST_FRA(this_Jpsi) );
  }

  // add (or overwrite) collision RHS to Vlasov RHS
  for (cdit.begin(); cdit.ok(); ++cdit)
  {
    if(m_cls_only){
      rhs_dfn[cdit].copy(rhs_cls[cdit]);}
    else{
      rhs_dfn[cdit].plus(rhs_cls[cdit]);}
  }

  // print parameters (e.g., stable dt) at the first time step
  if ((m_verbosity) && (m_first_call)) {printParameters(a_soln);}

   
  if (m_first_call && m_diagnostics) {
    diagnostics(rhs_cls, rhs_species, a_time);
    exit(1);
  }
  
  m_first_call = false;
}


void ConsDragDiff::computeSelfConsistFreq(LevelData<FArrayBox>&       a_cls_freq,
                                          const LevelData<FArrayBox>& a_density,
                                          const LevelData<FArrayBox>& a_temperature,
                                          const double                a_mass,
                                          const double                a_charge ) const

{
    //Get normalization parameters (units)
    double N, T, L;
    ParmParse ppunits( "units" );
    ppunits.get("number_density",N);  //[m^{-3}]
    ppunits.get("temperature",T);     //[eV]
    ppunits.get("length",L);          //[m]
  
    double pi = Constants::PI;
    double ech = Constants::ELEMENTARY_CHARGE;
    double eps0 = Constants::VACUUM_PERMITTIVITY;
    double coeff = 1/( 12.0 * eps0 * eps0 * pow(pi, 3.0/2.0) );

    // Compute normalized collision frequency, cls_norm = cls_phys * L/Sqrt[T/mp]
    // use Eq. (2.53) from E. Shi Thesis, but retain our expression for Coulomb logarithm
    double Coulomb_Lg = 23 - log( sqrt(2.0) * pow(a_charge,3) * sqrt(N)/1000.0 / pow(T, 3.0/2.0) );
    double cls_norm = coeff * N * pow(ech, 2) * pow(a_charge, 4) * L
                     / ( sqrt(a_mass) *  pow(T, 2)) * Coulomb_Lg;

    const DisjointBoxLayout& grids = a_cls_freq.disjointBoxLayout();
    DataIterator dit(grids.dataIterator());
    for (dit.begin(); dit.ok(); ++dit) {

      FArrayBox& this_cls_freq = a_cls_freq[dit];
      const FArrayBox& this_n = a_density[dit];
      const FArrayBox& this_T = a_temperature[dit];

      FORT_COMPUTE_SC_CLS_FREQ(CHF_BOX(this_cls_freq.box()),
                               CHF_FRA1(this_cls_freq,0),
                               CHF_CONST_FRA1(this_n,0),
                               CHF_CONST_FRA1(this_T,0));
  
      this_cls_freq.mult(cls_norm);

    }
}

void ConsDragDiff::ParseParameters(const ParmParse& a_pp) {

  a_pp.query("cls_only",m_cls_only);
  a_pp.query("verbosity",m_verbosity);
  a_pp.query("time_implicit", m_time_implicit);
  CH_assert( m_verbosity == 0 || m_verbosity == 1);
  
  CFG::GridFunctionLibrary* grid_library = CFG::GridFunctionLibrary::getInstance();
  std::string grid_function_name;
  
  if (a_pp.contains("cls_freq")) {
    a_pp.get("cls_freq", grid_function_name );
    m_cls_freq_func = grid_library->find( grid_function_name );
  }

  a_pp.query("update_cls_freq", m_update_freq );
}

inline void ConsDragDiff::printParameters( const KineticSpeciesPtrVect&  soln )
{
   // get stuff to calculate stability parameters on time step
   const KineticSpecies& soln_species( *(soln[0]) );
   const LevelData<FArrayBox>& soln_dfn = soln_species.distributionFunction();
   double mass = soln_species.mass();
   const DisjointBoxLayout& dbl = soln_dfn.getBoxes();
   const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();
   const ProblemDomain& phase_domain = phase_geom.domain();
   const Box& domain_box = phase_domain.domainBox();
   int num_mu_cells = domain_box.size(MU_DIR);
   int num_vp_cells = domain_box.size(VPARALLEL_DIR);

   // check to ensure that fourth-order one-sided calculations do not
   // reach into undefined physical ghosts
   if ( procID() == 0 && (num_vp_cells < 5 || num_mu_cells < 5) ){
     MayDay::Error( "total number of velocity cells is too few for physically meaningful results" );
   }

   // get magnetic field, density, and pressure, etc..
   const CFG::MagGeom & mag_geom = phase_geom.magGeom();
   const LevelData<FArrayBox>& inj_B = phase_geom.getBFieldMagnitude();
   IntVect ghost_phase = soln_dfn.ghostVect();
   CFG::IntVect ghost_cfg = phase_geom.config_restrict(ghost_phase);
   
   CFG::LevelData<CFG::FArrayBox> dens_cfg(mag_geom.grids(), 1, ghost_cfg);
   CFG::LevelData<CFG::FArrayBox> temp_cfg(mag_geom.grids(), 1, ghost_cfg);
   CFG::LevelData<CFG::FArrayBox> Upar_cfg(mag_geom.grids(), 1, ghost_cfg);

   soln_species.numberDensity(dens_cfg);
   soln_species.parallelVelocity(Upar_cfg);
   soln_species.pressure(temp_cfg,Upar_cfg);
   
   CFG::DataIterator cfg_dit = Upar_cfg.dataIterator();
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit)
   {
       temp_cfg[cfg_dit].divide(dens_cfg[cfg_dit]);
   }

   //Introduce a zero-valued vpar-shift, to make Justin old calculations consistent with the 
   //new pressure diagnostics, that now subtract v_par_shift from v_par --MD 02/28/14
   //CFG::LevelData<CFG::FArrayBox> vpar_zero(mag_geom.grids(), 3, CFG::IntVect::Zero);
   //CFG::DataIterator cfg_dit = vpar_zero.dataIterator();
   //for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit)
   //{
   //   vpar_zero[cfg_dit].setVal(0.0);
   //}


   LevelData<FArrayBox> density;
   LevelData<FArrayBox> Temp;
   LevelData<FArrayBox> Upar;
   phase_geom.injectConfigurationToPhase(dens_cfg, density);
   phase_geom.injectConfigurationToPhase(temp_cfg, Temp);
   phase_geom.injectConfigurationToPhase(Upar_cfg, Upar);

   // creat temperature/B
   LevelData<FArrayBox> TonB(dbl, 1, IntVect::Zero);
   DataIterator bdit= inj_B.dataIterator();
   for (bdit.begin(); bdit.ok(); ++bdit)
   {
     TonB[bdit].copy(Temp[bdit].divide(inj_B[bdit]));
   }

   // get velocity cell spacing
   Real dvp;
   Real dmu;
   for (bdit.begin(); bdit.ok(); ++bdit)
   {
     const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(dbl[bdit]);
     const RealVect& phase_dx =  block_coord_sys.dx();
     dvp = phase_dx[VPARALLEL_DIR];
     dmu = phase_dx[MU_DIR];
   }
   const Real mu_max = num_mu_cells*dmu;
   const Real vp_max = num_vp_cells*dvp/2.0;

   // get max values of T, T/B, and Upar
   const DisjointBoxLayout& grids = inj_B.getBoxes();
   Real local_maxTemp(-1000);     // -1000 is just a starting point
   Real local_maxTonB(-1000);     // -1000 is just a starting point
   Real local_maxNu(-1000);       // -1000 is just a starting point
   Real local_minTemp(1000);      // 1000 is just a starting point
   Real local_minTonB(1000);      // 1000 is just a starting point
   Real local_minUpar(1000);      // 1000 is just a starting point

   for (bdit.begin(); bdit.ok(); ++bdit)
   {
      Box box( grids[bdit] );
      Real box_maxTemp( Temp[bdit].max( box ) );
      Real box_maxTonB( TonB[bdit].max( box ) );
      Real box_maxNu( m_cls_freq[bdit].max( box ) );
      Real box_minTemp( Temp[bdit].min( box ) );
      Real box_minTonB( TonB[bdit].min( box ) );
      Real box_minUpar( Upar[bdit].min( box ) );
      
     
      local_maxTemp = Max( local_maxTemp, box_maxTemp );
      local_maxTonB = Max( local_maxTonB, box_maxTonB );
      local_maxNu = Max( local_maxNu, box_maxNu );
      local_minTemp = Min( local_minTemp, box_minTemp );
      local_minTonB = Min( local_minTonB, box_minTonB );
      local_minUpar = Min( local_minUpar, box_minUpar );

   }
   Real max_Temp( local_maxTemp );
   Real max_TonB( local_maxTonB );
   Real max_Nu( local_maxNu );
   Real min_Temp( local_minTemp );
   Real min_TonB( local_minTonB );
   Real min_Upar( local_minUpar );

#ifdef CH_MPI
   MPI_Allreduce( &local_maxTemp, &max_Temp, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
   MPI_Allreduce( &local_maxTonB, &max_TonB, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
   MPI_Allreduce( &local_maxNu,   &max_Nu,   1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
   MPI_Allreduce( &local_minTemp, &min_Temp, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD );
   MPI_Allreduce( &local_minTonB, &min_TonB, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD );
   MPI_Allreduce( &local_minUpar, &min_Upar, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD );
#endif

// cout << "max_Temp =" << max_Temp << endl;
// cout << "max_TonB =" << max_TonB << endl;

   // calculate min time step set by purely diffusive parts
   // dt*sum_i( D_i/dx_i^2 ) << 1/2;
   const Real alpha_vp = max_Nu*(max_Temp/mass)/(dvp*dvp);
   const Real alpha_mu = max_Nu*(2.0*mu_max*max_TonB)/(dmu*dmu);
   const Real dt_stable_alpha = 0.5/(alpha_vp+alpha_mu);

   // calculate min time step set by advective-diffusive
   // dt*sum_i( U_i^2/D_i ) << 2;
   const Real beta_vp = max_Nu*(vp_max-min_Upar)*(vp_max-min_Upar)/min_Temp/mass;
   const Real beta_mu = 2.0*max_Nu*mu_max/min_TonB;
   const Real dt_stable_beta = 2.0/(beta_vp+beta_mu);

   // take mininum of two time step restrictions
   Real dt_stable = Min(dt_stable_alpha, dt_stable_beta);

   // calculate Peclet number for mu and vp
   // (must be less than one for physically meaningful results)
   const Real Peclet_vp = dvp/sqrt(2.0*min_Temp/mass);
   const Real Peclet_mu = dmu/2.0/min_TonB;
   if ( procID() == 0 && !(Peclet_vp < 1) ){
     MayDay::Warning( "vpar grid spacing too large for physically meaningful results" );
   }
   if ( procID() == 0 && !(Peclet_mu < 1) ){
     MayDay::Warning( "mu grid spacing too large for physically meaningful results" );
   }

   std::string cls_only_str("false");
   if (m_cls_only){ cls_only_str = "true";}

   if (procID() == 0 ) {
     cout << "--------------------------------------------------------------------- " << endl;
     cout << "Conservative Drag Diffusion Collision Model being defined for species " << m_species_name <<endl;
     cout << "max_cls_freq  =  " << max_Nu << ", cls_only = " <<  cls_only_str << endl;
     cout << "Stable Diffusive time step           = " << dt_stable_alpha << endl;
     cout << "Stable Advective-Diffusive time step = " << dt_stable_beta << endl;
     cout << "Stable time step = " << dt_stable << endl;
     cout << "Peclet_vp = " << Peclet_vp << endl;
     cout << "Peclet_mu = " << Peclet_mu << endl;
     cout << "Implicit in time = " << m_time_implicit << endl;
     cout << "--------------------------------------------------------------------- " << endl;
   }
}

void ConsDragDiff::preTimeStep(const KineticSpeciesPtrVect& a_soln_mapped,
                               const int a_species,
                               const Real a_time,
                               const KineticSpeciesPtrVect& a_soln_physical )
{
   m_it_counter+=1;
}

Real ConsDragDiff::computeTimeScale( const KineticSpeciesPtrVect&  soln, const int a_idx )
{
   // replace with a good estimate of time scale
   return DBL_MAX;
}

Real ConsDragDiff::computeDtExplicitTI( const KineticSpeciesPtrVect&  soln, const int a_idx )
{
   // replace with a good estimate of dt
   return DBL_MAX;
}

Real ConsDragDiff::computeDtImExTI( const KineticSpeciesPtrVect&  soln, const int a_idx )
{
  if (m_time_implicit) {
    return DBL_MAX;
  } else {
   // replace with a good estimate of dt
   return DBL_MAX;
  }
}

void ConsDragDiff::diagnostics(const LevelData<FArrayBox>& a_rhs,
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
  moment_op.compute(parMom_src, a_rhs_species, a_rhs, ParallelVelKernel<FArrayBox>());
  phase_geom.plotConfigurationData( "parMom_src", parMom_src, a_time );
  
  //Plot energy source
  CFG::LevelData<CFG::FArrayBox> energy_src( mag_geom.grids(), 1, CFG::IntVect::Zero );
  moment_op.compute(energy_src, a_rhs_species, a_rhs, KineticEnergyKernel<FArrayBox>());
  phase_geom.plotConfigurationData( "energy_src", energy_src, a_time );
}

#include "NamespaceFooter.H"
