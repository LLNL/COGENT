#include <math.h>
#include "CONSTANTS.H"

#include "Directions.H"
#include "Anomalous.H"
#include "TransportF_F.H"

#include "MomentOp.H"
#include "Kernels.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FourthOrderUtil.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H" // has to be the last one

Anomalous::Anomalous( const string& a_species_name, ParmParse& a_pptpm, const int a_verbosity )
    : m_moment_op( MomentOp::instance() ),
      m_D_fluid(4,0),
      m_D_kinet(4,0),
      m_DN0(0),
      m_model_only(false),
      m_const_coeff(true),
      m_first_call(true),
      m_first_stage(true),
      m_arbitrary_grid(true),
      m_flux_aligned_grid(false),
      m_simple_diffusion(false),
      m_update_freq(-1),
      m_it_counter(0),
      m_verbosity(0)
{

   m_verbosity = a_verbosity;
   pptpm = a_pptpm;
   species_name = a_species_name;
   
   ParseParameters();
}

Anomalous::~Anomalous()
{
}

void Anomalous::evalTpmRHS( KineticSpecies&               a_rhs_species,
                            const KineticSpeciesPtrVect&  a_soln,
                            const int                     a_species,
                            const Real                    a_time )
{
  /*
    Evaluates the anomalous flux normal to magnetic surfaces as
    described in Ch. 5 of Justin Angus' PhD Thesis:
    df/dt_anomalous=-div(flux_anom)
    flux_anom = flux_psi*psi_hat
    flux_psi = U_psi*f - D_psi*nabla_psi(f)
    nabla_psi =  nabla_psi|v=const = nabla_psi|mu=const - mu/B*nabla_psi(B)*d/dmu
    U_psi = Dn0*nabla_psi(ln(n)) + DT0*nabla_psi(ln(T))
          + (vp^2/Vth^2-3/2)*(Dn2*nabla_psi(ln(n))+DT2*nabla_psi(ln(T)))
    D_psi = D0
  */
  
  CH_TIME("Anomalous::evalTpmRHS");
     
  // print parameters at the first time step
  // involves kinetic species to print stable dt, improve design later
  // move print parameters to the class constructor
  const KineticSpecies& soln_species( *(a_soln[a_species]) );
  if ((m_verbosity) && (m_first_call)) {printParameters(soln_species);}

  // get vlasov RHS for the current species
  LevelData<FArrayBox>& rhs_dfn = a_rhs_species.distributionFunction();

  // get solution distribution function (f*Bstarpar) for the current species
  const LevelData<FArrayBox>& soln_fB = soln_species.distributionFunction();
  double mass = soln_species.mass();
  const DisjointBoxLayout& dbl = soln_fB.getBoxes();

  // get coordinate system parameters
  const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();
  bool second_order = phase_geom.secondOrder();

  // set phase ghost width
  if (second_order) {
      m_ghostVect = IntVect::Unit;
  }
  else {
      m_ghostVect = 2 * IntVect::Unit;
  }

  // copy const soln_fB (along with the boundary values) to a temporary
  // we do it since m_fB can have smaller number of ghosts than soln_fB
  // and this can speed up moment calculations
  if ( !m_fB.isDefined()) m_fB.define( dbl, 1, m_ghostVect );
  DataIterator sdit = m_fB.dataIterator();
  for (sdit.begin(); sdit.ok(); ++sdit) {
    m_fB[sdit].copy(soln_fB[sdit]);
  }

  // compute configuration space moments
  computeCfgMoments(m_fB, soln_species);
  
  // initialize transport coefficients, Lame coefficients, and preconditioner
  if (m_first_call) {
    initializeData(phase_geom, a_time);
  }

  // create cell-centered dfB/dmu (presently we only have 2nd order implementation)
  if ( !m_dfBdmu_cc.isDefined()) m_dfBdmu_cc.define( dbl, 1, IntVect::Unit );
  CH_assert(soln_fB.ghostVect() > m_dfBdmu_cc.ghostVect());
  if (!m_simple_diffusion) {
    computeDfDmu(m_dfBdmu_cc, soln_fB, phase_geom);
  }
  
  // define member object to stor RHS
  if (!m_rhs_transport.isDefined()) m_rhs_transport.define(dbl, 1, IntVect::Zero);
      
  if (SpaceDim == 4) {
    // compute divergence for locally-orthogonal and flux-algined grid
    if (!m_arbitrary_grid) {
      computeRhsOrthogonalGrid(m_rhs_transport, phase_geom, mass);
    }
    // compute divergence for arbitrary grid
    else {
      computeRhsArbtryGrid(m_rhs_transport, phase_geom, mass);
    }
  }
  
  // compute simple diffusion for 3D
 if (SpaceDim == 5) {
   // NB: for some reasons BC classes do not properly fill codim2 ghosts
   // for the case of 3D SN field-aligned geometry, so use the
   // alternative functionality here 
    phase_geom.fillCoDim2BoundaryGhosts(m_fB);
    computeRhs3D(m_rhs_transport, phase_geom);
  }
   
  // add (or overwrite) transport RHS to Vlasov RHS
  DataIterator rdit = m_rhs_transport.dataIterator();
  for (rdit.begin(); rdit.ok(); ++rdit)
  {
    if(m_model_only) {
      rhs_dfn[rdit].copy(m_rhs_transport[rdit]);
    }
    else  {
      rhs_dfn[rdit].plus(m_rhs_transport[rdit]);
    }
  }

  m_first_call = false;
  m_first_stage = false;
}

void Anomalous::computeRhsOrthogonalGrid(LevelData<FArrayBox>&  a_rhs,
                                         const PhaseGeom&       a_phase_geom,
                                         const double           a_mass)
{
  
  const DisjointBoxLayout& dbl = a_rhs.getBoxes();
  if (!m_fluxNorm.isDefined()) m_fluxNorm.define(dbl, 1, IntVect::Zero);

  const LevelData<FArrayBox>& inj_B = a_phase_geom.getBFieldMagnitude();
  
  for (DataIterator dit( m_fluxNorm.dataIterator() ); dit.ok(); ++dit)
  {
    // get phase space dx on this patch
    const PhaseBlockCoordSys& block_coord_sys = a_phase_geom.getBlockCoordSys(dbl[dit]);
    const RealVect& phase_dx =  block_coord_sys.dx();

    const FArrayBox& fB_on_patch = m_fB[dit];
    const FArrayBox& dfBdmu_on_patch = m_dfBdmu_cc[dit];
    const FArrayBox& B_on_patch  = inj_B[dit];
    const FArrayBox& N_on_patch  = m_density[dit];
    const FArrayBox& U_on_patch  = m_Upar[dit];
    const FArrayBox& T_on_patch  = m_temperature[dit];
    const FArrayBox& C_on_patch  = m_fourth_coef[dit];
    const FArrayBox& P_on_patch  = m_perp_coef[dit];
         
    for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& thisfluxNorm = m_fluxNorm[dit][dir];
      const FArrayBox& lame_on_patch = m_lame_faces[dit][dir];
      const FArrayBox& D_kinet_on_patch = m_D_kinet_faces[dit][dir];
      
      const int simpleDiff = m_simple_diffusion? 1: 0;
      
      FORT_EVAL_ANOM_FLUX(CHF_CONST_INT(dir),
                          CHF_BOX(thisfluxNorm.box()),
                          CHF_CONST_REALVECT(phase_dx),
                          CHF_CONST_FRA(lame_on_patch),
                          CHF_CONST_FRA(D_kinet_on_patch),
                          CHF_CONST_REAL(a_mass),
                          CHF_CONST_INT(simpleDiff),
                          CHF_CONST_FRA1(fB_on_patch,0),
                          CHF_CONST_FRA1(dfBdmu_on_patch,0),
                          CHF_CONST_FRA1(B_on_patch,0),
                          CHF_CONST_FRA1(N_on_patch,0),
                          CHF_CONST_FRA1(U_on_patch,0),
                          CHF_CONST_FRA1(T_on_patch,0),
                          CHF_CONST_FRA1(C_on_patch,0),
                          CHF_CONST_FRA1(P_on_patch,0),
                          CHF_FRA1(thisfluxNorm,0) );
    }
  }
  
  a_phase_geom.averageAtBlockBoundaries(m_fluxNorm);
  a_phase_geom.mappedGridDivergenceFromFluxNormals(a_rhs, m_fluxNorm);
}

void Anomalous::computeRhsArbtryGrid(LevelData<FArrayBox>&  a_rhs,
                                     const PhaseGeom&       a_phase_geom,
                                     const double           a_mass)
{
  const DisjointBoxLayout& dbl = a_rhs.getBoxes();
  if (!m_flux.isDefined()) m_flux.define(dbl, SpaceDim, m_ghostVect-IntVect::Unit);
  
  const LevelData<FArrayBox>& inj_B = a_phase_geom.getBFieldMagnitude();
  
  for (DataIterator dit( m_flux.dataIterator() ); dit.ok(); ++dit)
  {
    // get phase space dx on this patch
    const PhaseBlockCoordSys& block_coord_sys = a_phase_geom.getBlockCoordSys(dbl[dit]);
    const RealVect& phase_dx =  block_coord_sys.dx();

    const FArrayBox& fB_on_patch = m_fB[dit];
    const FArrayBox& dfBdmu_on_patch = m_dfBdmu_cc[dit];
    const FArrayBox& B_on_patch  = inj_B[dit];
    const FArrayBox& N_on_patch  = m_density[dit];
    const FArrayBox& U_on_patch  = m_Upar[dit];
    const FArrayBox& T_on_patch  = m_temperature[dit];
    const FArrayBox& C_on_patch  = m_fourth_coef[dit];
    const FArrayBox& P_on_patch  = m_perp_coef[dit];
         
    for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& thisflux = m_flux[dit][dir];
      const FArrayBox& bunit_on_patch  = m_inj_bunit[dit][dir];
      const FArrayBox& NJinv_on_patch = m_inj_pointwiseNJinv[dit][dir];
      const FArrayBox& D_kinet_on_patch = m_D_kinet_faces[dit][dir];

      Box face_box( dbl[dit] );
      // for 4th order, we need an extra face in the mapped-grid,
      // transverse directions to handle 4th-order products
      bool second_order = a_phase_geom.secondOrder();
      if (!second_order) {
        for (int tdir(0); tdir<SpaceDim; tdir++) {
          if (tdir!=dir) {
            const int TRANSVERSE_GROW(1);
            face_box.grow( tdir, TRANSVERSE_GROW );
          }
        }
      }
   
      // This will turn cell-centered box into
      //face-centered box in the direction dir
      face_box.surroundingNodes( dir );
            
      const int simpleDiff = m_simple_diffusion? 1: 0;
      
      FORT_EVAL_ANOM_FLUX_ARBTRY_GRID(CHF_CONST_INT(dir),
                                      CHF_BOX(face_box),
                                      CHF_CONST_REALVECT(phase_dx),
                                      CHF_CONST_FRA(NJinv_on_patch),
                                      CHF_CONST_FRA(bunit_on_patch),
                                      CHF_CONST_FRA(D_kinet_on_patch),
                                      CHF_CONST_REAL(a_mass),
                                      CHF_CONST_INT(simpleDiff),
                                      CHF_CONST_FRA1(fB_on_patch,0),
                                      CHF_CONST_FRA1(dfBdmu_on_patch,0),
                                      CHF_CONST_FRA1(B_on_patch,0),
                                      CHF_CONST_FRA1(N_on_patch,0),
                                      CHF_CONST_FRA1(U_on_patch,0),
                                      CHF_CONST_FRA1(T_on_patch,0),
                                      CHF_CONST_FRA1(C_on_patch,0),
                                      CHF_CONST_FRA1(P_on_patch,0),
                                      CHF_FRA(thisflux) );
          
    }
  }

  a_phase_geom.applyAxisymmetricCorrection(m_flux);
  // The averaging at block boundaries is likely unnesessary (at least to 2nd order)
  // (as it is also done as part of mappedGridDivergence calculation)
  a_phase_geom.averageAtBlockBoundaries(m_flux);
  const bool OMIT_NT(false);
  a_phase_geom.mappedGridDivergence( a_rhs, m_flux, OMIT_NT );
  for (DataIterator dit( a_rhs.dataIterator() ); dit.ok(); ++dit) {
    const PhaseBlockCoordSys& block_coord_sys( a_phase_geom.getBlockCoordSys( dbl[dit] ) );
    double fac( 1.0 / block_coord_sys.getMappedCellVolume() );
    a_rhs[dit].mult( fac );
  }
}
  

void Anomalous::computeRhs3D(LevelData<FArrayBox>&  a_rhs,
                             const PhaseGeom&       a_phase_geom)
{
  const DisjointBoxLayout& dbl = a_rhs.getBoxes();
  if (!m_flux.isDefined()) m_flux.define(dbl, SpaceDim, m_ghostVect-IntVect::Unit);
  
  const LevelData<FluxBox>& inj_B_face = a_phase_geom.getBFieldFace();
  
  if (!m_inj_coords_fc.isDefined()) {
    const CFG::MagGeom& mag_geom = a_phase_geom.magGeom();
    const CFG::LevelData<CFG::FluxBox>& cfg_real_coords = mag_geom.getFaceCenteredRealCoords();
    a_phase_geom.injectConfigurationToPhase(cfg_real_coords, m_inj_coords_fc);
  }
  
  // Compute physical flux
  for (DataIterator dit( m_flux.dataIterator() ); dit.ok(); ++dit)
  {
    // get phase space dx on this patch
    const PhaseBlockCoordSys& block_coord_sys = a_phase_geom.getBlockCoordSys(dbl[dit]);
    const RealVect& phase_dx =  block_coord_sys.dx();
    
    for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& thisflux = m_flux[dit][dir];
      const FArrayBox& fB_on_patch = m_fB[dit];
      const FArrayBox& B_on_patch  = inj_B_face[dit][dir];
      const FArrayBox& coords_on_patch  = m_inj_coords_fc[dit][dir];
      const FArrayBox& NJinv_on_patch = m_inj_pointwiseNJinv[dit][dir];
      const FArrayBox& D_kinet_on_patch = m_D_kinet_faces[dit][dir];

      Box face_box( dbl[dit] );
      // for 4th order, we need an extra face in the mapped-grid,
      // transverse directions to handle 4th-order products
      bool second_order = a_phase_geom.secondOrder();
      if (!second_order) {
        for (int tdir(0); tdir<SpaceDim; tdir++) {
          if (tdir!=dir) {
            const int TRANSVERSE_GROW(1);
            face_box.grow( tdir, TRANSVERSE_GROW );
          }
        }
      }
   
      //This will turn cell-centered box into
      //face-centered box in the direction dir
      face_box.surroundingNodes( dir );
      
      const int flux_aligned = m_flux_aligned_grid? 1: 0;
      FORT_EVAL_ANOM_FLUX_3D(CHF_CONST_INT(dir),
                             CHF_BOX(face_box),
                             CHF_CONST_REALVECT(phase_dx),
                             CHF_CONST_FRA(NJinv_on_patch),
                             CHF_CONST_FRA(D_kinet_on_patch),
                             CHF_CONST_FRA(coords_on_patch),
                             CHF_CONST_FRA(B_on_patch),
                             CHF_CONST_INT(flux_aligned),
                             CHF_CONST_FRA1(fB_on_patch,0),
                             CHF_FRA(thisflux) );
    }
  }

  a_phase_geom.applyAxisymmetricCorrection(m_flux);
  const bool OMIT_NT(false);
  a_phase_geom.mappedGridDivergence( a_rhs, m_flux, OMIT_NT );
  for (DataIterator dit( a_rhs.dataIterator() ); dit.ok(); ++dit) {
    const PhaseBlockCoordSys& block_coord_sys( a_phase_geom.getBlockCoordSys( dbl[dit] ) );
    double fac( 1.0 / block_coord_sys.getMappedCellVolume() );
    a_rhs[dit].mult( fac );
  }
}


void Anomalous::initializeData(const PhaseGeom& a_phase_geom,
                               const Real       a_time)
{
  // get NJinverse and B-field data for arbitrary grid calculation
  const CFG::IntVect cfg_ghostVect(a_phase_geom.config_restrict(m_ghostVect));

  const CFG::MagGeom& mag_geom = a_phase_geom.magGeom();
  
  CFG::LevelData<CFG::FluxBox> pointwiseNJinv_cfg(mag_geom.grids(), CFG_DIM*CFG_DIM, cfg_ghostVect);
  mag_geom.getPointwiseNJInverse(pointwiseNJinv_cfg);
  a_phase_geom.injectConfigurationToPhase(pointwiseNJinv_cfg, m_inj_pointwiseNJinv);
  
  const CFG::LevelData<CFG::FluxBox>& bunit_cfg = mag_geom.getFCBFieldDir();
  a_phase_geom.injectConfigurationToPhase(bunit_cfg, m_inj_bunit);

  // get face centered lame coefficients (h_r, h_theta, and h_phi) for locally-orthogonal grids
  const DisjointBoxLayout& dbl = a_phase_geom.gridsFull();
  m_lame_faces.define(dbl, 3, m_ghostVect - IntVect::Unit);
  getFaceCenteredLameCoefficients(m_lame_faces, a_phase_geom, dbl);
         
  // set spatial depedence of transport coefficients
  // because we are using forthOrderCellToFaceCenters (update to 2nd order later)
  // we have to use a minimum of 2 ghost layers for the related quantities
  CFG::LevelData<CFG::FArrayBox> D_cfg( mag_geom.grids(), 5, 2*CFG::IntVect::Unit );
  CFG::DataIterator cfg_dit = D_cfg.dataIterator();
  if (m_const_coeff) {
    for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
      for (int iD=0; iD<4; iD++) {
        D_cfg[cfg_dit].setVal(m_D_kinet[iD],iD);
      }
      D_cfg[cfg_dit].setVal(m_DN0,4);
    }
  }
  else {
    setTransportCoeffSpatialDependence( D_cfg, mag_geom, a_time);
  }
  CFG::IntVect cfg_face_ghost = cfg_ghostVect - CFG::IntVect::Unit;
  CFG::LevelData<CFG::FluxBox> D_cfg_faces(mag_geom.grids(), D_cfg.nComp(), cfg_face_ghost);
  fourthOrderCellToFaceCenters(D_cfg_faces, D_cfg);
  a_phase_geom.injectConfigurationToPhase( D_cfg_faces, m_D_kinet_faces);
        
  // compute the preconditioner coeffficient for the implicit solver
  m_precond_D.define(mag_geom.grids(), CFG_DIM * CFG_DIM, CFG::IntVect::Unit);
  computePrecondCoefficient(m_precond_D, mag_geom, D_cfg_faces);
}



void Anomalous::computeDfDmu(LevelData<FArrayBox>&        a_dfBdmu_cc,
                             const LevelData<FArrayBox>&  a_soln_fB,
                             const PhaseGeom&             a_phase_geom)
{
  
  const DisjointBoxLayout& dbl = a_soln_fB.getBoxes();
  const ProblemDomain& phase_domain = a_phase_geom.domain();
  const Box& domain_box = phase_domain.domainBox();
  int num_mu_cells = domain_box.size(MU_DIR);
  
  DataIterator dit = a_dfBdmu_cc.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    // get phase space dx on this patch
    const PhaseBlockCoordSys& block_coord_sys = a_phase_geom.getBlockCoordSys(dbl[dit]);
    const RealVect& phase_dx =  block_coord_sys.dx();

    const FArrayBox& fB_on_patch = a_soln_fB[dit];
    FArrayBox& dfBdmu_on_patch = a_dfBdmu_cc[dit];

    FORT_DFBDMU_CELL_CENTER(CHF_BOX(dfBdmu_on_patch.box()),
                            CHF_CONST_REALVECT(phase_dx),
                            CHF_CONST_INT(num_mu_cells),
                            CHF_CONST_FRA1(fB_on_patch,0),
                            CHF_FRA1(dfBdmu_on_patch,0) );
  }
  
  bool second_order = a_phase_geom.secondOrder();
  if (!second_order) {
    a_dfBdmu_cc.exchange();
  }
}


void Anomalous::computeCfgMoments(const LevelData<FArrayBox>& a_fB,
                                  const KineticSpecies&       a_soln_species)
{
  
  CH_TIME("Anomalous::computeCfgMoments");
  
  const PhaseGeom& phase_geom = a_soln_species.phaseSpaceGeometry();
  const CFG::MagGeom& mag_geom = phase_geom.magGeom();
  const CFG::DisjointBoxLayout& mag_grids = mag_geom.grids();
  const CFG::IntVect cfg_ghostVect(phase_geom.config_restrict(m_ghostVect));

  // for simple_diffusion case we don't need these moments at all
  if (m_simple_diffusion && !m_density.isDefined()) {
    // Create unidentified density, mean parallel velocity, temperature, and fourth moment coefficient
    CFG::LevelData<CFG::FArrayBox> dens_cfg(mag_grids, 1, cfg_ghostVect);
    CFG::LevelData<CFG::FArrayBox> Upar_cfg(mag_grids, 1, cfg_ghostVect);
    CFG::LevelData<CFG::FArrayBox> temp_cfg(mag_grids, 1, cfg_ghostVect);
    CFG::LevelData<CFG::FArrayBox> four_cfg(mag_grids, 1, cfg_ghostVect);
    CFG::LevelData<CFG::FArrayBox> perp_cfg(mag_grids, 1, cfg_ghostVect);
         
    phase_geom.injectConfigurationToPhase(dens_cfg, m_density);
    phase_geom.injectConfigurationToPhase(Upar_cfg, m_Upar);
    phase_geom.injectConfigurationToPhase(temp_cfg, m_temperature);
    phase_geom.injectConfigurationToPhase(four_cfg, m_fourth_coef);
    phase_geom.injectConfigurationToPhase(perp_cfg, m_perp_coef);
  }

  else if ((m_update_freq < 0) || (m_it_counter % m_update_freq == 0 && m_first_stage) || m_first_call) {
    // Update density, mean parallel velocity, temperature, and fourth moment coefficient
    CFG::LevelData<CFG::FArrayBox> dens_cfg(mag_geom.grids(), 1, cfg_ghostVect);
    CFG::LevelData<CFG::FArrayBox> Upar_cfg(mag_geom.grids(), 1, cfg_ghostVect);
    CFG::LevelData<CFG::FArrayBox> temp_cfg(mag_geom.grids(), 1, cfg_ghostVect);
    CFG::LevelData<CFG::FArrayBox> four_cfg(mag_geom.grids(), 1, cfg_ghostVect);
    CFG::LevelData<CFG::FArrayBox> perp_cfg(mag_geom.grids(), 1, cfg_ghostVect);

    m_moment_op.compute(dens_cfg, a_soln_species, a_fB, DensityKernel<FArrayBox>());
    m_moment_op.compute(Upar_cfg, a_soln_species, a_fB, ParallelMomKernel<FArrayBox>());
    m_moment_op.compute(four_cfg, a_soln_species, a_fB, FourthMomentKernel<FArrayBox>());
    m_moment_op.compute(perp_cfg, a_soln_species, a_fB, PerpEnergyKernel<FArrayBox>());
    CFG::DataIterator cfg_dit = dens_cfg.dataIterator();
    for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
      Upar_cfg[cfg_dit].divide(dens_cfg[cfg_dit]);
    }
    m_moment_op.compute(temp_cfg, a_soln_species, a_fB, PressureKernel<FArrayBox>(Upar_cfg));

    for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
      four_cfg[cfg_dit].divide(temp_cfg[cfg_dit]); // dividing by pressure
      perp_cfg[cfg_dit].divide(temp_cfg[cfg_dit]); // dividing by pressure
    }
    for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
      temp_cfg[cfg_dit].divide(dens_cfg[cfg_dit]);
    }
    for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
      four_cfg[cfg_dit].divide(temp_cfg[cfg_dit]); // dividing by temperature
    }

    // Perform extrapolation into ghosts to avoid possible division by zero
    // (can occur for some BC choices, e.g., zero inflow BC)
    for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
      const CFG::MagBlockCoordSys& cfg_block_coord_sys = mag_geom.getBlockCoordSys(mag_grids[cfg_dit]);
      const CFG::ProblemDomain& cfg_domain = cfg_block_coord_sys.domain();
                     
      secondOrderCellExtrapAtDomainBdry(dens_cfg[cfg_dit], mag_grids[cfg_dit], cfg_domain);
      secondOrderCellExtrapAtDomainBdry(Upar_cfg[cfg_dit], mag_grids[cfg_dit], cfg_domain);
      secondOrderCellExtrapAtDomainBdry(temp_cfg[cfg_dit], mag_grids[cfg_dit], cfg_domain);
      secondOrderCellExtrapAtDomainBdry(four_cfg[cfg_dit], mag_grids[cfg_dit], cfg_domain);
      secondOrderCellExtrapAtDomainBdry(perp_cfg[cfg_dit], mag_grids[cfg_dit], cfg_domain);
    }
    mag_geom.fillInternalGhosts( dens_cfg );
    mag_geom.fillInternalGhosts( Upar_cfg );
    mag_geom.fillInternalGhosts( temp_cfg );
    mag_geom.fillInternalGhosts( four_cfg );
    mag_geom.fillInternalGhosts( perp_cfg );
         
    //Inject into phase-space
    phase_geom.injectConfigurationToPhase(dens_cfg, m_density);
    phase_geom.injectConfigurationToPhase(Upar_cfg, m_Upar);
    phase_geom.injectConfigurationToPhase(temp_cfg, m_temperature);
    phase_geom.injectConfigurationToPhase(four_cfg, m_fourth_coef);
    phase_geom.injectConfigurationToPhase(perp_cfg, m_perp_coef);
  }
}

void Anomalous::getFaceCenteredLameCoefficients( LevelData<FluxBox>&       a_lame_faces,
                                                 const PhaseGeom&          a_phase_geom,
                                                 const DisjointBoxLayout&  a_dbl )
{
   CH_TIME("Anomalous::getFaceCenteredLameCoefficients");
   
   // FluxBox to be filled needs 3 components on each face (hr, htheta, and hphi)
   CH_assert(a_lame_faces.nComp() == 3);
   DataIterator dit = a_lame_faces.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
     // get face centered N and dX/dq components on this patch
     const CFG::MagBlockCoordSys& mag_block_coord_sys = a_phase_geom.getMagBlockCoordSys(a_dbl[dit]);
     Box phase_box(a_dbl[dit]);
     CFG::Box box_config;
     a_phase_geom.projectPhaseToConfiguration(phase_box, box_config);
     box_config.grow(4);
     CFG::FluxBox N_cfg_face(box_config, CFG_DIM*CFG_DIM);  //4 components of N on each CFG_DIM face
     CFG::FluxBox dXdq_cfg_face(box_config, CFG_DIM*CFG_DIM);  //4 components on each CFG_DIM face
     mag_block_coord_sys.getPointwiseN(N_cfg_face);
     mag_block_coord_sys.getFaceCentereddXdxi(dXdq_cfg_face);
     FluxBox N_face, dXdq_face;
     a_phase_geom.injectConfigurationToPhase(N_cfg_face, N_face);
     a_phase_geom.injectConfigurationToPhase(dXdq_cfg_face, dXdq_face);
      
     for (int dir=0; dir<CFG_DIM; dir++)
       {
	 FArrayBox& lame_on_patch = a_lame_faces[dit][dir];
	 
	 FORT_LAME_COEFFICIENTS(CHF_BOX(lame_on_patch.box()),
				CHF_CONST_FRA(N_face[dir]),
				CHF_CONST_FRA(dXdq_face[dir]),
				CHF_FRA(lame_on_patch));
       }
   }  
}

void Anomalous::getCellCenteredLameCoefficients( FArrayBox&        a_lame_cells,
                                                 const PhaseGeom&  a_phase_geom,
                                                 const Box&        a_dbl )
{
   // get cell centered N and dX/dq components
   const CFG::MagBlockCoordSys& mag_block_coord_sys = a_phase_geom.getMagBlockCoordSys(a_dbl);
   Box phase_box(a_dbl);
   CFG::Box box_config;
   a_phase_geom.projectPhaseToConfiguration(phase_box, box_config);
   box_config.grow(2);
   CFG::FArrayBox N_cfg_cent(box_config, CFG_DIM*CFG_DIM);  //4 components of N on each cell center
   CFG::FArrayBox dXdq_cfg_cent(box_config, CFG_DIM*CFG_DIM);  //4 components on each cell center
   mag_block_coord_sys.getPointwiseN(N_cfg_cent);
   mag_block_coord_sys.getCellCentereddXdxi(dXdq_cfg_cent);
   FArrayBox N_cent, dXdq_cent;
   a_phase_geom.injectConfigurationToPhase(N_cfg_cent, N_cent);
   a_phase_geom.injectConfigurationToPhase(dXdq_cfg_cent, dXdq_cent);

   FORT_LAME_COEFFICIENTS( CHF_BOX(a_lame_cells.box()),
                           CHF_CONST_FRA(N_cent),
                           CHF_CONST_FRA(dXdq_cent),
                           CHF_FRA(a_lame_cells));

}

void Anomalous::setTransportCoeffSpatialDependence( CFG::LevelData<CFG::FArrayBox>& a_D_kinet_cfg,
                                                    const CFG::MagGeom&             a_mag_geom,
                                                    const Real                      a_time )
{
   CH_assert(a_D_kinet_cfg.nComp() == 5);
   CFG::DataIterator cfg_dit = a_D_kinet_cfg.dataIterator();
   const CFG::IntVect cfg_ghostVect(2*CFG::IntVect::Unit);
   CFG::LevelData<CFG::FArrayBox> D0_cfg, D1_cfg, D2_cfg, D3_cfg, DN_cfg;
   CFG::LevelData<CFG::FArrayBox> D0_fluid, D1_fluid, D2_fluid, D3_fluid;
   
   D0_cfg.define( a_mag_geom.grids(), 1, cfg_ghostVect );
   D1_cfg.define( a_mag_geom.grids(), 1, cfg_ghostVect );
   D2_cfg.define( a_mag_geom.grids(), 1, cfg_ghostVect );
   D3_cfg.define( a_mag_geom.grids(), 1, cfg_ghostVect );
   DN_cfg.define( a_mag_geom.grids(), 1, cfg_ghostVect );

   D0_fluid.define( a_mag_geom.grids(), 1, cfg_ghostVect );
   D1_fluid.define( a_mag_geom.grids(), 1, cfg_ghostVect );
   D2_fluid.define( a_mag_geom.grids(), 1, cfg_ghostVect );
   D3_fluid.define( a_mag_geom.grids(), 1, cfg_ghostVect );


   // assign shape functions to CFG vectors
   if( pptpm.contains("shape_function_D0") ) {
      m_shape_function_D0->assign( D0_cfg, a_mag_geom, a_time, false);
      m_shape_function_D0->assign( D0_fluid, a_mag_geom, a_time, false);
   } 
   else {
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
         D0_cfg[cfg_dit].setVal( 1.0 );
         D0_fluid[cfg_dit].setVal( 1.0 );
      }
   }
   if( pptpm.contains("shape_function_D1") ) {
      m_shape_function_D1->assign( D1_cfg, a_mag_geom, a_time, false);
      m_shape_function_D1->assign( D1_fluid, a_mag_geom, a_time, false);
   } 
   else {
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
         D1_cfg[cfg_dit].setVal( 1.0 );
         D1_fluid[cfg_dit].setVal( 1.0 );
      }
   }
   if( pptpm.contains("shape_function_D2") ) {
      m_shape_function_D2->assign( D2_cfg, a_mag_geom, a_time, false);
      m_shape_function_D2->assign( D2_fluid, a_mag_geom, a_time, false);
   } 
   else {
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
         D2_cfg[cfg_dit].setVal( 1.0 );
         D2_fluid[cfg_dit].setVal( 1.0 );
      }
   }
   if( pptpm.contains("shape_function_D3") ) {
      m_shape_function_D3->assign( D3_cfg, a_mag_geom, a_time, false);
      m_shape_function_D3->assign( D3_fluid, a_mag_geom, a_time, false);
   } 
   else {
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
         D3_cfg[cfg_dit].setVal( 1.0 );
         D3_fluid[cfg_dit].setVal( 1.0 );
      }
   }
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
      DN_cfg[cfg_dit].setVal( 0.0 );
   }


   // multiply by coefficients set in input file
   if (pptpm.contains("D_kinet")) { // spatial functions refer to kinetic transport coeffs
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
         D0_cfg[cfg_dit].mult( m_D_kinet[0] );
         D1_cfg[cfg_dit].mult( m_D_kinet[1] );
         D2_cfg[cfg_dit].mult( m_D_kinet[2] );
         D3_cfg[cfg_dit].mult( m_D_kinet[3] );
      }
   }


   // if spatial functions are for fluid coefficients, perform arithmetic
   if (pptpm.contains("D_fluid")) { // spatial functions refer to fluid transport coeffs
     
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
         D0_fluid[cfg_dit].mult( m_D_fluid[0] );
         D1_fluid[cfg_dit].mult( m_D_fluid[1] );
         D2_fluid[cfg_dit].mult( m_D_fluid[2] );
         D3_fluid[cfg_dit].mult( m_D_fluid[3] );
         if (m_D_fluid[0]==0) {
            D0_cfg[cfg_dit].copy( D3_fluid[cfg_dit].mult( 2.0/3.0 ) );
            D1_cfg[cfg_dit].copy( D1_fluid[cfg_dit].mult( -1.0 ) );
            D2_cfg[cfg_dit].copy( D2_fluid[cfg_dit].mult( -2.0/3.0 ) );
            D3_cfg[cfg_dit].copy( D1_fluid[cfg_dit] );
            DN_cfg[cfg_dit].copy( D0_fluid[cfg_dit] ); // this is free parameter used to insure D0 is not zero
         }
         else {
            D0_cfg[cfg_dit].copy( D0_fluid[cfg_dit] );               // DN0
            D1_cfg[cfg_dit].copy( D1_fluid[cfg_dit].mult( -1.0 ) );  // DT0
            D2_cfg[cfg_dit].copy( D2_fluid[cfg_dit].mult( -2.0/3.0 ) );
            D2_cfg[cfg_dit].plus( D0_fluid[cfg_dit] );            // DN2
            D3_cfg[cfg_dit].copy( D3_fluid[cfg_dit].mult( -2.0/3.0 ) );
            D3_cfg[cfg_dit].plus( D0_fluid[cfg_dit] );            
            D3_cfg[cfg_dit].plus( D1_fluid[cfg_dit] );            // DT2
         }
      }
   }
     

   // insert coefficients into appropriate components of Fluxbox
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
      a_D_kinet_cfg[cfg_dit].copy(D0_cfg[cfg_dit],0,0,1);
      a_D_kinet_cfg[cfg_dit].copy(D1_cfg[cfg_dit],0,1,1);
      a_D_kinet_cfg[cfg_dit].copy(D2_cfg[cfg_dit],0,2,1);
      a_D_kinet_cfg[cfg_dit].copy(D3_cfg[cfg_dit],0,3,1);
      a_D_kinet_cfg[cfg_dit].copy(DN_cfg[cfg_dit],0,4,1);
   }

}


void Anomalous::computePrecondCoefficient( CFG::LevelData<CFG::FluxBox>&         a_D,
                                          const CFG::MagGeom&                    a_mag_geom,
                                          const CFG::LevelData<CFG::FluxBox>&    a_D_faces )
{
   
   CFG::LevelData<CFG::FluxBox> NJinverse(a_mag_geom.grids(), CFG_DIM*CFG_DIM, CFG::IntVect::Unit);
   a_mag_geom.getPointwiseNJInverse(NJinverse);
   
   CFG::LevelData<CFG::FluxBox> N(a_mag_geom.grids(), CFG_DIM*CFG_DIM, CFG::IntVect::Unit);
   a_mag_geom.getPointwiseN(N);
   
   const CFG::LevelData<CFG::FluxBox>& bunit = a_mag_geom.getFCBFieldDir();
         
   CFG::LevelData<CFG::FluxBox> a_D0(a_mag_geom.grids(), 1, CFG::IntVect::Unit );

   for (CFG::DataIterator dit(a_D.dataIterator()); dit.ok(); ++dit) {
      a_D0[dit].copy(a_D_faces[dit],0,0,1);      

      for (int dir=0; dir<CFG_DIM; ++dir) {

         const CFG::FArrayBox& this_D0 = a_D0[dit][dir];
         CFG::FArrayBox& this_coeff = a_D[dit][dir];
         CFG::FArrayBox& this_N = N[dit][dir];
         CFG::FArrayBox& this_NJinverse = NJinverse[dit][dir];
         const CFG::FArrayBox& this_b = bunit[dit][dir];
         const CFG::Box& box = this_coeff.box();

         //Get coefficient in the cylindrical coordinate frame
         CFG::FArrayBox phys_coeff(box, CFG_DIM * CFG_DIM);
         CFG::BoxIterator bit(box);
         for (bit.begin(); bit.ok(); ++bit) {
            CFG::IntVect iv = bit();
            
            double bNmag = sqrt(pow(this_b(iv,0),2) + pow(this_b(iv,2),2));
            CFG::RealVect bN;
            bN[RADIAL_DIR] = this_b(iv,2) / bNmag;
            bN[POLOIDAL_DIR] = -this_b(iv,0) / bNmag;
            
            phys_coeff(iv,0) = this_D0(iv,0) * bN[RADIAL_DIR]   * bN[RADIAL_DIR];
            phys_coeff(iv,1) = this_D0(iv,0) * bN[POLOIDAL_DIR] * bN[RADIAL_DIR];
            phys_coeff(iv,2) = this_D0(iv,0) * bN[RADIAL_DIR]   * bN[POLOIDAL_DIR];
            phys_coeff(iv,3) = this_D0(iv,0) * bN[POLOIDAL_DIR] * bN[POLOIDAL_DIR];
         
         }
         
         CFG::FArrayBox tmp(box, 1);
         // Multiply the coefficient matrix times the NJInverse matrix
         CFG::FArrayBox phys_coeff_times_NJinv(box, CFG_DIM * CFG_DIM);
         phys_coeff_times_NJinv.setVal(0.);
         for (int row=0; row<CFG_DIM; ++row) {
            for (int col=0; col<CFG_DIM; ++col) {
               for (int m=0; m<CFG_DIM; ++m) {
                  tmp.copy(phys_coeff,CFG_DIM*row+m,0,1);
                  tmp.mult(this_NJinverse,CFG_DIM*m+col,0,1);
                  phys_coeff_times_NJinv.plus(tmp,0,CFG_DIM*row+col,1);
               }
            }
         }
         
         //Premultiply by the NTranspose matrix
         this_coeff.setVal(0.);
         for (int row=0; row<CFG_DIM; ++row) {
            for (int col=0; col<CFG_DIM; ++col) {
               for (int m=0; m<CFG_DIM; ++m) {
                  tmp.copy(this_N,CFG_DIM*m+row,0,1);
                  tmp.mult(phys_coeff_times_NJinv,CFG_DIM*m+col,0,1);
                  this_coeff.plus(tmp,0,CFG_DIM*row+col,1);
               }
            }
         }
      }
   }
}

void Anomalous::ParseParameters()

{

   pptpm.query("model_only",m_model_only);
   pptpm.query("simple_diffusion",m_simple_diffusion);
   pptpm.query("verbosity",m_verbosity);
   CH_assert( m_verbosity == 0 || m_verbosity == 1);
   pptpm.query("arbitrary_grid",m_arbitrary_grid);
   pptpm.query("flux_aligned_grid",m_flux_aligned_grid);
   pptpm.query("update_frequency",m_update_freq);
   
   if (pptpm.contains("const_coeff") ) {
     pptpm.get("const_coeff", m_const_coeff);

     if (!m_const_coeff) {
       
       CFG::GridFunctionLibrary* grid_library = CFG::GridFunctionLibrary::getInstance();
       std::string grid_function_name;
       //pptpm.query( "shape_function", grid_function_name );
       //m_shape_function = grid_library->find( grid_function_name );
      
       if( pptpm.query( "shape_function_D0", grid_function_name ) ) {
         m_shape_function_D0 = grid_library->find( grid_function_name );
       }
       if( pptpm.query( "shape_function_D1", grid_function_name ) ) {
         m_shape_function_D1 = grid_library->find( grid_function_name );
       }
       if ( pptpm.query( "shape_function_D2", grid_function_name ) ) {
         m_shape_function_D2 = grid_library->find( grid_function_name );
       }
       if( pptpm.query( "shape_function_D3", grid_function_name ) ) {
         m_shape_function_D3 = grid_library->find( grid_function_name );
       }
     }
   }

   // check to make sure coefficients properly set in input file
   if ( procID()==0 && !pptpm.contains("D_fluid") && !pptpm.contains("D_kinet") ) {
     cout << "NEITHER FLUID NOR KINETIC TRANSPORT COEFFICIENTS ARE SET FOR SPECIES " << species_name << endl;
     MayDay::Error( "TRANSPORT MODEL BEING USED, BUT COEFFICIENTS NOT PROPERLY DEFINED" );
   }
   if ( procID()==0 && pptpm.contains("D_fluid") && pptpm.contains("D_kinet") ) {
     cout << "BOTH FLUID AND KINETIC TRANSPORT COEFFICIENTS ARE SET FOR SPECIES " << species_name << endl;
     MayDay::Error( "TRANSPORT MODEL BEING USED, BUT COEFFICIENTS NOT PROPERLY DEFINED" );
   }
   if (pptpm.contains("D_fluid")) {
     pptpm.getarr("D_fluid",m_D_fluid,0,4);
     // convert fluid matrix components to kinetic coefficients
     if (m_D_fluid[0]==0) {
       m_D_kinet[0] = 2.0/3.0*m_D_fluid[3];
       m_D_kinet[1] = -m_D_fluid[1];
       m_D_kinet[2] = -2.0/3.0*m_D_fluid[2];
       m_D_kinet[3] = m_D_fluid[1];
       m_DN0 = m_D_kinet[0]; // this is free parameter used to insure D0 is not zero
     }
     else {
       m_D_kinet[0] = m_D_fluid[0];                                // D0
       m_D_kinet[1] = -m_D_fluid[1];                               // DT0
       m_D_kinet[2] = -2.0/3.0*m_D_fluid[2]+m_D_fluid[0];            // DN2
       m_D_kinet[3] = -2.0/3.0*m_D_fluid[3]+m_D_fluid[0]+m_D_fluid[1]; // DT2
     }
   }
   if (pptpm.contains("D_kinet")) {
     pptpm.getarr("D_kinet",m_D_kinet,0,4);
     // convert kinetic matrix components to fluid coefficients
     m_D_fluid[0] = m_D_kinet[0];                                  // DNN
     m_D_fluid[1] = -m_D_kinet[1];                                 // DNT
     m_D_fluid[2] = 3.0/2.0*(m_D_kinet[0]-m_D_kinet[2]);             // DTN
     m_D_fluid[3] = 3.0/2.0*(m_D_kinet[0]-m_D_kinet[1]-m_D_kinet[3]);  // DTT
   }

   // make sure that Anomalous Diff. coeff. D0 is positive
   if ( procID()==0 && (m_D_kinet[0] <= 0) ) {
     cout << "DIFFUSION COEFFICIENT D0 FOR SPECIES " << species_name << " IS <= ZERO " << endl;
     if (pptpm.contains("D_fluid")) {
       MayDay::Error( "FLUID COEFFICIENT D11 = D_fluid[0] MUST BE >= ZERO " );
     }
     else {
       MayDay::Error( "KINETIC DIFFUSION COEFFICIENT D0 = D_kinet[0] MUST BE > ZERO " );
     }
   }

}

inline void Anomalous::printParameters( const KineticSpecies& soln )
{

   Real dt_stable = computeStableDt(soln); 

   std::string tpt_only_str("false");
   std::string tpt_arbitraryGrid_str("true");
   std::string tpt_simpleDiffusion_str("false");
   if (m_model_only){ tpt_only_str = "true";}
   if (!m_arbitrary_grid){ tpt_arbitraryGrid_str = "false";}
   if (m_simple_diffusion){ tpt_simpleDiffusion_str = "true";}

   if (procID() == 0 ) {
     cout << "--------------------------------------------------------------------- " << endl;
     cout << "Anomalous Transport Model being defined for species " << species_name <<endl;
     cout << "tpt_only = " <<  tpt_only_str << endl;
     cout << "Arbitrary grid calculation  = " <<  tpt_arbitraryGrid_str << endl;
     cout << "simple diffusion model  = " <<  tpt_simpleDiffusion_str << endl;
     cout << "Fluid Tranpsort Matrix:  " << " D11 = " << m_D_fluid[0] << ", D12 = " << m_D_fluid[1] << endl;
     cout << "                         " << " D21 = " << m_D_fluid[2] << ", D22 = " << m_D_fluid[3] << endl;
     cout << "Kinetic Diffusion Coef:  " << " D0  = " << m_D_kinet[0] << endl;
     cout << "Kinetic Advection Coefs: " << " DN0 = " << m_DN0        << ", DT0  = " << m_D_kinet[1] << endl;
     cout << "                         " << " DN2 = " << m_D_kinet[2] << ", DT2  = " << m_D_kinet[3] << endl;
     //cout << "Stable Diffusive time step           = " << dt_stable_diff << endl;
     //cout << "Stable Advective-Diffusive time step = " << dt_stable_adv_diff << endl;
     cout << "Stable time step < " << dt_stable << endl;
     cout << "--------------------------------------------------------------------- " << endl;
   }
}

Real Anomalous::computeStableDt( const KineticSpecies& a_soln_species )
{
   CH_TIMERS("Anomalous::computeStableDt");
   CH_TIMER("MPI_AllReduce", t_MPI_AllReduce);
   CH_TIMER("find_minMax", t_findMinMax);

   // get stuff to calculate stability parameters on time step and number of mu cells
   //const KineticSpecies& a_soln_species( *(soln[0]) );
   const LevelData<FArrayBox>& soln_fB = a_soln_species.distributionFunction();
   const DisjointBoxLayout& dbl = soln_fB.getBoxes();
   const PhaseGeom& phase_geom = a_soln_species.phaseSpaceGeometry();

   // calculate cell-centered beta=U_psi^2/D_psi
   if (!m_beta.isDefined())  m_beta.define(dbl, 1, IntVect::Zero);
   
   if (!m_simple_diffusion) {
      computeBeta(m_beta, a_soln_species);
   }
   else {
     for (DataIterator dit(dbl); dit.ok(); ++dit) {
       m_beta[dit].setVal(1.0e-10);
     }
   }
   
   // get mapped dr (current mappings assume dr to be the same in all blocks)
   const PhaseCoordSys& coord_sys = phase_geom.phaseCoordSys();
   const PhaseBlockCoordSys& block_coord_sys = (const PhaseBlockCoordSys&)(*coord_sys.getCoordSys(0));
   const RealVect& phase_dx =  block_coord_sys.dx();
   Real dr = phase_dx[RADIAL_DIR];
   
   // get r lame coefficient
   if (!m_hr.isDefined()) {
     m_hr.define(dbl, 1, IntVect::Zero);
     for (DataIterator dit(dbl); dit.ok(); ++dit) {
       // get cell-centered h_r, h_theta, and h_phi on this patch
       Box phase_box(dbl[dit]);
       FArrayBox lame_cells(phase_box, 3);
       getCellCenteredLameCoefficients(lame_cells, phase_geom, dbl[dit]);
       m_hr[dit].copy(lame_cells,0,0);
     }
   }
 
   // convert fluid matrix components to Anomalous coefficients
   // const Real D2 = 2.0/3.0*D_fluid[2]-D_fluid[0];
   //const Real D0 = 2.0/3.0*D_fluid[3]-D_fluid[1] - 7.0/2.0*D2;
   //const Real D1 = 2.0/3.0*D_fluid[3]-D_fluid[1] + D_fluid[0] - four/three*D_fluid[2];
   //const Real D3 = 3.0/2.0*D2 - D_fluid[1];
   //const Real D4 = -2.0*D2;
   //const Real D5 = -4.0/3.0*( 2.0/3.0*D_fluid[3]-D_fluid[1] );

   // get the min value of hr and Dpsi and max values of Dpsi and Upsi
   Real local_min_hr(1000);  // 1000 is just a starting point
   Real local_max_beta(-1000);

   CH_START(t_findMinMax);
   for (DataIterator dit(dbl); dit.ok(); ++dit) {
     Box box( dbl[dit] );
     Real box_min_hr( m_hr[dit].min (box) );
     Real box_max_beta( m_beta[dit].max (box) );
     local_min_hr = Min( local_min_hr, box_min_hr );
     local_max_beta = Max( local_max_beta, box_max_beta );
   }
   CH_STOP(t_findMinMax);
   
   CH_START(t_MPI_AllReduce);
   Real min_hr( local_min_hr );
   Real max_beta( local_max_beta );
#ifdef CH_MPI
   MPI_Allreduce( &local_min_hr, &min_hr, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD );
   MPI_Allreduce( &local_max_beta, &max_beta, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
#endif
   CH_STOP(t_MPI_AllReduce);

   // calculate min time step set by purely diffusive parts
   // dt*sum_i( D_i/dx_i^2 ) << 1/2
   const Real alpha_r = m_D_kinet[0]/(min_hr*min_hr*dr*dr);
   const Real dt_stable_diff = 0.5/alpha_r;

   // calculate min time step set by advective-diffusive
   // dt*sum_i( U_i^2/D_i ) << 2
   //const Real beta_r = max_Upsi*max_Upsi/min_Dpsi;
   const Real dt_stable_adv_diff = 2.0/max_beta;
   
   // take mininum stable time step
   Real stableDt = Min(dt_stable_diff, dt_stable_adv_diff);

   // JRA prevent -nan return that occurs before first step and causes code not to run   
   if(stableDt != stableDt) stableDt = 123456.7;

   return stableDt;


}

void
Anomalous::computeBeta(LevelData<FArrayBox>& a_beta,
                       const KineticSpecies& a_soln_species )
{
   /*
    Computes cell-centered beta=U_psi^2/D_psi
    */
   
   // get stuff to calculate stability parameters on time step and number of mu cells
   //const KineticSpecies& a_soln_species( *(soln[0]) );
   const LevelData<FArrayBox>& soln_fB = a_soln_species.distributionFunction();
   const DisjointBoxLayout& dbl = soln_fB.getBoxes();
   const PhaseGeom& phase_geom = a_soln_species.phaseSpaceGeometry();
   const CFG::MagGeom & mag_geom = phase_geom.magGeom();
   const CFG::DisjointBoxLayout& mag_grids = mag_geom.grids();
   const LevelData<FArrayBox>& inj_B = phase_geom.getBFieldMagnitude();
   double mass = a_soln_species.mass();
   const ProblemDomain& phase_domain = phase_geom.domain();
   const Box& domain_box = phase_domain.domainBox();
   int num_r_cells = domain_box.size(0);

   // copy const soln_fB (along wiht the boundary values) to a temporary
   if ( !m_fB_beta.isDefined()) m_fB_beta.define( dbl, 1, IntVect::Unit );
   DataIterator sdit = m_fB_beta.dataIterator();
   for (sdit.begin(); sdit.ok(); ++sdit) {
     m_fB_beta[sdit].copy(soln_fB[sdit]);
   }

   // get density, mean parallel velocity, temperature, and fourth moment coefficient
   if ((m_update_freq < 0) || (m_it_counter % m_update_freq == 0) || m_first_call) {
      
      CFG::IntVect	cfg_ghostVect =	CFG::IntVect::Unit;
      CFG::LevelData<CFG::FArrayBox> dens_cfg(mag_grids, 1, cfg_ghostVect);
      CFG::LevelData<CFG::FArrayBox> Upar_cfg(mag_grids, 1, cfg_ghostVect);
      CFG::LevelData<CFG::FArrayBox> temp_cfg(mag_grids, 1, cfg_ghostVect);
      CFG::LevelData<CFG::FArrayBox> four_cfg(mag_grids, 1, cfg_ghostVect);
      CFG::LevelData<CFG::FArrayBox> perp_cfg(mag_grids, 1, cfg_ghostVect);

      m_moment_op.compute(dens_cfg, a_soln_species, m_fB_beta, DensityKernel<FArrayBox>());
      m_moment_op.compute(Upar_cfg, a_soln_species, m_fB_beta, ParallelMomKernel<FArrayBox>());
      m_moment_op.compute(four_cfg, a_soln_species, m_fB_beta, FourthMomentKernel<FArrayBox>());
      m_moment_op.compute(perp_cfg, a_soln_species, m_fB_beta, PerpEnergyKernel<FArrayBox>());
      CFG::DataIterator cfg_dit = dens_cfg.dataIterator();
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
         Upar_cfg[cfg_dit].divide(dens_cfg[cfg_dit]);
      }
      m_moment_op.compute(temp_cfg, a_soln_species, m_fB_beta, PressureKernel<FArrayBox>(Upar_cfg));

      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
         four_cfg[cfg_dit].divide(temp_cfg[cfg_dit]); // dividing by pressure
         perp_cfg[cfg_dit].divide(temp_cfg[cfg_dit]); // dividing by pressure
      }
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
         temp_cfg[cfg_dit].divide(dens_cfg[cfg_dit]);
      }
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
         four_cfg[cfg_dit].divide(temp_cfg[cfg_dit]); // dividing by temperature
      }

      // Perform extrapolation into ghosts to avoid division by zero
      // (can occur for some BC choises)
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
         const CFG::MagBlockCoordSys& cfg_block_coord_sys = mag_geom.getBlockCoordSys(mag_grids[cfg_dit]);
         const CFG::ProblemDomain& cfg_domain = cfg_block_coord_sys.domain();
                  
         secondOrderCellExtrapAtDomainBdry(dens_cfg[cfg_dit], mag_grids[cfg_dit], cfg_domain);
         secondOrderCellExtrapAtDomainBdry(Upar_cfg[cfg_dit], mag_grids[cfg_dit], cfg_domain);
         secondOrderCellExtrapAtDomainBdry(temp_cfg[cfg_dit], mag_grids[cfg_dit], cfg_domain);
         secondOrderCellExtrapAtDomainBdry(four_cfg[cfg_dit], mag_grids[cfg_dit], cfg_domain);
         secondOrderCellExtrapAtDomainBdry(perp_cfg[cfg_dit], mag_grids[cfg_dit], cfg_domain);
      }
      mag_geom.fillInternalGhosts( dens_cfg );
      mag_geom.fillInternalGhosts( Upar_cfg );
      mag_geom.fillInternalGhosts( temp_cfg );
      mag_geom.fillInternalGhosts( four_cfg );
      mag_geom.fillInternalGhosts( perp_cfg );
      
      phase_geom.injectConfigurationToPhase(dens_cfg, m_density_beta);
      phase_geom.injectConfigurationToPhase(Upar_cfg, m_Upar_beta);
      phase_geom.injectConfigurationToPhase(temp_cfg, m_temperature_beta);
      phase_geom.injectConfigurationToPhase(four_cfg, m_fourth_coef_beta);
      phase_geom.injectConfigurationToPhase(perp_cfg, m_perp_coef_beta);
   }
   
   DataIterator bdit = a_beta.dataIterator();
   for (bdit.begin(); bdit.ok(); ++bdit)
   {
     // get phase space dx on this patch
     const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(dbl[bdit]);
     const RealVect& phase_dx =  block_coord_sys.dx();

     const FArrayBox& B_on_patch  = inj_B[bdit];
     const FArrayBox& N_on_patch  = m_density_beta[bdit];
     const FArrayBox& U_on_patch  = m_Upar_beta[bdit];
     const FArrayBox& T_on_patch  = m_temperature_beta[bdit];
     const FArrayBox& C_on_patch  = m_fourth_coef_beta[bdit];
     const FArrayBox& P_on_patch  = m_perp_coef_beta[bdit];

     FORT_EVAL_BETA( CHF_BOX(a_beta[bdit].box()),
                     CHF_CONST_REALVECT(phase_dx),
                     CHF_CONST_VR(m_D_kinet),
                     CHF_CONST_REAL(mass),
                     CHF_CONST_INT(num_r_cells),
                     CHF_CONST_FRA1(B_on_patch,0),
                     CHF_CONST_FRA1(N_on_patch,0),
                     CHF_CONST_FRA1(U_on_patch,0),
                     CHF_CONST_FRA1(T_on_patch,0),
                     CHF_CONST_FRA1(C_on_patch,0),
                     CHF_CONST_FRA1(P_on_patch,0),
                     CHF_FRA1(a_beta[bdit],0) );
   }

}

void Anomalous::preTimeStep(const KineticSpeciesPtrVect& a_soln_mapped,
                            const int a_species,
                            const Real a_time,
                            const KineticSpeciesPtrVect& a_soln_physical )
{
   m_it_counter+=1;
   m_first_stage = true;
}

#include "NamespaceFooter.H"
