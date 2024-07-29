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
      m_update_cls_freq(1),
      m_update_operator_coeff(1),
      m_nbands(5),
      m_my_pc_idx(-1),
      m_skip_stage_update(false),
      m_time_implicit(true),
      m_diagnostics(false),
      m_verbosity(a_verbosity),
      m_cls_freq_func(NULL)
{
   
  m_opt_string = a_ppcls_str;
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
    df/dt_coll=div(flux_coll), where flux_coll = coll_freq*[(v-U)*f+vth^2*dfdv].
    Note that this operator exactly conservs density, momentum,
    and thermal energy while also relaxing the species to a Maxwellian state.
    Details of the implementation can be found in Justin's thesis (Chapter 4).
  */

  CH_TIME("ConsDragDiff::evalClsRHS");
  
  // get vlasov RHS for the current species
  KineticSpecies& rhs_species( *(a_rhs[a_species]) );
  LevelData<FArrayBox>& rhs_dfn = rhs_species.distributionFunction();

  // get solution distribution function (f*J*Bstarpar) for the current species
  const KineticSpecies& soln_species( *(a_soln[a_species]) );
  const LevelData<FArrayBox>& soln_fBJ = soln_species.distributionFunction();
  double mass = soln_species.mass();
  double charge = soln_species.charge();

  // get coordinate system parameters
  const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();
  const CFG::MagGeom & mag_geom = phase_geom.magGeom();
  const DisjointBoxLayout& dbl = soln_fBJ.getBoxes();
  
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
  
  // compute the divergence of individual velocity space fluxes
  // i.e., Psi quantities from Justin's thesis Chapter 4 (multiplied by J)
  // these terms are then combined with self-consistenly computed
  // Upar and Temp that provide discrete conservation of momentum and energy
  if (!m_Jpsi.isDefined()) m_Jpsi.define(dbl, 3, IntVect::Zero);
  computeVelFluxesDiv(m_Jpsi, m_fBJ_vel_ghost, phase_geom, mass);
  
  if (!m_skip_stage_update) {
    // compute conservative values of Upar and Temp
    if ((m_it_counter % m_update_operator_coeff == 0) || m_first_call) {
      computeConsUparTemp(m_Jpsi, soln_species, phase_geom);
    }
    
    // compute conservative values of Upar and Temp
    if ((m_it_counter % m_update_cls_freq == 0) || m_first_call) {
      computeClsFreq(soln_species, phase_geom, dbl, a_time);
    }
  }
    
  // create conservative drag-diffusion collisional RHS
  if (!m_rhs_cls.isDefined()) m_rhs_cls.define(dbl, 1, IntVect::Zero);
  DataIterator cdit = m_rhs_cls.dataIterator();
  for (cdit.begin(); cdit.ok(); ++cdit)
  {
    const FArrayBox& this_Upar = m_Upar[cdit];
    const FArrayBox& this_Temp = m_Temp[cdit];
    const FArrayBox& this_Jpsi = m_Jpsi[cdit];

    FORT_EVAL_CONSDRAGDIFF_RHS( CHF_BOX(m_rhs_cls[cdit].box()),
                                CHF_FRA1(m_rhs_cls[cdit],0),
                                CHF_CONST_FRA1(m_cls_freq[cdit],0),
                                CHF_CONST_FRA1(this_Upar,0),
                                CHF_CONST_FRA1(this_Temp,0),
                                CHF_CONST_FRA(this_Jpsi) );
  }

  // add (or overwrite) collision RHS to Vlasov RHS
  for (cdit.begin(); cdit.ok(); ++cdit)
  {
    if(m_cls_only){
      rhs_dfn[cdit].copy(m_rhs_cls[cdit]);}
    else{
      rhs_dfn[cdit].plus(m_rhs_cls[cdit]);}
  }

  // print parameters (e.g., stable dt) at the first time step
  if ((m_verbosity) && (m_first_call)) {printParameters(a_soln, a_species);}

  // check conservation properties
  if (m_first_call && m_diagnostics) {
    diagnostics(m_rhs_cls, rhs_species, a_time);
    exit(1);
  }
  
  m_first_call = false;
}


void ConsDragDiff::computeVelFluxesDiv(LevelData<FArrayBox>&        a_Jpsi,
                                       const LevelData<FArrayBox>&  a_fBJ,
                                       const PhaseGeom&             a_phase_geom,
                                       const double                 a_mass)
{
  /*This computes Jpsi divergence quantities.
   See Justin's thesis Chapter 4*/
  
  // get injected magnetic field
  const LevelData<FArrayBox>& inj_B = a_phase_geom.getBFieldMagnitude();
  
  // get injected cell-centered velocity spatial normalization
  const int use_spatial_vel_norm = a_phase_geom.spatialVelNorm();
  const LevelData<FArrayBox>& inj_vel_norm = a_phase_geom.getVelNorm();
      
  // get problem domain and number of vpar and mu cells
  const ProblemDomain& phase_domain = a_phase_geom.domain();
  const Box& domain_box = phase_domain.domainBox();
  int num_vp_cells = domain_box.size(VPARALLEL_DIR);
  int num_mu_cells = domain_box.size(MU_DIR);

  // create the temporary fluxes needed to compute the conservative
  // mean velocity and conservative temperature
  const DisjointBoxLayout& dbl = a_fBJ.getBoxes();
  if (!m_fluxes.isDefined())  m_fluxes.define(dbl, 3, IntVect::Zero);
  DataIterator dit = m_fluxes.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    // get phase space dx
    const PhaseBlockCoordSys& block_coord_sys = a_phase_geom.getBlockCoordSys(dbl[dit]);
    const RealVect& phase_dx =  block_coord_sys.dx();
    
    const FArrayBox& fBJ_on_patch = a_fBJ[dit];
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
      FORT_EVAL_CONSDRAGDIFF_FLUX(CHF_BOX(m_fluxes[dit][dir].box()),
                                  CHF_FRA(m_fluxes[dit][dir]),
                                  CHF_CONST_FRA1(fBJ_on_patch,0),
                                  CHF_CONST_FRA1(B_on_patch,0),
                                  CHF_CONST_FRA((*this_velnormptr)),
                                  CHF_CONST_REAL(a_mass),
                                  CHF_CONST_REALVECT(phase_dx),
                                  CHF_CONST_INT(dir),
                                  CHF_CONST_INT(num_vp_cells),
                                  CHF_CONST_INT(num_mu_cells),
                                  CHF_CONST_INT(use_spatial_vel_norm));

    }
  }

  // compute the divergence of the velocity space fluxes
  a_phase_geom.mappedGridDivergenceFromFluxNormals(a_Jpsi, m_fluxes);
  
}


void ConsDragDiff::computeConsUparTemp(LevelData<FArrayBox>&  a_Jpsi,
                                       const KineticSpecies&  a_soln_species,
                                       const PhaseGeom&       a_phase_geom )
{
  
  const CFG::MagGeom & mag_geom = a_phase_geom.magGeom();
  
  // evaluate moments needed to compute conservative Upar and Temperature
  CFG::LevelData<CFG::FArrayBox> vpar_moms(mag_geom.grids(), 3, CFG::IntVect::Zero);
  CFG::LevelData<CFG::FArrayBox> pres_moms(mag_geom.grids(), 3, CFG::IntVect::Zero);

  //need actual Upar for pressure moment
  CFG::LevelData<CFG::FArrayBox> ushift(mag_geom.grids(), 1, CFG::IntVect::Zero);
  a_soln_species.parallelVelocity(ushift);

  moment_op.compute(vpar_moms, a_soln_species, a_Jpsi, ParallelVelKernel<FArrayBox>());
  moment_op.compute(pres_moms, a_soln_species, a_Jpsi, PressureKernel<FArrayBox>(ushift));
  
  // this will also define m_inj_vpar_moms and m_inj_pres_moms during the first call
  a_phase_geom.injectConfigurationToPhase(vpar_moms, m_inj_vpar_moms);
  a_phase_geom.injectConfigurationToPhase(pres_moms, m_inj_pres_moms);

  // create conservative mean parallel velocity and Temperature
  if (!m_Upar.isDefined()) m_Upar.define(m_inj_vpar_moms.getBoxes(), 1, IntVect::Zero);
  if (!m_Temp.isDefined()) m_Temp.define(m_inj_vpar_moms.getBoxes(), 1, IntVect::Zero);
  DataIterator mdit = m_Upar.dataIterator();
  for (mdit.begin(); mdit.ok(); ++mdit)
  {
    const FArrayBox& this_vpar_moms = m_inj_vpar_moms[mdit];
    const FArrayBox& this_pres_moms = m_inj_pres_moms[mdit];

    FORT_EVAL_CONS_UPAR_TEMP( CHF_BOX(m_Upar[mdit].box()),
                              CHF_FRA1(m_Upar[mdit],0),
                              CHF_FRA1(m_Temp[mdit],0),
                              CHF_CONST_FRA(this_vpar_moms),
                              CHF_CONST_FRA(this_pres_moms) );
  }
}


void ConsDragDiff::computeClsFreq(const KineticSpecies&     a_soln_species,
                                  const PhaseGeom&          a_phase_geom,
                                  const DisjointBoxLayout&  a_grids,
                                  const Real                a_time)
{
  
  // define collision frequency
  if (!m_cls_freq.isDefined()) m_cls_freq.define(a_grids, 1, IntVect::Zero);
  
  const CFG::MagGeom & mag_geom = a_phase_geom.magGeom();
  
  // set prescribed collision frequency
  if (m_cls_freq_func != NULL && m_first_call) {
    CFG::LevelData<CFG::FArrayBox> cls_freq_cfg(mag_geom.grids(), 1, CFG::IntVect::Zero);
    m_cls_freq_func->assign( cls_freq_cfg, mag_geom, a_time);
    a_phase_geom.injectAndExpandConfigurationToPhase(cls_freq_cfg, m_cls_freq);
  }
  
  // compute self-consistent collision frequency
  if (m_cls_freq_func == NULL) {
    
    CFG::LevelData<CFG::FArrayBox> density(mag_geom.grids(), 1, CFG::IntVect::Zero);
    a_soln_species.numberDensity(density);
    mag_geom.divideJonValid(density);
    a_phase_geom.injectConfigurationToPhase(density, m_dens_inj);

    CFG::LevelData<CFG::FArrayBox> temperature(mag_geom.grids(), 1, CFG::IntVect::Zero);
    a_soln_species.temperature(temperature);
    a_phase_geom.injectConfigurationToPhase(temperature, m_temp_inj);
          
    double mass = a_soln_species.mass();
    double charge = a_soln_species.charge();
    
    computeSelfConsistFreq(m_cls_freq, m_dens_inj, m_temp_inj, mass, charge);
  }
  
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


void ConsDragDiff::defineBlockPC( std::vector<Preconditioner<ODEVector,AppCtxt>*>&  a_pc,
                                  std::vector<DOFList>&                             a_dof_list,
                                  const ODEVector&                                  a_X,
                                  void*                                             a_ops,
                                  const std::string&                                a_out_string,
                                  const std::string&                                a_opt_string,
                                  bool                                              a_im,
                                  const KineticSpecies&                             a_species,
                                  const GlobalDOFKineticSpecies&                    a_global_dofs,
                                  const int                                         a_species_index,
                                  const int                                         a_id )
{
  if (a_im && m_time_implicit) {
 
    CH_assert(a_pc.size() == a_dof_list.size());

    if (!procID()) {
      std::cout << "  Kinetic Species " << a_species_index
                << " - ConsDragDiff collision term: "
                << " creating " << _BASIC_GK_PC_ << " preconditioner "
                << " (index = " << a_pc.size() << ").\n";
    }

    Preconditioner<ODEVector,AppCtxt> *pc;
    pc = new BasicGKPreconditioner<ODEVector,AppCtxt>;
    pc->setSysID(a_id);
    DOFList dof_list(0);

    const LevelData<FArrayBox>& soln_dfn    (a_species.distributionFunction());
    const DisjointBoxLayout&    grids       (soln_dfn.disjointBoxLayout());
    const int                   n_comp      (soln_dfn.nComp());
    const LevelData<FArrayBox>& pMapping    (a_global_dofs.data());
  
    for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const Box& grid = grids[dit];
      const FArrayBox& pMap = pMapping[dit];
      for (BoxIterator bit(grid); bit.ok(); ++bit) {
        IntVect ic = bit();
        for (int n(0); n < n_comp; n++) {
          dof_list.push_back((int) pMap.get(ic ,n) - a_global_dofs.mpiOffset());
        }
      }
    }

    m_my_pc_idx = a_pc.size();

    {
      BasicGKPreconditioner<ODEVector,AppCtxt>* this_pc
        = dynamic_cast< BasicGKPreconditioner<ODEVector,AppCtxt>* > (pc);
      this_pc->define(a_X, a_ops, m_opt_string, m_opt_string, a_im, dof_list);
      
      const GlobalDOF* global_dof = a_X.getGlobalDOF();
      int n_local = a_X.getVectorSize();
      int n_total = n_local;
#ifdef CH_MPI
      MPI_Allreduce(  &n_local,
                      &n_total,
                      1,
                      MPI_INT,
                      MPI_SUM,
                      MPI_COMM_WORLD );
#endif
      if (!procID()) {
        std::cout << "    Setting up banded matrix with " << n_total << " rows ";
        std::cout << "and " << m_nbands << " bands for the preconditioner.\n";
      }
      BandedMatrix& pc_mat(this_pc->getBandedMatrix());
      pc_mat.define(n_local, m_nbands, global_dof->mpiOffset());
    }

    a_pc.push_back(pc);
    a_dof_list.push_back(dof_list);

  }

  return;
}

void ConsDragDiff::updateBlockPC( std::vector<Preconditioner<ODEVector,AppCtxt>*>&  a_pc,
                                  const KineticSpecies&                             a_species,
                                  const GlobalDOFKineticSpecies&                    a_global_dofs,
                                  const Real                                        a_time,
                                  const int                                         a_step,
                                  const int                                         a_stage,
                                  const Real                                        a_shift,
                                  const bool                                        a_im,
                                  const int                                         a_species_index )
{
  if (a_im && m_time_implicit) {
    CH_assert(m_my_pc_idx >= 0);
    CH_assert(a_pc.size() > m_my_pc_idx);

    if (!procID()) {
      std::cout << "    ==> Updating " <<_BASIC_GK_PC_ << " preconditioner "
                << " for FP collision term of kinetic species " << a_species_index << ".\n";
    }

    BasicGKPreconditioner<ODEVector,AppCtxt> *pc
      = dynamic_cast<BasicGKPreconditioner<ODEVector,AppCtxt>*>
        (a_pc[m_my_pc_idx]);
    CH_assert(pc != NULL);

    BandedMatrix& pc_mat(pc->getBandedMatrix());
    pc_mat.setToIdentityMatrix();
    assemblePrecondMatrix(  pc_mat,
                            a_species,
                            a_global_dofs,
                            a_step,
                            a_stage,
                            a_shift);
    pc_mat.finalAssembly();
  }
  return;
}

void ConsDragDiff::assemblePrecondMatrix( BandedMatrix&                   a_P,
                                          const KineticSpecies&           a_species,
                                          const GlobalDOFKineticSpecies&  a_global_dofs,
                                          const int                       a_step,
                                          const int                       a_stage,
                                          const Real                      a_shift )
{

  if (m_time_implicit) {

    const LevelData<FArrayBox>& soln_dfn    (a_species.distributionFunction());
    const DisjointBoxLayout&    grids       (soln_dfn.disjointBoxLayout());
    const PhaseGeom&            phase_geom  (a_species.phaseSpaceGeometry());
    const int                   n_comp      (soln_dfn.nComp());
    const VEL::VelCoordSys&     vel_coords  (phase_geom.velSpaceCoordSys());
    const VEL::RealVect&        vel_dx      (vel_coords.dx());
    const double                mass        (a_species.mass());
    const LevelData<FArrayBox>& pMapping    (a_global_dofs.data());
    
    Real  dv = 1.0/vel_dx[0], dmu = 1.0/vel_dx[1], dv_sq = dv*dv, dmu_sq = dmu*dmu;

    // get injected B-field
    const LevelData<FArrayBox>& inj_B = phase_geom.getBFieldMagnitude();
    
    // get injected cell-centered velocity spatial normalization
    const int use_spatial_vel_norm = phase_geom.spatialVelNorm();
    //const LevelData<FArrayBox>& inj_vel_norm = phase_geom.getVelNorm();

    // iterate over box to set PC matrix coefficients
    DataIterator dit = grids.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
      const Box& grid = grids[dit];
      const FArrayBox& pMap = pMapping[dit];
  
      const Box& Bbox = inj_B[dit].box();
      int vp_inj = Bbox.smallEnd(VPARALLEL_DIR);
      int mu_inj = Bbox.smallEnd(MU_DIR);
              
      BoxIterator bit(grid);
      for (bit.begin(); bit.ok(); ++bit) {
        /* this point */
        IntVect ic = bit();
        /* neighboring points */
        IntVect ie(ic);
        IntVect iw(ic);
        IntVect in(ic);
        IntVect is(ic);
        /* north-south is along mu; east-west is along v|| */
        ie[VPARALLEL_DIR]++;                     /* east  */
        iw[VPARALLEL_DIR]--;                     /* west  */
        in[MU_DIR]++;                            /* north */
        is[MU_DIR]--;                            /* south */
  
        /* iv index for injected quantities*/
        IntVect iv_inj(ic);
        iv_inj[VPARALLEL_DIR] = vp_inj;
        iv_inj[MU_DIR] = mu_inj;
        
        /* collision frequency */
        Real nu = m_cls_freq[dit](iv_inj);
        //Real spatial_vel_norm = inj_vel_norm[dit](iv_inj);
        
        Real DvL, DvR, DmuL, DmuR, DvvL, DvvR, DmumuL, DmumuR;

        /* Dv = nu * (vpar - m_Upar) */
        Real u_par = m_Upar[dit](iv_inj);
        Real vparL = ic[VPARALLEL_DIR]*vel_dx[0];
        Real vparR = ie[VPARALLEL_DIR]*vel_dx[0];
        DvL = nu * (vparL - u_par);
        DvR = nu * (vparR - u_par);
        
        /* Dmu = 2.0 * nu * mu*/
        Real muL = ic[MU_DIR]*vel_dx[1];
        Real muR = in[MU_DIR]*vel_dx[1];
        DmuL = 2.0 * nu * muL;
        DmuR = 2.0 * nu * muR;

        /* Dvv = nu * m_Temp/mass */
        Real temperature = m_Temp[dit](iv_inj);
        DvvL = nu * temperature/mass;
        DvvR = nu * temperature/mass;
        
        /* Dmumu = 4 * nu * mu * m_Temp/Bmag */
        Real Bmag = inj_B[dit](iv_inj);
        DmumuL = 4.0 * nu * muL * temperature / Bmag;
        DmumuR = 4.0 * nu * muR * temperature / Bmag;
  
        for (int n(0); n < n_comp; n++) {
          /* global row/column numbers */
          int pc, pe, pw, pn, ps, pne, pnw, pse, psw;//, pee, pww, pnn, pss;
          pc  = (int) pMap.get(ic ,n);
          pn  = (int) pMap.get(in ,n);
          ps  = (int) pMap.get(is ,n);
          pe  = (int) pMap.get(ie ,n);
          pw  = (int) pMap.get(iw ,n);
  
          /* coefficients */
          
          Real ac = 0;
          Real an = 0;
          Real as = 0;
          Real ae = 0;
          Real aw = 0;

          /* Advection terms along vpar */
          ae += 0.5 * DvR * dv;
          aw -= 0.5 * DvL * dv;
          ac += 0.5 * (-DvL + DvR) * dv;
          
          /* Advection terms along mu */
          an += 0.5 * DmuR * dmu;
          as -= 0.5 * DmuL * dmu;
          ac += 0.5 * (-DmuL + DmuR) * dmu;

          /* Laplacian term along vpar */
          ae += DvvR * dv_sq;
          aw += DvvL * dv_sq;
          ac -= (DvvR + DvvL) * dv_sq;

          /* Laplacian term along mu */
          an += DmumuR * dmu_sq;
          as += DmumuL * dmu_sq;
          ac -= (DmumuR + DmumuL) * dmu_sq;

          int  ncols = m_nbands, ix = 0;
          int  *icols = (int*)  calloc (ncols,sizeof(int));
          Real *data  = (Real*) calloc (ncols,sizeof(Real));
  
          /* center element */
          icols[ix] = pc;
          data[ix] = a_shift - ac;
          ix++;
  
          /* east element */
          if (pe >= 0) {
            icols[ix] = pe;
            data[ix] = -ae;
            ix++;
          }
          /* west element */
          if (pw >= 0) {
            icols[ix] = pw;
            data[ix] = -aw;
            ix++;
          }
          /* north element */
          if (pn >= 0) {
            icols[ix] = pn;
            data[ix] = -an;
            ix++;
          }
          /* south element */
          if (ps >= 0) {
            icols[ix] = ps;
            data[ix] = -as;
            ix++;
          }
  
          CH_assert(ix <= m_nbands);
          CH_assert(ix <= a_P.getNBands());
          a_P.setRowValues(pc,ix,icols,data);
          free(data);
          free(icols);
        }
      }
    }
  }
  return;
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

  a_pp.query("update_cls_freq", m_update_cls_freq );
  a_pp.query("update_operator_coeff", m_update_operator_coeff );
  a_pp.query("skip_stage_update", m_skip_stage_update );
  
}

inline void ConsDragDiff::printParameters(const KineticSpeciesPtrVect&  a_soln_mapped,
                                          const int                     a_species)
{
   // get stuff to calculate stability parameters on time step
   const KineticSpecies& soln_species( *(a_soln_mapped[a_species]) );
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
   const Real beta_vp = max_Nu*(vp_max-min_Upar)*(vp_max-min_Upar)/(min_Temp/mass);
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
  
  /*
   The collision operator is nonlinear: both {Upar, T} coefficients and
   the self-consisten collision frequency depends of dfn. The skip_stage_update
   option pre-computes Upar, T, and cls_freq in preTimeStep() and freezes them
   during a time step, thereby making the operator "linear" for implicit
   solve puropses. This, however, may degrade conservative and time integration
   properties.
   
   If skip_stage_update is not used that may cause issues during the implicit solve.
   Indeed, at some iterations we can have a largely negative dfn, which can result in negative
   values of the temeprature. In this case the evaluation of a self_consistent frequency
   that involves T^(3/2) will cause the code to crash. N.B. we can think of introducing
   abs(T^(3/2)) in the self-consisistent frequency evaluation to possibly mitigate this issue.
   */
  
  if (m_skip_stage_update) {
  
    // get solution distribution function (f*J*Bstarpar) for the current species
    const KineticSpecies& soln_species( *(a_soln_mapped[a_species]) );
    const LevelData<FArrayBox>& soln_fBJ = soln_species.distributionFunction();
    double mass = soln_species.mass();

    // get coordinate system parameters
    const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();
    const DisjointBoxLayout& dbl = soln_fBJ.getBoxes();
    
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
    
    // compute the divergence of individual velocity space fluxes
    // i.e., Psi quantities from Justin's thesis Chapter 4 (multiplied by J)
    // these terms are then combined with self-consistenly computed
    // Upar and Temp that provide discrete conservation of momentum and energy
    if (!m_Jpsi.isDefined()) m_Jpsi.define(dbl, 3, IntVect::Zero);
    computeVelFluxesDiv(m_Jpsi, m_fBJ_vel_ghost, phase_geom, mass);
    
    // compute conservative values of Upar and Temp
    if ((m_it_counter % m_update_operator_coeff == 0) || m_first_call) {
      computeConsUparTemp(m_Jpsi, soln_species, phase_geom);
    }

    // compute conservative values of Upar and Temp
    if ((m_it_counter % m_update_cls_freq == 0) || m_first_call) {
      computeClsFreq(soln_species, phase_geom, dbl, a_time);
    }

    
    // print parameters (e.g., stable dt) at the first time step
    if ((m_verbosity) && (m_first_call)) {printParameters(a_soln_mapped, a_species);}
    
    m_first_call = false;
  }
  
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
