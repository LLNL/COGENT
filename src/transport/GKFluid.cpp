#include <math.h>
#include "CONSTANTS.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "MagBlockCoordSys.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "PhaseGeom.H"
#include "PhaseGeomF_F.H"
#include "PhaseBlockCoordSys.H"
#include "GKFluid.H"
#include "TransportF_F.H"

#include "MomentOp.H"
#include "Kernels.H"

#include "NamespaceHeader.H" // has to be the last one

GKFluid::GKFluid( const string& a_species_name, ParmParse& a_pptpm, const int a_verbosity )
    : D_fluid(4,0),
      D_kinet(4,0),
      model_only(false),
      moment_op( MomentOp::instance() ),
      m_first_step(true),
      verbosity(0)
{

   verbosity = a_verbosity;
   pptpm = a_pptpm;
   species_name = a_species_name;

}

GKFluid::~GKFluid()
{
}

void GKFluid::evalTpmRHS( KineticSpeciesPtrVect&        rhs,
                            const KineticSpeciesPtrVect&  soln,
                            const int                     species,
                            const Real                    time )
{
  /*
    Evaluates the anomalous radial flux across magnetic surfaces:
    df/dt_anomalous=-div(flux_anom). Currently, this computes
    second order diffusion of density with constant Diffusion
    coefficient (i.e. dn/dt = -div( -D*grad_r(n)) ) by computing
    the following eq. for fB: dfB/dt = -div( -D*grad_r(fB) )
  */
      // parse the pptpm database for initial condition items for this species
  if ( m_first_step ) {ParseParameters();}

      // print parameters at the first time step
      if ((verbosity) && (m_first_step)) {printParameters(soln);}

      // get vlasov RHS for the current species
      KineticSpecies& rhs_species( *(rhs[species]) );
      LevelData<FArrayBox>& rhs_dfn = rhs_species.distributionFunction();

      // get solution distribution function (f*Bstarpar) for the current species
      const KineticSpecies& soln_species( *(soln[species]) );
      const LevelData<FArrayBox>& soln_fB = soln_species.distributionFunction();
      double mass = soln_species.mass();
      const DisjointBoxLayout& dbl = soln_fB.getBoxes();

      // get coordinate system parameters
      const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();
      const CFG::MultiBlockLevelGeom & mag_geom = phase_geom.magGeom();
      const LevelData<FArrayBox>& inj_B = phase_geom.getBFieldMagnitude();
      const ProblemDomain& phase_domain = phase_geom.domain();
      const Box& domain_box = phase_domain.domainBox();
      int num_r_cells = domain_box.size(0);

      // copy const soln_fB to a temporary 
      // and put boundary values in ghost cells
      const IntVect ghostVect(IntVect::Unit);
      LevelData<FArrayBox> fB(dbl, 1, ghostVect); 
      DataIterator sdit = fB.dataIterator();
      for (sdit.begin(); sdit.ok(); ++sdit) {
        fB[sdit].copy(soln_fB[sdit]);
      }
      fB.exchange();

      // get density, mean parallel velocity, and temperature
      CFG::LevelData<CFG::FArrayBox> dens_cfg(mag_geom.grids(), 1, CFG::IntVect::Unit);
      CFG::LevelData<CFG::FArrayBox> Upar_cfg(mag_geom.grids(), 1, CFG::IntVect::Unit);
      CFG::LevelData<CFG::FArrayBox> temp_cfg(mag_geom.grids(), 1, CFG::IntVect::Unit);
      CFG::LevelData<CFG::FArrayBox> four_cfg(mag_geom.grids(), 1, CFG::IntVect::Unit);
      moment_op.compute(dens_cfg, soln_species, fB, DensityKernel());
      moment_op.compute(Upar_cfg, soln_species, fB, ParallelMomKernel());
      CFG::DataIterator cfg_dit = dens_cfg.dataIterator();
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
        Upar_cfg[cfg_dit].divide(dens_cfg[cfg_dit]);
      }
      moment_op.compute(temp_cfg, soln_species, fB, PressureKernel(Upar_cfg));
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
        temp_cfg[cfg_dit].divide(dens_cfg[cfg_dit]);
      }
      LevelData<FArrayBox> density;
      LevelData<FArrayBox> Upar;
      LevelData<FArrayBox> temperature;
      LevelData<FArrayBox> fourth_coef;
      phase_geom.injectConfigurationToPhase(dens_cfg, density);
      phase_geom.injectConfigurationToPhase(Upar_cfg, Upar);
      phase_geom.injectConfigurationToPhase(temp_cfg, temperature);

      // get face centered metrics h_r, h_theta, and h_phi on each CFG_DIM face
      if (m_first_step) {
        metrics_faces.define(dbl, 3, IntVect::Zero);
        getFaceCenteredMetrics(metrics_faces, phase_geom, dbl);
      }

      // calculate face-averaged D*dfB/dr
      LevelData<FluxBox> fluxA(dbl, 1, IntVect::Zero);
      DataIterator dit = fluxA.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
      {
        // get phase space dx on this patch
        const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(dbl[dit]);
        const RealVect& phase_dx =  block_coord_sys.dx();

        const FArrayBox& fB_on_patch = fB[dit];
        const FArrayBox& B_on_patch  = inj_B[dit];
        const FArrayBox& N_on_patch  = density[dit];
        const FArrayBox& U_on_patch  = Upar[dit];
        const FArrayBox& T_on_patch  = temperature[dit];

        for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& thisfluxA = fluxA[dit][dir];
          const FArrayBox& metrics_on_patch = metrics_faces[dit][dir];

          FORT_EVAL_GKFLUID_FLUX(CHF_CONST_INT(dir),
                          CHF_BOX(thisfluxA.box()),
                          CHF_CONST_REALVECT(phase_dx),
                          CHF_CONST_FRA(metrics_on_patch),
                          CHF_CONST_VR(D_kinet),
                          CHF_CONST_REAL(mass),
                          CHF_CONST_INT(num_r_cells),
                          CHF_CONST_FRA1(fB_on_patch,0),
                          CHF_CONST_FRA1(B_on_patch,0),
                          CHF_CONST_FRA1(N_on_patch,0),
                          CHF_CONST_FRA1(U_on_patch,0),
                          CHF_CONST_FRA1(T_on_patch,0),
                          CHF_FRA1(thisfluxA,0));
        }
      }

      // calculate div(flux)
      phase_geom.averageAtBlockBoundaries(fluxA);
      LevelData<FArrayBox> rhs_transport;
      rhs_transport.define(rhs_dfn);
      phase_geom.mappedGridDivergenceFromFluxNormals(rhs_transport, fluxA);
      DataIterator rdit = rhs_transport.dataIterator();

      // add (or overwrite) transport RHS to Vlasov RHS
      for (rdit.begin(); rdit.ok(); ++rdit)
      {
        if(model_only){
        rhs_dfn[rdit].copy(rhs_transport[rdit]);}
        else{
        rhs_dfn[rdit].plus(rhs_transport[rdit]);}
      }

  m_first_step = false;  

}

void GKFluid::getFaceCenteredMetrics( LevelData<FluxBox>&       metrics_faces,
                                        const PhaseGeom&          a_phase_geom,
                                        const DisjointBoxLayout&  a_dbl )
{
   // FluxBox to be filled needs 3 components on each face (hr, htheta, and hphi)
   CH_assert(metrics_faces.nComp() == 3);
   DataIterator dit = metrics_faces.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {

   // get face centered N components on this patch to formulate face centered metrics h_r and h_theta
   const CFG::MagBlockCoordSys& mag_block_coord_sys = a_phase_geom.getMagBlockCoordSys(a_dbl[dit]);
   const CFG::RealVect& real_dx = mag_block_coord_sys.dx();
   Box phase_box(a_dbl[dit]);
   CFG::Box box_config;
   a_phase_geom.projectPhaseToConfiguration(phase_box, box_config);
   box_config.grow(2);
   CFG::FluxBox N_cfg_face(box_config, CFG_DIM*CFG_DIM);  //4 components of N on each CFG_DIM face
   mag_block_coord_sys.getPointwiseN(N_cfg_face);
   FluxBox N_face;
   a_phase_geom.injectConfigurationToPhase(N_cfg_face, N_face);

   // get face-centered 2*pi*Rmaj on this patch (note that 2piRmaj=h_toroidalangle)
   CFG::FluxBox TwoPiR_cfg_face(box_config, 1);
   for (int dir=0; dir<CFG_DIM; ++dir)
   {
     CFG::RealVect offset = real_dx;
     offset *= 0.5;
     offset[dir] = 0.0;
     CFG::FArrayBox& this_fab = TwoPiR_cfg_face[dir];
     CFG::BoxIterator bit(this_fab.box());
     for (bit.begin(); bit.ok(); ++bit)
     {
       CFG::IntVect iv = bit();
       CFG::RealVect mapped_loc = iv*real_dx + offset;
       double TwoPiRmaj = 2. * Pi * mag_block_coord_sys.majorRadius(mapped_loc);
       this_fab(iv) = TwoPiRmaj;
     }
   }
   FluxBox TwoPiR_face;
   a_phase_geom.injectConfigurationToPhase(TwoPiR_cfg_face,TwoPiR_face);

   for (int dir=0; dir<CFG_DIM; dir++)
   {
      FArrayBox& metrics_on_patch = metrics_faces[dit][dir];

      FORT_METRICS_FACES( CHF_BOX(metrics_on_patch.box()),
                         CHF_CONST_INT(dir),
                         CHF_CONST_FRA1(TwoPiR_face[dir],0),
                         CHF_CONST_FRA(N_face[dir]),
                         CHF_FRA(metrics_on_patch));
   }
   }

}

void GKFluid::getCellCenteredMetrics( FArrayBox&        metrics_cells,
                                        const PhaseGeom&  a_phase_geom,
                                        const Box&        a_dbl )
{
   // get cell centered N components on this patch to formulate cell centered metrics h_r and h_theta
   const CFG::MagBlockCoordSys& mag_block_coord_sys = a_phase_geom.getMagBlockCoordSys(a_dbl);
   const CFG::RealVect& real_dx = mag_block_coord_sys.dx();
   Box phase_box(a_dbl);
   CFG::Box box_config;
   a_phase_geom.projectPhaseToConfiguration(phase_box, box_config);
   box_config.grow(2);
   CFG::FArrayBox N_cfg_cent(box_config, CFG_DIM*CFG_DIM);  //4 components of N on each cell center
   mag_block_coord_sys.getPointwiseN(N_cfg_cent);
   FArrayBox N_cent;
   a_phase_geom.injectConfigurationToPhase(N_cfg_cent, N_cent);

   // get cell centered 2*pi*Rmaj on this patch (note that 2piRmaj=h_toroidalangle)
   CFG::FArrayBox TwoPiR_cfg_cent(box_config, 1);
   CFG::BoxIterator bit(TwoPiR_cfg_cent.box());
   CFG::RealVect offset = real_dx;
   offset *= 0.5;
   for (bit.begin(); bit.ok(); ++bit)
   {
     CFG::IntVect iv = bit();
     CFG::RealVect mapped_loc = iv*real_dx + offset;
     double TwoPiRmaj = 2. * Pi * mag_block_coord_sys.majorRadius(mapped_loc);
     TwoPiR_cfg_cent(iv) = TwoPiRmaj;
   }
   FArrayBox TwoPiR_cent;
   a_phase_geom.injectConfigurationToPhase(TwoPiR_cfg_cent,TwoPiR_cent);

   FORT_METRICS_CELLS( CHF_BOX(metrics_cells.box()),
                       CHF_CONST_FRA1(TwoPiR_cent,0),
                       CHF_CONST_FRA(N_cent),
                       CHF_FRA(metrics_cells));

}

void GKFluid::ParseParameters()

{

   if ( procID()==0 && pptpm.contains("verbose") ) {
     pptpm.get("verbose", verbosity);
   }
   if (pptpm.contains("model_only") ) {
     pptpm.get("model_only",model_only);
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
     pptpm.getarr("D_fluid",D_fluid,0,4);
     // convert fluid matrix components to kinetic coefficients
     D_kinet[2] = 2.0/3.0*D_fluid[2]-D_fluid[0];
     D_kinet[0] = 2.0/3.0*D_fluid[3]-D_fluid[1] - 7.0/2.0*D_kinet[2];
     D_kinet[1] = 2.0/3.0*D_fluid[3]-D_fluid[1] + D_fluid[0] - four/three*D_fluid[2];
     D_kinet[3] = 3.0/2.0*D_kinet[2]-D_fluid[1];

     // set D0 and D2 to zero exactly if they should be
     if ( abs(D_kinet[0]) <= 1e-8 ) {
       D_kinet[0] = 0.0;
     }
     if ( abs(D_kinet[2]) <= 1e-8 ) {
       D_kinet[2] = 0.0;
     }
   }
   if (pptpm.contains("D_kinet")) {
     pptpm.getarr("D_kinet",D_kinet,0,4);
     // convert kinetic matrix components to fluid coefficients
     D_fluid[0] = D_kinet[0]-D_kinet[1] + 3.0/2.0*D_kinet[2];
     D_fluid[1] = 3.0/2.0*D_kinet[2] - 1.0*D_kinet[3];
     D_fluid[2] = 3.0/2.0*(D_kinet[0]-D_kinet[1]) + 15.0/4.0*D_kinet[2];
     D_fluid[3] = 3.0/2.0*(D_kinet[0]-D_kinet[3]) + 15.0/2.0*D_kinet[2];
   }

   // make sure that transport parameters are such that GKFluid
   // Diff. coeff: D_a = D0 +D2*v^2/v_th^2 is globally positive
   if ( procID()==0 && !(D_kinet[0] >= 0) ) {
     cout << "DIFFUSION COEFFICIENT D0 FOR SPECIES " << species_name << " IS < ZERO " << endl;
     if (pptpm.contains("D_fluid")) {
       MayDay::Error( "KINETIC DIFFUSION COEFFICIENT D0 = (2/3*D_fluid[3]-D_fluid[1]) - 7/2*(2/3*D_fluid[2]-D_fluid[0]) MUST BE >= ZERO " );
     }
     else {
       MayDay::Error( "KINETIC DIFFUSION COEFFICIENT D0 = D_kinet[0] MUST BE >= ZERO " );
     }
   }
   if ( procID()==0 && (D_kinet[0] >=0) && !(D_kinet[2] >= 0)) {
     cout << "DIFFUSION COEFFICIENT D2 FOR SPECIES " << species_name << " IS < ZERO " << endl;
     if (pptpm.contains("D_fluid")) {
       MayDay::Error( "KINETIC DIFFUSION COEFFICIENT D2 = 2/3*D_fluid[2]-D_fluid[0] MUST BE >= ZERO " );
     }
     else {
       MayDay::Error( "KINETIC DIFFUSION COEFFICIENT D2 = D_kinet[2] MUST BE >= ZERO " );
     }
   }
   if ( procID()==0 && (D_kinet[0]+D_kinet[2] == 0)) {
     cout << "GLOBAL DIFFUSION COEFFICIENT D0+D2 FOR " << species_name << " IS  ZERO " << endl;
     MayDay::Error( "PURELY ADVECTIVE SCHEME  WITH FINITE CENTRAL DIFFERENCE IS UNSTABLE " );
   }

}

inline void GKFluid::printParameters( const KineticSpeciesPtrVect& soln )
{
   // get stuff to calculate stability parameters on time step and number of mu cells
   const KineticSpecies& soln_species( *(soln[0]) );
   const LevelData<FArrayBox>& soln_fB = soln_species.distributionFunction();
   const DisjointBoxLayout& dbl = soln_fB.getBoxes();
   const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();
   const CFG::MultiBlockLevelGeom & mag_geom = phase_geom.magGeom();
   const LevelData<FArrayBox>& inj_B = phase_geom.getBFieldMagnitude();
   double mass = soln_species.mass();
   const ProblemDomain& phase_domain = phase_geom.domain();
   const Box& domain_box = phase_domain.domainBox();
   int num_r_cells = domain_box.size(0);

   // copy const soln_fB to a temporary
   const IntVect ghostVect(IntVect::Unit);
   LevelData<FArrayBox> fB(dbl, 1, ghostVect); 
   DataIterator sdit = fB.dataIterator();
   for (sdit.begin(); sdit.ok(); ++sdit) {
     fB[sdit].copy(soln_fB[sdit]);
   }


   // get density, mean parallel velocity, and temperature
   CFG::LevelData<CFG::FArrayBox> dens_cfg(mag_geom.grids(), 1, CFG::IntVect::Unit);
   CFG::LevelData<CFG::FArrayBox> Upar_cfg(mag_geom.grids(), 1, CFG::IntVect::Unit);
   CFG::LevelData<CFG::FArrayBox> temp_cfg(mag_geom.grids(), 1, CFG::IntVect::Unit);
   CFG::LevelData<CFG::FArrayBox> four_cfg(mag_geom.grids(), 1, CFG::IntVect::Unit);
   moment_op.compute(dens_cfg, soln_species, fB, DensityKernel());
   moment_op.compute(Upar_cfg, soln_species, fB, ParallelMomKernel());
   CFG::DataIterator cfg_dit = dens_cfg.dataIterator();
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
     Upar_cfg[cfg_dit].divide(dens_cfg[cfg_dit]);
   }
   moment_op.compute(temp_cfg, soln_species, fB, PressureKernel(Upar_cfg));
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
     temp_cfg[cfg_dit].divide(dens_cfg[cfg_dit]);
   }

   LevelData<FArrayBox> density;
   LevelData<FArrayBox> Upar;
   LevelData<FArrayBox> temperature;
   phase_geom.injectConfigurationToPhase(dens_cfg, density);
   phase_geom.injectConfigurationToPhase(Upar_cfg, Upar);
   phase_geom.injectConfigurationToPhase(temp_cfg, temperature);

   // calculate cell-centered D_psi and beta=U_psi^2/D_psi
   LevelData<FArrayBox> Dpsi(dbl, 1, IntVect::Zero);
   LevelData<FArrayBox> beta(dbl, 1, IntVect::Zero);
   DataIterator D_dit = Dpsi.dataIterator();
   for (D_dit.begin(); D_dit.ok(); ++D_dit)
   {
     // get phase space dx on this patch
     const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(dbl[D_dit]);
     const RealVect& phase_dx =  block_coord_sys.dx();

     const FArrayBox& B_on_patch  = inj_B[D_dit];
     const FArrayBox& N_on_patch  = density[D_dit];
     const FArrayBox& U_on_patch  = Upar[D_dit];
     const FArrayBox& T_on_patch  = temperature[D_dit];

     FORT_EVAL_DPSI_AND_BETA( CHF_BOX(Dpsi[D_dit].box()),
                        CHF_CONST_REALVECT(phase_dx),
                        CHF_CONST_VR(D_kinet),
                        CHF_CONST_REAL(mass),
                        CHF_CONST_INT(num_r_cells),
                        CHF_CONST_FRA1(B_on_patch,0),
                        CHF_CONST_FRA1(N_on_patch,0),
                        CHF_CONST_FRA1(U_on_patch,0),
                        CHF_CONST_FRA1(T_on_patch,0),
                        CHF_FRA1(Dpsi[D_dit],0),
                        CHF_FRA1(beta[D_dit],0) );
   }

   // get r metric hr
   LevelData<FArrayBox> hr(dbl, 1, IntVect::Zero);
   DataIterator dit= hr.dataIterator();
   Real dr;
   for (dit.begin(); dit.ok(); ++dit)
   {
     const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(dbl[dit]);
     const RealVect& phase_dx =  block_coord_sys.dx();
     dr = phase_dx[0];

     // get face centered metrics h_r, h_theta, and h_phi on this patch
     Box phase_box(dbl[dit]);
     FArrayBox metrics_cells(phase_box, 3);
     getCellCenteredMetrics(metrics_cells, phase_geom, dbl[dit]);
     hr[dit].copy(metrics_cells,0,0);

   }
   hr.exchange();

   // convert fluid matrix components to GKFluid coefficients
   // const Real D2 = 2.0/3.0*D_fluid[2]-D_fluid[0];
   //const Real D0 = 2.0/3.0*D_fluid[3]-D_fluid[1] - 7.0/2.0*D2;
   //const Real D1 = 2.0/3.0*D_fluid[3]-D_fluid[1] + D_fluid[0] - four/three*D_fluid[2];
   //const Real D3 = 3.0/2.0*D2 - D_fluid[1];
   //const Real D4 = -2.0*D2;
   //const Real D5 = -4.0/3.0*( 2.0/3.0*D_fluid[3]-D_fluid[1] );

   // get the min value of hr and Dpsi and max values of Dpsi and Upsi
   Real local_min_hr(1000);  // 1000 is just a starting point
   Real local_max_Dpsi(0);
   Real local_max_beta(-1000);

   for (dit.begin(); dit.ok(); ++dit)
   {
     Box box( dbl[dit] );
     Real box_min_hr( hr[dit].min (box) );
     Real box_max_Dpsi( Dpsi[dit].max (box) );
     Real box_max_beta( beta[dit].max (box) );
     local_min_hr = Min( local_min_hr, box_min_hr );
     local_max_Dpsi = Max( local_max_Dpsi, box_max_Dpsi );
     local_max_beta = Max( local_max_beta, box_max_beta );
   }

   Real min_hr( local_min_hr );
   Real max_Dpsi( local_max_Dpsi );
   Real max_beta( local_max_beta );
#ifdef CH_MPI
   MPI_Allreduce( &local_min_hr, &min_hr, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD );
   MPI_Allreduce( &local_max_Dpsi, &max_Dpsi, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
   MPI_Allreduce( &local_max_beta, &max_beta, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
#endif

// cout << "min_Dpsi =" << min_Dpsi << endl;
// cout << "max_Dpsi =" << max_Dpsi << endl;
// cout << "max_beta =" << max_beta << endl;

   // calculate min time step set by purely diffusive parts
   // dt*sum_i( D_i/dx_i^2 ) << 1/2
   const Real alpha_r = max_Dpsi/(min_hr*min_hr*dr*dr);
   const Real dt_stable_diff = 0.5/alpha_r;

   // calculate min time step set by advective-diffusive
   // dt*sum_i( U_i^2/D_i ) << 2
   //const Real beta_r = max_Upsi*max_Upsi/min_Dpsi;
   const Real dt_stable_adv_diff = 2.0/max_beta;

   // take mininum stable time step
   Real dt_stable = Min(dt_stable_diff, dt_stable_adv_diff);

   std::string tpt_only_str("false");
   if (model_only){ tpt_only_str = "true";}

   if (procID() == 0 ) {
     cout << "--------------------------------------------------------------------- " << endl;
     cout << "GKFluid Transport Model being defined for species " << species_name <<endl;
     cout << "tpt_only = " <<  tpt_only_str << endl;
     cout << "Fluid Tranpsort Matrix:  " << " D11 = " << D_fluid[0] << ", D12 = " << D_fluid[1] << endl;
     cout << "                         " << " D21 = " << D_fluid[2] << ", D22 = " << D_fluid[3] << endl;
     cout << "Kinetic Diffusion Coef: " << " D0 = " << D_kinet[0] << ", D2  = " << D_kinet[2] << endl;
     cout << "Kinetic Advection Coef: " << " DN = " << D_kinet[1] << ", DT  = " << D_kinet[3] << endl;
     cout << "Stable Diffusive time step           = " << dt_stable_diff << endl;
     cout << "Stable Advective-Diffusive time step = " << dt_stable_adv_diff << endl;
     cout << "Stable time step < " << dt_stable << endl;
     cout << "--------------------------------------------------------------------- " << endl;
   }
}


#include "NamespaceFooter.H"
