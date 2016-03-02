#include <math.h>

#include "PhaseGeom.H"
#include "PhaseBlockCoordSys.H"
#include "ConsDragDiff.H"
#include "CollisionsF_F.H"

#include "MomentOp.H"
#include "Kernels.H"

#include "NamespaceHeader.H" //Should be the last one


ConsDragDiff::ConsDragDiff( const string& a_species_name, ParmParse& a_ppcls, const int a_verbosity )
    : cls_freq(0),
      cls_only(false),
      moment_op( MomentOp::instance() ),
      verbosity(0)
{
   verbosity = a_verbosity;
   ppcls = a_ppcls;
   species_name = a_species_name;
}

ConsDragDiff::~ConsDragDiff()
{
}

void ConsDragDiff::evalClsRHS( KineticSpeciesPtrVect&        rhs,
                               const KineticSpeciesPtrVect&  soln,
                               const int                     species,
                               const Real                    time,
                               const int                     flag )
{
  /*
    Evaluates the like-like species fully conservative drag-diffusion
    velocity space collision operator: df/dt_coll=div(flux_coll). Where
    flux_coll = coll_freq*[(v-U)*f+vth*dfdv]. Note that this operator
    exactly conservs density, momentum, and thermal energy while also
    relaxing the species to a Maxwellian state.
  */
      // parse the pptpm database for initial condition items for this species
      ParseParameters();

      // print parameters at the first time step
      if ((verbosity) && (time==0)) {printParameters(soln);}

      // get vlasov RHS for the current species
      KineticSpecies& rhs_species( *(rhs[species]) );
      LevelData<FArrayBox>& rhs_dfn = rhs_species.distributionFunction();

      // get solution distribution function (f*J*Bstarpar) for the current species
      const KineticSpecies& soln_species( *(soln[species]) );
      const LevelData<FArrayBox>& soln_fBJ = soln_species.distributionFunction();
      double mass = soln_species.mass();
      const DisjointBoxLayout& dbl = soln_fBJ.getBoxes();

      // copy soln_fBJ so can perform exchange to ghost cells
      LevelData<FArrayBox> copy_fBJ;
      copy_fBJ.define(soln_fBJ);
      copy_fBJ.exchange();

      // get coordinate system parameters
      const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();
      const CFG::MultiBlockLevelGeom & mag_geom = phase_geom.magGeom();
      const LevelData<FArrayBox>& inj_B = phase_geom.getBFieldMagnitude();
      const ProblemDomain& phase_domain = phase_geom.domain();
      const Box& domain_box = phase_domain.domainBox();
      int num_vp_cells = domain_box.size(2);
      int num_mu_cells = domain_box.size(3);

      // create the temporary fluxes needed to compute the conservative
      // mean velocity and conservative temperature
      LevelData<FluxBox> fluxes(dbl, 3, IntVect::Zero);
      DataIterator dit = fluxes.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
      {
        // get phase space dx
        const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(dbl[dit]);
        const RealVect& phase_dx =  block_coord_sys.dx();

        const FArrayBox& fBJ_on_patch = copy_fBJ[dit];
        const FArrayBox& B_on_patch   = inj_B[dit];

        for (int dir=0; dir<SpaceDim; dir++)
        {
          FORT_EVAL_CONSDRAGDIFF_FLUX( CHF_BOX(fluxes[dit][dir].box()),
                                       CHF_FRA(fluxes[dit][dir]),
                                       CHF_CONST_FRA1(fBJ_on_patch,0),
                                       CHF_CONST_FRA1(B_on_patch,0),
                                       CHF_CONST_REAL(mass),
                                       CHF_CONST_REALVECT(phase_dx),
                                       CHF_CONST_INT(dir),
                                       CHF_CONST_INT(num_vp_cells),
                                       CHF_CONST_INT(num_mu_cells) );

        }
      }

      // compute the divergence of the velocity space fluxes
      LevelData<FArrayBox> Jpsi(dbl, 3, IntVect::Zero);
      phase_geom.mappedGridDivergenceFromFluxNormals(Jpsi, fluxes);

      // take momentum needed to compute conservative Upar and Temperature
      CFG::LevelData<CFG::FArrayBox> vpar_moms(mag_geom.grids(), 3, CFG::IntVect::Zero);
      CFG::LevelData<CFG::FArrayBox> pres_moms(mag_geom.grids(), 3, CFG::IntVect::Zero);

      //Introduce a zero-valued vpar-shift, to make Justin old calculations consistent with the 
      //new pressure diagnostics, that now subtract v_par_shift from v_par --MD 02/28/14
      CFG::LevelData<CFG::FArrayBox> vpar_zero(mag_geom.grids(), 3, CFG::IntVect::Zero);
      CFG::DataIterator cfg_dit = vpar_zero.dataIterator();
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit)
      {
        vpar_zero[cfg_dit].setVal(0.0);
      }

      moment_op.compute(vpar_moms, soln_species, Jpsi, ParallelMomKernel());
      moment_op.compute(pres_moms, soln_species, Jpsi, PressureKernel(vpar_zero));
      LevelData<FArrayBox> inj_vpar_moms;
      LevelData<FArrayBox> inj_pres_moms;
      phase_geom.injectConfigurationToPhase(vpar_moms, inj_vpar_moms);
      phase_geom.injectConfigurationToPhase(pres_moms, inj_pres_moms);

      // create conseravetive mean parallel velocity and Temperature
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
                                    CHF_CONST_REAL(cls_freq),
                                    CHF_CONST_FRA1(this_Upar,0),
                                    CHF_CONST_FRA1(this_Temp,0),
                                    CHF_CONST_FRA(this_Jpsi) );
      }

      // add (or overwrite) collision RHS to Vlasov RHS
      for (cdit.begin(); cdit.ok(); ++cdit)
      {
        if(cls_only){
        rhs_dfn[cdit].copy(rhs_cls[cdit]);}
        else{
        rhs_dfn[cdit].plus(rhs_cls[cdit]);}
      }

}

void ConsDragDiff::ParseParameters( ){

   ppcls.query("cls_freq",cls_freq);
   CH_assert(cls_freq >=0 );
   ppcls.query("cls_only",cls_only);

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
   int num_mu_cells = domain_box.size(3);
   int num_vp_cells = domain_box.size(2);

   // get magnetic field, density, and pressure, etc..
   const CFG::MultiBlockLevelGeom & mag_geom = phase_geom.magGeom();
   const LevelData<FArrayBox>& inj_B = phase_geom.getBFieldMagnitude();
   CFG::LevelData<CFG::FArrayBox> dens_cfg(mag_geom.grids(), 1, CFG::IntVect::Zero);
   CFG::LevelData<CFG::FArrayBox> pres_cfg(mag_geom.grids(), 1, CFG::IntVect::Zero);
   CFG::LevelData<CFG::FArrayBox> Upar_cfg(mag_geom.grids(), 1, CFG::IntVect::Zero);

   soln_species.chargeDensity(dens_cfg);

   //Introduce a zero-valued vpar-shift, to make Justin old calculations consistent with the 
   //new pressure diagnostics, that now subtract v_par_shift from v_par --MD 02/28/14
   CFG::LevelData<CFG::FArrayBox> vpar_zero(mag_geom.grids(), 3, CFG::IntVect::Zero);
   CFG::DataIterator cfg_dit = vpar_zero.dataIterator();
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit)
   {
      vpar_zero[cfg_dit].setVal(0.0);
   }

   soln_species.pressureMoment(pres_cfg,vpar_zero);
   soln_species.ParallelMomentum(Upar_cfg);

   CFG::LevelData<CFG::FArrayBox> temp_cfg;
   temp_cfg.define(pres_cfg);

   LevelData<FArrayBox> density;
   LevelData<FArrayBox> pressure;
   LevelData<FArrayBox> Upar;
   phase_geom.injectConfigurationToPhase(dens_cfg, density);
   phase_geom.injectConfigurationToPhase(pres_cfg, pressure);
   phase_geom.injectConfigurationToPhase(Upar_cfg, Upar);

   // create temperature, temperature/B, and Upar
   LevelData<FArrayBox> TonB(dbl, 1, IntVect::Zero);
   LevelData<FArrayBox> Temp(dbl, 1, IntVect::Zero);
   DataIterator bdit= inj_B.dataIterator();
   for (bdit.begin(); bdit.ok(); ++bdit)
   {
     Temp[bdit].copy(pressure[bdit].divide(density[bdit]));
     TonB[bdit].copy(Temp[bdit]);
     TonB[bdit].divide(inj_B[bdit]);
     Upar[bdit].divide(density[bdit]);
   }

   // get velocity cell spacing
   Real dvp;
   Real dmu;
   for (bdit.begin(); bdit.ok(); ++bdit)
   {
     const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(dbl[bdit]);
     const RealVect& phase_dx =  block_coord_sys.dx();
     dvp = phase_dx[2];
     dmu = phase_dx[3];
   }
   const Real mu_max = num_mu_cells*dmu;
   const Real vp_max = num_vp_cells*dvp/2.0;

   // get max values of T, T/B, and Upar
   const DisjointBoxLayout& grids = inj_B.getBoxes();
   Real local_maxTemp(-1000);   // -1000 is just a starting point
   Real local_maxTonB(-1000);   // -1000 is just a starting point
   Real local_minTemp(1000);    // 1000 is just a starting point
   Real local_minTonB(1000);    // 1000 is just a starting point
   Real local_minUpar(1000);    // 1000 is just a starting point
   for (bdit.begin(); bdit.ok(); ++bdit)
   {
      Box box( grids[bdit] );
      Real box_maxTemp( Temp[bdit].max( box ) );
      Real box_maxTonB( TonB[bdit].max( box ) );
      Real box_minTemp( Temp[bdit].min( box ) );
      Real box_minTonB( TonB[bdit].min( box ) );
      Real box_minUpar( Upar[bdit].min( box ) );
      local_maxTemp = Max( local_maxTemp, box_maxTemp );
      local_maxTonB = Max( local_maxTonB, box_maxTonB );
      local_minTemp = Min( local_minTemp, box_minTemp );
      local_minTonB = Min( local_minTonB, box_minTonB );
      local_minUpar = Min( local_minUpar, box_minUpar );
   }
   Real max_Temp( local_maxTemp );
   Real max_TonB( local_maxTonB );
   Real min_Temp( local_minTemp );
   Real min_TonB( local_minTonB );
   Real min_Upar( local_minUpar );

#ifdef CH_MPI
   MPI_Allreduce( &local_maxTemp, &max_Temp, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
   MPI_Allreduce( &local_maxTonB, &max_TonB, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
   MPI_Allreduce( &local_minTemp, &min_Temp, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD );
   MPI_Allreduce( &local_minTonB, &min_TonB, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD );
   MPI_Allreduce( &local_minUpar, &min_Upar, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD );
#endif

// cout << "max_Temp =" << max_Temp << endl;
// cout << "max_TonB =" << max_TonB << endl;

   // calculate min time step set by purely diffusive parts
   // dt*sum_i( D_i/dx_i^2 ) << 1/2;
   const Real alpha_vp = cls_freq*(max_Temp/mass)/(dvp*dvp);
   const Real alpha_mu = cls_freq*(2.0*mu_max*max_TonB)/(dmu*dmu);
   const Real dt_stable_alpha = 0.5/(alpha_vp+alpha_mu);

   // calculate min time step set by advective-diffusive
   // dt*sum_i( U_i^2/D_i ) << 2;
   const Real beta_vp = cls_freq*(vp_max-min_Upar)*(vp_max-min_Upar)/min_Temp/mass;
   const Real beta_mu = 2.0*cls_freq*mu_max/min_TonB;
   const Real dt_stable_beta = 2.0/(beta_vp+beta_mu);

   // take mininum of two time step restrictions
   Real dt_stable = Min(dt_stable_alpha, dt_stable_beta);

   // calculate Peclet number for mu and vp
   // (must be less than one for physically meaningful results)
   const Real Peclet_vp = dvp/sqrt(2.0*min_Temp/mass);
   const Real Peclet_mu = dmu/2.0/min_TonB;
   if ( procID() == 0 && !(Peclet_vp < 1) ){
     MayDay::Warning( "vpar grid spacing too large for physically meaninful results" );
   }
   if ( procID() == 0 && !(Peclet_mu < 1) ){
     MayDay::Warning( "mu grid spacing too large for physically meaninful results" );
   }

   std::string cls_only_str("false");
   if (cls_only){ cls_only_str = "true";}

   if (procID() == 0 ) {
     cout << "--------------------------------------------------------------------- " << endl;
     cout << "Conservative Drag Diffusion Collision Model being defined for species " << species_name <<endl;
     cout << "cls_freq  =  " << cls_freq << ", cls_only = " <<  cls_only_str << endl;
     cout << "Stable Diffusive time step           = " << dt_stable_alpha << endl;
     cout << "Stable Advective-Diffusive time step = " << dt_stable_beta << endl;
     cout << "Stable time step = " << dt_stable << endl;
     cout << "Peclet_vp = " << Peclet_vp << endl;
     cout << "Peclet_mu = " << Peclet_mu << endl;
     cout << "--------------------------------------------------------------------- " << endl;
   }
}

Real ConsDragDiff::computeDt( const KineticSpeciesPtrVect&  soln )
{
   return DBL_MAX;
}

#include "NamespaceFooter.H"
