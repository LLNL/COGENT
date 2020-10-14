#include <math.h>

#include "FourthOrderUtil.H"
#include "Directions.H"
#include "EdgeToCell.H"
#include "Linearized.H"
#include "CollisionsF_F.H"
#include "KineticFunctionLibrary.H"
#include "ConstFact.H"
#include "MomentOp.H"
#include "Kernels.H"


#include "NamespaceHeader.H" //Should be the last one



Linearized::Linearized( const std::string& a_ppcls_str, const int a_verbosity )
    : m_verbosity(a_verbosity),
      m_time_implicit(true),
      m_cls_freq(-1.0),
      m_conserve_momentum(true),
      m_conserve_energy(true),
      m_second_order(false),
      m_first_step(true),
      m_moment_op( MomentOp::instance() )
{
   ParmParse ppcls(a_ppcls_str.c_str());
   parseParameters( ppcls );
   if (m_verbosity>0) {
      printParameters();
   }
   
   m_fixed_cls_freq = (m_cls_freq<0.0) ? false : true ;
}

Linearized::~Linearized()
{
}

void Linearized::evalClsRHS( KineticSpeciesPtrVect&       a_rhs,
                             const KineticSpeciesPtrVect& a_soln,
                             const int                    a_species,
                             const int                    a_species_bkgr,
                             const Real                   a_time )
{
   CH_TIME("Linearized::evalClsRHS");
  
   // Get solution distribution function (J*Bstar_par*dfn) for the current species
   const KineticSpecies& soln_species( *(a_soln[a_species]) );
   const LevelData<FArrayBox>& soln_dfn( soln_species.distributionFunction() );
   
   // Get phase geometry
   const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();
   const DisjointBoxLayout& grids( soln_dfn.getBoxes() );

   if (m_first_step) {
      // Get reference (J*Bstar_par*dfn_init) distribution
      KineticSpeciesPtr ref_species( soln_species.clone( IntVect::Unit, false ) );
      m_ref_func->assign( *ref_species, a_time );
      const LevelData<FArrayBox>& ref_dfn( ref_species->distributionFunction() );
      
      if ( !m_F0.isDefined()) m_F0.define( grids, 1, IntVect::Zero );
      for (DataIterator dit(m_F0.dataIterator()); dit.ok(); ++dit) {
         m_F0[dit].copy(ref_dfn[dit]);
      }
      
      // Get reference temeprature
      const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
      CFG::LevelData<CFG::FArrayBox> ref_temperature( mag_geom.grids(), 1, CFG::IntVect::Zero );
      m_ref_temp->assign( ref_temperature, mag_geom, a_time);
      phase_geom.injectConfigurationToPhase( ref_temperature, m_temperature);

      //Create reference density distribution distribution and self-consistent cls_freq
      if (!m_fixed_cls_freq){
         CFG::LevelData<CFG::FArrayBox> ref_density( mag_geom.grids(), 1, CFG::IntVect::Zero );
         m_ref_dens->assign( ref_density, mag_geom, a_time);

         LevelData<FArrayBox> inj_density;
         phase_geom.injectConfigurationToPhase( ref_density, inj_density);

         m_sc_cls_freq.define(soln_dfn.disjointBoxLayout(), 1, IntVect::Zero);
         computeSelfConsistFreq(m_sc_cls_freq, inj_density, soln_species.mass(), soln_species.charge());
      }
   }

   // Compute the difference from the reference (or initial) solution
   const int n_comp( soln_dfn.nComp() );
   if ( !m_delta_dfn.isDefined()) m_delta_dfn.define( grids, n_comp, IntVect::Zero );

   for (DataIterator sdit(soln_dfn.dataIterator()); sdit.ok(); ++sdit) {
      m_delta_dfn[sdit].copy( soln_dfn[sdit] );
      m_delta_dfn[sdit].minus( m_F0[sdit] );
   }

   //Covert delta_dfn to cell centers
   if (!m_second_order) convertToCellCenters(phase_geom, m_delta_dfn);
   
   // Calculate test-particle (TP) collisional RHS
   const double mass = soln_species.mass();
   if ( !m_tp_rhs_coll.isDefined()) m_tp_rhs_coll.define( grids, n_comp, IntVect::Zero );
   testPartCollRHS(m_tp_rhs_coll, m_delta_dfn, phase_geom, mass);

   // Add conservative terms (field-particle terms)
   KineticSpecies& rhs_species( *(a_rhs[a_species]) );
   if (enforceConservation()) {
      addConservativeTerms( rhs_species,
                            m_tp_rhs_coll,
                            m_delta_dfn,
                            a_time );
   }

   // Add test-particle collisional RHS to Vlasov RHS (rhs_dfn has no ghost cells)
   LevelData<FArrayBox>& rhs_dfn( rhs_species.distributionFunction() );
   for (DataIterator rdit(rhs_dfn.dataIterator()); rdit.ok(); ++rdit) {
      if (m_fixed_cls_freq) {m_tp_rhs_coll[rdit].mult( m_cls_freq );}
      else {m_tp_rhs_coll[rdit].mult( m_sc_cls_freq[rdit] );}
      rhs_dfn[rdit].plus( m_tp_rhs_coll[rdit] );
   }  

   m_first_step = false;
}

void Linearized::addConservativeTerms(KineticSpecies& a_rhs_species,
                                      const LevelData<FArrayBox>& a_tp_rhs_coll,
                                      const LevelData<FArrayBox>& a_delta_F,
                                      const Real a_time )
{
   CH_TIME("Linearized::addConservativeTerms");
   
   const PhaseGeom& phase_geom = a_rhs_species.phaseSpaceGeometry();
   const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
   const DisjointBoxLayout& grids( a_tp_rhs_coll.getBoxes() );

   // Compute energy restoring factor
   CFG::LevelData<CFG::FArrayBox> zero_cfg( mag_geom.grids(), 1, CFG::IntVect::Zero);
   CFG::DataIterator cfg_dit = zero_cfg.dataIterator();
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
      zero_cfg[cfg_dit].setVal(0.);
   }
   CFG::LevelData<CFG::FArrayBox> energ_rest( mag_geom.grids(), 1, CFG::IntVect::Zero );
   m_moment_op.compute(energ_rest, a_rhs_species, a_tp_rhs_coll, PressureKernel(zero_cfg));
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
      energ_rest[cfg_dit].mult(3.0);
      energ_rest[cfg_dit].divide(a_rhs_species.mass());
   }

   // Compute momentum restoring factor
   CFG::LevelData<CFG::FArrayBox> moment_rest( mag_geom.grids(), 1, CFG::IntVect::Zero );
   m_moment_op.compute(moment_rest, a_rhs_species, a_tp_rhs_coll, ParallelMomKernel());

   //Calculate normalizing factors
   if (m_first_step) {

      LevelData<FArrayBox> kern_moment_norm( grids, 1, IntVect::Zero );
      LevelData<FArrayBox> kern_energ_norm( grids, 1, IntVect::Zero );
      evaluateNormKern( kern_moment_norm, kern_energ_norm, phase_geom, a_rhs_species.mass() );
      
      m_norm_momentum.define( mag_geom.grids(), 1, CFG::IntVect::Zero );
      m_norm_energy.define( mag_geom.grids(), 1, CFG::IntVect::Zero );
      
      m_moment_op.compute(m_norm_energy, a_rhs_species, kern_energ_norm, PressureKernel(zero_cfg));
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
         m_norm_energy[cfg_dit].mult(3.0);
         m_norm_energy[cfg_dit].divide(a_rhs_species.mass());
      }
      
      m_moment_op.compute(m_norm_momentum, a_rhs_species, kern_moment_norm, ParallelMomKernel());
   }

   // Calculate field-particle (FP) collisional RHS
   if (!m_fp_rhs_coll.isDefined()) m_fp_rhs_coll.define( grids, 1, IntVect::Zero );
   fieldPartCollRHS( m_fp_rhs_coll,
                     a_delta_F,
                     phase_geom,
                     energ_rest,
                     moment_rest,
                     a_rhs_species.mass() );

   // Add field-particle collisional RHS to Vlasov RHS
   LevelData<FArrayBox>& rhs_dfn( a_rhs_species.distributionFunction() );
   for (DataIterator rdit( rhs_dfn.dataIterator() ); rdit.ok(); ++rdit) {
      if (m_fixed_cls_freq) {m_fp_rhs_coll[rdit].mult( m_cls_freq );}
      else {m_fp_rhs_coll[rdit].mult( m_sc_cls_freq[rdit] );}
      rhs_dfn[rdit].plus( m_fp_rhs_coll[rdit] );
   }
}

void Linearized::testPartCollRHS(LevelData<FArrayBox>& a_rhs_coll,
                                 const LevelData<FArrayBox>& a_delta_dfn,
                                 const PhaseGeom& a_phase_geom,
                                 const double a_mass ) const
{

   /*
      Computes collisions with a Maxwellian background
    */

   CH_TIMERS("Linearized::testPartCollRHS");
   CH_TIMER("create_tmp_dfn",t_create_tmp_dfn);
   
   // Get coordinate system parameters and injected B
   const DisjointBoxLayout& grids( a_rhs_coll.getBoxes() );
   const LevelData<FArrayBox>& injected_B = a_phase_geom.getBFieldMagnitude();

   // Get velocity coordinate system parameters
   const VEL::VelCoordSys& vel_coords = a_phase_geom.velSpaceCoordSys();
   const VEL::ProblemDomain& vel_domain = vel_coords.domain();
   const VEL::Box& domain_box = vel_domain.domainBox();
   const VEL::RealVect& vel_dx = vel_coords.dx();
   const int num_vpar_cells = domain_box.size(0);
   const int num_mu_cells = domain_box.size(1);

   //Currently we only have 4th order implementation,
   //which require at least two layers of ghosts
   //see if we can optimize by growing only in v-space
   CH_START(t_create_tmp_dfn);
   IntVect ghostVect = IntVect::Zero;
   for (int dir=CFG_DIM; dir<SpaceDim; dir++) {
     ghostVect[dir] = 3;
   }
   if ( !m_tmp_dfn_withGhosts.isDefined()) m_tmp_dfn_withGhosts.define( grids, 1, ghostVect );

   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      m_tmp_dfn_withGhosts[dit].setVal(0.);
      m_tmp_dfn_withGhosts[dit].copy(a_delta_dfn[dit]);
   }
   m_tmp_dfn_withGhosts.exchange();
   CH_STOP(t_create_tmp_dfn);

   //Compute cell-centered collision fluxes (0 comp - vpar_dir, 1 comp - mu_dir)
   ghostVect = IntVect::Zero;
   for (int dir=CFG_DIM; dir<SpaceDim; dir++) {
     // Because we use different approach to calcuation of face-centered fluxes
     // for second-order, we need to have extra ghost layer here
     if (m_second_order) ghostVect[dir] = 1;
   }
   if ( !m_flux_vmu_cc.isDefined()) m_flux_vmu_cc.define( grids, 2, ghostVect );
   
   for (DataIterator dit( a_rhs_coll.dataIterator() ); dit.ok(); ++dit) {
      FArrayBox& this_flux_cc = m_flux_vmu_cc[dit];
      const FArrayBox& this_delta_dfn_tmp = m_tmp_dfn_withGhosts[dit];
      const FArrayBox& this_b = injected_B[dit];
      
      //Compute Lorentz pitch-angle scattering term
      const FArrayBox& this_temperature = m_temperature[dit];
      FORT_EVALUATE_TP_LORENTZ(CHF_FRA(this_flux_cc),
                               CHF_CONST_FRA1(this_delta_dfn_tmp,0),
                               CHF_CONST_FRA1(this_b,0),
                               CHF_CONST_FRA1(this_temperature,0),
                               CHF_BOX(this_flux_cc.box()),
                               CHF_CONST_REALVECT(vel_dx),
                               CHF_CONST_REAL(a_mass));

      //Add energy-diffusion part of the TP collisions
      FORT_EVALUATE_TP_ENERG_DIFF(CHF_FRA(this_flux_cc),
                               CHF_CONST_FRA1(this_delta_dfn_tmp,0),
                               CHF_CONST_FRA1(this_b,0),
                               CHF_CONST_FRA1(this_temperature,0),
                               CHF_BOX(this_flux_cc.box()),
                               CHF_CONST_REALVECT(vel_dx),
                               CHF_CONST_REAL(a_mass));
   }

   //Calculate collision flux cell-average
   if (!m_second_order) convertToCellAverage(a_phase_geom, m_flux_vmu_cc);

   //Calculate collision fluxes on cell faces
   if ( !m_flux_vmu.isDefined()) m_flux_vmu.define( grids, 2, IntVect::Zero );
   convertToCellFaces(m_flux_vmu, m_flux_vmu_cc);

   //Create final (combined r , theta, mu, and v_par) rhs flux
   if ( !m_flux.isDefined()) m_flux.define( grids, SpaceDim, IntVect::Zero );
   for (DataIterator dit(m_flux.dataIterator() ); dit.ok(); ++dit) {
      const FArrayBox& this_flux_cc = m_flux_vmu_cc[dit];

      FArrayBox tmp_flux_cell(grow(this_flux_cc.box(),1),2);
      tmp_flux_cell.setVal(0.);
      tmp_flux_cell.copy(this_flux_cc);

      for (int dir=0; dir<SpaceDim; dir++) {
          FArrayBox& this_flux_rhs = m_flux[dit][dir];
          FArrayBox& this_flux_face = m_flux_vmu[dit][dir];
          FORT_EVALUATE_COLL_FLUX_COMBINED(CHF_FRA(this_flux_rhs),
                                           CHF_CONST_INT(dir),
                                           CHF_BOX(this_flux_rhs.box()),
                                           CHF_CONST_FRA(this_flux_face),
                                           CHF_CONST_FRA(tmp_flux_cell),
                                           CHF_CONST_INT(num_vpar_cells),
                                           CHF_CONST_INT(num_mu_cells));
      }
   }

   a_phase_geom.mappedGridDivergence(a_rhs_coll, m_flux, true);
   for (DataIterator dit( a_rhs_coll.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
      double fac = 1. / block_coord_sys.getMappedCellVolume();
      a_rhs_coll[dit].mult(fac);
    }
}

void Linearized::fieldPartCollRHS(LevelData<FArrayBox>& a_rhs_coll,
                                  const LevelData<FArrayBox>& a_delta_F,
                                  const PhaseGeom& a_phase_geom,
                                  const CFG::LevelData<CFG::FArrayBox>& a_rest_energy,
                                  const CFG::LevelData<CFG::FArrayBox>& a_rest_momentum,
                                  const double a_mass ) const
{
   CH_TIME("Linearized::fieldPartCollRHS");
   
   const DisjointBoxLayout& grids( a_rhs_coll.getBoxes() );
   //const int n_comp( a_rhs_coll.nComp() );
   const LevelData<FArrayBox>& B_injected( a_phase_geom.getBFieldMagnitude() );
   const VEL::VelCoordSys& vel_coords = a_phase_geom.velSpaceCoordSys();
   const VEL::ProblemDomain& vel_domain = vel_coords.domain();
   const VEL::Box& domain_box = vel_domain.domainBox();
   const VEL::RealVect& vel_dx = vel_coords.dx();
   int num_vpar_cells = domain_box.size(0);
   int num_mu_cells = domain_box.size(1);


   if (m_conserve_momentum) {

     LevelData<FArrayBox> inj_MRest;
     a_phase_geom.injectConfigurationToPhase( a_rest_momentum, inj_MRest );

     LevelData<FArrayBox> inj_MNorm;
     a_phase_geom.injectConfigurationToPhase( m_norm_momentum, inj_MNorm );

     for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
        //const PhaseBlockCoordSys& block_coord_sys( a_phase_geom.getBlockCoordSys(grids[dit]) );

        FArrayBox& this_RHS( a_rhs_coll[dit] );
        const FArrayBox& this_MRest( inj_MRest[dit] );
        const FArrayBox& this_MNorm( inj_MNorm[dit] );
      
        const FArrayBox& this_B( B_injected[dit] );
        const FArrayBox& this_Temp_Distr( m_temperature[dit] );
        
        FORT_EVALUATE_FULL_MR(CHF_BOX(this_RHS.box()),
                              CHF_FRA1(this_RHS,0),
                              CHF_CONST_FRA1(this_MRest,0),
                              CHF_CONST_FRA1(this_MNorm,0),
                              CHF_CONST_FRA1(this_Temp_Distr,0),
                              CHF_CONST_FRA1(this_B,0),
                              CHF_CONST_REALVECT(vel_dx),
                              CHF_CONST_REAL(a_mass));
     }
   }

   //Add the energy-conservative term in the divergent form
   if (m_conserve_energy) {

     LevelData<FArrayBox> inj_ERest;
     a_phase_geom.injectConfigurationToPhase( a_rest_energy, inj_ERest );

     LevelData<FArrayBox> inj_ENorm;
     a_phase_geom.injectConfigurationToPhase( m_norm_energy, inj_ENorm );

     if (!m_flux_full_ERest.isDefined()) m_flux_full_ERest.define(grids, SpaceDim, IntVect::Zero);
     if (!m_rhs_ERest.isDefined()) m_rhs_ERest.define(grids, 1, IntVect::Zero);

     for (DataIterator dit(m_flux_full_ERest.dataIterator() ); dit.ok(); ++dit) {
        FArrayBox& this_ERest = inj_ERest[dit];
        FArrayBox& this_ENorm = inj_ENorm[dit];
        const FArrayBox& this_TempDistr = m_temperature[dit];
        const FArrayBox& this_B = B_injected[dit];
        
        m_flux_full_ERest[dit].setVal(0.);
        for (int dir=CFG_DIM; dir<SpaceDim; dir++) {
           FArrayBox& thisFlux_FullERest = m_flux_full_ERest[dit][dir];
           FORT_EVALUATE_FULL_ER_FLUX(CHF_CONST_INT(dir),
                                      CHF_BOX(thisFlux_FullERest.box()),
                                      CHF_FRA(thisFlux_FullERest),
                                      CHF_CONST_FRA1(this_ERest,0),
                                      CHF_CONST_FRA1(this_ENorm,0),
                                      CHF_CONST_FRA1(this_TempDistr,0),
                                      CHF_CONST_FRA1(this_B,0),
                                      CHF_CONST_REALVECT(vel_dx),
                                      CHF_CONST_INT(num_vpar_cells),
                                      CHF_CONST_INT(num_mu_cells),
                                      CHF_CONST_REAL(a_mass));
        }
     }

     // Calculate energy restoring term
     a_phase_geom.mappedGridDivergence(m_rhs_ERest, m_flux_full_ERest, true);
     for (DataIterator dit(m_rhs_ERest.dataIterator() ); dit.ok(); ++dit) {
        const PhaseBlockCoordSys& block_coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
        double fac = 1. / block_coord_sys.getMappedCellVolume();
        m_rhs_ERest[dit].mult(fac);
        if (m_conserve_momentum) { 
            a_rhs_coll[dit].plus(m_rhs_ERest[dit]);
        }
        else { 
            a_rhs_coll[dit].copy(m_rhs_ERest[dit]);
        }
     }
   }
}

void Linearized::evaluateNormKern( LevelData<FArrayBox>& a_kern_moment_norm,
                                   LevelData<FArrayBox>& a_kern_energ_norm,
                                   const PhaseGeom& a_phase_geom,
                                   const double a_mass ) const
{
   CH_TIME("Linearized::evaluateNormKern");
   
   const LevelData<FArrayBox>& B_injected = a_phase_geom.getBFieldMagnitude();
   const DisjointBoxLayout& grids ( a_kern_moment_norm.disjointBoxLayout() );

   const VEL::VelCoordSys& vel_coords = a_phase_geom.velSpaceCoordSys();
   const VEL::ProblemDomain& vel_domain = vel_coords.domain();
   const VEL::Box& domain_box = vel_domain.domainBox();
   const VEL::RealVect& vel_dx = vel_coords.dx();
   int num_vpar_cells = domain_box.size(0);
   int num_mu_cells = domain_box.size(1);

   //Calculate "nu_E" using its divergence representation 
   //and storing the result in a_kern_energ_norm
   LevelData<FluxBox> flux_norm_ERest(grids, SpaceDim, IntVect::Zero);
   for (DataIterator dit(flux_norm_ERest.dataIterator() ); dit.ok(); ++dit) {
      const FArrayBox& this_TempDistr = m_temperature[dit];
      const FArrayBox& this_B = B_injected[dit];

      flux_norm_ERest[dit].setVal(0.);
      for (int dir=0; dir<SpaceDim; dir++) {
         FArrayBox& this_flux_NormERest = flux_norm_ERest[dit][dir];

         FORT_EVALUATE_NORM_ER_FLUX(CHF_CONST_INT(dir),
                                    CHF_BOX(this_flux_NormERest.box()),
                                    CHF_FRA(this_flux_NormERest),
                                    CHF_CONST_FRA1(this_TempDistr,0),
                                    CHF_CONST_FRA1(this_B,0),
                                    CHF_CONST_REALVECT(vel_dx),
                                    CHF_CONST_INT(num_vpar_cells),
                                    CHF_CONST_INT(num_mu_cells),
                                    CHF_CONST_REAL(a_mass));
      }
   }

   a_phase_geom.mappedGridDivergence(a_kern_energ_norm, flux_norm_ERest, true);
   for (DataIterator dit(a_kern_energ_norm.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
      double fac = 1. / block_coord_sys.getMappedCellVolume();
      a_kern_energ_norm[dit].mult(fac);
   }

   //Calculate kernels for the conserving terms normalization factors
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      //const PhaseBlockCoordSys& block_coord_sys( a_phase_geom.getBlockCoordSys(grids[dit]) );

      FArrayBox& this_MN( a_kern_moment_norm[dit] );
      const FArrayBox& this_TempDistr = m_temperature[dit];
      const FArrayBox& this_B = B_injected[dit];
      FORT_EVALUATE_NORM_MR(CHF_BOX(this_MN.box()),
                            CHF_FRA1(this_MN,0),
                            CHF_CONST_FRA1(this_TempDistr,0),
                            CHF_CONST_FRA1(this_B,0),
                            CHF_CONST_REALVECT(vel_dx),
                            CHF_CONST_REAL(a_mass));
   }
}

double Linearized::nu_D(double x) const
{
   double pi=3.14159265358979;
   double G=(erf(x)-x*2.0/sqrt(pi)*exp(-x*x))/(2.0*x*x);
   double F=erf(x) ;

   return 1/(x*x*x)*(F-G);

}

double Linearized::nu_s(double x) const
{
   double pi=3.14159265358979;
   double G=(erf(x)-x*2.0/sqrt(pi)*exp(-x*x))/(2.0*x*x);
   return 4.0*G/(x);
}


void Linearized::computeSelfConsistFreq(LevelData<FArrayBox>& a_cls_freq,
                                  const LevelData<FArrayBox>& a_density,
                                  const double                a_mass,
                                  const double                a_charge      ) const

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
    double coeff = 1/( 4.0 * pi * eps0 * eps0 * pow(2, 3.0/2.0) ); 

    //Compute normalized collision frequency
    double Coulomb_Lg = 23 - log( sqrt(2.0) * pow(a_charge,3) * sqrt(N)/1000.0 / pow(T, 3.0/2.0) ); 
    double cls_norm = coeff * N * pow(ech, 2) * pow(a_charge, 4) * L 
                     / ( sqrt(a_mass) *  pow(T, 2)) * Coulomb_Lg;

    const DisjointBoxLayout& grids = a_cls_freq.disjointBoxLayout();
    DataIterator dit(grids.dataIterator());
    for (dit.begin(); dit.ok(); ++dit) { 

      FArrayBox& this_cls_freq = a_cls_freq[dit];
      const FArrayBox& this_n = a_density[dit];
      const FArrayBox& this_T = m_temperature[dit];

      FORT_COMPUTE_SC_CLS_FREQ(CHF_BOX(this_cls_freq.box()),
                               CHF_FRA1(this_cls_freq,0),
                               CHF_CONST_FRA1(this_n,0),
                               CHF_CONST_FRA1(this_T,0));
  
    this_cls_freq.mult(cls_norm); 

    }

}

void Linearized::convertToCellAverage(const PhaseGeom&  a_phase_geom,
                                      LevelData<FArrayBox>&      a_dfn) const
{
   // Convert from cell centered values to cell averages.  We need to use
   // a temporary with an extra ghost cell for the interpolation

   LevelData<FArrayBox> dfn_tmp(a_dfn.disjointBoxLayout(),
                                a_dfn.nComp(),
                                a_dfn.ghostVect()+IntVect::Unit);

   const DisjointBoxLayout& grids = a_dfn.disjointBoxLayout();

   DataIterator dit(grids.dataIterator());
   for (dit.begin(); dit.ok(); ++dit) {
      dfn_tmp[dit].copy(a_dfn[dit]);
   }
   dfn_tmp.exchange();

   for (dit.begin(); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
      const ProblemDomain& domain = coord_sys.domain();

      fourthOrderAverageCell(dfn_tmp[dit], domain, grids[dit]);
   }
   dfn_tmp.exchange();

   for (dit.begin(); dit.ok(); ++dit) {
      a_dfn[dit].copy(dfn_tmp[dit]);
   }
}

void Linearized::convertToCellCenters(const PhaseGeom&  a_phase_geom,
                                      LevelData<FArrayBox>&      a_dfn) const
{
   // Convert from cell centered values to cell averages.  We need to use
   // a temporary with an extra ghost cell for the interpolation

   LevelData<FArrayBox> dfn_tmp(a_dfn.disjointBoxLayout(),
                                a_dfn.nComp(),
                                a_dfn.ghostVect()+IntVect::Unit);

   const DisjointBoxLayout& grids = a_dfn.disjointBoxLayout();

   DataIterator dit(grids.dataIterator());
   for (dit.begin(); dit.ok(); ++dit) {
      dfn_tmp[dit].copy(a_dfn[dit]);
   }
   dfn_tmp.exchange();

   for (dit.begin(); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
      const ProblemDomain& domain = coord_sys.domain();

      fourthOrderAverageCell(dfn_tmp[dit], domain, grids[dit], -1);
   }
   dfn_tmp.exchange();

   for (dit.begin(); dit.ok(); ++dit) {
      a_dfn[dit].copy(dfn_tmp[dit]);
   }
}

void Linearized::convertToCellFaces(LevelData<FluxBox>& a_faceData,
                                    const LevelData<FArrayBox>& a_cellData) const
{
    // Convert from cell centered values to cell averages.  We need to use
    // a temporary with two extra ghost cells for the interpolation. The ghost cells at the domain
    // boundaries are filled with zeros for our purposes. Since, the boundary fluxes in the velocity directions
    // will be overwrited later, and the fluxes in the configuration domansions do not contribute to the RHS.

    if (!m_second_order) {
      LevelData<FArrayBox> dfn_tmp(a_cellData.disjointBoxLayout(),
                                   a_cellData.nComp(),
                                   a_cellData.ghostVect()+2*IntVect::Unit);
      
      const DisjointBoxLayout& grids = a_cellData.disjointBoxLayout();

      DataIterator dit(grids.dataIterator());
      for (dit.begin(); dit.ok(); ++dit) {
         dfn_tmp[dit].setVal(0.0);
         dfn_tmp[dit].copy(a_cellData[dit]);
      }
      dfn_tmp.exchange();
      
      fourthOrderCellToFace(a_faceData, dfn_tmp);
    }
   
    else {
       const DisjointBoxLayout& grids = a_cellData.disjointBoxLayout();
       DataIterator dit(grids.dataIterator());
       for (dit.begin(); dit.ok(); ++dit) {
          a_faceData[dit].setVal(0.);
          for (int dir=CFG_DIM; dir<SpaceDim; dir++) {
             FORT_SECOND_ORDER_CELL_TO_FACE(CHF_CONST_FRA(a_cellData[dit]),
                                            CHF_FRA(a_faceData[dit][dir]),
                                            CHF_BOX(a_faceData[dit][dir].box()),
                                            CHF_CONST_INT(dir));
          }
       }
    }
}

inline
void Linearized::parseParameters( ParmParse& a_ppcls )
{
   a_ppcls.query("time_implicit", m_time_implicit);
   a_ppcls.query("second_order", m_second_order);
   a_ppcls.query("cls_freq",m_cls_freq);
   a_ppcls.query("conserve_momentum",m_conserve_momentum);
   a_ppcls.query("conserve_energy",m_conserve_energy);

   KineticFunctionLibrary* library = KineticFunctionLibrary::getInstance();
   std::string function_name;
   a_ppcls.query( "ref_function", function_name );
   m_ref_func = library->find( function_name );

   CFG::GridFunctionLibrary* grid_library = CFG::GridFunctionLibrary::getInstance();
   std::string grid_function_name;

   a_ppcls.query( "ref_temperature", grid_function_name );
   m_ref_temp = grid_library->find( grid_function_name );

   a_ppcls.query( "ref_density", grid_function_name );
   m_ref_dens = grid_library->find( grid_function_name );
}

inline
void Linearized::printParameters()
{
   if (procID()==0) {
      std::cout << "Linearized collisions parameters:" << std::endl;
      std::cout << "  cls_freq  =  " << m_cls_freq
                << ", conserve_momentum = " << m_conserve_momentum 
                << ", conserve_energy = " << m_conserve_energy 
                << ", implicit in time = " << m_time_implicit
                << ", Second order = " << m_second_order
                << std::endl;
      std::cout << "  Reference Function:" << std::endl;
      m_ref_func->printParameters();
      std::cout << "  Reference Temperature:" << std::endl;
      m_ref_temp->printParameters();
      std::cout << "  Reference Density:" << std::endl;
      m_ref_dens->printParameters();

   }
}



Real Linearized::computeTimeScale(const KineticSpeciesPtrVect& soln, const int a_idx)
{
   // Replace with an estimate of time scale
   return DBL_MAX;
}

Real Linearized::computeDtExplicitTI(const KineticSpeciesPtrVect& soln, const int a_idx)
{
   // Replace with an estimate of dt
   return DBL_MAX;
}

Real Linearized::computeDtImExTI(const KineticSpeciesPtrVect& soln, const int a_idx)
{
  if (m_time_implicit) {
    return DBL_MAX;
  } else {
    return computeDtExplicitTI(soln, a_idx);
  }
}

#include "NamespaceFooter.H"
