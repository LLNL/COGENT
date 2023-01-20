#include <math.h>

#include "FourthOrderUtil.H"
#include "Directions.H"
#include "EdgeToCell.H"
#include "LinearizedUnlike.H"
#include "LinearizedUnlikeF_F.H"
#include "CollisionsF_F.H"
#include "KineticFunctionLibrary.H"
#include "ConstFact.H"

#include "MomentOp.H"
#include "Kernels.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FluxSurface.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H" //Should be the last one



LinearizedUnlike::LinearizedUnlike( const std::string& a_ppcls_str, const int a_verbosity )
    : m_verbosity(a_verbosity),
      m_time_implicit(false),
      m_update_freq(-1),
      m_it_counter(0),
      m_second_order(true),
      m_conserve_momentum(true),
      m_conserve_energy(true),
      m_fixed_bkgr(false),
      m_first_call(true),
      m_first_stage(true),
      m_include_bkgr_cls(true),
      m_moment_op( MomentOp::instance() )
{
   ParmParse ppcls(a_ppcls_str.c_str());
   parseParameters( ppcls );
   if (m_verbosity>0) {
      printParameters();
   }
}

LinearizedUnlike::~LinearizedUnlike()
{
}

void LinearizedUnlike::evalClsRHS(KineticSpeciesPtrVect&       a_rhs,
                                  const KineticSpeciesPtrVect& a_soln,
                                  const int                    a_species,
                                  const int                    a_species_bkgr,
                                  const Real                   a_time )
{
   CH_TIME("LinearizedUnlike::evalClsRHS");
   
   // Get solution distribution function (J*Bstar_par*dfn) for the colliding species (a)
   const KineticSpecies& soln_species_a( *(a_soln[a_species]) );
   const LevelData<FArrayBox>& soln_dfn_a( soln_species_a.distributionFunction() );

   // Get solution distribution function (J*Bstar_par*dfn) for the background species (b)
   const KineticSpecies& soln_species_b( *(a_soln[a_species_bkgr]) );
   const LevelData<FArrayBox>& soln_dfn_b( soln_species_b.distributionFunction() );

   // Get geometry parameters
   const PhaseGeom& phase_geom = soln_species_a.phaseSpaceGeometry();
   const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
   const DisjointBoxLayout& grids( soln_dfn_a.getBoxes() );
   const CFG::DisjointBoxLayout& grids_cfg = mag_geom.grids();
   
   //Get species parameters
   m_mass_a = soln_species_a.mass();
   m_mass_b = soln_species_b.mass();
   m_charge_a = soln_species_a.charge();
   m_charge_b = soln_species_b.charge();

   // Get reference data
   if ( !m_F0a.isDefined()) m_F0a.define( grids, 1, IntVect::Zero );
   if ( !m_F0b.isDefined()) m_F0b.define( grids, 1, IntVect::Zero );

   if ( !m_Ta.isDefined()) m_Ta.define( grids_cfg, 1, CFG::IntVect::Zero );
   if ( !m_Tb.isDefined()) m_Tb.define( grids_cfg, 1, CFG::IntVect::Zero );

   if ( !m_Na.isDefined()) m_Na.define( grids_cfg, 1, CFG::IntVect::Zero );
   if ( !m_Nb.isDefined()) m_Nb.define( grids_cfg, 1, CFG::IntVect::Zero );
   
   if ( !m_Rab.isDefined()) m_Rab.define( grids, 1, IntVect::Zero );
   if ( !m_Qab.isDefined()) m_Qab.define( grids, 1, IntVect::Zero );
   
   if (m_fixed_bkgr && m_first_call) {
      // Get prescribed background data
      getReferenceData(m_F0a, m_Ta, m_Na, m_Ta_inj, m_Na_inj, soln_species_a, m_ref_func_a, m_ref_Ta, m_ref_Na, a_time);
      getReferenceData(m_F0b, m_Tb, m_Nb, m_Tb_inj, m_Nb_inj, soln_species_b, m_ref_func_b, m_ref_Tb, m_ref_Nb, a_time);
      // Compute cls frequency normalization and FP terms
      getClsFreqNorm(m_cls_norm, m_Ta, m_Tb, m_Na, m_Nb, phase_geom);
      computeRandQTerms(m_Rab, m_Qab, soln_species_a);
   }
   
   if ( !m_fixed_bkgr && ((m_update_freq < 0) || (m_it_counter % m_update_freq == 0 && m_first_stage) || m_first_call)) {
      // Get background data by fitting Maxwellian
      computeMaxwellianFit(m_F0a, m_Ta, m_Na,  m_Ta_inj, m_Na_inj, soln_species_a, soln_dfn_a, a_time );
      computeMaxwellianFit(m_F0b, m_Tb, m_Nb,  m_Tb_inj, m_Nb_inj, soln_species_b, soln_dfn_b, a_time );
      
      // Compute cls frequency normalization and FP terms
      getClsFreqNorm(m_cls_norm, m_Ta, m_Tb, m_Na, m_Nb, phase_geom);
      computeRandQTerms(m_Rab, m_Qab, soln_species_a);
   }
      
   // Get delta-f distributions; presently some calculations are hard-wired 4th order
   // therefore we need at least two layers of ghosts; calculations however
   // do not use any information from physical ghosts
   // (currently delta_dfn at large vpar and mu is set to 0, perhaps change to extrapolation later)
   if ( !m_delta_dfn_a.isDefined()) m_delta_dfn_a.define( grids, 1, IntVect::Zero );
   if ( !m_delta_dfn_b.isDefined()) m_delta_dfn_b.define( grids, 1, IntVect::Zero );
   
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      m_delta_dfn_a[dit].copy(soln_dfn_a[dit]);
      m_delta_dfn_a[dit].minus(m_F0a[dit]);
      
      m_delta_dfn_b[dit].copy(soln_dfn_b[dit]);
      m_delta_dfn_b[dit].minus(m_F0b[dit]);
   }

   if (!m_second_order) {
     convertToCellCenters(phase_geom, m_delta_dfn_a);
     convertToCellCenters(phase_geom, m_delta_dfn_b);
   }
   
   //Compute RHS for the colliding speices (a)
   KineticSpecies& rhs_species( *(a_rhs[a_species]) );
   LevelData<FArrayBox>& rhs_dfn( rhs_species.distributionFunction() );

   if ( !m_rhs_tmp.isDefined()) m_rhs_tmp.define( grids, 1, IntVect::Zero );

   // Compute CTP_ab[Fa0b,F0b] and add to RHS
   if (m_include_bkgr_cls) {
     testPartCollRHS(m_rhs_tmp, m_F0a, m_Nb_inj, m_Tb_inj, soln_species_a, m_mass_a, m_mass_b);
     for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
       rhs_dfn[dit].plus(m_rhs_tmp[dit]);
     }
   }
   
   // Compute CTP_ab[delta_fa,F0b] and add to RHS
   testPartCollRHS(m_rhs_tmp, m_delta_dfn_a, m_Nb_inj, m_Tb_inj, soln_species_a, m_mass_a, m_mass_b);
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
     rhs_dfn[dit].plus(m_rhs_tmp[dit]);
   }

   // Compute CFP_ab and add to RHS
   fieldPartCollRHS(m_rhs_tmp, m_delta_dfn_b, soln_species_b);
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
     rhs_dfn[dit].plus(m_rhs_tmp[dit]);
   }

   m_first_call = false;
   m_first_stage = false;
}

void LinearizedUnlike::testPartCollRHS(LevelData<FArrayBox>& a_rhs_coll,
                                       const LevelData<FArrayBox>& a_delta_dfn_a,
                                       const LevelData<FArrayBox>& a_Nb_inj,
                                       const LevelData<FArrayBox>& a_Tb_inj,
                                       const KineticSpecies&   a_species,
                                       const double a_ma,
                                       const double a_mb) const
{
   /*
      Computes collisions of dfn (species a) with a Maxwellian background (species b)
    */

   CH_TIMERS("LinearizedUnlike::testPartCollRHS");
   CH_TIMER("create_tmp_dfn",t_create_tmp_dfn); 
   
   // Get coordinate system parameters and injected B
   const PhaseGeom& phase_geom = a_species.phaseSpaceGeometry();
   const DisjointBoxLayout& grids( a_rhs_coll.getBoxes() );
   const LevelData<FArrayBox>& injected_B = phase_geom.getBFieldMagnitude();

   // Get velocity coordinate system parameters
   const VEL::VelCoordSys& vel_coords = phase_geom.velSpaceCoordSys();
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
      m_tmp_dfn_withGhosts[dit].copy(a_delta_dfn_a[dit]);
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
      const FArrayBox& this_delta_dfn = m_tmp_dfn_withGhosts[dit];
      const FArrayBox& this_b = injected_B[dit];
      
      //Compute Lorentz pitch-angle scattering term
      const FArrayBox& this_temperature = a_Tb_inj[dit];
      FORT_EVALUATE_TP_LORENTZ_UNLIKE(CHF_FRA(this_flux_cc),
                                      CHF_CONST_FRA1(this_delta_dfn,0),
                                      CHF_CONST_FRA1(this_b,0),
                                      CHF_CONST_FRA1(this_temperature,0),
                                      CHF_BOX(this_flux_cc.box()),
                                      CHF_CONST_REALVECT(vel_dx),
                                      CHF_CONST_REAL(a_ma),
                                      CHF_CONST_REAL(a_mb));

      //Add energy-diffusion part of the TP collisions
      FORT_EVALUATE_TP_ENERG_DIFF_UNLIKE(CHF_FRA(this_flux_cc),
                                         CHF_CONST_FRA1(this_delta_dfn,0),
                                         CHF_CONST_FRA1(this_b,0),
                                         CHF_CONST_FRA1(this_temperature,0),
                                         CHF_BOX(this_flux_cc.box()),
                                         CHF_CONST_REALVECT(vel_dx),
                                         CHF_CONST_REAL(a_ma),
                                         CHF_CONST_REAL(a_mb));
   }

   //Calculate collision flux cell-average
   if (!m_second_order) convertToCellAverage(phase_geom, m_flux_vmu_cc);

   //Calculate collision fluxes on cell faces
   if ( !m_flux_vmu.isDefined()) m_flux_vmu.define( grids, 2, IntVect::Zero );
   convertToCellFaces(m_flux_vmu, m_flux_vmu_cc);

   //Create final (combined r , theta, mu, and v_par) rhs flux
   if ( !m_flux.isDefined()) m_flux.define( grids, SpaceDim, IntVect::Zero );
   for (DataIterator dit(m_flux.dataIterator() ); dit.ok(); ++dit) {
      const FArrayBox& this_flux_cell = m_flux_vmu_cc[dit];

      FArrayBox tmp_flux_cell(grow(this_flux_cell.box(),1),2);
      tmp_flux_cell.setVal(0.);
      tmp_flux_cell.copy(this_flux_cell);

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
   
   // The following works for the identity-mapping the velocity space
   phase_geom.mappedGridDivergence(a_rhs_coll, m_flux, true);
   
   for (DataIterator dit( a_rhs_coll.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(grids[dit]);
      double fac = 1. / block_coord_sys.getMappedCellVolume();
      a_rhs_coll[dit].mult(fac);
      
      FORT_CLS_FREQ_SCALE(CHF_BOX(a_rhs_coll[dit].box()),
                          CHF_FRA1(a_rhs_coll[dit],0),
                          CHF_CONST_FRA1(a_Nb_inj[dit],0),
                          CHF_CONST_FRA1(m_cls_norm[dit],0),
                          CHF_CONST_REAL(a_ma));
    }
}


void LinearizedUnlike::fieldPartCollRHS(LevelData<FArrayBox>& a_rhs_coll,
                                        const LevelData<FArrayBox>& a_delta_dfn_b,
                                        const KineticSpecies&       a_species_b) const
{
   /*
    This function computes model field-particle operator C_FP_ab
    given in Eqs. (24) of Kolesnikov et al, 2010.
    */

   CH_TIME("LinearizedUnlike::fieldPartCollRHS");
   
   //Get geometry parameters
   const PhaseGeom& phase_geom = a_species_b.phaseSpaceGeometry();
   const CFG::MultiBlockLevelGeom& mag_geom( phase_geom.magGeom() );
   const DisjointBoxLayout& grids( m_F0a.getBoxes() );

   if ( !m_tmp_dfn_noGhosts.isDefined()) m_tmp_dfn_noGhosts.define( grids, 1, IntVect::Zero );
   
   /*
     Add momentum conserving term Rab * vpar * delta_P (see Eqs. (24)-(25))
     */
   
   //Compute CTP_ba(delta_f_b) and place result in m_tmp_dfn_noGhosts
   testPartCollRHS(m_tmp_dfn_noGhosts, m_delta_dfn_b, m_Na_inj, m_Ta_inj, a_species_b, m_mass_b, m_mass_a);
   
   //Compute int(vpar*CTP_ba(delta_f_b),d3v)
   CFG::LevelData<CFG::FArrayBox> tmp_cfg( mag_geom.grids(), 1, CFG::IntVect::Zero);
   m_moment_op.compute(tmp_cfg, a_species_b, m_tmp_dfn_noGhosts, ParallelVelKernel<FArrayBox>());
   
   LevelData<FArrayBox> injected_data_R;
   phase_geom.injectConfigurationToPhase(tmp_cfg, injected_data_R);
   
   for (DataIterator dit( a_rhs_coll.dataIterator() ); dit.ok(); ++dit) {
      
      if (m_conserve_momentum) {
         a_rhs_coll[dit].copy(m_Rab[dit]);
      }
      else {
         a_rhs_coll[dit].setVal(0.);
      }
      a_rhs_coll[dit].mult(m_mass_b);
      a_rhs_coll[dit].negate();
      
      FORT_MULT_CFG_DATA(CHF_BOX(a_rhs_coll[dit].box()),
                         CHF_FRA1(a_rhs_coll[dit],0),
                         CHF_CONST_FRA1(injected_data_R[dit],0));
   }
   
   
   /*
     Add energy conserving term Qab * delta_EP (see Eqs. (24), (26))
    */

   //Compute int(mb*v^2/3*CTP_ba(delta_f_b),d3v)
   CFG::LevelData<CFG::FArrayBox> zero_cfg( mag_geom.grids(), 1, CFG::IntVect::Zero);
   CFG::DataIterator cfg_dit = zero_cfg.dataIterator();
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
      zero_cfg[cfg_dit].setVal(0.);
   }
   m_moment_op.compute(tmp_cfg, a_species_b, m_tmp_dfn_noGhosts, PressureKernel<FArrayBox>(zero_cfg));
   
   LevelData<FArrayBox> injected_data_Q;
   phase_geom.injectConfigurationToPhase(tmp_cfg, injected_data_Q);
   
   for (DataIterator dit( a_rhs_coll.dataIterator() ); dit.ok(); ++dit) {
      
      m_tmp_dfn_noGhosts[dit].copy(m_Qab[dit]);
      m_tmp_dfn_noGhosts[dit].mult(3.0/2.0);
      m_tmp_dfn_noGhosts[dit].negate();
      
      FORT_MULT_CFG_DATA(CHF_BOX(m_tmp_dfn_noGhosts[dit].box()),
                         CHF_FRA1(m_tmp_dfn_noGhosts[dit],0),
                         CHF_CONST_FRA1(injected_data_Q[dit],0));
      
      if (m_conserve_energy) {
         a_rhs_coll[dit].plus(m_tmp_dfn_noGhosts[dit]);
      }
   }
   
}

void LinearizedUnlike::computeRandQTerms(LevelData<FArrayBox>&       a_R,
                                         LevelData<FArrayBox>&       a_Q,
                                         const KineticSpecies&       a_species) const
{
 
   /*
    This function computes RabV_par (=a_R) and Qab (=a_Q) terms given in
    Eqs. (30)-(31) of Kolesnikov et al, 2010.
    */

   CH_TIME("LinearizedUnlike::computeRandQTerms");
   
   //Get geometry parameters
   const PhaseGeom& phase_geom = a_species.phaseSpaceGeometry();
   const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
   //const CFG::MultiBlockLevelGeom& mag_geom( phase_geom.magGeom() );
   const DisjointBoxLayout& grids( m_F0a.getBoxes() );

   if ( !m_tmp_dfn_noGhosts.isDefined()) m_tmp_dfn_noGhosts.define( grids, 1, IntVect::Zero );
   
   /*
     Compute RabVpar term
     */
   
   //Multiply F0a by vpar * mass_a
   for (DataIterator dit( m_F0a.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(grids[dit]);
      FArrayBox velocityRealCoords(m_tmp_dfn_noGhosts[dit].box(), VEL_DIM);
      block_coord_sys.getVelocityRealCoords(velocityRealCoords);
      
      m_tmp_dfn_noGhosts[dit].copy(m_F0a[dit]);
      
      FORT_MULT_MVPAR(CHF_BOX(m_tmp_dfn_noGhosts[dit].box()),
                      CHF_FRA1(m_tmp_dfn_noGhosts[dit],0),
                      CHF_CONST_FRA(velocityRealCoords),
                      CHF_CONST_REAL(m_mass_a));
   }
         
   //Compute CTP_ab(ma*vpar*F0a)
   testPartCollRHS(a_R, m_tmp_dfn_noGhosts, m_Nb_inj, m_Tb_inj, a_species, m_mass_a, m_mass_b);
      
   //Compute int(vpar*CTP_ab(ma*vpar*F0a),d3v)
   CFG::LevelData<CFG::FArrayBox> tmp_cfg(mag_geom.grids(), 1, CFG::IntVect::Zero);
   m_moment_op.compute(tmp_cfg, a_species, a_R, ParallelVelKernel<FArrayBox>());
      
   LevelData<FArrayBox> injected_data_R;
   phase_geom.injectConfigurationToPhase(tmp_cfg, injected_data_R);
   
   for (DataIterator dit( a_R.dataIterator() ); dit.ok(); ++dit) {
      
      a_R[dit].divide(m_mass_a);
      
      FORT_DIVIDE_CFG_DATA(CHF_BOX(a_R[dit].box()),
                           CHF_FRA1(a_R[dit],0),
                           CHF_CONST_FRA1(injected_data_R[dit],0));
   }
   
   /*
     Compute Qab term
    */

   //Multiply F0a by ma*v^2
   const LevelData<FArrayBox>& B_injected = phase_geom.getBFieldMagnitude();
   for (DataIterator dit( m_F0a.dataIterator() ); dit.ok(); ++dit) {
      
      const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(grids[dit]);
      FArrayBox velocityRealCoords(m_tmp_dfn_noGhosts[dit].box(), VEL_DIM);
      block_coord_sys.getVelocityRealCoords(velocityRealCoords);

      m_tmp_dfn_noGhosts[dit].copy(m_F0a[dit]);
      
      FORT_MULT_MV2(CHF_BOX(m_tmp_dfn_noGhosts[dit].box()),
                    CHF_FRA1(m_tmp_dfn_noGhosts[dit],0),
                    CHF_CONST_FRA(velocityRealCoords),
                    CHF_CONST_FRA1(B_injected[dit],0),
                    CHF_CONST_REAL(m_mass_a));
   }
   
   //Compute CTP_ab(ma*v^2*F0a)
   testPartCollRHS(a_Q, m_tmp_dfn_noGhosts, m_Nb_inj, m_Tb_inj, a_species, m_mass_a, m_mass_b);
   
   //Compute int(ma*v^2/3*CTP_ab(ma*vpar*F0a),d3v)
   CFG::LevelData<CFG::FArrayBox> zero_cfg( mag_geom.grids(), 1, CFG::IntVect::Zero);
   CFG::DataIterator cfg_dit = zero_cfg.dataIterator();
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
      zero_cfg[cfg_dit].setVal(0.);
   }

   m_moment_op.compute(tmp_cfg, a_species, a_Q, PressureKernel<FArrayBox>(zero_cfg));
   
   LevelData<FArrayBox> injected_data_Q;
   phase_geom.injectConfigurationToPhase(tmp_cfg, injected_data_Q);
   
   for (DataIterator dit( a_Q.dataIterator() ); dit.ok(); ++dit) {
      
      a_Q[dit].divide(3.0/2.0);
      
      FORT_DIVIDE_CFG_DATA(CHF_BOX(a_Q[dit].box()),
                           CHF_FRA1(a_Q[dit],0),
                           CHF_CONST_FRA1(injected_data_Q[dit],0));
   }
}

void LinearizedUnlike::computeMaxwellianFit(LevelData<FArrayBox>&             a_F0,
                                            CFG::LevelData<CFG::FArrayBox>&   a_temperature,
                                            CFG::LevelData<CFG::FArrayBox>&   a_density,
                                            LevelData<FArrayBox>&             a_temperature_inj,
                                            LevelData<FArrayBox>&             a_density_inj,
                                            const KineticSpecies&             a_species,
                                            const LevelData<FArrayBox>&       a_dfn,
                                            const Real                        a_time) const
{
    CH_TIME("LinearizedUnlike::computeMaxwellianFit");

    const PhaseGeom& phase_geom = a_species.phaseSpaceGeometry();
    const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
  
    a_species.numberDensity( a_density );
    mag_geom.divideJonValid(a_density);
   
    CFG::LevelData<CFG::FArrayBox> ParallelMom( mag_geom.grids(), 1, CFG::IntVect::Zero );
    a_species.parallelParticleFlux( ParallelMom );
    mag_geom.divideJonValid(ParallelMom);
    
    for (CFG::DataIterator dit(a_density.dataIterator()); dit.ok(); ++dit) {
      ParallelMom[dit].divide(a_density[dit]);
    }

    a_species.pressure(a_temperature, ParallelMom);
    mag_geom.divideJonValid(a_temperature);
   
    for (CFG::DataIterator dit(a_density.dataIterator()); dit.ok(); ++dit) {
      a_temperature[dit].divide(a_density[dit]);
    }
   
    // We reconstruct zero-order Maxwellian, with zero velocity shift
    // We assume that parallel velocity is drift-ordered and therefore is part of delta-f correction
    CFG::LevelData<CFG::FArrayBox> zero_cfg( mag_geom.grids(), 1, CFG::IntVect::Zero);
    CFG::DataIterator cfg_dit = zero_cfg.dataIterator();
    for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
       zero_cfg[cfg_dit].setVal(0.);
    }
   
    MaxwellianKernel<FArrayBox> maxwellian(a_density, a_temperature, zero_cfg);
    maxwellian.eval(a_F0,a_species);
    phase_geom.multJonValid(a_F0);
   
    // Compute injected density and temperature
    phase_geom.injectConfigurationToPhase(a_temperature, a_temperature_inj);
    phase_geom.injectConfigurationToPhase(a_density, a_density_inj);
}

void LinearizedUnlike::getReferenceData(LevelData<FArrayBox>& a_F0,
                                        CFG::LevelData<CFG::FArrayBox>& a_T,
                                        CFG::LevelData<CFG::FArrayBox>& a_N,
					LevelData<FArrayBox>& a_T_inj,
                                        LevelData<FArrayBox>& a_N_inj,
                                        const KineticSpecies&  a_species,
					const RefCountedPtr<KineticFunction>& a_ref_func,
					const RefCountedPtr<CFG::GridFunction>& a_ref_T,
					const RefCountedPtr<CFG::GridFunction>& a_ref_N,
                                        const Real a_time) const
{
      
   // Get geometry parameters
   const PhaseGeom& phase_geom = a_species.phaseSpaceGeometry();
   const CFG::MultiBlockLevelGeom& mag_geom( phase_geom.magGeom() );
   
   // Create reference distributions
   KineticSpeciesPtr ref_species( a_species.clone( IntVect::Unit, false ) );

   a_ref_func->assign( *ref_species, a_time );
   const LevelData<FArrayBox>& dfn( ref_species->distributionFunction() );
   for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
     a_F0[dit].copy(dfn[dit]);
   }

   // Create reference temperature and density profiles
   a_ref_T->assign( a_T, mag_geom, a_time);
   a_ref_N->assign( a_N, mag_geom, a_time);

   // Inject to phase space
   phase_geom.injectConfigurationToPhase(a_T, a_T_inj);
   phase_geom.injectConfigurationToPhase(a_N, a_N_inj);
}
                                              
void LinearizedUnlike::getClsFreqNorm(LevelData<FArrayBox>&             a_cls_norm,
                                      CFG::LevelData<CFG::FArrayBox>&   a_Ta,
                                      CFG::LevelData<CFG::FArrayBox>&   a_Tb,
                                      CFG::LevelData<CFG::FArrayBox>&   a_Na,
                                      CFG::LevelData<CFG::FArrayBox>&   a_Nb,
				      const PhaseGeom&                  a_phase_geom) const

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

    CFG::LevelData<CFG::FArrayBox> cls_norm_cfg(a_Ta.getBoxes(), 1, CFG::IntVect::Zero);
    CFG::DataIterator dit = cls_norm_cfg.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
       CFG::BoxIterator bit(cls_norm_cfg[dit].box());
       for (bit.begin(); bit.ok(); ++bit) {
         CFG::IntVect iv = bit();
         double fac1 = abs(m_charge_a) * abs(m_charge_b) * (m_mass_a + m_mass_b) /
                           (m_mass_a * a_Tb[dit](iv,0) + m_mass_b * a_Ta[dit](iv,0));
         double fac2 = sqrt(a_Na[dit](iv,0)/a_Ta[dit](iv,0)*pow(m_charge_a,2) +
                            a_Nb[dit](iv,0)/a_Tb[dit](iv,0)*pow(m_charge_b,2));
          
         double Coulomb_Lg = 23 - log( fac1 * fac2 * sqrt(N)/1000.0 / pow(T, 3.0/2.0) );
         cls_norm_cfg[dit](iv,0) = coeff * N * pow(ech, 2) * pow(m_charge_a, 2) * pow(m_charge_b, 2) * L
                                 / ( pow(T, 2)) * Coulomb_Lg;
       }
    }
    a_phase_geom.injectConfigurationToPhase(cls_norm_cfg, a_cls_norm);
}

void LinearizedUnlike::convertToCellAverage(const PhaseGeom&  a_phase_geom,
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

void LinearizedUnlike::convertToCellCenters(const PhaseGeom&  a_phase_geom,
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

void LinearizedUnlike::convertToCellFaces(LevelData<FluxBox>& a_faceData,
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
void LinearizedUnlike::parseParameters( ParmParse& a_ppcls )
{
   a_ppcls.query("time_implicit", m_time_implicit);
   a_ppcls.query("fixed_background", m_fixed_bkgr);
   a_ppcls.query("update_frequency", m_update_freq);
   a_ppcls.query("second_order", m_second_order);
   a_ppcls.query("conserve_momentum", m_conserve_momentum);
   a_ppcls.query("conserve_energy", m_conserve_energy);
   a_ppcls.query("include_background_collisions", m_include_bkgr_cls); 
   a_ppcls.query( "verbose", m_verbosity);

   if (m_fixed_bkgr) {
      // Parse reference distribution functions
      KineticFunctionLibrary* library = KineticFunctionLibrary::getInstance();
      std::string function_name;

      a_ppcls.query( "ref_function_a", function_name );
      m_ref_func_a = library->find( function_name );

      a_ppcls.query( "ref_function_b", function_name );
      m_ref_func_b = library->find( function_name );
   
      // Parse reference temeprature and density profiles
      CFG::GridFunctionLibrary* grid_library = CFG::GridFunctionLibrary::getInstance();
      std::string grid_function_name;

      a_ppcls.query( "ref_Ta", grid_function_name );
      m_ref_Ta = grid_library->find( grid_function_name );

      a_ppcls.query( "ref_Tb", grid_function_name );
      m_ref_Tb = grid_library->find( grid_function_name );

      a_ppcls.query( "ref_Na", grid_function_name );
      m_ref_Na = grid_library->find( grid_function_name );

      a_ppcls.query( "ref_Nb", grid_function_name );
      m_ref_Nb = grid_library->find( grid_function_name );
   }
}

inline
void LinearizedUnlike::printParameters()
{
   if (procID()==0) {
      std::cout << "LinearizedUnlike collisions parameters:" << std::endl;
      std::cout << "  fixed_background  =  " << m_fixed_bkgr
                << ", conserve_momentum = " << m_conserve_momentum 
                << ", conserve_energy = " << m_conserve_energy 
                << ", implicit in time = " << m_time_implicit
                << ", update frequency = " << m_update_freq
		<< ", include background collisions  = " << m_include_bkgr_cls
                << std::endl;
      
      if (m_fixed_bkgr) {
         std::cout << "  Reference Function species a:" << std::endl;
         m_ref_func_a->printParameters();
         std::cout << "  Reference Temperature species a:" << std::endl;
         m_ref_Ta->printParameters();
         std::cout << "  Reference Density species a:" << std::endl;
         m_ref_Na->printParameters();

         std::cout << "  Reference Function species b:" << std::endl;
         m_ref_func_b->printParameters();
         std::cout << "  Reference Temperature species b:" << std::endl;
         m_ref_Tb->printParameters();
         std::cout << "  Reference Density species b:" << std::endl;
         m_ref_Nb->printParameters();
      }
   }
}

void LinearizedUnlike::preTimeStep(const KineticSpeciesPtrVect& a_soln_mapped,
                                   const int a_species,
                                   const Real a_time,
                                   const KineticSpeciesPtrVect& a_soln_physical )
{
   m_it_counter+=1;
   m_first_stage = true;
}
                       
Real LinearizedUnlike::computeTimeScale(const KineticSpeciesPtrVect& soln, const int a_idx)
{
   // Replace with an estimate of time scale
   return DBL_MAX;
}

Real LinearizedUnlike::computeDtExplicitTI(const KineticSpeciesPtrVect& soln, const int a_idx)
{
   // Replace with an estimate of dt
   return DBL_MAX;
}

Real LinearizedUnlike::computeDtImExTI(const KineticSpeciesPtrVect& soln, const int a_idx)
{
  if (m_time_implicit) {
    return DBL_MAX;
  } else {
    return computeDtExplicitTI(soln, a_idx);
  }
}

void LinearizedUnlike::diagnostics(const KineticSpecies&  a_rhs_species,
				   const int              a_species) const
{

   // Get geometry parameters
   const PhaseGeom& phase_geom = a_rhs_species.phaseSpaceGeometry();
   const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
   
   //Diagnostic for conservation properties
   std::string file_name;
   file_name = "RHS_data_" + to_string(a_species);
   const LevelData<FArrayBox>& rhs_dfn( a_rhs_species.distributionFunction() );
   phase_geom.plotAtConfigurationIndex(file_name, CFG::IntVect::Unit, rhs_dfn, 0.);
   
   CFG::LevelData<CFG::FArrayBox> cfg_data( mag_geom.grids(), 1, CFG::IntVect::Zero);
   a_rhs_species.numberDensity( cfg_data );
   file_name = "ParticleSink_" + to_string(a_species);
   mag_geom.plotCellData(file_name, cfg_data, 0.);
   
   a_rhs_species.parallelParticleFlux( cfg_data );
   CFG::DataIterator cfg_dit = cfg_data.dataIterator();
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
      cfg_data[cfg_dit].mult(m_mass_a);
   }
   mag_geom.divideJonValid(cfg_data);
   
   file_name = "MomentumSink_" + to_string(a_species);
   mag_geom.plotCellData(file_name, cfg_data, 0.);

   const CFG::LevelData<CFG::FArrayBox>& Bmag = mag_geom.getCCBFieldMag();
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
     cfg_data[cfg_dit].divide(Bmag[cfg_dit]);
   }

   CFG::FluxSurface flux_surface(mag_geom, false);
   CFG::LevelData<CFG::FArrayBox> RadialFlux;
   RadialFlux.define(cfg_data);
   flux_surface.averageAndSpread(cfg_data, RadialFlux);
   file_name = "RadialFlux_" + to_string(a_species);
   mag_geom.plotCellData(file_name, RadialFlux, 0.);
  
   CFG::LevelData<CFG::FArrayBox> zero_cfg( mag_geom.grids(), 1, CFG::IntVect::Zero);
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
      zero_cfg[cfg_dit].setVal(0.);
   }
   a_rhs_species.pressure(cfg_data, zero_cfg);
   mag_geom.divideJonValid(cfg_data);
 
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
      //cfg_data[cfg_dit].mult(m_mass_a);
   }
   file_name = "EnergySink_" + to_string(a_species);
   mag_geom.plotCellData(file_name, cfg_data, 0.);
   
   

}

#include "NamespaceFooter.H"
