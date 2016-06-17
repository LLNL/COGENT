#include <math.h>

#include "FourthOrderUtil.H"
#include "Directions.H"
#include "EdgeToCell.H"
#include "Lorentz.H"
#include "CollisionsF_F.H"
#include "KineticFunctionLibrary.H"
#include "ConstFact.H"

#include "NamespaceHeader.H" //Should be the last one



Lorentz::Lorentz( ParmParse& a_ppcls, const int a_verbosity )
    : m_verbosity(a_verbosity),
      m_cls_freq(-1.0),
      m_constant_freq(false),
      m_conserve_momentum(true),
      m_first_step(true)
{
   parseParameters( a_ppcls );
   if (m_verbosity>0) {
      printParameters();
   }

   m_fixed_cls_freq = (m_cls_freq<0.0) ? false : true ;
}

Lorentz::~Lorentz()
{
}

void Lorentz::evalClsRHS( KineticSpeciesPtrVect&       a_rhs,
                          const KineticSpeciesPtrVect& a_soln,
                          const int                    a_species,
                          const Real                   a_time,
                          const int                    a_flag)

// NB: a_soln is on the computational grid and has 4 ghost cells (passed here as Nans)
{
   // Get solution distribution function (J*Bstar_par*dfn) for the current species
   const KineticSpecies& soln_species( *(a_soln[a_species]) );
   const LevelData<FArrayBox>& soln_dfn( soln_species.distributionFunction() );
   
   // Create reference (J*Bstar_par*dfn_init) distribution
   KineticSpeciesPtr ref_species( soln_species.clone( IntVect::Unit, false ) );
   m_ref_func->assign( *ref_species, a_time );
   const LevelData<FArrayBox>& init_dfn( ref_species->distributionFunction() );

   //Create reference temperature distribution
   const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();
   if (m_first_step) {
    const CFG::MultiBlockLevelGeom& mag_geom( phase_geom.magGeom() );
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
   const DisjointBoxLayout& grids( soln_dfn.getBoxes() );
   const int n_comp( soln_dfn.nComp() );
   LevelData<FArrayBox> delta_dfn( grids, n_comp, IntVect::Zero );
   
   for (DataIterator sdit(soln_dfn.dataIterator()); sdit.ok(); ++sdit) {
      delta_dfn[sdit].copy( soln_dfn[sdit] );
      delta_dfn[sdit].minus( init_dfn[sdit] );
   }

   //Covert delta_dfn to cell centers
   convertToCellCenters(phase_geom, delta_dfn);
   

   // Calculate test-particle (TP) collisional RHS
   const double mass = soln_species.mass();
   LevelData<FArrayBox> tp_rhs_coll(grids, n_comp, IntVect::Zero);
   testPartCollRHS(tp_rhs_coll, delta_dfn, phase_geom, mass);

   // Add conservative terms (field-particle terms)
   KineticSpecies& rhs_species( *(a_rhs[a_species]) );
   if (enforceConservation()) {
      addConservativeTerms( rhs_species,
                            soln_species,
                            tp_rhs_coll,
                            delta_dfn,
                            a_time );
   }

   // Add test-particle collisional RHS to Vlasov RHS (rhs_dfn has no ghost cells)
   LevelData<FArrayBox>& rhs_dfn( rhs_species.distributionFunction() );
   for (DataIterator rdit(rhs_dfn.dataIterator()); rdit.ok(); ++rdit) {
      if (m_fixed_cls_freq) {tp_rhs_coll[rdit].mult( m_cls_freq );}
      else {tp_rhs_coll[rdit].mult( m_sc_cls_freq[rdit] );}
      rhs_dfn[rdit].plus( tp_rhs_coll[rdit] );
   }  

   m_first_step = false;
}


void Lorentz::addConservativeTerms( KineticSpecies& a_rhs_species,
                                  const KineticSpecies& a_soln_species,
                                  const LevelData<FArrayBox>& a_tp_rhs_coll,
                                  const LevelData<FArrayBox>& a_delta_dfn,
                                  const Real a_time )
{

   const LevelData<FArrayBox>& soln_dfn( a_soln_species.distributionFunction() );
   const DisjointBoxLayout& grids( soln_dfn.getBoxes() );
   const int n_comp( soln_dfn.nComp() );

   //Create kernels for momentum conervative terms
   LevelData<FArrayBox> kern_moment( grids, n_comp, IntVect::Zero );
   const PhaseGeom& phase_geom( a_soln_species.phaseSpaceGeometry() );
   const double mass( a_soln_species.mass() );
   collKernels( kern_moment, a_tp_rhs_coll, phase_geom, mass );
   
   // Copy const rhs_dfn to a temporary
   LevelData<FArrayBox>& rhs_dfn( a_rhs_species.distributionFunction() );
   LevelData<FArrayBox> rhs_copy( grids, n_comp, IntVect::Zero );
   for (DataIterator rhsdit(rhs_copy.dataIterator()); rhsdit.ok(); ++rhsdit) {
      rhs_copy[rhsdit].copy(rhs_dfn[rhsdit]);
   }
     
   //Calculate restoring terms using DensityOperator, which calculates dvpar*dmu integrals
   const CFG::MultiBlockLevelGeom& mag_geom( phase_geom.magGeom() );
   CFG::LevelData<CFG::FArrayBox> moment_rest( mag_geom.grids(), 1, CFG::IntVect::Zero );
   for (DataIterator rdit( rhs_dfn.dataIterator() ); rdit.ok(); ++rdit) {
      rhs_dfn[rdit].copy( kern_moment[rdit] );
   }
   a_rhs_species.numberDensity( moment_rest );
   
   //Calculate normalizing factors
   if (m_first_step) {
      
      LevelData<FArrayBox> kern_moment_norm( grids, n_comp, IntVect::Zero );
      evaluateNormKern( kern_moment_norm, phase_geom, mass);
      
      m_norm_momentum.define( mag_geom.grids(), 1, CFG::IntVect::Zero );
      
      for (DataIterator rdit( rhs_dfn.dataIterator() ); rdit.ok(); ++rdit) {
         rhs_dfn[rdit].copy( kern_moment_norm[rdit] );
      }
      a_rhs_species.numberDensity( m_norm_momentum );
   }
   
   // Calculate field-particle (FP) collisional RHS
   LevelData<FArrayBox> fp_rhs_coll( grids, n_comp, IntVect::Zero );
   fieldPartCollRHS( fp_rhs_coll,
                     a_delta_dfn,
                     phase_geom,
                     moment_rest,
                     mass );
   
   // Add field-particle collisional RHS to Vlasov RHS
   for (DataIterator rdit( rhs_dfn.dataIterator() ); rdit.ok(); ++rdit) {
      rhs_dfn[rdit].copy( rhs_copy[rdit] );
      if (m_fixed_cls_freq) {fp_rhs_coll[rdit].mult( m_cls_freq );}
      else {fp_rhs_coll[rdit].mult( m_sc_cls_freq[rdit] );}
      rhs_dfn[rdit].plus( fp_rhs_coll[rdit] );
   }
}

void Lorentz::testPartCollRHS(LevelData<FArrayBox>& a_rhs_coll,
                        const LevelData<FArrayBox>& a_delta_dfn,
                        const PhaseGeom& a_phase_geom,
                        const double a_mass ) const
{

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

   //Create temporary delta_dfn with two extra layers of ghost cells
   LevelData<FArrayBox> delta_dfn_tmp(a_delta_dfn.disjointBoxLayout(),
                                a_delta_dfn.nComp(),
                                a_delta_dfn.ghostVect()+2*IntVect::Unit);

   for (DataIterator dit( a_rhs_coll.dataIterator() ); dit.ok(); ++dit) {
      delta_dfn_tmp[dit].setVal(0.0);
      delta_dfn_tmp[dit].copy(a_delta_dfn[dit],grids[dit]);
   }
   delta_dfn_tmp.exchange();

   //Create cell-centered collsion fluxes (0 comp - vpar_dir, 1 comp - mu_dir)
   LevelData<FArrayBox> flux_cell(grids, 2, IntVect::Zero);
   for (DataIterator dit( a_rhs_coll.dataIterator() ); dit.ok(); ++dit) {
      FArrayBox& this_flux_cell = flux_cell[dit];
      const FArrayBox& this_delta_dfn_tmp = delta_dfn_tmp[dit];
      const FArrayBox& this_b = injected_B[dit];

      if(m_constant_freq) {
         FORT_EVALUATE_TP_LORENTZ_CONST_NUD(CHF_FRA(this_flux_cell),
                                            CHF_CONST_FRA1(this_delta_dfn_tmp,0),
                                            CHF_CONST_FRA1(this_b,0),
                                            CHF_BOX(this_flux_cell.box()),
                                            CHF_CONST_REALVECT(vel_dx),
                                            CHF_CONST_REAL(a_mass));
       }

       else {
         const FArrayBox& this_temperature = m_temperature[dit];
         FORT_EVALUATE_TP_LORENTZ(CHF_FRA(this_flux_cell),
                                  CHF_CONST_FRA1(this_delta_dfn_tmp,0),
                                  CHF_CONST_FRA1(this_b,0),
                                  CHF_CONST_FRA1(this_temperature,0),
                                  CHF_BOX(this_flux_cell.box()),
                                  CHF_CONST_REALVECT(vel_dx),
                                  CHF_CONST_REAL(a_mass));
      }
   }

   //Calculate collision flux cell-average
   convertToCellAverage(a_phase_geom, flux_cell);


   //Calculate collision fluxes on cell faces 
   LevelData<FluxBox> flux_face(grids, 2, IntVect::Zero);
   convertToCellFaces(flux_face, flux_cell);

   //Create final (combined r , theta, mu, and v_par) rhs flux
   LevelData<FluxBox> flux_rhs(grids, SpaceDim, IntVect::Zero);
   for (DataIterator dit(flux_rhs.dataIterator() ); dit.ok(); ++dit) {
      const FArrayBox& this_flux_cell = flux_cell[dit];

      FArrayBox tmp_flux_cell(grow(this_flux_cell.box(),1),2);
      tmp_flux_cell.setVal(0.);
      tmp_flux_cell.copy(this_flux_cell);

      for (int dir=0; dir<SpaceDim; dir++) {
          FArrayBox& this_flux_rhs = flux_rhs[dit][dir];
          FArrayBox& this_flux_face = flux_face[dit][dir];
          FORT_EVALUATE_COLL_FLUX_COMBINED(CHF_FRA(this_flux_rhs),
                                           CHF_CONST_INT(dir),
                                           CHF_BOX(this_flux_rhs.box()),
                                           CHF_CONST_FRA(this_flux_face),
                                           CHF_CONST_FRA(tmp_flux_cell),
                                           CHF_CONST_INT(num_vpar_cells),
                                           CHF_CONST_INT(num_mu_cells));
      }
   }

   a_phase_geom.mappedGridDivergence(a_rhs_coll, flux_rhs, true);
   for (DataIterator dit( a_rhs_coll.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
      double fac = 1. / block_coord_sys.getMappedCellVolume();
      a_rhs_coll[dit].mult(fac);
    }

}

void Lorentz::fieldPartCollRHS(LevelData<FArrayBox>& a_rhs_coll,
                               const LevelData<FArrayBox>& a_delta_dfn,
                               const PhaseGeom& a_phase_geom,
                               const CFG::LevelData<CFG::FArrayBox>& a_moment_rest,
                               const double a_mass ) const
{
   const DisjointBoxLayout& grids( a_rhs_coll.getBoxes() );
   const int n_comp( a_rhs_coll.nComp() );

   LevelData<FArrayBox> full_MR( grids, n_comp, IntVect::Zero );

   LevelData<FArrayBox> inj_MRest;
   a_phase_geom.injectConfigurationToPhase(a_moment_rest, inj_MRest);
   
   LevelData<FArrayBox> inj_MNorm;
   a_phase_geom.injectConfigurationToPhase(m_norm_momentum, inj_MNorm);

   const LevelData<FArrayBox>& B_injected = a_phase_geom.getBFieldMagnitude();

   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys( a_phase_geom.getBlockCoordSys(grids[dit]) );

      FArrayBox& this_RHS( a_rhs_coll[dit] );
      FArrayBox& this_full_MR( full_MR[dit] );

      const FArrayBox& this_MRest( inj_MRest[dit] );
      const FArrayBox& this_MNorm( inj_MNorm[dit] );

      const FArrayBox& this_B( B_injected[dit] );
      const FArrayBox& this_Temp_Distr( m_temperature[dit] );
      
      const Box& box_B( this_B.box() );
      const int vpB( box_B.smallEnd( VPARALLEL_DIR ) );
      const int muB( box_B.smallEnd( MU_DIR ) );

      const Box& box_PR( this_MRest.box() );
      const int vpR( box_PR.smallEnd( VPARALLEL_DIR ) );
      const int muR( box_PR.smallEnd( MU_DIR ) );

      // Get the physical velocity coordinates (Vpar, mu) for this part
      // of phase space; Vpar and Mu are defined at the cell centers
      FArrayBox velocityRealCoords( this_RHS.box(), VEL_DIM );
      block_coord_sys.getVelocityRealCoords( velocityRealCoords );

      //Iterate over the points inside the box, iv is a 4-component int vect
      for (BoxIterator bit( this_full_MR.box() ); bit.ok(); ++bit) {
         IntVect iv( bit() );
         IntVect ivB( iv );
         ivB[VPARALLEL_DIR] = vpB;
         ivB[MU_DIR] = muB;
         
         double v_parallel( velocityRealCoords(iv,0) );
         double mu( velocityRealCoords(iv,1) );
         double v_perp_sq( mu * this_B(ivB) / a_mass );
         double v_sq( v_parallel * v_parallel + v_perp_sq );
         double v_th_sq( 2.0 * this_Temp_Distr(ivB,0) / a_mass );
         double x_sq( v_sq / v_th_sq );
         double x( sqrt(v_sq / v_th_sq) );
         double exp_m_x_sq( exp( -x_sq ) );

         //Calculate momentum restoring term
         IntVect ivR( iv );
         ivR[VPARALLEL_DIR] = vpR;
         ivR[MU_DIR] = muR;

         this_full_MR(iv,0) = -( this_MRest(ivR) / this_MNorm(ivR) );
         this_full_MR(iv,0) *= v_parallel * exp_m_x_sq;

         this_RHS(iv,0) = this_full_MR(iv,0);
         this_RHS(iv,0) *= (m_constant_freq) ? 1.0 : nu_D(x);         
      }
   }

}

void Lorentz::collKernels(LevelData<FArrayBox>& a_kern_moment,
                         const LevelData<FArrayBox>& a_test_part_coll_RHS,
                         const PhaseGeom& a_phase_geom,
                         const double a_mass ) const


{
   const LevelData<FArrayBox>& B_injected( a_phase_geom.getBFieldMagnitude() );

   const DisjointBoxLayout& grids( a_test_part_coll_RHS.getBoxes() );
   const int n_comp( a_test_part_coll_RHS.nComp() );
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys( a_phase_geom.getBlockCoordSys(grids[dit]) );

      FArrayBox& this_KM( a_kern_moment[dit] );
      const FArrayBox& this_RHS( a_test_part_coll_RHS[dit] );
      const FArrayBox& this_B( B_injected[dit] );

      // Get the physical velocity coordinates for this part of phase space
      const Box& box_kern( this_KM.box() );
      FArrayBox velocityRealCoords( box_kern, VEL_DIM );
      block_coord_sys.getVelocityRealCoords( velocityRealCoords );

      const Box& box_B( this_B.box() );
      const int vp_index( box_B.smallEnd( VPARALLEL_DIR ) );
      const int mu_index( box_B.smallEnd( MU_DIR ) );
      for (BoxIterator bit( box_kern ); bit.ok(); ++bit) {
         IntVect iv( bit() );
         IntVect ivB( iv );
         ivB[VPARALLEL_DIR] = vp_index;
         ivB[MU_DIR] = mu_index;
         const double v_parallel( velocityRealCoords(iv,0) );

         // Calculate the contributions to the parallel velocity and
         // number of particle moments created by the test-particle operator
         for (int dfn_comp(0); dfn_comp<n_comp; ++dfn_comp) {
            this_KM(iv,dfn_comp) = this_RHS(iv,dfn_comp) * v_parallel;
         }
      }
   }
}


void Lorentz::evaluateNormKern( LevelData<FArrayBox>& a_kern_moment_norm,
                              const PhaseGeom& a_phase_geom,
                              const double a_mass ) const
{
   const LevelData<FArrayBox>& B_injected( a_phase_geom.getBFieldMagnitude() );

   const DisjointBoxLayout& grids( a_kern_moment_norm.disjointBoxLayout() );
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys( a_phase_geom.getBlockCoordSys(grids[dit]) );

      FArrayBox& this_MN( a_kern_moment_norm[dit] );
      const FArrayBox& this_temp_distr( m_temperature[dit] );
      const FArrayBox& this_B( B_injected[dit] );

      const Box& box_B( this_B.box() );
      const Box& box_k( this_MN.box() );

      const int vp_index( box_B.smallEnd( VPARALLEL_DIR ) );
      const int mu_index( box_B.smallEnd( MU_DIR ) );

      // Get the physical velocity coordinates for this part of phase space
      FArrayBox velocityRealCoords( box_k, VEL_DIM );
      block_coord_sys.getVelocityRealCoords( velocityRealCoords );

      for (BoxIterator bit( velocityRealCoords.box() ); bit.ok(); ++bit) {
          IntVect iv( bit() );
          IntVect ivB( iv );
          ivB[VPARALLEL_DIR] = vp_index;
          ivB[MU_DIR] = mu_index;
          double v_parallel( velocityRealCoords(iv,0) );
          double v_par_sq( v_parallel * v_parallel );
          double mu( velocityRealCoords(iv,1) );
          double v_perp_sq( mu * this_B(ivB) / a_mass );
          double v_sq( v_par_sq + v_perp_sq );
          double v_th_sq( 2.0 * this_temp_distr(ivB,0) / a_mass );
          double x_sq( v_sq / v_th_sq );
          double x( sqrt(v_sq / v_th_sq) );

          this_MN(iv,0) = v_par_sq * exp( -x_sq );
          this_MN(iv,0) *= (m_constant_freq) ? 1.0 : nu_D(x);
      }
   }
}

void Lorentz::computeSelfConsistFreq(LevelData<FArrayBox>& a_cls_freq,
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


void Lorentz::convertToCellAverage(const PhaseGeom&       a_phase_geom,
                                    LevelData<FArrayBox>& a_dfn) const
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

void Lorentz::convertToCellCenters(const PhaseGeom&       a_phase_geom,
                                    LevelData<FArrayBox>& a_dfn) const
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



void Lorentz::convertToCellFaces(LevelData<FluxBox>& a_faceData,
                                 const LevelData<FArrayBox>& a_cellData) const
{
   // Convert from cell centered values to cell averages.  We need to use
   // a temporary with two extra ghost cells for the interpolation. The ghost cells at the domain 
   // boundaries are filled with zeros for our purposes. Since, the boundary fluxes in the velocity directions 
   // will be overwrited later, and the fluxes in the configuration domansions do not contribute to the RHS. 

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



double Lorentz::nu_D(double x) const 
{
      double pi=3.14159265358979;
      double G=(erf(x)-x*2.0/sqrt(pi)*exp(-x*x))/(2.0*x*x);
      double F=erf(x) ;

      return 1/(x*x*x)*(F-G);

}



inline
void Lorentz::parseParameters( ParmParse& a_ppcls )
{
   a_ppcls.query("cls_freq",m_cls_freq);
   a_ppcls.query("const_freq",m_constant_freq);
   a_ppcls.query("conserve_momentum",m_conserve_momentum);

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
void Lorentz::printParameters()
{
   if (procID()==0) {
      std::cout << "Lorentz collisions parameters:" << std::endl;
      std::cout << "  cls_freq  =  " << m_cls_freq
                << ", const_freq = " << m_constant_freq
                << ", conserve_momentum = " << m_conserve_momentum << std::endl;
      std::cout << "  Reference Function:" << std::endl;
      m_ref_func->printParameters();
      std::cout << "  Reference Temperature:" << std::endl;
      m_ref_temp->printParameters();
      std::cout << "  Reference Density:" << std::endl;
      m_ref_dens->printParameters();

   }
}


Real Lorentz::computeDt(const KineticSpeciesPtrVect& soln)
{
   return DBL_MAX;
}

#include "NamespaceFooter.H"
