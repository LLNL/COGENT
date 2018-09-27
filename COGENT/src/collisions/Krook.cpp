#include <math.h>

#include "Krook.H"

#include "FourthOrderUtil.H"
#include "Directions.H"
#include "EdgeToCell.H"
#include "PhaseGeom.H"
#include "PhaseBlockCoordSys.H"
#include "CollisionsF_F.H"
#include "KineticFunctionLibrary.H"
#include "ConstFact.H"
#include "inspect.H"

#include "NamespaceHeader.H" 


Krook::Krook( const std::string& a_ppcls_str, const int a_verbosity )
   : m_verbosity(a_verbosity),
     m_time_implicit(true),
     m_cls_freq(-1.0),
     m_conserve_particle(false),
     m_conserve_momentum(false),
     m_first_step(true)
{
   ParmParse ppcls( a_ppcls_str.c_str() );
   parseParameters( ppcls );
   if (m_verbosity>0) {
      printParameters();
   }

   m_fixed_cls_freq = (m_cls_freq<0.0) ? false : true ;
}


Krook::~Krook()
{
}


void Krook::evalClsRHS( KineticSpeciesPtrVect&       a_rhs,
                        const KineticSpeciesPtrVect& a_soln,
                        const int                    a_species,
                        const Real                   a_time )
// NB: a_soln is on the computational grid
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
   LevelData<FArrayBox> delta_F( grids, n_comp, IntVect::Zero );
   
   for (DataIterator sdit(soln_dfn.dataIterator()); sdit.ok(); ++sdit) {
      delta_F[sdit].copy( soln_dfn[sdit] );
      delta_F[sdit].minus( init_dfn[sdit] );
   }
   
   // Calculate test-particle (TP) collisional RHS
   LevelData<FArrayBox> tp_rhs_coll( grids, n_comp, IntVect::Zero );
   testPartCollRHS( tp_rhs_coll, delta_F );
   
   // Get right-hand side distribution function for the current species
   KineticSpecies& rhs_species( *(a_rhs[a_species]) );
   if (enforceConservation()) {
      addConservativeTerms( rhs_species,
                            soln_species,
                            tp_rhs_coll,
                            delta_F,
                            a_time );
   }
   
   // Add test-particle collisional RHS to Vlasov RHS
   LevelData<FArrayBox>& rhs_dfn( rhs_species.distributionFunction() );
   for (DataIterator rdit(soln_dfn.dataIterator()); rdit.ok(); ++rdit) {
      if (m_fixed_cls_freq) {tp_rhs_coll[rdit].mult( m_cls_freq );}
      else {tp_rhs_coll[rdit].mult( m_sc_cls_freq[rdit] );}
      rhs_dfn[rdit].plus( tp_rhs_coll[rdit] );
   }  

   m_first_step = false;
}


void Krook::addConservativeTerms( KineticSpecies& a_rhs_species,
                                  const KineticSpecies& a_soln_species,
                                  const LevelData<FArrayBox>& a_tp_rhs_coll,
                                  const LevelData<FArrayBox>& a_delta_F,
                                  const Real a_time )
{
   const LevelData<FArrayBox>& soln_dfn( a_soln_species.distributionFunction() );
   const DisjointBoxLayout& grids( soln_dfn.getBoxes() );
   const int n_comp( soln_dfn.nComp() );

   //Create kernels for particle and momentum conervative terms
   LevelData<FArrayBox> kern_part( grids, n_comp, IntVect::Zero );
   LevelData<FArrayBox> kern_moment( grids, n_comp, IntVect::Zero );
   const PhaseGeom& phase_geom( a_soln_species.phaseSpaceGeometry() );
   const double mass( a_soln_species.mass() );
   collKernels( kern_part, kern_moment, a_tp_rhs_coll, phase_geom, mass );
   
   // Copy const rhs_dfn to a temporary
   LevelData<FArrayBox>& rhs_dfn( a_rhs_species.distributionFunction() );
   LevelData<FArrayBox> rhs_copy( grids, n_comp, IntVect::Zero );
   for (DataIterator rhsdit(rhs_copy.dataIterator()); rhsdit.ok(); ++rhsdit) {
      rhs_copy[rhsdit].copy(rhs_dfn[rhsdit]);
   }
     
   //Calculate restoring terms using DensityOperator, which calculates dvpar*dmu integrals
   const CFG::MultiBlockLevelGeom& mag_geom( phase_geom.magGeom() );
   CFG::LevelData<CFG::FArrayBox> part_rest( mag_geom.grids(), 1, CFG::IntVect::Zero );
   for (DataIterator rdit( rhs_dfn.dataIterator() ); rdit.ok(); ++rdit) {
      rhs_dfn[rdit].copy( kern_part[rdit] );
   }
   a_rhs_species.numberDensity( part_rest );
   
   CFG::LevelData<CFG::FArrayBox> moment_rest( mag_geom.grids(), 1, CFG::IntVect::Zero );
   for (DataIterator rdit( rhs_dfn.dataIterator() ); rdit.ok(); ++rdit) {
      rhs_dfn[rdit].copy( kern_moment[rdit] );
   }
   a_rhs_species.numberDensity( moment_rest );
   
   //Calculate normalizing factors
   if (m_first_step) {
      
      LevelData<FArrayBox> kern_moment_norm( grids, n_comp, IntVect::Zero );
      LevelData<FArrayBox> kern_part_norm( grids, n_comp, IntVect::Zero );
      evaluateNormKern( kern_moment_norm, kern_part_norm, phase_geom, mass);
      
      m_norm_momentum.define( mag_geom.grids(), 1, CFG::IntVect::Zero );
      m_norm_particle.define( mag_geom.grids(), 1, CFG::IntVect::Zero );
      
      for (DataIterator rdit( rhs_dfn.dataIterator() ); rdit.ok(); ++rdit) {
         rhs_dfn[rdit].copy( kern_moment_norm[rdit] );
      }
      a_rhs_species.numberDensity( m_norm_momentum );
      
      for (DataIterator rdit( rhs_dfn.dataIterator() ); rdit.ok(); ++rdit) {
         rhs_dfn[rdit].copy( kern_part_norm[rdit] );
      }
      a_rhs_species.numberDensity( m_norm_particle );
   }
   
   // Calculate field-particle (FP) collisional RHS
   LevelData<FArrayBox> fp_rhs_coll( grids, n_comp, IntVect::Zero );
   fieldPartCollRHS( fp_rhs_coll,
                     a_delta_F,
                     phase_geom,
                     part_rest,
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


void Krook::testPartCollRHS( LevelData<FArrayBox>& a_rhs_coll,
                             const LevelData<FArrayBox>& a_delta_F ) const
{
   const Real fact( -1.0 );
   for (DataIterator rdit( a_rhs_coll.dataIterator() ); rdit.ok(); ++rdit) {
      a_rhs_coll[rdit].copy(a_delta_F[rdit]);
      a_rhs_coll[rdit].mult(fact);
   }
}


void Krook::fieldPartCollRHS( LevelData<FArrayBox>& a_rhs_coll,
                              const LevelData<FArrayBox>& a_delta_F,
                              const PhaseGeom& a_phase_geom,
                              const CFG::LevelData<CFG::FArrayBox>& a_rest_particle,
                              const CFG::LevelData<CFG::FArrayBox>& a_rest_momentum,
                              const double a_mass ) const
{
   const DisjointBoxLayout& grids( a_rhs_coll.getBoxes() );
   const int n_comp( a_rhs_coll.nComp() );

   LevelData<FArrayBox> full_PR( grids, n_comp, IntVect::Zero );
   LevelData<FArrayBox> full_MR( grids, n_comp, IntVect::Zero );

   LevelData<FArrayBox> inj_PRest;
   a_phase_geom.injectConfigurationToPhase( a_rest_particle, inj_PRest );

   LevelData<FArrayBox> inj_MRest;
   a_phase_geom.injectConfigurationToPhase( a_rest_momentum, inj_MRest );

   LevelData<FArrayBox> inj_PNorm;
   a_phase_geom.injectConfigurationToPhase( m_norm_particle, inj_PNorm );

   LevelData<FArrayBox> inj_MNorm;
   a_phase_geom.injectConfigurationToPhase( m_norm_momentum, inj_MNorm );
   
   const LevelData<FArrayBox>& B_injected( a_phase_geom.getBFieldMagnitude() );

   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys( a_phase_geom.getBlockCoordSys(grids[dit]) );

      FArrayBox& this_RHS( a_rhs_coll[dit] );
      FArrayBox& this_full_PR( full_PR[dit] );
      FArrayBox& this_full_MR( full_MR[dit] );

      const FArrayBox& this_PRest( inj_PRest[dit] );
      const FArrayBox& this_MRest( inj_MRest[dit] );
      const FArrayBox& this_PNorm( inj_PNorm[dit] );
      const FArrayBox& this_MNorm( inj_MNorm[dit] );
      
      const FArrayBox& this_B( B_injected[dit] );
      const FArrayBox& this_Temp_Distr( m_temperature[dit] );
      
      const Box& box_B( this_B.box() );
      const int vpB( box_B.smallEnd( VPARALLEL_DIR ) );
      const int muB( box_B.smallEnd( MU_DIR ) );

      const Box& box_PR( this_PRest.box() );
      const int vpR( box_PR.smallEnd( VPARALLEL_DIR ) );
      const int muR( box_PR.smallEnd( MU_DIR ) );

      // Get the physical velocity coordinates (Vpar, mu) for this part
      // of phase space; Vpar and Mu are defined at the cell centers
      FArrayBox velocityRealCoords( this_RHS.box(), VEL_DIM );
      block_coord_sys.getVelocityRealCoords( velocityRealCoords );
      
      //Iterate over the points inside the box, iv is a 4-component int vect
      for (BoxIterator bit( this_full_PR.box() ); bit.ok(); ++bit) {
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
         double exp_m_x_sq( exp( -x_sq ) );

         //Calculate restoring terms Fm0*(Vpar*U+P)
         IntVect ivR( iv );
         ivR[VPARALLEL_DIR] = vpR;
         ivR[MU_DIR] = muR;

         this_full_MR(iv,0) = -( this_MRest(ivR) / this_MNorm(ivR) );
         this_full_MR(iv,0) *= v_parallel * exp_m_x_sq;

         this_full_PR(iv,0) = -( this_PRest(ivR) / this_PNorm(ivR) );
         this_full_PR(iv,0) *= exp_m_x_sq;
         
         this_RHS(iv,0) = 0.0;
         this_RHS(iv,0) += (m_conserve_particle) ? this_full_PR(iv,0) : 0;
         this_RHS(iv,0) += (m_conserve_momentum) ? this_full_MR(iv,0) : 0;
      }
   }
}


void Krook::collKernels( LevelData<FArrayBox>& a_kern_part,
                         LevelData<FArrayBox>& a_kern_moment,
                         const LevelData<FArrayBox>& a_test_part_coll_RHS,
                         const PhaseGeom& a_phase_geom,
                         const double a_mass ) const

{
   const LevelData<FArrayBox>& B_injected( a_phase_geom.getBFieldMagnitude() );

   const DisjointBoxLayout& grids( a_test_part_coll_RHS.getBoxes() );
   const int n_comp( a_test_part_coll_RHS.nComp() );
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys( a_phase_geom.getBlockCoordSys(grids[dit]) );

      FArrayBox& this_KP( a_kern_part[dit] );
      FArrayBox& this_KM( a_kern_moment[dit] );
      const FArrayBox& this_RHS( a_test_part_coll_RHS[dit] );
      const FArrayBox& this_B( B_injected[dit] );

      // Get the physical velocity coordinates for this part of phase space
      const Box& box_kern( this_KP.box() );
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
            this_KP(iv,dfn_comp) = this_RHS(iv,dfn_comp);
            this_KM(iv,dfn_comp) = this_RHS(iv,dfn_comp) * v_parallel;
         }
      }
   }
}



void Krook::evaluateNormKern( LevelData<FArrayBox>& a_kern_moment_norm,
                              LevelData<FArrayBox>& a_kern_part_norm,
                              const PhaseGeom& a_phase_geom,
                              const double a_mass ) const
{
   const LevelData<FArrayBox>& B_injected( a_phase_geom.getBFieldMagnitude() );

   const DisjointBoxLayout& grids( a_kern_moment_norm.disjointBoxLayout() );
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys( a_phase_geom.getBlockCoordSys(grids[dit]) );

      FArrayBox& this_PN( a_kern_part_norm[dit] );
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
         
         this_PN(iv,0) = exp( -x_sq );
         this_MN(iv,0) = v_par_sq * this_PN(iv,0);
      }
   }
}


void Krook::computeSelfConsistFreq(LevelData<FArrayBox>& a_cls_freq,
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

inline
void Krook::parseParameters( ParmParse& a_ppcls )
{
   a_ppcls.query( "time_implicit", m_time_implicit);
   a_ppcls.query( "cls_freq", m_cls_freq );
   a_ppcls.query( "conserve_particle", m_conserve_particle );
   a_ppcls.query( "conserve_momentum", m_conserve_momentum );

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
void Krook::printParameters()
{
   if (procID()==0) {
      std::cout << "Krook collisions parameters:" << std::endl;
      std::cout << "  cls_freq  =  " << m_cls_freq
                << ", conserve_particle = " << m_conserve_particle
                << ", conserve_momentum = " << m_conserve_momentum 
                << ", implicit in time = " << m_time_implicit
                << std::endl;
      std::cout << "  Reference Function:" << std::endl;
      m_ref_func->printParameters();
      std::cout << "  Reference Temperature:" << std::endl;
      m_ref_temp->printParameters();
      std::cout << "  Reference Density:" << std::endl;
      m_ref_dens->printParameters();

   }
}

Real Krook::collisionFrequency()
{
  /* implemented for fixed collision frequency for now */
  CH_assert(m_fixed_cls_freq);
  return m_cls_freq;
}

void Krook::addReferenceDfn( KineticSpecies& a_result,
                             const Real      a_time, 
                             const Real      a_scale )
{
   LevelData<FArrayBox>& result_dfn( a_result.distributionFunction() );
         
   KineticSpeciesPtr ref_species( a_result.clone( IntVect::Unit, false ) );
   m_ref_func->assign( *ref_species, a_time );
   const LevelData<FArrayBox>& ref_dfn( ref_species->distributionFunction() );
   for (DataIterator dit(result_dfn.dataIterator()); dit.ok(); ++dit) {
      result_dfn[dit].plus( ref_dfn[dit], a_scale );
   }
}

Real Krook::computeDt(const KineticSpeciesPtrVect& soln)
{
   if (m_fixed_cls_freq)  return 1.0/m_cls_freq;
   else {
      Real max_freq = 0.0;
      const DisjointBoxLayout & grids = m_sc_cls_freq.disjointBoxLayout();
      for (DataIterator dit(m_sc_cls_freq.dataIterator()); dit.ok(); ++dit) {
         FArrayBox& this_data = m_sc_cls_freq[dit];
         Box box(grids[dit]);
         if (this_data.max(box) > max_freq) max_freq = this_data.max(box);
      }
      return 1.0/max_freq;
   } 
}

Real Krook::TimeScale(const KineticSpeciesPtrVect& soln)
{
   if (m_fixed_cls_freq)  return 1.0/m_cls_freq;
   else {
      Real max_freq = 0.0;
      const DisjointBoxLayout & grids = m_sc_cls_freq.disjointBoxLayout();
      for (DataIterator dit(m_sc_cls_freq.dataIterator()); dit.ok(); ++dit) {
         FArrayBox& this_data = m_sc_cls_freq[dit];
         Box box(grids[dit]);
         if (this_data.max(box) > max_freq) max_freq = this_data.max(box);
      }
      return 1.0/max_freq;
   } 
}

#include "NamespaceFooter.H"
