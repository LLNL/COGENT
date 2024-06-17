#include <math.h>

#include "BGK.H"

#include "FourthOrderUtil.H"
#include "Directions.H"
#include "EdgeToCell.H"
#include "PhaseGeom.H"
#include "PhaseBlockCoordSys.H"
#include "CollisionsF_F.H"
#include "KineticFunctionLibrary.H"
#include "Kernels.H"
#include "MomentOp.H"
#include "ConstFact.H"
#include "inspect.H"

#include "NamespaceHeader.H"


BGK::BGK( const std::string& a_ppcls_str, const int a_verbosity )
   : m_verbosity(a_verbosity),
     m_time_implicit(true),
     m_cls_freq(-1.0),
     m_conserve_two_moments(false),
     m_conserve_three_moments(false),
     m_update_freq(-1),
     m_it_counter(0),
     m_first_call(true),
     m_first_stage(true)
{
   ParmParse ppcls( a_ppcls_str.c_str() );
   parseParameters( ppcls );
   if (m_verbosity>0) {
      printParameters();
   }

   m_fixed_cls_freq = (m_cls_freq<0.0) ? false : true ;
   
}


BGK::~BGK()
{
}


void BGK::evalClsRHS(KineticSpeciesPtrVect&       a_rhs,
                     const KineticSpeciesPtrVect& a_soln,
                     const int                    a_species,
                     const int                    a_species_bkgr,
                     const Real                   a_time )
{
   
   /*
      Computes the BGK collision operator model:
      C[JBF] = -nu * (JBF  - coeff[0] * JBFref
                           - coeff[1] * vpar * JBFref
                           - coeff[2] * energy * JBFref )
  
      energy = (mass*v_par^2 + mu*B)/3.0
      (N.B. we use 1/3 factor, so that we could use exisitng fast calculation of pressure moment)
    
      Fref is a maxwellian function with the same n,Vpar, and T
      as the current solution. Coefficients (coeff[0,1,2] are computed from
      the conditions the the operator numerically conserves two moments
      (i.e.,density and parallel momentum) or three moments (i.e.,density
      and parallel momentum and energy).
    */
   
   // Get solution distribution function (J*Bstar_par*dfn) for the current species
   const KineticSpecies& soln_species( *(a_soln[a_species]) );
   const LevelData<FArrayBox>& soln_dfn( soln_species.distributionFunction() );

   const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();
   const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
   
   const DisjointBoxLayout& grids = soln_dfn.disjointBoxLayout();
   
   // Create Maxwellian reference function
   if (m_first_call) {
      m_ref_species = soln_species.clone( IntVect::Zero, false );
      m_coeff_inj.define(grids, 3, IntVect::Zero);
   }
   
   // Compute reference maxwellian dfn, and its moments (e.g., m_density)
   computeReferenceSolution(*m_ref_species, soln_species, a_time);
      
   if (!m_fixed_cls_freq) {
         
      if (!m_sc_cls_freq.isDefined()) {
         m_sc_cls_freq.define(grids, 1, IntVect::Zero);
         m_density_inj.define(grids, 1, IntVect::Zero);
         m_temperature_inj.define(grids, 1, IntVect::Zero);
      }
         
      phase_geom.injectConfigurationToPhase( m_density, m_density_inj );
      phase_geom.injectConfigurationToPhase( m_temperature, m_temperature_inj );
         
      computeSelfConsistFreq(m_sc_cls_freq,
                             m_density_inj,
                             m_temperature_inj,
                             soln_species.mass(),
                             soln_species.charge());
   }
      
   CFG::LevelData<CFG::FArrayBox> coeff_cfg( mag_geom.grids(), 3, CFG::IntVect::Zero );
      
   if (m_conserve_two_moments) {
      computeTwoConservationCoeff(coeff_cfg);
   }
   else if (m_conserve_three_moments) {
      computeThreeConservationCoeff(coeff_cfg, soln_species);
   }
   else {
      for (CFG::DataIterator dit(coeff_cfg.dataIterator()); dit.ok(); ++dit) {
         coeff_cfg[dit].setVal(1.,0);
         coeff_cfg[dit].setVal(0.,1);
         coeff_cfg[dit].setVal(0.,2);
      }
   }

   // Inject cfg coefficient into the phase space
   phase_geom.injectConfigurationToPhase( coeff_cfg, m_coeff_inj);
   
   // Compute collision operator
   KineticSpecies& rhs_species( *(a_rhs[a_species]) );
   LevelData<FArrayBox>& rhs_dfn( rhs_species.distributionFunction() );
   
   LevelData<FArrayBox>& ref_dfn( m_ref_species->distributionFunction() );
   const LevelData<FArrayBox>& B_injected( phase_geom.getBFieldMagnitude() );
   
   double mass = soln_species.mass();
   
   DataIterator dit(grids.dataIterator());
   for (dit.begin(); dit.ok(); ++dit) {

      FArrayBox collision_op(grids[dit],1);
      const FArrayBox& this_ref_dfn = ref_dfn[dit];
      const FArrayBox& this_soln_dfn = soln_dfn[dit];
      const FArrayBox& this_coeff = m_coeff_inj[dit];
      const FArrayBox& this_B = B_injected[dit];
      const FArrayBox& this_velocity = phase_geom.getVelocityRealCoords()[dit];

      FORT_COMPUTE_BGK_OP(CHF_BOX(collision_op.box()),
                          CHF_FRA1(collision_op,0),
                          CHF_CONST_FRA1(this_ref_dfn,0),
                          CHF_CONST_FRA1(this_soln_dfn,0),
                          CHF_CONST_FRA1(this_B,0),
                          CHF_CONST_FRA(this_coeff),
                          CHF_CONST_FRA(this_velocity),
                          CHF_CONST_REAL(mass));
      
      if (m_fixed_cls_freq) collision_op.mult( m_cls_freq );
      else  collision_op.mult( m_sc_cls_freq[dit] );
      
      rhs_dfn[dit].plus( collision_op );
   }
      
   m_first_call = false;
   m_first_stage = false;
}


void BGK::computeReferenceSolution(KineticSpecies& a_ref_species,
                                   const KineticSpecies& a_soln_species,
                                   const Real a_time)
{
   const PhaseGeom& phase_geom = a_soln_species.phaseSpaceGeometry();
   const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
   
   LevelData<FArrayBox>& ref_maxw( a_ref_species.distributionFunction() );
   
   if (!m_density.isDefined()) {
      m_density.define( mag_geom.grids(), 1, CFG::IntVect::Zero );
      m_parallel_part_flux.define( mag_geom.grids(), 1, CFG::IntVect::Zero );
      m_temperature.define( mag_geom.grids(), 1, CFG::IntVect::Zero );
   }
   
   // Compute solution moments (n, nVpar, T)
   a_soln_species.numberDensity( m_density );
   a_soln_species.parallelParticleFlux( m_parallel_part_flux );
  
   CFG::LevelData<CFG::FArrayBox> parallel_velocity( mag_geom.grids(), 1, CFG::IntVect::Zero );
   for (CFG::DataIterator dit(m_density.dataIterator()); dit.ok(); ++dit) {
      parallel_velocity[dit].copy(m_parallel_part_flux[dit]);
      parallel_velocity[dit].divide(m_density[dit]);
   }
   CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.grids(), 1, CFG::IntVect::Zero );
   a_soln_species.pressure(pressure, parallel_velocity);
  
   for (CFG::DataIterator dit(m_density.dataIterator()); dit.ok(); ++dit) {
      m_temperature[dit].copy(pressure[dit]);
      m_temperature[dit].divide(m_density[dit]);
   }
  
   // Construct reference maxwellian (we don't need to multiply by J, because it is already contained
   // in the density moment we evalulated on the soln_comp
   MaxwellianKernel<FArrayBox> maxwellian(m_density,m_temperature,parallel_velocity);
   maxwellian.eval(ref_maxw,a_soln_species);
   
}


void BGK::computeTwoConservationCoeff(CFG::LevelData<CFG::FArrayBox>& a_coeff)
{
   /*
      a_coeff are obtained from the following system:
    
      n - coeff[0] * n_ref - coeff[1] * (nVpar)_ref = 0
      nVpar - acoeff[0] * (nVpar)_ref - coeff[1] * Ppar_ref = 0
    
      coeff[0] = (n - coeff[1] * (nVpar)_ref)/n_ref
      coeff[1] = -(nVpar - n/n_ref * (nVpar)_ref)/((nVpar)_ref^2/n_ref - Ppar_ref)
      
    */
   
   const PhaseGeom& phase_geom =  m_ref_species->phaseSpaceGeometry();

   const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
   const CFG::DisjointBoxLayout& grids_cfg( mag_geom.grids() );
   
   // Compute moments of the reference distribution
   CFG::LevelData<CFG::FArrayBox> ref_density( mag_geom.grids(), 1, CFG::IntVect::Zero );
   m_ref_species->numberDensity( ref_density );
   
   CFG::LevelData<CFG::FArrayBox> ref_parallel_part_flux( mag_geom.grids(), 1, CFG::IntVect::Zero );
   m_ref_species->parallelParticleFlux( ref_parallel_part_flux );
   
   CFG::LevelData<CFG::FArrayBox> ref_parallel_pressure( mag_geom.grids(), 1, CFG::IntVect::Zero );
   CFG::LevelData<CFG::FArrayBox> zero( mag_geom.grids(), 1, CFG::IntVect::Zero );
   for (CFG::DataIterator dit(zero.dataIterator()); dit.ok(); ++dit) {
      zero[dit].setVal(0.);
   }
   m_ref_species->parallelPressure(ref_parallel_pressure, zero);
   const double mass(m_ref_species->mass() );
   for (CFG::DataIterator dit(ref_parallel_pressure.dataIterator()); dit.ok(); ++dit) {
      ref_parallel_pressure[dit].divide(mass);
   }
   
   // Compute coefficients
   for (CFG::DataIterator dit(a_coeff.dataIterator()); dit.ok(); ++dit) {

      const CFG::Box& box = grids_cfg[dit];

      CFG::FArrayBox numerator(box, 1);
      CFG::FArrayBox denominator(box, 1);
      
      //Compute coeff[1]
      numerator.copy(ref_parallel_part_flux[dit]);
      numerator.mult(m_density[dit]);
      numerator.divide(ref_density[dit]);
      numerator.mult(-1.0);
      numerator.plus(m_parallel_part_flux[dit]);
      
      denominator.copy(ref_parallel_part_flux[dit]);
      denominator.mult(ref_parallel_part_flux[dit]);
      denominator.divide(ref_density[dit]);
      denominator.minus(ref_parallel_pressure[dit]);
      
      a_coeff[dit].copy(numerator, 0, 1, 1);
      a_coeff[dit].divide(denominator, 0, 1, 1);
      a_coeff[dit].mult(-1.0, 1);
      
      
      //Compute coeff[0]
      CFG::FArrayBox tmp(box, 1);
      tmp.copy(a_coeff[dit], 1, 0, 1);
      tmp.mult(-1.0);
      
      a_coeff[dit].copy(ref_parallel_part_flux[dit], 0, 0, 1);
      a_coeff[dit].mult(tmp,0,0,1);
      a_coeff[dit].plus(m_density[dit], 0, 0, 1);
      a_coeff[dit].divide(ref_density[dit], 0, 0, 1);
      
      //Set coeff[2] to zero
      a_coeff[dit].setVal(0, 2);
   }
}

void BGK::computeThreeConservationCoeff(CFG::LevelData<CFG::FArrayBox>& a_coeff,
                                        const KineticSpecies& a_soln_species)
{
   /*
      C[JBF] = -nu * (JBF - a_coeff[0] * JBFref
                          - a_coeff[1] * vpar * JBFref
                          - a_coeff[2] * energy * JBFref )
    
      energy = (mass*v_par^2 + mu*B)/3.0
      (N.B. we use 1/3 factor, so that we could use exisitng fast calculation of pressure moment)
    
      a_coefficients are obtained from the following system
      (N.B. the numerical system is multipled by J)  :
    
      a_0 * int(Fref,dv) + a_1*int(vpar*Fref,dv) + a_2*int(energy * Fref) = n
      a_0 * int(vpar*Fref,dv) + a_1*int(vpar^2*Fref,dv) + a_2*int(vpar*energy * Fref) = nVpar
      a_0 * int(energy*Fref,dv) + a_1*int(vpar*energy*Fref,dv) + a_2*int(energy^2 * Fref) = PressureMoment(Vshift=0)
        
    */

   CH_assert(a_coeff.nComp()==3);
   
   const LevelData<FArrayBox>& ref_dfn( m_ref_species->distributionFunction() );
   const DisjointBoxLayout& grids = ref_dfn.disjointBoxLayout();
   
   const PhaseGeom& phase_geom =  m_ref_species->phaseSpaceGeometry();
   const LevelData<FArrayBox>& B_injected( phase_geom.getBFieldMagnitude() );
   
   const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
   const CFG::DisjointBoxLayout& grids_cfg( mag_geom.grids() );

   double mass = a_soln_species.mass();
   
   if (!m_kernels.isDefined()) m_kernels.define(grids, 9, IntVect::Zero);

   DataIterator dit(grids.dataIterator());
   for (dit.begin(); dit.ok(); ++dit) {

      FArrayBox& this_kernels = m_kernels[dit];
      const FArrayBox& this_ref_dfn = ref_dfn[dit];
      const FArrayBox& this_B = B_injected[dit];
      const FArrayBox& this_velocity = phase_geom.getVelocityRealCoords()[dit];

      FORT_COMPUTE_BGK_CONS_KERNELS(CHF_BOX(this_kernels.box()),
                                    CHF_FRA(this_kernels),
                                    CHF_CONST_FRA1(this_ref_dfn,0),
                                    CHF_CONST_FRA1(this_B,0),
                                    CHF_CONST_FRA(this_velocity),
                                    CHF_CONST_REAL(mass));
   }
   
      
   CFG::LevelData<CFG::FArrayBox> matrix_coeff( mag_geom.grids(), 9, CFG::IntVect::Zero );

   MomentOp& moment_op = MomentOp::instance();
   moment_op.compute(matrix_coeff, *m_ref_species, m_kernels, DensityKernel<FArrayBox>());
   
   CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.grids(), 1, CFG::IntVect::Zero );
   CFG::LevelData<CFG::FArrayBox> zero( mag_geom.grids(), 1, CFG::IntVect::Zero );
   for (CFG::DataIterator dit(zero.dataIterator()); dit.ok(); ++dit) {
      zero[dit].setVal(0.);
   }
   a_soln_species.pressure(pressure, zero);
   
   // Solve trilinear system
   CFG::LevelData<CFG::FArrayBox> rhs( mag_geom.grids(), 3, CFG::IntVect::Zero );
   for (CFG::DataIterator dit(rhs.dataIterator()); dit.ok(); ++dit) {
      rhs[dit].copy(m_density[dit],0,0,1);
      rhs[dit].copy(m_parallel_part_flux[dit],0,1,1);
      rhs[dit].copy(pressure[dit],0,2,1);
   }
   
   CFG::SpaceUtils::solveTrilinearSystem(a_coeff,
                                         matrix_coeff,
                                         rhs);

}

void BGK::computeSelfConsistFreq(LevelData<FArrayBox>&         a_cls_freq,
                                 const LevelData<FArrayBox>&   a_density,
                                 const LevelData<FArrayBox>&   a_temperature,
                                 const double                  a_mass,
                                 const double                  a_charge      ) const

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
       const FArrayBox& this_T = a_temperature[dit];

       FORT_COMPUTE_SC_CLS_FREQ(CHF_BOX(this_cls_freq.box()),
                                CHF_FRA1(this_cls_freq,0),
                                CHF_CONST_FRA1(this_n,0),
                                CHF_CONST_FRA1(this_T,0));
  
       this_cls_freq.mult(cls_norm);

    }

}

inline
void BGK::parseParameters( ParmParse& a_ppcls )
{
   a_ppcls.query( "time_implicit", m_time_implicit);
   a_ppcls.query( "cls_freq", m_cls_freq );
   a_ppcls.query( "conserve_two_moments", m_conserve_two_moments );
   a_ppcls.query( "conserve_three_moments", m_conserve_three_moments );
   
   if (m_conserve_two_moments && m_conserve_three_moments) {
      MayDay::Error("BGK:: choose whether to conserve two or three moments ");
   }
   
   a_ppcls.query( "update_frequency", m_update_freq);

}


inline
void BGK::printParameters()
{
   if (procID()==0) {
      std::cout << "BGK collisions parameters:" << std::endl;
      std::cout << "  cls_freq  =  " << m_cls_freq
                << ", conserve_two_moments = " << m_conserve_two_moments
                << ", conserve_three_moments = " << m_conserve_three_moments
                << ", update_frequency = " << m_update_freq
                << ", implicit in time = " << m_time_implicit
                << std::endl;
   }
}

Real BGK::collisionFrequency()
{
  /* implemented for fixed collision frequency for now */
  CH_assert(m_fixed_cls_freq);
  return m_cls_freq;
}

void BGK::preTimeStep(const KineticSpeciesPtrVect& a_soln_mapped,
                      const int a_species,
                      const Real a_time,
                      const KineticSpeciesPtrVect& a_soln_physical )
{
   m_it_counter+=1;
   m_first_stage = true;
}

Real BGK::computeDtExplicitTI(const KineticSpeciesPtrVect& soln, const int a_idx)
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

Real BGK::computeDtImExTI(const KineticSpeciesPtrVect& soln, const int a_idx)
{
  if (m_time_implicit) {
    return DBL_MAX;
  } else {
    return computeDtExplicitTI(soln, a_idx);
  }
}

Real BGK::computeTimeScale(const KineticSpeciesPtrVect& soln, const int a_idx)
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
