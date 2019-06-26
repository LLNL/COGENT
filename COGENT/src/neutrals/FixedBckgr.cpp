#include <math.h>

#include "FixedBckgr.H"

#include "Directions.H"
#include "PhaseGeom.H"
#include "PhaseBlockCoordSys.H"
#include "NeutralsF_F.H"
#include "KineticFunctionLibrary.H"
#include "ConstFact.H"
#include "Kernels.H"
#include "inspect.H"

#include "NamespaceHeader.H" 


FixedBckgr::FixedBckgr( ParmParse& a_ppntr, const int a_verbosity )
   : m_verbosity(a_verbosity),
     m_fixed_source_dfn(false),
     m_neutr_vel(NULL),
     m_neutr_temp(NULL)
{
   parseParameters( a_ppntr );

   if (m_verbosity>0) {
      printParameters();
   }
   m_first_call = true;
}


FixedBckgr::~FixedBckgr()
{
}


void FixedBckgr::evalNtrRHS(KineticSpecies&              a_rhs_species,
                            const KineticSpeciesPtrVect& a_soln,
                            const int                    a_species,
                            const Real                   a_time )
{
   /*
    Evaluates the ionization and charge-exchange contribution to the RHS as
    if m_iz_fixed_source = false:
       d(J*fB)/dt = J * (<sigmaV_ionization> * neutr_dens * fB
    if m_iz_fixed_source = true:
      d(J*fB)/dt = J * (<sigmaV_ionization> * neutr_dens * fB_neutral
 
    if m_chx_model_friction = false:
                  + <sigmaV_chx> (n_ion * fB_neutral - n_neutral * fB)
    if m_chx_model_friction = true:
                  + <sigmaV_chx> * n_neutral * vpar * upar * mass * fB_unshifted_Maxw / Ti)
    */
   
   // Get physical solution distribution function (Bstar_par*dfn) for the current species
   const KineticSpecies& soln_species( *(a_soln[a_species]) );
   const LevelData<FArrayBox>& soln_dfn( soln_species.distributionFunction() );

   //Get geometry
   const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();
   const CFG::MultiBlockLevelGeom& mag_geom( phase_geom.magGeom() );

   //Get neutral profiles
   
   if (m_first_call) {
      
      double dens_norm, vel_norm, temp_norm;
      computeChxNormalization(m_ionization_norm, m_chx_norm, dens_norm, vel_norm, temp_norm, m_SI_input);
      
      m_neutral_density_cfg.define( mag_geom.grids(), 1, CFG::IntVect::Zero );
      m_neutr_dens->assign( m_neutral_density_cfg, mag_geom, a_time);
      for (CFG::DataIterator dit(mag_geom.grids()); dit.ok(); ++dit) {
         m_neutral_density_cfg[dit].divide(dens_norm);
      }
      phase_geom.injectConfigurationToPhase( m_neutral_density_cfg, m_neutral_density);

      CFG::LevelData<CFG::FArrayBox> ionization_rate( mag_geom.grids(), 1, CFG::IntVect::Zero );
      m_ioniz_rate->assign( ionization_rate, mag_geom, a_time);
      phase_geom.injectConfigurationToPhase( ionization_rate, m_ionization_rate);

      if (m_include_chx) {
         m_neutral_velocity_cfg.define( mag_geom.grids(), 1, CFG::IntVect::Zero );
         m_neutr_vel->assign( m_neutral_velocity_cfg, mag_geom, a_time);
         for (CFG::DataIterator dit(mag_geom.grids()); dit.ok(); ++dit) {
            m_neutral_velocity_cfg[dit].divide(vel_norm);
         }
      
        m_neutral_dfn.define( phase_geom.gridsFull(), 1, IntVect::Zero );
      
        if (m_neutr_temp != NULL) {
           CFG::LevelData<CFG::FArrayBox> neutral_temperature_cfg( mag_geom.grids(), 1, CFG::IntVect::Zero );
           m_neutr_temp->assign( neutral_temperature_cfg, mag_geom, a_time);

           for (CFG::DataIterator dit(mag_geom.grids()); dit.ok(); ++dit) {
              neutral_temperature_cfg[dit].divide(temp_norm);
           }
         
           MaxwellianKernel maxwellian(m_neutral_density_cfg, neutral_temperature_cfg, m_neutral_velocity_cfg);
           maxwellian.eval(m_neutral_dfn,soln_species);

        }
      
      }      
   }
  
 
   //Calculate ionization RHS
   LevelData<FArrayBox> ionization_rhs( soln_dfn.getBoxes(), soln_dfn.nComp(), IntVect::Zero );

   if (m_fixed_source_dfn) {
      KineticSpeciesPtr ref_species( soln_species.clone( IntVect::Unit, false ) );
      m_neutr_func->assign( *ref_species, a_time );
      LevelData<FArrayBox>& neutr_dfn( ref_species->distributionFunction() );
      phase_geom.divideJonValid(neutr_dfn);    
      computeIonizationRHS( ionization_rhs, neutr_dfn, m_neutral_density );
    }   
    else
      computeIonizationRHS( ionization_rhs,  soln_dfn, m_neutral_density );

   //Calculate charge-exchange RHS
   LevelData<FArrayBox> chx_rhs( soln_dfn.getBoxes(), soln_dfn.nComp(), IntVect::Zero );
   if (m_include_chx) {
      
      //Get the size of ghost vector
      IntVect ghost = soln_dfn.ghostVect();
      CFG::IntVect ghost_cfg;
      ghost_cfg[0] = ghost[0];
      ghost_cfg[1] = ghost[1];
     
      //Get ion species density
      CFG::LevelData<CFG::FArrayBox> ion_density_cfg( mag_geom.grids(), 1, ghost_cfg );
      soln_species.numberDensity( ion_density_cfg );
      LevelData<FArrayBox> ion_density;
      phase_geom.injectConfigurationToPhase(ion_density_cfg, ion_density);

      //Get ion species temperature
      CFG::LevelData<CFG::FArrayBox> ion_parallelVel_cfg( mag_geom.grids(), 1, ghost_cfg );
      soln_species.ParallelMomentum( ion_parallelVel_cfg );
      for (CFG::DataIterator dit(ion_density_cfg.dataIterator()); dit.ok(); ++dit) {
         ion_parallelVel_cfg[dit].divide(ion_density_cfg[dit]);
      }
      CFG::LevelData<CFG::FArrayBox> ion_temperature_cfg( mag_geom.grids(), 1, ghost_cfg );
      soln_species.pressureMoment(ion_temperature_cfg, ion_parallelVel_cfg);
      for (CFG::DataIterator dit(ion_density_cfg.dataIterator()); dit.ok(); ++dit) {
        ion_temperature_cfg[dit].divide(ion_density_cfg[dit]);
      }
      LevelData<FArrayBox> ion_temperature;
      phase_geom.injectConfigurationToPhase(ion_temperature_cfg, ion_temperature);

      if (!m_chx_model_friction) {
         //Compute neutral dfn using the ion temperature
         if (m_neutr_temp == NULL) {
            MaxwellianKernel maxwellian(m_neutral_density_cfg, ion_temperature_cfg, m_neutral_velocity_cfg);
            maxwellian.eval(m_neutral_dfn, soln_species);
         }

         computeChargeExchangeRHS( chx_rhs, soln_dfn, ion_density, ion_temperature, m_neutral_density, m_neutral_dfn );
      }

      if (m_chx_model_friction) {
         //Create unshifted ion Maxwellian and compute the parallel velocity difference
         CFG::LevelData<CFG::FArrayBox> zero_vel_cfg( mag_geom.grids(), 1, ghost_cfg );
         CFG::LevelData<CFG::FArrayBox> u_par_cfg( mag_geom.grids(), 1, ghost_cfg );
         for (CFG::DataIterator dit(ion_density_cfg.dataIterator()); dit.ok(); ++dit) {
            zero_vel_cfg[dit].setVal(0.0);
            u_par_cfg[dit].copy( m_neutral_velocity_cfg[dit]);
            u_par_cfg[dit].minus( ion_parallelVel_cfg[dit]);
         }
      
         LevelData<FArrayBox> ion_unshifted_maxw( phase_geom.gridsFull(), 1, IntVect::Zero );
         MaxwellianKernel maxwellian(ion_density_cfg, ion_temperature_cfg, zero_vel_cfg);
         maxwellian.eval(ion_unshifted_maxw, soln_species);
         
         LevelData<FArrayBox> u_par;
         phase_geom.injectConfigurationToPhase(u_par_cfg, u_par);
         
         double mass = soln_species.mass();
         computeModelChargeExchangeRHS( chx_rhs, m_neutral_density, ion_temperature, u_par, ion_unshifted_maxw, mass, phase_geom );

      }
   }
   
   //Multiply by J to convert to the computational space
   phase_geom.multJonValid(ionization_rhs);
   phase_geom.multJonValid(chx_rhs);

   //Add neutral-model RHS to the total RHS
   LevelData<FArrayBox>& rhs_dfn( a_rhs_species.distributionFunction() );
   for (DataIterator rdit(soln_dfn.dataIterator()); rdit.ok(); ++rdit) {
       rhs_dfn[rdit].plus( ionization_rhs[rdit] );
//       rhs_dfn[rdit].copy( ionization_rhs[rdit] ); // debug check to see effect of ionization term alone
       if (m_include_chx) {
          rhs_dfn[rdit].plus( chx_rhs[rdit] );
       }
   }  

   m_first_call = false;
}



void FixedBckgr::computeIonizationRHS( LevelData<FArrayBox>& a_rhs,
                                       const LevelData<FArrayBox>& a_soln_dfn,
                                       const LevelData<FArrayBox>& a_neutral_density) const
{
   for (DataIterator dit( a_rhs.dataIterator() ); dit.ok(); ++dit) {

      FORT_COMPUTE_IONIZATION(CHF_BOX(a_rhs[dit].box()),
                              CHF_FRA1(a_rhs[dit],0),
                              CHF_CONST_FRA1(a_soln_dfn[dit],0),
                              CHF_CONST_FRA1(a_neutral_density[dit],0),
                              CHF_CONST_FRA1(m_ionization_rate[dit],0));
      
      a_rhs[dit].mult(m_ionization_norm);
   }
}


void FixedBckgr::computeChargeExchangeRHS(LevelData<FArrayBox>& a_rhs,
                                          const LevelData<FArrayBox>& a_soln_dfn,
                                          const LevelData<FArrayBox>& a_ion_density,
                                          const LevelData<FArrayBox>& a_ion_temperature,
                                          const LevelData<FArrayBox>& a_neutral_density,
                                          const LevelData<FArrayBox>& a_neutral_dfn) const
{

   for (DataIterator dit( a_rhs.dataIterator() ); dit.ok(); ++dit) {
      
      FORT_COMPUTE_CHARGE_EXCHANGE(CHF_BOX(a_rhs[dit].box()),
                                   CHF_FRA1(a_rhs[dit],0),
                                   CHF_CONST_FRA1(a_soln_dfn[dit],0),
                                   CHF_CONST_FRA1(a_ion_density[dit],0),
                                   CHF_CONST_FRA1(a_ion_temperature[dit],0),
                                   CHF_CONST_FRA1(a_neutral_density[dit],0),
                                   CHF_CONST_FRA1(a_neutral_dfn[dit],0));
      
      
      a_rhs[dit].mult(m_chx_norm);
   }
}


void FixedBckgr::computeModelChargeExchangeRHS(LevelData<FArrayBox>&       a_rhs,
                                               const LevelData<FArrayBox>& a_neutral_density,
                                               const LevelData<FArrayBox>& a_ion_temperature,
                                               const LevelData<FArrayBox>& a_par_vel_shift,
                                               const LevelData<FArrayBox>& a_ion_unshifted_maxw,
                                               const double                a_mass,
                                               const PhaseGeom&            a_phase_geom) const
{
   
   const DisjointBoxLayout& grids = a_ion_unshifted_maxw.disjointBoxLayout();
   
   for (DataIterator dit( a_rhs.dataIterator() ); dit.ok(); ++dit) {
      
      // Get the physical velocity coordinates for this part of phase space
      const PhaseBlockCoordSys& block_coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
      FArrayBox velocityRealCoords(a_rhs[dit].box(), VEL_DIM);
      block_coord_sys.getVelocityRealCoords(velocityRealCoords);
      
      FORT_COMPUTE_MODEL_CHARGE_EXCHANGE(CHF_BOX(a_rhs[dit].box()),
                                         CHF_FRA1(a_rhs[dit],0),
                                         CHF_CONST_FRA1(a_neutral_density[dit],0),
                                         CHF_CONST_FRA1(a_ion_temperature[dit],0),
                                         CHF_CONST_FRA1(a_par_vel_shift[dit],0),
                                         CHF_CONST_FRA1(a_ion_unshifted_maxw[dit],0),
                                         CHF_CONST_FRA1(velocityRealCoords,0),
                                         CHF_CONST_REAL(a_mass));
      
      
      a_rhs[dit].mult(m_chx_norm);
   }
}


void FixedBckgr::computeChxNormalization(double&       a_ioniz_norm,
                                         double&       a_chx_norm,
                                         double&       a_dens_norm,
                                         double&       a_vel_norm,
                                         double&       a_temp_norm,
                                         const bool    a_SI_input) const
{

   //If a_SI_input = true: expect that the physical parameters are provided as
   //neutral_density[1/m^3], <sigmaV>[m^3/s], neutral_Vpar[m/s], neutral_Temp[eV]
   //If a_SI_input = false: expect all input quantities are provided in COGENT units
   
   //Universal constants (in CGS)
   double mp = 1.6726e-24;

   //Get normalization parameters (units)
   double N, T, L;
   ParmParse ppunits( "units" );
   ppunits.get("number_density",N);  //[m^{-3}]
   ppunits.get("temperature",T);     //[eV]
   ppunits.get("length",L);          //[m]

   double Tcgs = 1.602e-12 * T; //[erg]
   double Lcgs  = 1.0e2 * L;   //[cm]
   
   double time_norm = Lcgs / sqrt(Tcgs/mp); //[s]

   if (a_SI_input) {
      a_ioniz_norm = time_norm * N;
      a_chx_norm = time_norm * 1.98e-14 * sqrt(T) * N;
      a_dens_norm = N;
      a_vel_norm = sqrt(Tcgs/mp);
      a_temp_norm = T;
   }
   
   else {
      a_ioniz_norm = 1.0;
      a_chx_norm = time_norm * 1.98e-14 * sqrt(T) * N;
      a_dens_norm = 1.0;
      a_vel_norm = 1.0;
      a_temp_norm = 1.0;

   }
   
}

inline
void FixedBckgr::parseParameters( ParmParse& a_ppntr )
{

   if (a_ppntr.contains("include_charge_exchange")) {
      a_ppntr.get( "include_charge_exchange", m_include_chx );
   }
   else{
      m_include_chx = false;
   }

   if (a_ppntr.contains("charge_exchange_model_friction")) {
      a_ppntr.get( "charge_exchange_model_friction", m_chx_model_friction );
   }
   else{
      m_chx_model_friction = false;
   }

   if (a_ppntr.contains("SI_input")) {
      a_ppntr.get( "SI_input", m_SI_input );
   }
   else{
      m_SI_input = false;
   }

   KineticFunctionLibrary* library = KineticFunctionLibrary::getInstance();
   std::string function_name;
   if (a_ppntr.contains("neutral_phase_func")) {
       a_ppntr.get( "neutral_phase_func", function_name );
       m_neutr_func = library->find( function_name );
       m_fixed_source_dfn = true;
   }
   
   CFG::GridFunctionLibrary* grid_library = CFG::GridFunctionLibrary::getInstance();
   std::string grid_function_name;

   if (a_ppntr.contains("density")) {
      a_ppntr.get( "density", grid_function_name );
      m_neutr_dens = grid_library->find( grid_function_name );
   }
   else{
      MayDay::Error("FixedBckg:: neutral density must be specified ");
   }

   if (a_ppntr.contains("ionization_rate")) {
      a_ppntr.get( "ionization_rate", grid_function_name );
      m_ioniz_rate = grid_library->find( grid_function_name );
   }
   else{
      MayDay::Error("FixedBckg:: ionization rate must be specified ");
   }

   if (m_include_chx) {
   
      if (a_ppntr.contains("parallel_velocity")) {
         a_ppntr.get( "parallel_velocity", grid_function_name );
         m_neutr_vel = grid_library->find( grid_function_name );
      }
      else{
         MayDay::Error("FixedBckg:: neutral parallel velocity must be specified ");
      }

   
      if (a_ppntr.contains("temperature")) {
         a_ppntr.get( "temperature", grid_function_name );
         m_neutr_temp = grid_library->find( grid_function_name );
      }
      else{
         MayDay::Warning("FixedBckg:: ion temperature is used for neutral temperature");
      }
   }
  
}


inline
void FixedBckgr::printParameters()
{
   if (procID()==0) {
      std::cout << "FixedBckgr neutral parameters:" << std::endl;
      std::cout << "  SI_input  =  " << m_SI_input;
      std::cout << "  include_charge_exchange  =  " << m_include_chx;
      std::cout << "  fixed neutral distribution function = " << m_fixed_source_dfn;
      if (m_fixed_source_dfn) {
         std::cout << "  electron density:" << std::endl;
         m_neutr_dens->printParameters();
      }
      else {
         std::cout << "  Neutral Density:" << std::endl;
         m_neutr_dens->printParameters();

         if (m_include_chx) {

            std::cout << "  Neutral parallel velocity:" << std::endl;
            m_neutr_vel->printParameters();
         
            if (m_neutr_temp != NULL) {
               std::cout << "  Neutral Temperature:" << std::endl;
               m_neutr_temp->printParameters();
            }
         }
      }
   }
}

Real FixedBckgr::computeDtExplicitTI(const KineticSpeciesPtrVect& a_soln, const int a_species)
{
  return TimeScale(a_soln, a_species);
}

Real FixedBckgr::computeDtImExTI(const KineticSpeciesPtrVect& a_soln, const int a_species)
{
  return TimeScale(a_soln, a_species);
}

Real FixedBckgr::TimeScale(const KineticSpeciesPtrVect& a_soln, const int a_species)
{
  //Simple calculation that assumes normalized neutral density the order of unity,
  //and that charge-exchange process is the stiffest process  


  static bool first_call_tscale = true;
  if (first_call_tscale) {
    double dens_norm, vel_norm, temp_norm;
    computeChxNormalization(m_ionization_norm, m_chx_norm, dens_norm, vel_norm, temp_norm, m_SI_input);
  }

  Real time_scale;
  if (m_include_chx) {
    time_scale = 1.0/m_chx_norm;
  }
  else {
    time_scale = 1.0/m_ionization_norm;
  }

  first_call_tscale = false;

  return time_scale;
}

void FixedBckgr::diagnostics(const LevelData<FArrayBox>& a_rhs,
                             const KineticSpecies&       a_rhs_species,
			     const double                a_time) const
{

  //Get geometry                                                                                                                                                            
  const PhaseGeom& phase_geom = a_rhs_species.phaseSpaceGeometry();
  const CFG::MultiBlockLevelGeom& mag_geom( phase_geom.magGeom() );

  //Plot particle source
  CFG::LevelData<CFG::FArrayBox> particle_src( mag_geom.grids(), 1, CFG::IntVect::Zero );                                                                                 
  a_rhs_species.numberDensity( particle_src );                                                                                                                              
  phase_geom.plotConfigurationData( "particle_src", particle_src, a_time );                                                                                                 

  //Plot parallel momentum source                                                                                                                                         
  CFG::LevelData<CFG::FArrayBox> parMom_src( mag_geom.grids(), 1, CFG::IntVect::Zero );
  a_rhs_species.ParallelMomentum( parMom_src );
  phase_geom.plotConfigurationData( "parMom_src", parMom_src, a_time );

}

#include "NamespaceFooter.H"
