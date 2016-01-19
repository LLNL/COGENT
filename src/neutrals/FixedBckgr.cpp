#include <math.h>

#include "FixedBckgr.H"

#include "Directions.H"
#include "PhaseGeom.H"
#include "PhaseBlockCoordSys.H"
#include "NeutralsF_F.H"
#include "KineticFunctionLibrary.H"
#include "ConstFact.H"
#include "inspect.H"

#include "NamespaceHeader.H" 


FixedBckgr::FixedBckgr( ParmParse& a_ppntr, const int a_verbosity )
   : m_verbosity(a_verbosity),
     m_ntr_freq(-1.0),
     m_first_step(true)
{
   parseParameters( a_ppntr );
   if (m_verbosity>0) {
      printParameters();
   }

   m_fixed_ntr_freq = (m_ntr_freq<0.0) ? false : true ;
}


FixedBckgr::~FixedBckgr()
{
}


void FixedBckgr::evalNtrRHS( KineticSpeciesPtrVect&  a_rhs,
                        const KineticSpeciesPtrVect& a_soln,
                        const int                    a_species,
                        const Real                   a_time )
// NB: a_soln is on the computational grid
{
   // Get solution distribution function (J*Bstar_par*dfn) for the current species
   const KineticSpecies& soln_species( *(a_soln[a_species]) );
   const LevelData<FArrayBox>& soln_dfn( soln_species.distributionFunction() );

   // Create reference (J*Bstar_par*dfn_init) distribution
   KineticSpeciesPtr ref_species( soln_species.clone( IntVect::Zero, false ) );
   m_ref_func->assign( *ref_species, a_time );
   const LevelData<FArrayBox>& init_dfn( ref_species->distributionFunction() );

   //Create reference temperature distribution
   const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();

   if (m_first_step) {
      const CFG::MultiBlockLevelGeom& mag_geom( phase_geom.magGeom() );
      CFG::LevelData<CFG::FArrayBox> ref_temperature( mag_geom.grids(), 1, CFG::IntVect::Zero );
      m_ref_temp->assign( ref_temperature, mag_geom, a_time);
      phase_geom.injectConfigurationToPhase( ref_temperature, m_temperature);
      
      //Create reference density distribution distribution and self-consistent ntr_freq
      if (!m_fixed_ntr_freq){
         CFG::LevelData<CFG::FArrayBox> ref_density( mag_geom.grids(), 1, CFG::IntVect::Zero );
         m_ref_dens->assign( ref_density, mag_geom, a_time);
         
         LevelData<FArrayBox> inj_density;
         phase_geom.injectConfigurationToPhase( ref_density, inj_density);
         
         m_sc_ntr_freq.define(soln_dfn.disjointBoxLayout(), 1, IntVect::Zero);
         computeSelfConsistFreq(m_sc_ntr_freq, inj_density, soln_species.mass(), soln_species.charge());
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
   
   // Calculate test-particle (TP) neutral-model RHS
   LevelData<FArrayBox> tp_rhs_coll( grids, n_comp, IntVect::Zero );
   testPartCollRHS( tp_rhs_coll, delta_F );
    
   // Add test-particle neutral-model RHS to Vlasov RHS
   KineticSpecies& rhs_species( *(a_rhs[a_species]) );
   LevelData<FArrayBox>& rhs_dfn( rhs_species.distributionFunction() );
   for (DataIterator rdit(soln_dfn.dataIterator()); rdit.ok(); ++rdit) {
      if (m_fixed_ntr_freq) {tp_rhs_coll[rdit].mult( m_ntr_freq );}
      else {tp_rhs_coll[rdit].mult( m_sc_ntr_freq[rdit] );}
      rhs_dfn[rdit].plus( tp_rhs_coll[rdit] );
   }  

   m_first_step = false;
}



void FixedBckgr::testPartCollRHS( LevelData<FArrayBox>& a_rhs_coll,
                             const LevelData<FArrayBox>& a_delta_F ) const
{
   const Real fact( -1.0 );
   for (DataIterator rdit( a_rhs_coll.dataIterator() ); rdit.ok(); ++rdit) {
      a_rhs_coll[rdit].copy(a_delta_F[rdit]);
      a_rhs_coll[rdit].mult(fact);
   }
}



void FixedBckgr::computeSelfConsistFreq(LevelData<FArrayBox>& a_ntr_freq,
                                  const LevelData<FArrayBox>& a_density,
                                  const double                a_mass,
                                  const double                a_charge      ) const

{
    //This is just an example from the collision models
    
    
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
    double ntr_norm = coeff * N * pow(ech, 2) * pow(a_charge, 4) * L 
                     / ( sqrt(a_mass) *  pow(T, 2)) * Coulomb_Lg;

    const DisjointBoxLayout& grids = a_ntr_freq.disjointBoxLayout();
    DataIterator dit(grids.dataIterator());
    for (dit.begin(); dit.ok(); ++dit) { 

      FArrayBox& this_ntr_freq = a_ntr_freq[dit];
      const FArrayBox& this_n = a_density[dit];
      const FArrayBox& this_T = m_temperature[dit];

      FORT_COMPUTE_SC_NTR_FREQ(CHF_BOX(this_ntr_freq.box()),
                               CHF_FRA1(this_ntr_freq,0),
                               CHF_CONST_FRA1(this_n,0),
                               CHF_CONST_FRA1(this_T,0));
  
    this_ntr_freq.mult(ntr_norm); 

    }

}

inline
void FixedBckgr::parseParameters( ParmParse& a_ppntr )
{
   a_ppntr.query( "ntr_freq", m_ntr_freq );

   KineticFunctionLibrary* library = KineticFunctionLibrary::getInstance();
   std::string function_name;
   a_ppntr.query( "ref_function", function_name );
   m_ref_func = library->find( function_name );

   CFG::GridFunctionLibrary* grid_library = CFG::GridFunctionLibrary::getInstance();
   std::string grid_function_name;

   a_ppntr.query( "ref_temperature", grid_function_name );
   m_ref_temp = grid_library->find( grid_function_name );

   a_ppntr.query( "ref_density", grid_function_name );
   m_ref_dens = grid_library->find( grid_function_name );


}


inline
void FixedBckgr::printParameters()
{
   if (procID()==0) {
      std::cout << "FixedBckgr neutrals parameters:" << std::endl;
      std::cout << "  ntr_freq  =  " << m_ntr_freq;
      std::cout << "  Reference Function:" << std::endl;
      m_ref_func->printParameters();
      std::cout << "  Reference Temperature:" << std::endl;
      m_ref_temp->printParameters();
      std::cout << "  Reference Density:" << std::endl;
      m_ref_dens->printParameters();

   }
}


Real FixedBckgr::computeDt(const KineticSpeciesPtrVect& soln)
{
   return DBL_MAX;
}

#include "NamespaceFooter.H"
