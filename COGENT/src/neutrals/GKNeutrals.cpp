#include "GKNeutrals.H"

#include "NTRInterface.H"
#include "FixedBckgr.H"
#include "PrescribedSources.H"
#include "FluidNeutrals.H"
#include "NullNTR.H"

#include <float.h>
#include <sstream>

#include "NamespaceHeader.H"

GKNeutrals::GKNeutrals( const int a_verbose )
   : m_verbose(a_verbose)
{
   bool more_kinetic_species(true);
   int count(0);
   while (more_kinetic_species) {

      // look for data specifying another kinetic species
      std::stringstream s;
      s << "kinetic_species." << count+1;
      ParmParse ppspecies( s.str().c_str() );

      std::string species_name("Invalid");
      if (ppspecies.contains("name")) {
         ppspecies.get("name", species_name);
         
         std::string ntr_type("None");
         NTRInterface* ntr(NULL);
         
         if (ppspecies.contains( "ntr" )) {
            ppspecies.get( "ntr", ntr_type );
            const std::string prefix( "NTR." + species_name );
            ParmParse ppntr( prefix.c_str() );
            
            if (ntr_type == "FixedBckgr") {
               ntr = new FixedBckgr( ppntr, m_verbose );
            }

           if (ntr_type == "PrescribedSources") {
              ntr = new PrescribedSources( ppntr, m_verbose );
           }
           
            if (ntr_type == "FluidNeutrals") {
               ntr = new FluidNeutrals( ppntr, m_verbose );
            }
         }
         else {
            if (procID()==0) {
              cout << "Unrecognized neutral model; setting model to NULL" << endl;
            }
            ntr = new NullNTR();
         }
         
         m_neutral_model.push_back( ntr );
         typedef std::map<std::string,int>::value_type valType;
         m_species_map.insert( valType( species_name, count ) );
	 m_neutral_model_name.insert( valType( ntr_type, count ) );
         count++;
      }
      else {
         more_kinetic_species = false;
      }
   }
}


GKNeutrals::~GKNeutrals()
{
   for (int i(0); i<m_neutral_model.size(); i++ ) {
      delete m_neutral_model[i];
   }
}


NTRInterface& GKNeutrals::neutralModel( const std::string& a_name )
{
   typedef std::map<std::string,int>::iterator mapIterator;
   const mapIterator it( m_species_map.find( a_name ) );
   CH_assert(it!=m_species_map.end());
   const int index((*it).second);
   return *(m_neutral_model[index]);
}


void GKNeutrals::accumulateRHS( KineticSpeciesPtrVect&            a_rhs,
                                const KineticSpeciesPtrVect&      a_kinetic_species_phys,
                                const CFG::FluidSpeciesPtrVect&   a_fluid_species_phys,
                                const Real                        a_time )
{
   for (int species(0); species<a_rhs.size(); species++) {
      KineticSpecies& rhs_species( *(a_rhs[species]) );
      const std::string species_name( rhs_species.name() );
      NTRInterface& NTR( neutralModel( species_name ) );
      NTR.evalNtrRHS( rhs_species,
                      a_kinetic_species_phys,
                      a_fluid_species_phys,
                      species,
                      a_time );
   }
}


Real GKNeutrals::computeDtExplicitTI( const KineticSpeciesPtrVect& soln )
{
  Real dt(DBL_MAX);
  int count(0);
  std::map<std::string,int>::iterator it;
  for (it=m_neutral_model_name.begin(); it!=m_neutral_model_name.end(); ++it) {
    Real tmp = m_neutral_model[it->second]->computeDtExplicitTI(soln,it->second);
    dt = (tmp < dt ? tmp : dt);
    count++;
  }
  return (count ? dt : -1);

}

Real GKNeutrals::computeDtImExTI( const KineticSpeciesPtrVect& soln )
{
  Real dt(DBL_MAX);
  int count(0);
  std::map<std::string,int>::iterator it;
  for (it=m_neutral_model_name.begin(); it!=m_neutral_model_name.end(); ++it) {
    Real tmp = m_neutral_model[it->second]->computeDtImExTI(soln,it->second);
    dt = (tmp < dt ? tmp : dt);
    count++;
  }
  return (count ? dt : -1);

}

Real GKNeutrals::computeTimeScale( const KineticSpeciesPtrVect& soln )
{
  std::map<std::string,int>::iterator it;
  Real scale = DBL_MAX;
  int count = 0;
  for (it=m_neutral_model_name.begin(); it!=m_neutral_model_name.end(); ++it) {
    Real tmp = m_neutral_model[it->second]->TimeScale(soln,it->second);
    scale = (tmp < scale ? tmp : scale);
    count++;
  }
  return (count ? scale : -1);

}



#include "NamespaceFooter.H"
