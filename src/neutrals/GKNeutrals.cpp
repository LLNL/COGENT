#include "GKNeutrals.H"

#include "NTRInterface.H"
#include "FixedBckgr.H"
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


void GKNeutrals::accumulateRHS( KineticSpeciesPtrVect&       a_rhs,
                                  const KineticSpeciesPtrVect& a_soln,
                                  const Real                   a_time )
{
   for (int species(0); species<a_rhs.size(); species++) {
      KineticSpecies& rhs_species( *(a_rhs[species]) );
      const std::string species_name( rhs_species.name() );
      NTRInterface& NTR( neutralModel( species_name ) );
      NTR.evalNtrRHS( a_rhs, a_soln, species, a_time );
   }
}


Real GKNeutrals::computeDt( const KineticSpeciesPtrVect& soln )
{
   return DBL_MAX;
}



#include "NamespaceFooter.H"
