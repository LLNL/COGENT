#include "GKCollisions.H"

#include "CLSInterface.H"
#include "Krook.H"
#include "MyKrook.H"
#include "Lorentz.H"
#include "Linearized.H"
#include "FokkerPlanck.H"
#include "ConsDragDiff.H"
#include "NullCLS.H"

#include <float.h>
#include <sstream>

#include "NamespaceHeader.H"

GKCollisions::GKCollisions( const int a_verbose )
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
         
         std::string cls_type("None");
         CLSInterface* cls(NULL);
         
         if (ppspecies.contains( "cls" )) {
            ppspecies.get( "cls", cls_type );
            const std::string prefix( "CLS." + species_name );
            ParmParse ppcls( prefix.c_str() );
            
            if (cls_type == "Krook") {
               cls = new Krook( ppcls, m_verbose );
            }
            else if (cls_type == "MyKrook") {
               cls = new MyKrook( ppcls, m_verbose );
            }
            else if (cls_type == "Lorentz") {
               cls = new Lorentz( ppcls, m_verbose );
            }
            else if (cls_type == "Linearized") {
               cls = new Linearized( ppcls, m_verbose );
            }
            else if (cls_type == "FokkerPlanck") {
               cls = new FokkerPlanck( ppcls, m_verbose );
            }
            else if (cls_type == "ConsDragDiff") {
               cls = new ConsDragDiff( species_name, ppcls, m_verbose );
            }
         }
         else {
            if (procID()==0) {
              cout << "Unrecognized collision model; setting model to NULL" << endl;
            }
            cls = new NullCLS();
         }
         
         m_collision_model.push_back( cls );
         typedef std::map<std::string,int>::value_type valType;
         m_species_map.insert( valType( species_name, count ) );
         m_collision_model_name.insert( valType( cls_type, count ) );
         count++;
      }
      else {
         more_kinetic_species = false;
      }
   }
}


GKCollisions::~GKCollisions()
{
   for (int i(0); i<m_collision_model.size(); i++ ) {
      delete m_collision_model[i];
   }
}


CLSInterface& GKCollisions::collisionModel( const std::string& a_name )
{
   typedef std::map<std::string,int>::iterator mapIterator;
   const mapIterator it( m_species_map.find( a_name ) );
   CH_assert(it!=m_species_map.end());
   const int index((*it).second);
   return *(m_collision_model[index]);
}


void GKCollisions::accumulateRHS( KineticSpeciesPtrVect&       a_rhs,
                                  const KineticSpeciesPtrVect& a_soln,
                                  const Real                   a_time )
{
   for (int species(0); species<a_rhs.size(); species++) {
      KineticSpecies& rhs_species( *(a_rhs[species]) );
      const std::string species_name( rhs_species.name() );
      CLSInterface& CLS( collisionModel( species_name ) );
      CLS.evalClsRHS( a_rhs, a_soln, species, a_time );
   }
}


Real GKCollisions::computeDt( const KineticSpeciesPtrVect& soln )
{
   std::map<std::string,int>::iterator it;
   for (it=m_collision_model_name.begin(); it!=m_collision_model_name.end(); ++it) {
      if (it->first == "Krook") return ((Krook*)m_collision_model[it->second])->computeDt(soln);
      else                      return DBL_MAX;
   }
}

Real GKCollisions::computeTimeScale( const KineticSpeciesPtrVect& soln )
{
   std::map<std::string,int>::iterator it;
   for (it=m_collision_model_name.begin(); it!=m_collision_model_name.end(); ++it) {
      if (it->first == "Krook") return ((Krook*)m_collision_model[it->second])->TimeScale(soln);
      else                      return DBL_MAX;
   }
}

#include "NamespaceFooter.H"
