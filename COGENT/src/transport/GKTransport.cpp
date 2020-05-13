#include "GKTransport.H"

#include "TPMInterface.H"
#include "CONSTANTS.H"
#include "Anomalous.H"
#include "NullTPM.H"
#include "TransportF_F.H"

#include <float.h>
#include <sstream>

#include "NamespaceHeader.H"

GKTransport::GKTransport( const int a_verbose )
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

         std::string tpm_type("None");
         TPMInterface* tpm(NULL);

         if (ppspecies.contains( "tpm" )) {
            ppspecies.get( "tpm", tpm_type );
            const std::string prefix( "TPM." + species_name );
            ParmParse pptpm( prefix.c_str() );
            if (tpm_type == "Anomalous") {
               tpm = new Anomalous(species_name, pptpm, m_verbose );
            }
            else {
              if (procID()==0) {
                cout << "Unrecognized transport model for species " << species_name  << endl;
                MayDay::Error("Unrecognized transport model defined");
              }
              tpm = new NullTPM();
            }
         }
         else {
           if (procID()==0) {
                cout << "Transport model not set for species " << species_name << "; setting model to NULL" << endl;
              }
           tpm = new NullTPM();
         }

         m_transport_model.push_back( tpm );
         typedef std::map<std::string,int>::value_type valType;
         m_species_map.insert( valType( species_name, count ) );
         m_transport_model_name.insert( valType( tpm_type, count ) );
         if (!procID()) {
           cout << "Transport model for " << count << "\t" << species_name << ":\t";
           cout << tpm_type << "\n" ;
         }
         count++;
      }
      else {
         more_kinetic_species = false;
      }
   }
}


GKTransport::~GKTransport()
{
   for (int i(0); i<m_transport_model.size(); i++ ) {
      delete m_transport_model[i];
   }
}


TPMInterface& GKTransport::transportModel( const std::string& a_name )
{
   typedef std::map<std::string,int>::iterator mapIterator;
   const mapIterator it( m_species_map.find( a_name ) );
   CH_assert(it!=m_species_map.end());
   const int index((*it).second);
   return *(m_transport_model[index]);
}


void GKTransport::accumulateRHS( KineticSpeciesPtrVect&       a_rhs,
                                 const KineticSpeciesPtrVect& a_soln,
                                 const Real                   a_time )
{
   for (int species(0); species<a_rhs.size(); species++) {
      KineticSpecies& rhs_species( *(a_rhs[species]) );
      const std::string species_name( rhs_species.name() );
      TPMInterface& TPM( transportModel( species_name ) );
      TPM.evalTpmRHS( rhs_species, a_soln, species, a_time );
   }
}

Real GKTransport::computeTimeScale( const KineticSpeciesPtrVect& soln )
{
   Real dt_stable(DBL_MAX);
   for (int species(0); species<soln.size(); species++) {
      KineticSpecies& this_species( *(soln[species]) );
      const std::string species_name( this_species.name() );
      TPMInterface& TPM( transportModel( species_name ) );
      dt_stable = TPM.computeStableDt(this_species);
   } 
   return dt_stable;
}

Real GKTransport::computeDtExplicitTI( const KineticSpeciesPtrVect& soln )
{
   // This is copied from computeDt in GKCollision
   /*
   Real dt(DBL_MAX);
   int count(0);
   std::map<std::string,int>::iterator it;
   for (it=m_transport_model_name.begin(); it!=m_transport_model_name.end(); ++it) {
      Real tmp = m_transport_model[it->second]->computeDtExplicitTI(soln,it->second);
      dt = (tmp < dt ? tmp : dt);
      count++;
   }
   return (count ? dt : -1);
   */   

   Real dt(DBL_MAX);
   
   return dt;
}

Real GKTransport::computeDtImExTI( const KineticSpeciesPtrVect& soln )
{
   // This is copied from computeDt in GKCollision
   /*
   Real dt(DBL_MAX);
   int count(0);
   std::map<std::string,int>::iterator it;
   for (it=m_transport_model_name.begin(); it!=m_transport_model_name.end(); ++it) {
      Real tmp = m_transport_model[it->second]->computeDtImExTI(soln,it->second);
      dt = (tmp < dt ? tmp : dt);
      count++;
   }
   return (count ? dt : -1);
   */   

   Real dt(DBL_MAX);
   
   return dt;
}

void GKTransport::preTimeStep( const KineticSpeciesPtrVect& a_soln,
                               const Real a_time,
                               const KineticSpeciesPtrVect& a_soln_physical )

{
  for (int species(0); species<a_soln.size(); species++) {
    KineticSpecies&           soln_species(*(a_soln[species]));
    const std::string         species_name(soln_species.name());
    TPMInterface& TPM( transportModel( species_name ) );
    TPM.preTimeStep(a_soln, species, a_time, a_soln_physical);
  }
}

#include "NamespaceFooter.H"


