#include "GKSystemIC.H"

#include "Box.H"
//#include "FourthOrderUtil.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
//#include "FourthOrderUtil.H"
#include "GridFunctionLibrary.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "KineticFunctionLibrary.H"

#include "NamespaceHeader.H"

GKSystemIC::GKSystemIC( ParmParse& a_pp,
                        const PhaseGeom& a_phase_geometry,
                        const KineticSpeciesPtrVect& a_species )
   : m_phase_geometry( a_phase_geometry ),
     m_mag_geometry( a_phase_geometry.magGeom() ),
     m_verbosity(0)
{
   a_pp.query( "gksystem.verbosity", m_verbosity );

   parsePotential( a_pp );
   parseSpecies( a_pp, a_species );
}


GKSystemIC::~GKSystemIC()
{
}


void GKSystemIC::initializeSystem( KineticSpeciesPtrVect& a_species,
                                   CFG::LevelData<CFG::FArrayBox>& a_potential,
                                   const Real& a_time ) const
{
   initializePotentialField( a_potential, a_time );
   initializeDistributionFunctions( a_species, a_time );
}


void GKSystemIC::initializePotentialField(
   CFG::LevelData<CFG::FArrayBox>& a_potential,
   const Real& a_time  ) const
{
   const CFG::GridFunction& icpot( configSpaceIC( "potential" ) );
   icpot.assign( a_potential, m_mag_geometry, a_time );
}


bool GKSystemIC::hasInitialCondition( std::string& a_name ) const
{
   bool ic_found(false);
   SpeciesICMap::const_iterator sit = m_phase_ics.find( a_name );
   if ( sit != m_phase_ics.end() ) {
      ic_found = true;
   }

   ConfigICMap::const_iterator cit = m_config_ics.find( a_name );
   if ( cit != m_config_ics.end() ) {
      ic_found = true;
   }

   return ic_found;
}


const KineticFunction& GKSystemIC::phaseSpaceIC( const std::string& a_name ) const
{
   SpeciesICMap::const_iterator it = m_phase_ics.find( a_name );
   if ( it == m_phase_ics.end() ) {
      const std::string msg( "GKSystemIC: Initial condition " + a_name + " not found!" );
      MayDay::Error( msg.c_str() );
   }
   return *(it->second);
}


const CFG::GridFunction& GKSystemIC::configSpaceIC( const std::string& a_name ) const
{
   ConfigICMap::const_iterator it = m_config_ics.find( a_name );
   if ( it == m_config_ics.end() ) {
      const std::string msg( "GKSystemIC: Initial condition " + a_name + " not found!" );
      MayDay::Error( msg.c_str() );
   }
   return *(it->second);
}


void GKSystemIC::parsePotential( ParmParse& a_pp )
{
   CFG::GridFunctionLibrary* library = CFG::GridFunctionLibrary::getInstance();
   ParmParse pppot( "IC.potential" );
   std::string function_name;
   pppot.query( "function", function_name );
   RefCountedPtr<CFG::GridFunction> ic( library->find( function_name ) );
   m_config_ics.insert( ConfigICMap::value_type( "potential", ic ) );
}


void GKSystemIC::parseSpecies( ParmParse& a_pp,
                               const KineticSpeciesPtrVect& a_species )
{
   KineticFunctionLibrary* library = KineticFunctionLibrary::getInstance();
   for (int s(0); s<a_species.size(); s++) {
      const std::string& species_name( a_species[s]->name() );
      const std::string prefix( "IC." + species_name );
      ParmParse ppsp( prefix.c_str() );
      std::string function_name;
      ppsp.query( "function", function_name );
      RefCountedPtr<KineticFunction> ic( library->find( function_name ) );
      m_phase_ics.insert( SpeciesICMap::value_type( species_name, ic ) );
   }
}


inline
void GKSystemIC::initializeDistributionFunctions(
   KineticSpeciesPtrVect& a_species,
   const Real& a_time ) const
{
   for (int species(0); species<a_species.size(); species++) {
      KineticSpecies& soln_species( *(a_species[species]) );
      const KineticFunction& ksic( phaseSpaceIC( soln_species.name() ) );
      ksic.assign( soln_species, a_time );
   }
}



#include "NamespaceFooter.H"
