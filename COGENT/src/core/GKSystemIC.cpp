#include "GKSystemIC.H"

#include "Box.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "GridFunction.H"
#include "GridFunctionLibrary.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "KineticFunctionLibrary.H"

#include "NamespaceHeader.H"

GKSystemIC::GKSystemIC( ParmParse& a_pp,
                        const GKState& a_state )
   : m_phase_geometry( *(a_state.geometry()) ),
     m_mag_geometry( m_phase_geometry.magGeom() ),
     m_verbosity(0)
{
   a_pp.query( "gksystem.verbosity", m_verbosity );

   parsePotential( a_pp );
   parseKineticSpecies( a_pp, a_state.dataKinetic() );
}


GKSystemIC::~GKSystemIC()
{
}


void GKSystemIC::initializeSystem( GKState& a_state,
                                   CFG::LevelData<CFG::FArrayBox>& a_potential,
                                   const Real& a_time ) const
{
   initializePotentialField( a_potential, a_time );
   initializeKineticSpecies( a_state.dataKinetic(), a_time );
   initializeScalars( a_state.dataScalar(), a_time );
}

inline
void GKSystemIC::initializeKineticSpecies( KineticSpeciesPtrVect& a_species,
                                           const Real& a_time ) const
{
   for (int species(0); species<a_species.size(); species++) {
      KineticSpecies& soln_species( *(a_species[species]) );
      const KineticFunction& ksic( kineticSpeciesIC( soln_species.name() ) );
      ksic.assign( soln_species, a_time );
   }
}

inline
void GKSystemIC::initializeScalars( ScalarPtrVect&      a_scalars,
                                    const Real&         a_time ) const
{
   // Until we invent another mechanism, all scalars are initialized to zero
   for (int scalar(0); scalar<a_scalars.size(); scalar++) {
      Vector<Real>& scalar_data = a_scalars[scalar]->data();
      for (int i=0; i<scalar_data.size(); ++i) {
         scalar_data[i] = 0.;
      }
   }
}


/// deprecated ////////////////////
void GKSystemIC::initializeSystem( KineticSpeciesPtrVect& a_species,
                                   CFG::LevelData<CFG::FArrayBox>& a_potential,
                                   const Real& a_time ) const
{
   initializePotentialField( a_potential, a_time );
   initializeKineticSpecies( a_species, a_time );
}
/// deprecated ////////////////////


void GKSystemIC::initializePotentialField(
   CFG::LevelData<CFG::FArrayBox>& a_potential,
   const Real& a_time  ) const
{
   const CFG::GridFunction& icpot( fluidSpeciesIC( "potential" ) );
   icpot.assign( a_potential, m_mag_geometry, a_time, false );
}


bool GKSystemIC::hasInitialCondition( std::string& a_name ) const
{
   bool ic_found(false);
   KineticSpeciesICMap::const_iterator ksit = m_kinetic_ics.find( a_name );
   if ( ksit != m_kinetic_ics.end() ) {
      ic_found = true;
   }

   FluidSpeciesICMap::const_iterator fsit = m_fluid_ics.find( a_name );
   if ( fsit != m_fluid_ics.end() ) {
      ic_found = true;
   }

   return ic_found;
}


const KineticFunction& GKSystemIC::kineticSpeciesIC( const std::string& a_name ) const
{
   KineticSpeciesICMap::const_iterator it = m_kinetic_ics.find( a_name );
   if ( it == m_kinetic_ics.end() ) {
      const std::string msg( "GKSystemIC: Initial condition " + a_name + " not found!" );
      MayDay::Error( msg.c_str() );
   }
   return *(it->second);
}


const CFG::GridFunction& GKSystemIC::fluidSpeciesIC( const std::string& a_name ) const
{
   FluidSpeciesICMap::const_iterator it = m_fluid_ics.find( a_name );
   if ( it == m_fluid_ics.end() ) {
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
   m_fluid_ics.insert( FluidSpeciesICMap::value_type( "potential", ic ) );
}


void GKSystemIC::parseKineticSpecies( ParmParse& a_pp,
                                      const KineticSpeciesPtrVect& a_species )
{
   if ( a_species.size() > 0 ) {
      KineticFunctionLibrary* library = KineticFunctionLibrary::getInstance();
      for (int s(0); s<a_species.size(); s++) {
         const std::string& species_name( a_species[s]->name() );
         const std::string prefix( "IC." + species_name );
         ParmParse ppsp( prefix.c_str() );
         std::string function_name;
         ppsp.query( "function", function_name );
         RefCountedPtr<KineticFunction> ic( library->find( function_name ) );
         m_kinetic_ics.insert( KineticSpeciesICMap::value_type( species_name, ic ) );
      }
   }
}


#include "NamespaceFooter.H"
