#include "GKSystemIC.H"

#include "Box.H"
//#include "FourthOrderUtil.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
//#include "FourthOrderUtil.H"
#include "GridFunction.H"
#include "GridFunctionLibrary.H"
#include "FluidFunction.H"
#include "FluidFunctionLibrary.H"
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
   parseFluidSpecies( a_pp, a_state.dataFluid() );
   parseFields( a_pp, a_state.dataField() );
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
   initializeFluidSpecies( a_state.dataFluid(), a_time );
   initializeFields( a_state.dataField(), a_time );
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
void GKSystemIC::initializeFluidSpecies( CFG::FluidSpeciesPtrVect& a_species,
                                         const Real& a_time ) const
{
   for (int species(0); species<a_species.size(); species++) {
      CFG::FluidSpecies& soln_species( *(a_species[species]) );
      const CFG::FluidFunction& fsic( fluidSpeciesIC( soln_species.name() ) );
      fsic.assign( soln_species, a_time );
   }
}

inline
void GKSystemIC::initializeFields( CFG::FieldPtrVect& a_fields,
                                   const Real& a_time ) const
{
   for (int comp(0); comp<a_fields.size(); comp++) {
      CFG::Field& field( *(a_fields[comp]) );
      const CFG::GridFunction& fsic( fieldIC( field.name() ) );
      fsic.assign( field.data(), field.configurationSpaceGeometry(), a_time );
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
   const CFG::GridFunction& icpot( fieldIC( "potential" ) );
   icpot.assign( a_potential, m_mag_geometry, a_time );
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

   FieldICMap::const_iterator cit = m_field_ics.find( a_name );
   if ( cit != m_field_ics.end() ) {
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


const CFG::FluidFunction& GKSystemIC::fluidSpeciesIC( const std::string& a_name ) const
{
   FluidSpeciesICMap::const_iterator it = m_fluid_ics.find( a_name );
   if ( it == m_fluid_ics.end() ) {
      const std::string msg( "GKSystemIC: Initial condition " + a_name + " not found!" );
      MayDay::Error( msg.c_str() );
   }
   return *(it->second);
}


const CFG::GridFunction& GKSystemIC::fieldIC( const std::string& a_name ) const
{
   FieldICMap::const_iterator it = m_field_ics.find( a_name );
   if ( it == m_field_ics.end() ) {
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
   m_field_ics.insert( FieldICMap::value_type( "potential", ic ) );
}


void GKSystemIC::parseKineticSpecies( ParmParse& a_pp,
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
      m_kinetic_ics.insert( KineticSpeciesICMap::value_type( species_name, ic ) );
   }
}


void GKSystemIC::parseFluidSpecies( ParmParse& a_pp,
                                    const CFG::FluidSpeciesPtrVect& a_species )
{
   CFG::GridFunctionLibrary* library = CFG::GridFunctionLibrary::getInstance();
   for (int s(0); s<a_species.size(); s++) {
      const std::string& species_name( a_species[s]->name() );
      const std::string prefix( "IC." + species_name );
      ParmParse ppfsp( prefix.c_str() );
      std::string function_name;
      ppfsp.query( "function", function_name );
      RefCountedPtr<CFG::GridFunction> ic( library->find( function_name ) );
      m_fluid_ics.insert( FluidSpeciesICMap::value_type( species_name, ic ) );
   }
}


void GKSystemIC::parseFields( ParmParse& a_pp,
                              const CFG::FieldPtrVect& a_field )
{
   CFG::GridFunctionLibrary* library = CFG::GridFunctionLibrary::getInstance();
   for (int s(0); s<a_field.size(); s++) {
      const std::string& field_name( a_field[s]->name() );
      const std::string prefix( "IC." + field_name );
      ParmParse ppfield( prefix.c_str() );
      std::string function_name;
      ppfield.query( "function", function_name );
      RefCountedPtr<CFG::GridFunction> ic( library->find( function_name ) );
      m_field_ics.insert( FieldICMap::value_type( field_name, ic ) );
   }
}


#include "NamespaceFooter.H"
