#include "GKSystemBC.H"

#include "KineticSpeciesBCFactory.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "MagGeom.H"
#include "PotentialBCFactory.H"
#include "FluidSpeciesBC.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"
namespace CFG = CFG_NAMESPACE;


inline
std::string determineCoordSysType( const PhaseGeom& a_phase_geometry )
{
   const PhaseCoordSys& coord_sys( a_phase_geometry.phaseCoordSys() );
   return coord_sys.type();
}


GKSystemBC::GKSystemBC( ParmParse& a_pp,
                        const GKState& a_state )
   : m_verbosity(0),
     m_phase_geometry( *(a_state.geometry()) )

{
   a_pp.query( "gksystem.verbosity", m_verbosity );

   const std::string coord_sys_type( determineCoordSysType( m_phase_geometry ) );

   parsePotential( a_pp, coord_sys_type );
   parseKineticSpecies( a_pp, coord_sys_type, a_state.dataKinetic() );
   parseFluidSpecies( a_pp, coord_sys_type, a_state.dataFluid() );
   parseFields( a_pp, coord_sys_type, a_state.dataField() );

   const CFG::MagGeom& mag_geom = m_phase_geometry.magGeom();
   mag_geom.definePotentialBC( *m_potential_bcs );
}


GKSystemBC::~GKSystemBC()
{
   for (int i(0); i<m_kinetic_bcs.size(); i++ ) {
      delete m_kinetic_bcs[i];
   }
}


inline
void GKSystemBC::executeInternalExchanges( KineticSpecies& a_species ) const
{
   LevelData<FArrayBox>& dfn( a_species.distributionFunction() );
   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
   // Fill ghost cells except for those on physical boundaries
   geometry.fillInternalGhosts(dfn);
}


inline
void GKSystemBC::executeInternalExchanges( CFG::FluidSpecies& a_fluid ) const
{
   CFG::LevelData<CFG::FArrayBox>& fld( a_fluid.data() );
   const CFG::MagGeom& geometry( a_fluid.configurationSpaceGeometry() );
   // Fill ghost cells except for those on physical boundaries
   geometry.fillInternalGhosts( fld );
}

inline
void GKSystemBC::executeInternalExchanges( CFG::Field& a_field ) const
{
   CFG::LevelData<CFG::FArrayBox>& fld( a_field.data() );
   const CFG::MagGeom& geometry( a_field.configurationSpaceGeometry() );
   // Fill ghost cells except for those on physical boundaries
   geometry.fillInternalGhosts( fld );
}


void GKSystemBC::fillGhostCells(
   GKState& a_state,
   const CFG::LevelData<CFG::FArrayBox>& a_phi,
   const LevelData<FluxBox>& a_E_field,
   const Real& a_time ) const
{
   fillKineticSpeciesGhostCells( a_state.dataKinetic(), a_phi, a_E_field, a_time );
   fillFluidSpeciesGhostCells( a_state.dataFluid(), a_phi, a_E_field, a_time );
   fillFieldsGhostCells( a_state.dataField(), a_phi, a_E_field, a_time );
}


inline
void GKSystemBC::fillKineticSpeciesGhostCells(
   KineticSpeciesPtrVect& a_species,
   const CFG::LevelData<CFG::FArrayBox>& a_phi,
   const LevelData<FluxBox>& a_E_field,
   const Real& a_time ) const
{
   for (int s(0); s<a_species.size(); s++) {

      KineticSpecies& species_physical( *(a_species[s]) );
      executeInternalExchanges( species_physical );

      LevelData<FluxBox> velocity;
      species_physical.computeMappedVelocity( velocity, a_E_field );

      KineticSpeciesBC& ksbc( kineticSpeciesBC( species_physical.name() ) );
      ksbc.apply( species_physical, a_phi, velocity, a_time );
   }
}

inline
void GKSystemBC::fillFluidSpeciesGhostCells(
   CFG::FluidSpeciesPtrVect& a_species,
   const CFG::LevelData<CFG::FArrayBox>& a_phi,
   const LevelData<FluxBox>& a_E_field,
   const Real& a_time ) const
{
   for (int s(0); s<a_species.size(); s++) {

      CFG::FluidSpecies& species_physical( *(a_species[s]) );
      executeInternalExchanges( species_physical );

      CFG::LevelData<CFG::FluxBox> velocity;
//      species_physical.computeMappedVelocity( velocity, a_E_field );

      CFG::FluidSpeciesBC& fsbc( fluidSpeciesBC( species_physical.name() ) );
      fsbc.apply( species_physical, a_phi, velocity, a_time );
   }
}

inline
void GKSystemBC::fillFieldsGhostCells(
   CFG::FieldPtrVect& a_fields,
   const CFG::LevelData<CFG::FArrayBox>& a_phi,
   const LevelData<FluxBox>& a_E_field,
   const Real& a_time ) const
{
   for (int comp(0); comp<a_fields.size(); comp++) {

      CFG::Field& field_physical( *(a_fields[comp]) );
      executeInternalExchanges( field_physical );

      CFG::LevelData<CFG::FluxBox> velocity;
//      field_physical.computeMappedVelocity( velocity, a_E_field );

      CFG::GridFunction& fsbc( fieldBC( field_physical.name() ) );
      fsbc.assign( field_physical.data(),
                   field_physical.configurationSpaceGeometry(),
                   a_time );
   }
}


/// deprecated ////////////////////
void GKSystemBC::fillGhostCells(
   KineticSpeciesPtrVect& a_species, 
   const CFG::LevelData<CFG::FArrayBox>& a_phi,
   const LevelData<FluxBox>& a_E_field,
   const Real& a_time ) const
{
   fillKineticSpeciesGhostCells( a_species, a_phi, a_E_field, a_time );
}
/// deprecated ////////////////////


bool GKSystemBC::hasBoundaryCondition( std::string& a_name ) const
{
   bool bc_found(false);
   for (int i(0); i<m_kinetic_bcs.size(); i++ ) {
      if (m_kinetic_bcs[i]->isForVariable( a_name )) {
         bc_found = true;
         break;
      }
   }

   for (int i(0); i<m_fluid_bcs.size(); i++ ) {
      if (m_fluid_bcs[i]->isForVariable( a_name )) {
         bc_found = true;
         break;
      }
   }

   FieldBCMap::const_iterator cit = m_field_bcs.find( a_name );
   if ( cit != m_field_bcs.end() ) {
      bc_found = true;
   }

   return bc_found;
}


KineticSpeciesBC& GKSystemBC::kineticSpeciesBC( const std::string& a_name ) const
{
   int index(-1);
   for (int i(0); i<m_kinetic_bcs.size(); i++ ) {
      if (m_kinetic_bcs[i]->isForVariable( a_name )) {
         index = i;
         break;
      }
   }
   if (index<0) {
      const std::string msg( "GKSystemBC: Boundary condition for species " + a_name + " not found!" );
      MayDay::Error( msg.c_str() );
   }
   return *(m_kinetic_bcs[index]);
}


CFG::FluidSpeciesBC& GKSystemBC::fluidSpeciesBC( const std::string& a_name ) const
{
   int index(-1);
   for (int i(0); i<m_fluid_bcs.size(); i++ ) {
      if (m_fluid_bcs[i]->isForVariable( a_name )) {
         index = i;
         break;
      }
   }
   if (index<0) {
      const std::string msg( "GKSystemBC: Boundary condition for fluid species" + a_name + " not found!" );
      MayDay::Error( msg.c_str() );
   }
   return *(m_fluid_bcs[index]);
}


CFG::GridFunction& GKSystemBC::fieldBC( const std::string& a_name ) const
{
   FieldBCMap::const_iterator it = m_field_bcs.find( a_name );
   if ( it == m_field_bcs.end() ) {
      const std::string msg( "GKSystemIC: Boundary condition for field " + a_name + " not found!" );
      MayDay::Error( msg.c_str() );
   }
   return *(it->second);
}


void GKSystemBC::parsePotential( ParmParse& a_pp,
                                 const std::string& a_coord_sys_type )
{
   const std::string name("potential");
   const std::string prefix( "BC." + name );
   ParmParse ppsp( prefix.c_str() );
   CFG::PotentialBCFactory potential_bc_factory;
   m_potential_bcs = potential_bc_factory.create( name,
                                                  ppsp,
                                                  a_coord_sys_type,
                                                  m_verbosity );
}


void GKSystemBC::parseKineticSpecies( ParmParse& a_pp,
                                      const std::string& a_coord_sys_type,
                                      const KineticSpeciesPtrVect& a_species )
{
   for (int s(0); s<a_species.size(); s++) {
      const std::string& name( a_species[s]->name() );
      const std::string prefix( "BC." + name );
      ParmParse ppsp( prefix.c_str() );
      KineticSpeciesBCFactory ksp_factory;
      KineticSpeciesBC* bc( ksp_factory.create( name,
                                                ppsp,
                                                a_coord_sys_type,
                                                m_verbosity ) );
      m_kinetic_bcs.push_back( bc );
   }
}


void GKSystemBC::parseFluidSpecies( ParmParse& a_pp,
                                    const std::string& a_coord_sys_type,
                                    const CFG::FluidSpeciesPtrVect& a_species )
{
}


void GKSystemBC::parseFields( ParmParse& a_pp,
                              const std::string& a_coord_sys_type,
                              const CFG::FieldPtrVect& a_fields )
{
}


#include "NamespaceFooter.H"
