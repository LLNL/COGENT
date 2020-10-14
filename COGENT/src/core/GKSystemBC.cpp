#include "GKSystemBC.H"
#include "SingleNullPhaseCoordSys.H"
#include "KineticSpeciesBCFactory.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "MagGeom.H"
#include "EllipticOpBCFactory.H"
#include "FluidSpecies.H"
#include "FluidVarBC.H"
#include "FluidVarBCFactory.H"
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
     m_phase_geometry( *(a_state.geometry())),
     m_potential_bcs(NULL)
{
   a_pp.query( "gksystem.verbosity", m_verbosity );

   const std::string coord_sys_type( determineCoordSysType( m_phase_geometry ) );

   parsePotential( a_pp, coord_sys_type );
   parseKineticSpecies( a_pp, coord_sys_type, a_state.dataKinetic() );
   //   parseFluidSpecies( a_pp, coord_sys_type, a_state.dataFluid() );

   const CFG::MagGeom& mag_geom = m_phase_geometry.magGeom();
   mag_geom.defineEllipticOpBC( *m_potential_bcs );
   
   int ghost = (m_phase_geometry.secondOrder()) ? 0 : 1;
   const DisjointBoxLayout& dbl = m_phase_geometry.gridsFull();
   m_mapped_velocity.define(dbl, SpaceDim, ghost * IntVect::Unit);
}


GKSystemBC::~GKSystemBC()
{
   for (int i(0); i<m_kinetic_bcs.size(); i++ ) {
      delete m_kinetic_bcs[i];
   }

   if (m_potential_bcs) delete m_potential_bcs;
}


inline
void GKSystemBC::executeInternalExchanges( KineticSpecies& a_species ) const
{
   LevelData<FArrayBox>& dfn( a_species.distributionFunction() );
   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
   // Fill ghost cells except for those on physical boundaries
   geometry.fillInternalGhosts(dfn, true );
}


inline
void GKSystemBC::executeInternalExchanges( CFG::FluidSpecies& a_fluid ) const
{
   const CFG::MagGeom& geometry( a_fluid.configurationSpaceGeometry() );

   for (int n=0; n<a_fluid.num_cell_vars(); ++n) {
      CFG::LevelData<CFG::FArrayBox>& fld( a_fluid.cell_var(n) );
      // Fill ghost cells except for those on physical boundaries
      geometry.fillInternalGhosts( fld );
   }
}


void GKSystemBC::fillGhostCells( GKState&                               a_state,
                                 const CFG::LevelData<CFG::FArrayBox>&  a_phi,
                                 const LevelData<FluxBox>&              a_E_field,
                                 const Real&                            a_time ) const
{
   fillKineticSpeciesGhostCells( a_state.dataKinetic(), a_phi, a_E_field, a_time );
   fillFluidSpeciesGhostCells( a_state.dataFluid(), a_E_field, a_time );
}


inline
void GKSystemBC::fillKineticSpeciesGhostCells( KineticSpeciesPtrVect&                 a_species,
                                               const CFG::LevelData<CFG::FArrayBox>&  a_phi,
                                               const LevelData<FluxBox>&              a_E_field,
                                               const Real&                            a_time ) const
{
   CH_TIMERS("GKSystemBC::fillKineticSpeciesGhostCells");
   CH_TIMER("computeMappedVelocity", t_compute_mapped_velocity);
   CH_TIMER("executeExchanges", t_execute_exchanges);
   CH_TIMER("apply_bc", t_apply_bc);
   for (int s_index(0); s_index<a_species.size(); s_index++) {

      // It is important to fill internal ghosts prior to physical
      // as they amay be reqired for co-dim ghost evaluation 
      CH_START(t_execute_exchanges);
      KineticSpecies& species_physical( *(a_species[s_index]) );
      executeInternalExchanges( species_physical );
      CH_STOP(t_execute_exchanges);

      CH_START(t_compute_mapped_velocity);
      species_physical.computeMappedVelocity( m_mapped_velocity, a_E_field, false, a_time );
      CH_STOP(t_compute_mapped_velocity);

      CH_START(t_apply_bc);
      KineticSpeciesBC& ksbc( kineticSpeciesBC( species_physical.name() ) );
      ksbc.apply( a_species, s_index, a_phi, m_mapped_velocity, a_time );
      CH_STOP(t_apply_bc);

      // For the case of a SN sheared geoemtry need to
      // refill internal ghosts again to handle the saw-tooth BCs
      //const MultiBlockCoordSys& coord_sys( *(m_phase_geometry.coordSysPtr()) );
      const CFG::MagGeom& mag_geom = m_phase_geometry.magGeom();
      if (mag_geom.mixedBoundaries()) {
	executeInternalExchanges( species_physical );
      }
      
   }
}

inline
void GKSystemBC::fillFluidSpeciesGhostCells( CFG::FluidSpeciesPtrVect&              a_species,
                                             const LevelData<FluxBox>&              a_E_field,
                                             const Real&                            a_time ) const
{
   for (int s(0); s<a_species.size(); s++) {

      // Here, we only fill ghost cells for FluidSpecies objects, assuming that
      // operators on other fluid species (e.g., the potential in VorticityOp) provide
      // their own methods for handling boundary conditions

      try {
         CFG::FluidSpecies& species_physical( dynamic_cast<CFG::FluidSpecies&>(*(a_species[s])) );
         executeInternalExchanges( species_physical );
         CFG::FluidVarBC& fsbc( fluidVarBC( species_physical.name() ) );
         fsbc.apply( species_physical, a_time );
      }
      catch (std::bad_cast& bc) {}
   }
}


/// deprecated ////////////////////

void GKSystemBC::fillGhostCells( KineticSpeciesPtrVect&                 a_species, 
                                 const CFG::LevelData<CFG::FArrayBox>&  a_phi,
                                 const LevelData<FluxBox>&              a_E_field,
                                 const Real&                            a_time ) const
{
   fillKineticSpeciesGhostCells( a_species, a_phi, a_E_field, a_time );
}

void GKSystemBC::fillGhostCells(CFG::FluidSpeciesPtrVect&  a_species,
				const LevelData<FluxBox>&  a_E_field,
				const Real&                a_time ) const
{
  fillFluidSpeciesGhostCells( a_species, a_E_field, a_time );
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


CFG::FluidVarBC& GKSystemBC::fluidVarBC( const std::string& a_name ) const
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
   CFG::EllipticOpBCFactory elliptic_op_bc_factory;
   m_potential_bcs = elliptic_op_bc_factory.create( name,
                                                    ppsp,
                                                    *(m_phase_geometry.magGeom().getCoordSys()),
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


void GKSystemBC::parseFluidSpecies( ParmParse&                       a_pp,
                                    const std::string&               a_coord_sys_type,
                                    const CFG::FluidSpeciesPtrVect&  a_species )
{
#if 0
   for (int s(0); s<a_species.size(); s++) {
      const std::string& name( a_species[s]->name() );
      const std::string prefix( "BC." + name );
      ParmParse ppsp( prefix.c_str() );
      CFG::FluidVarBCFactory ksp_factory;
      CFG::FluidVarBC* bc( ksp_factory.create( name,
                                               ppsp,
                                               a_coord_sys_type,
                                               m_verbosity ) );
      m_fluid_bcs.push_back( bc );
   }
#endif
}


#include "NamespaceFooter.H"
