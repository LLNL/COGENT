#include "GKSystemBC.H"

#include "KineticSpeciesBCFactory.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "PotentialBCFactory.H"
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
                        const PhaseGeom& a_phase_geometry,
                        const KineticSpeciesPtrVect& a_species )
   : m_verbosity(0),
     m_phase_geometry( a_phase_geometry )

{
   a_pp.query( "gksystem.verbosity", m_verbosity );

   const std::string coord_sys_type( determineCoordSysType( a_phase_geometry ) );

   parsePotential( a_pp, coord_sys_type );
   parseSpecies( a_pp, coord_sys_type, a_species );

   const CFG::MagGeom& mag_geom = a_phase_geometry.magGeom();
   mag_geom.definePotentialBC( *m_potential_bcs );
}


GKSystemBC::~GKSystemBC()
{
   for (int i(0); i<m_phase_bcs.size(); i++ ) {
      delete m_phase_bcs[i];
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


void GKSystemBC::fillGhostCells(
   KineticSpeciesPtrVect& a_species_phys,
   const CFG::LevelData<CFG::FArrayBox>& a_phi,
   const LevelData<FluxBox>& a_E_field,
   const Real& a_time ) const
{
   for (int s(0); s<a_species_phys.size(); s++) {

      KineticSpecies& species_phys( *(a_species_phys[s]) );
      executeInternalExchanges( species_phys );

      LevelData<FluxBox> velocity;
      species_phys.computeMappedVelocity( velocity, a_E_field );

      KineticSpeciesBC& ksbc( phaseSpaceBC( species_phys.name() ) );
      ksbc.apply( species_phys, a_phi, velocity, a_time );
   }
}


bool GKSystemBC::hasBoundaryCondition( std::string& a_name ) const
{
   bool bc_found(false);
   for (int i(0); i<m_phase_bcs.size(); i++ ) {
      if (m_phase_bcs[i]->isForVariable( a_name )) {
         bc_found = true;
         break;
      }
   }
   return bc_found;
}


KineticSpeciesBC& GKSystemBC::phaseSpaceBC( const std::string& a_name ) const
{
   int index(-1);
   for (int i(0); i<m_phase_bcs.size(); i++ ) {
      if (m_phase_bcs[i]->isForVariable( a_name )) {
         index = i;
         break;
      }
   }
   if (index<0) {
      const std::string msg( "GKSystemBC: Boundary condition for species " + a_name + " not found!" );
      MayDay::Error( msg.c_str() );
   }
   return *(m_phase_bcs[index]);
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


void GKSystemBC::parseSpecies( ParmParse& a_pp,
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
      m_phase_bcs.push_back( bc );
   }
}


#include "NamespaceFooter.H"
