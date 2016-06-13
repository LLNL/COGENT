#ifndef _GK_SYSTEM_BC_
#define _GK_SYSTEM_BC_

#include "KineticSpecies.H"
#include "KineticSpeciesBC.H"
#include "GridFunction.H"
#include "GKTimeIntegration.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FluidSpeciesBC.H"
#include "GridFunction.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

class GKSystemBC
{
   public:

      GKSystemBC( ParmParse& pp, const GKState& a_state );

      ~GKSystemBC();

      // fills ghost cells on the physical domain
      void fillGhostCells(
         GKState& a_state,
         const CFG::LevelData<CFG::FArrayBox>& a_phi,
         const LevelData<FluxBox>& a_E_field,
         const Real& a_time ) const;

/// deprecated ////////////////////
      void fillGhostCells(
         KineticSpeciesPtrVect& a_species, 
         const CFG::LevelData<CFG::FArrayBox>& a_phi,
         const LevelData<FluxBox>& a_E_field,
         const Real& a_time ) const;
/// deprecated ////////////////////

      CFG::PotentialBC& getPotentialBC() const {return *m_potential_bcs;}

   private:

      // prevent copying
      GKSystemBC( const GKSystemBC& );
      const GKSystemBC& operator=( const GKSystemBC& );

#if 0
      void fillBlockBoundaryCellsInCompSpace(
         KineticSpecies& species_comp,
         const KineticSpecies& species_phys ) const;
#endif

      bool hasBoundaryCondition( std::string& name ) const;

      KineticSpeciesBC& kineticSpeciesBC( const std::string& name ) const;

      CFG::FluidSpeciesBC& fluidSpeciesBC( const std::string& a_name ) const;

      CFG::GridFunction& fieldBC( const std::string& a_name ) const;

      inline
      void executeInternalExchanges( KineticSpecies& a_species ) const;

      inline
      void executeInternalExchanges( CFG::FluidSpecies& a_fluid ) const;

      inline
      void executeInternalExchanges( CFG::Field& a_field ) const;

      inline
      void fillKineticSpeciesGhostCells(
         KineticSpeciesPtrVect& a_species,
         const CFG::LevelData<CFG::FArrayBox>& a_phi,
         const LevelData<FluxBox>& a_E_field,
         const Real& a_time ) const;
   
      inline
      void fillFluidSpeciesGhostCells(
         CFG::FluidSpeciesPtrVect& a_species,
         const CFG::LevelData<CFG::FArrayBox>& a_phi,
         const LevelData<FluxBox>& a_E_field,
         const Real& a_time ) const;

      inline
      void fillFieldsGhostCells( CFG::FieldPtrVect& a_fields,
                                 const CFG::LevelData<CFG::FArrayBox>& a_phi,
                                 const LevelData<FluxBox>& a_E_field,
                                 const Real& a_time ) const;

      void parsePotential( ParmParse& a_pp,
                           const std::string& a_coord_sys_type );

      void parseKineticSpecies( ParmParse& a_pp,
                                const std::string& a_coord_sys_type,
                                const KineticSpeciesPtrVect& a_species  );

      void parseFluidSpecies( ParmParse& a_pp,
                              const std::string& a_coord_sys_type,
                              const CFG::FluidSpeciesPtrVect& a_species );
  
      void parseFields( ParmParse& a_pp,
                        const std::string& a_coord_sys_type,
                        const CFG::FieldPtrVect& a_fields );

      int                                 m_verbosity;
      const PhaseGeom&                    m_phase_geometry;
      std::vector<KineticSpeciesBC*>      m_kinetic_bcs;
      std::vector<CFG::FluidSpeciesBC*>   m_fluid_bcs;
      CFG::PotentialBC*                   m_potential_bcs;

      typedef std::map<std::string,RefCountedPtr<CFG::GridFunction> > FieldBCMap;
      FieldBCMap m_field_bcs;
};

#include "NamespaceFooter.H"

#endif
