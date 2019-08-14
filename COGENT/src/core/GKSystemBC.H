#ifndef _GK_SYSTEM_BC_
#define _GK_SYSTEM_BC_

#include "KineticSpecies.H"
#include "KineticSpeciesBC.H"
#include "GridFunction.H"
#include "GKTimeIntegration.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FluidVarBC.H"
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
      void fillGhostCells( GKState&                               state,
                           const CFG::LevelData<CFG::FArrayBox>&  phi,
                           const LevelData<FluxBox>&              E_field,
                           const Real&                            time ) const;

/// deprecated ////////////////////

      void fillGhostCells( KineticSpeciesPtrVect&                 species, 
                           const CFG::LevelData<CFG::FArrayBox>&  phi,
                           const LevelData<FluxBox>&              E_field,
                           const Real&                            time ) const;

      void fillGhostCells( CFG::FluidSpeciesPtrVect&              species,
                           const LevelData<FluxBox>&              E_field,
                           const Real&                            time ) const;

/// deprecated ////////////////////

      CFG::EllipticOpBC& getEllipticOpBC() const {return *m_potential_bcs;}

   private:

      // prevent copying
      GKSystemBC( const GKSystemBC& );
      const GKSystemBC& operator=( const GKSystemBC& );

      bool hasBoundaryCondition( std::string& name ) const;

      KineticSpeciesBC& kineticSpeciesBC( const std::string& name ) const;

      CFG::FluidVarBC& fluidVarBC( const std::string& a_name ) const;

      CFG::GridFunction& fieldBC( const std::string& a_name ) const;

      inline
      void executeInternalExchanges( KineticSpecies& a_species ) const;

      inline
      void executeInternalExchanges( CFG::FluidSpecies& a_fluid ) const;

      inline
      void fillKineticSpeciesGhostCells( KineticSpeciesPtrVect&                 species,
                                         const CFG::LevelData<CFG::FArrayBox>&  phi,
                                         const LevelData<FluxBox>&              E_field,
                                         const Real&                            time ) const;
   
      inline
      void fillFluidSpeciesGhostCells( CFG::FluidSpeciesPtrVect&  species,
                                       const LevelData<FluxBox>&  E_field,
                                       const Real&                time ) const;

      void parsePotential( ParmParse& a_pp,
                           const std::string& a_coord_sys_type );

      void parseKineticSpecies( ParmParse& a_pp,
                                const std::string& a_coord_sys_type,
                                const KineticSpeciesPtrVect& a_species  );

      void parseFluidSpecies( ParmParse& a_pp,
                              const std::string& a_coord_sys_type,
                              const CFG::FluidSpeciesPtrVect& a_species );

      int                                 m_verbosity;
      const PhaseGeom&                    m_phase_geometry;
      std::vector<KineticSpeciesBC*>      m_kinetic_bcs;
      std::vector<CFG::FluidVarBC*>       m_fluid_bcs;
      CFG::EllipticOpBC*                  m_potential_bcs;
      mutable LevelData<FluxBox>          m_mapped_velocity;

      typedef std::map<std::string,RefCountedPtr<CFG::GridFunction> > FieldBCMap;
      FieldBCMap m_field_bcs;
};

#include "NamespaceFooter.H"

#endif
