#ifndef _GK_SYSTEM_BC_
#define _GK_SYSTEM_BC_

#include "KineticSpecies.H"
#include "KineticSpeciesBC.H"

#include "NamespaceHeader.H"

class GKSystemBC
{
   public:

      GKSystemBC( ParmParse& pp,
                  const PhaseGeom& geometry,
                  const KineticSpeciesPtrVect& species );

      ~GKSystemBC();

      // fills ghost cells on the physical domain
      void fillGhostCells(
         KineticSpeciesPtrVect& a_species_phys,
         const CFG::LevelData<CFG::FArrayBox>& a_phi,
         const LevelData<FluxBox>& a_E_field,
         const Real& a_time ) const;

      CFG::PotentialBC& getPotentialBC() const {return *m_potential_bcs;}

   private:

      // prevent copying
      GKSystemBC( const GKSystemBC& );
      const GKSystemBC& operator=( const GKSystemBC& );

      void executeInternalExchanges( KineticSpecies& a_species ) const;

#if 0
      void fillBlockBoundaryCellsInCompSpace(
         KineticSpecies& species_comp,
         const KineticSpecies& species_phys ) const;
#endif

      bool hasBoundaryCondition( std::string& name ) const;

      KineticSpeciesBC& phaseSpaceBC( const std::string& name ) const;

      void parsePotential( ParmParse& a_pp,
                           const std::string& a_coord_sys_type );

      void parseSpecies( ParmParse& a_pp,
                         const std::string& a_coord_sys_type,
                         const KineticSpeciesPtrVect& a_species  );

      int                                 m_verbosity;
      const PhaseGeom&                    m_phase_geometry;
      std::vector<KineticSpeciesBC*>      m_phase_bcs;
      CFG::PotentialBC*                   m_potential_bcs;
};

#include "NamespaceFooter.H"

#endif
