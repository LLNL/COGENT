#ifndef _GK_SYSTEM_IC_
#define _GK_SYSTEM_IC_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "ParmParse.H"
#include "MultiBlockLevelGeom.H"
#include "PhaseGeom.H"
#include "KineticSpecies.H"
#include "KineticFunction.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "IntVect.H"
#include "GKPoisson.H"
#include "GridFunction.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include <vector>
#include <map>

#include "NamespaceHeader.H"

class GKSystemIC
{
   public:

      GKSystemIC( ParmParse& pp,
                  const PhaseGeom& geometry,
                  const KineticSpeciesPtrVect& species );

      ~GKSystemIC();

      void initializeSystem( KineticSpeciesPtrVect& a_species,
                             CFG::LevelData<CFG::FArrayBox>& a_potential,
                             const Real& a_time ) const;

      void initializePotentialField(
         CFG::LevelData<CFG::FArrayBox>& a_potential,
         const Real& a_time ) const;

   private:

      bool hasInitialCondition( std::string& name ) const;

      const KineticFunction& phaseSpaceIC( const std::string& name ) const;

      const CFG::GridFunction& configSpaceIC( const std::string& name ) const;

      void parsePotential( ParmParse& a_pp );

      void parseSpecies( ParmParse& a_pp, const KineticSpeciesPtrVect& a_species  );

      inline
      void initializeDistributionFunctions(
         KineticSpeciesPtrVect& a_species,
         const Real& a_time ) const;
      const PhaseGeom&                m_phase_geometry;
      const CFG::MultiBlockLevelGeom& m_mag_geometry;

      typedef std::map<std::string,RefCountedPtr<CFG::GridFunction> > ConfigICMap;
      ConfigICMap m_config_ics;

      typedef std::map<std::string,RefCountedPtr<KineticFunction> > SpeciesICMap;
      SpeciesICMap m_phase_ics;

      const CFG::IntVect              m_cfg_nghosts;
      int                             m_verbosity;


};

#include "NamespaceFooter.H"

#endif
