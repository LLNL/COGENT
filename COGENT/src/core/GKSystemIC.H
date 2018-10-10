#ifndef _GK_SYSTEM_IC_
#define _GK_SYSTEM_IC_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "ParmParse.H"
#include "MultiBlockLevelGeom.H"
#include "PhaseGeom.H"
#include "KineticSpecies.H"
#include "KineticSpeciesBC.H"
#include "KineticFunction.H"
#include "GKTimeIntegration.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "IntVect.H"
#include "GKPoisson.H"
#include "GridFunction.H"
#include "FluidSpecies.H"
#include "FluidVarBC.H"
//#include "FluidFunction.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include <vector>
#include <map>

#include "NamespaceHeader.H"
namespace CFG = CFG_NAMESPACE;

class GKSystemIC
{
   public:

      GKSystemIC( ParmParse& a_pp, const GKState& a_state );

      ~GKSystemIC();

      void initializeSystem( GKState& a_state,
                             CFG::LevelData<CFG::FArrayBox>& a_potential,
                             const Real& a_time ) const;

/// deprecated ////////////////////
      void initializeSystem( KineticSpeciesPtrVect& a_species,
                             CFG::LevelData<CFG::FArrayBox>& a_potential,
                             const Real& a_time ) const;
/// deprecated ////////////////////

      void initializePotentialField(
         CFG::LevelData<CFG::FArrayBox>& a_potential,
         const Real& a_time ) const;

   private:

      bool hasInitialCondition( std::string& name ) const;

      const KineticFunction& kineticSpeciesIC( const std::string& name ) const;

      const CFG::GridFunction& fluidSpeciesIC( const std::string& a_name ) const;

      void parsePotential( ParmParse& a_pp );

      void parseKineticSpecies( ParmParse& a_pp,
                                const KineticSpeciesPtrVect& a_species );

      inline
      void initializeKineticSpecies( KineticSpeciesPtrVect& a_species,
                                     const Real& a_time ) const;

      inline
      void  initializeScalars( ScalarPtrVect& a_scalars,
                               const Real&    a_time ) const;

      const PhaseGeom&                m_phase_geometry;
      const CFG::MultiBlockLevelGeom& m_mag_geometry;

      typedef std::map<std::string,RefCountedPtr<CFG::GridFunction> > FluidSpeciesICMap;
      FluidSpeciesICMap m_fluid_ics;

      typedef std::map<std::string,RefCountedPtr<KineticFunction> > KineticSpeciesICMap;
      KineticSpeciesICMap m_kinetic_ics;

      const CFG::IntVect              m_cfg_nghosts;
      int                             m_verbosity;


};

#include "NamespaceFooter.H"

#endif
