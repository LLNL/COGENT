#ifndef _SN_SYSTEM_IC_
#define _SN_SYSTEM_IC_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "ParmParse.H"
#include "MagCoordSys.H"
#include "MagGeom.H"
#include "MultiBlockLevelGeom.H"
#include "MultiBlockLevelExchangeAverage.H"
#include "IntVect.H"
#include "FluidFunction.H"
#include "SNTimeIntegration.H"

#include <string>
#include <map>

#include "NamespaceHeader.H"

class SNSystemIC
{
   public:

      SNSystemIC( ParmParse& pp, const SNState& solution );

      ~SNSystemIC();

      void initializeSystem( SNState& solution,  const Real& a_time ) const;

   private:

      const FluidFunction& initialCondition( const std::string& name ) const;

      void parseFluidSpecies( ParmParse& a_pp, const FluidSpeciesPtrVect& a_species );

      typedef std::map<std::string,RefCountedPtr<FluidFunction> > FluidSpeciesICMap;
      FluidSpeciesICMap m_ics;

      const IntVect m_nghosts;
      int m_verbosity;
};

#include "NamespaceFooter.H"

#endif
