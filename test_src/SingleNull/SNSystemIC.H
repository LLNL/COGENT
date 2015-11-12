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
#include "GridFunction.H"

#include <string>
#include <map>

#include "NamespaceHeader.H"

class SNSystemIC
{
   public:

      SNSystemIC( ParmParse& pp,
                  const MagGeom& geometry,
                  const LevelData<FArrayBox>& solution );

      ~SNSystemIC();

      void initializeSystem( LevelData<FArrayBox>& solution,
                             const Real& a_time ) const;

   private:

      const GridFunction& initialCondition( const std::string& name ) const;

      void parseInput( ParmParse& a_pp, const LevelData<FArrayBox>& a_solution  );

      const MagGeom& m_mag_geometry;

      typedef std::map<std::string,RefCountedPtr<GridFunction> > ICMap;
      ICMap                      m_ics;

      const IntVect              m_nghosts;
      int                        m_verbosity;
};

#include "NamespaceFooter.H"

#endif
