#ifndef _SN_SYSTEM_BC_
#define _SN_SYSTEM_BC_

#include "FluxBox.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "MultiBlockLevelExchangeAverage.H"
#include "MagGeom.H"
#include "ParmParse.H"
#include "SNBC.H"
#include "REAL.H"

#include <string>
#include <vector>

#include "NamespaceHeader.H"

class SNSystemBC
{
   public:

      SNSystemBC( ParmParse& pp,
                  const MagGeom& geometry,
                  const LevelData<FArrayBox>& solution );

      ~SNSystemBC();

      // fills ghost cells on the physical domain
      void fillGhostCells(
         LevelData<FArrayBox>& a_solution_phys,
         const LevelData<FArrayBox>& a_solution_comp,
         const LevelData<FluxBox>& a_velocity,
         const Real& a_time ) const;

   private:

      // prevent copying
      SNSystemBC( const SNSystemBC& );
      const SNSystemBC& operator=( const SNSystemBC& );

      void executeInternalExchanges( LevelData<FArrayBox>& a_solution ) const;

      void fillBlockBoundaryCellsInCompSpace(
         LevelData<FArrayBox>& solution_comp,
         const LevelData<FArrayBox>& solution_phys ) const;

      void fillDomainGhostCellsInPhysSpace(
         LevelData<FArrayBox>& a_solution_phys,
         const LevelData<FArrayBox>& a_solution_comp ) const;

      void parseInput( ParmParse& a_pp,
                       const std::string& a_coord_sys_type,
                       const LevelData<FArrayBox>& a_solution  );
   
      int                        m_verbosity;
      const MagGeom&             m_geometry;
      SNBC*                      m_bcs;
      MultiBlockLevelExchange*   m_mblexPtr;
};

#include "NamespaceFooter.H"

#endif
