#ifndef _SN_SYSTEM_BC_
#define _SN_SYSTEM_BC_

#include "FluxBox.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "MultiBlockLevelExchangeAverage.H"
#include "MagGeom.H"
#include "ParmParse.H"
//#include "SNBC.H"
#include "REAL.H"
#include "FluidSpecies.H"
#include "FluidSpeciesBC.H"
#include "SNTimeIntegration.H"

#include <string>
#include <vector>

#include "NamespaceHeader.H"

class SNSystemBC
{
   public:

      SNSystemBC( ParmParse& pp, const SNState& a_state );

      ~SNSystemBC();

      // fills ghost cells on the physical domain
      void fillFluidSpeciesGhostCells( FluidSpeciesPtrVect& a_species,
                                       const Real& a_time ) const;
   
   private:

      // prevent copying
      SNSystemBC( const SNSystemBC& );
      const SNSystemBC& operator=( const SNSystemBC& );

      void executeInternalExchanges( FluidSpecies& a_fluid ) const;

      void parseFluidSpecies( ParmParse& a_pp,
                              const std::string& a_coord_sys_type,
                              const FluidSpeciesPtrVect& a_solution  );
 
      FluidSpeciesBC& fluidSpeciesBC( const std::string& a_name ) const;
  
      int                        m_verbosity;
      const MagGeom&             m_geometry;
      std::vector<FluidSpeciesBC*> m_bcs;
      MultiBlockLevelExchange*   m_mblexPtr;
};

#include "NamespaceFooter.H"

#endif
