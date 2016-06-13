#ifndef _FLUID_SPECIES_BC_H_
#define _FLUID_SPECIES_BC_H_

#include "Dimensions.H"
#include "REAL.H"
#include "LevelData.H"
#include "FluxBox.H"
#include "FArrayBox.H"
#include "FluidSpecies.H"

#include "NamespaceHeader.H"

class FluidSpeciesBC
{
   public:

      virtual void apply( FluidSpecies& a_species_phys,
                          const LevelData<FArrayBox>& a_phi,
                          const LevelData<FluxBox>& a_velocity,
                          const Real& a_time ) = 0;

      virtual bool isForVariable( const std::string& name ) const = 0;
};

#include "NamespaceFooter.H"

#endif
