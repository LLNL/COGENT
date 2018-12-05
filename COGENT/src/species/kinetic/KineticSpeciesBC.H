#ifndef _KINETIC_SPECIES_BC_H_
#define _KINETIC_SPECIES_BC_H_

#define CH_SPACEDIM PDIM

#include "Dimensions.H"
#include "REAL.H"
#include "LevelData.H"
#include "FluxBox.H"
#include "KineticSpecies.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "LevelData.H"
#include "FArrayBox.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"
namespace CFG = CFG_NAMESPACE;

class KineticSpeciesBC
{
   public:

      virtual ~KineticSpeciesBC() {;}

      virtual void apply( KineticSpecies& a_species_phys,
                          const CFG::LevelData<CFG::FArrayBox>& a_phi,
                          const LevelData<FluxBox>& a_velocity,
                          const Real& a_time ) = 0;

      virtual bool isForVariable( const std::string& name ) const = 0;
};

#include "NamespaceFooter.H"

#endif
