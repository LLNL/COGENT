#ifndef _SN_BC_H_
#define _SN_BC_H_

#include "REAL.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "MagGeom.H"
#include "FluxBox.H"

#include "NamespaceHeader.H"

class SNBC
{
   public:

      virtual void apply( LevelData<FArrayBox>&     a_species_phys,
                          const MagGeom&            a_geometry,
                          const LevelData<FluxBox>& a_velocity,
                          const Real&               a_time ) = 0;
};

#include "NamespaceFooter.H"

#endif
