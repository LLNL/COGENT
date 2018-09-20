#ifndef _KINETIC_FUNCTION_UTILS_H_
#define _KINETIC_FUNCTION_UTILS_H_

#include "FourthOrderUtil.H"

#include "NamespaceHeader.H"

namespace KineticFunctionUtils {
   
inline
void convertToCellAverage( const PhaseGeom&       a_geometry,
                           LevelData<FArrayBox>&  a_dfn ) 
{
   if ( !(a_dfn.ghostVect() >= IntVect::Unit) ) {
      MayDay::Error("convertToCellAverage(): Insufficient number of ghost cells");
   }

   fourthOrderAverage(a_dfn);
}

}

#include "NamespaceFooter.H"

#endif
