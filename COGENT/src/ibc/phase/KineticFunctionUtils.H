#ifndef _KINETIC_FUNCTION_UTILS_H_
#define _KINETIC_FUNCTION_UTILS_H_

#include "FourthOrderUtil.H"
#include "FourthOrderUtilsClass.H.multidim"

#include "NamespaceHeader.H"

namespace KineticFunctionUtils {
   
inline
void convertToCellAverage( const PhaseGeom&       a_geometry,
                           LevelData<FArrayBox>&  a_dfn, 
						   const bool useSG=false ) 
{
   if ( !(a_dfn.ghostVect() >= IntVect::Unit) ) {
      MayDay::Error("convertToCellAverage(): Insufficient number of ghost cells");
   }
   FourthOrderUtil FourthOrderOperators;
   FourthOrderOperators.setSG(useSG);

   //fourthOrderAverage(a_dfn); // old Chombo-only version
   FourthOrderOperators.fourthOrderAverageGen(a_dfn);
}

}

#include "NamespaceFooter.H"

#endif
