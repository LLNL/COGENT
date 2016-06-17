#ifndef _KINETIC_FUNCTION_UTILS_H_
#define _KINETIC_FUNCTION_UTILS_H_

#include <iostream>
#include <typeinfo>

#include "Box.H"
#include "DataIterator.H"
#include "Directions.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "FourthOrderUtil.H"
#include "LevelData.H"
#include "MultiBlockCoordSys.H"
#include "PhaseBlockCoordSys.H"

#include "NamespaceHeader.H"

namespace KineticFunctionUtils {
   
inline
void convertToCellAverage( const MultiBlockCoordSys&  a_coord_sys,
                           LevelData<FArrayBox>&      a_dfn ) 
{
   const DisjointBoxLayout& grids( a_dfn.disjointBoxLayout() );

   a_dfn.exchange();
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const Box& box( grids[dit] );
      const int block_number( a_coord_sys.whichBlock( box ) );
      const PhaseBlockCoordSys* coord_sys
         = dynamic_cast<const PhaseBlockCoordSys*>( a_coord_sys.getCoordSys( block_number ) );

      fourthOrderAverageCell( a_dfn[dit], coord_sys->domain(), box );
   }
   a_dfn.exchange();
}

};

#include "NamespaceFooter.H"

#endif
