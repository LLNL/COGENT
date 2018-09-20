//#ifndef _BC_UTILS_H_
//#define _BC_UTILS_H_

#include "Box.H"
#include "LoHiSide.H"
#include "BlockBoundary.H"
#include "Tuple.H"
#include "Vector.H"
#include "MultiBlockCoordSys.H"

#include "NamespaceHeader.H"

namespace BCUtils {

   inline
   Box getGhostBox( const Box& a_domain_box,
                    const Box& a_box,
                    const int& a_dir,
                    const Side::LoHiSide& a_side,
                    const int& a_num_ghosts )
   {
      const Box gBox( adjCellBox( a_box, a_dir, a_side, a_num_ghosts ) );
      const Box gDomainBox( adjCellBox( a_domain_box, a_dir, a_side, a_num_ghosts ) );
      return (gBox & gDomainBox);
   }


   inline
   const Box& blockDomain( const Box& a_box,
                           const MultiBlockCoordSys& a_coord_sys )
   {
      const Vector<Box>& blocks( a_coord_sys.mappingBlocks() );
      const int block( a_coord_sys.whichBlock( a_box ) );
      return blocks[block];
   }


   inline
   bool isPhysicalBoundary( const MultiBlockCoordSys& a_coord_sys,
                            const Box&                a_box,
                            const int&                a_dir,
                            const Side::LoHiSide&     a_side )
   {
      const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries( a_coord_sys.boundaries() );
      const int block( a_coord_sys.whichBlock( a_box ) );
      const int index( a_side * SpaceDim + a_dir );
      return boundaries[block][index].isDomainBoundary();
   }

   inline
   bool touchesPhysicalBoundary( const MultiBlockCoordSys& a_coord_sys,
                                 const Box&                a_box,
                                 const int&                a_dir,
                                 const Side::LoHiSide&     a_side )
   {
      if (BCUtils::isPhysicalBoundary( a_coord_sys, a_box, a_dir, a_side ) ) {
         const Box& block_box( BCUtils::blockDomain( a_box, a_coord_sys ) );
         const Box g_box( adjCellBox( a_box, a_dir, a_side, 1 ) );
         const Box g_block_box( adjCellBox( block_box, a_dir, a_side, 1 ) );
         return g_box.intersects( g_block_box );
      }
      return false;
   }

}

#include "NamespaceFooter.H"

//#endif
