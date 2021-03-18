#include "BoundaryBoxLayout.H.multidim"

#include <map>

#include "LayoutIterator.H"
#include "LoHiSide.H"
#include "Vector.H"

#include "BCUtils.H"

#include "NamespaceHeader.H"

inline
bool buildBoxLayout( DisjointBoxLayout&            a_bdry_grids,
                     std::map<Box,LayoutIterator>& a_box_map,
                     const DisjointBoxLayout&      a_grids,
                     const MultiBlockCoordSys&     a_coord_sys,
                     const Box&                    a_domain_box,
                     const int&                    a_dir,
                     const Side::LoHiSide&         a_side,
                     const IntVect&                a_nghosts )
{
   Vector<Box> boxes;
   Vector<int> proc_ids;

   for (LayoutIterator lit( a_grids.layoutIterator() ); lit.ok(); ++lit) {

      if (BCUtils::isPhysicalBoundary( a_coord_sys, a_grids[lit], a_dir, a_side )) {

         Box boundary_box( BCUtils::getGhostBox( a_domain_box,
                                                 a_grids[lit],
                                                 a_dir,
                                                 a_side,
                                                 a_nghosts[a_dir]) );

         if (!boundary_box.isEmpty()) {
            boxes.push_back( boundary_box );
            proc_ids.push_back( a_grids.procID( lit() ) );

            typedef std::map<Box,LayoutIterator>::value_type bValType;
            a_box_map.insert( bValType( boundary_box, lit ) );
         }
      }
   }

   Box boundary_domain_box = adjCellBox(a_domain_box,
                                        a_dir,
                                        a_side,
                                        a_nghosts[a_dir]);
   
   ProblemDomain boundary_domain(boundary_domain_box);
   
   if (boxes.size()>0) {
     a_bdry_grids.define( boxes, proc_ids, boundary_domain );
      return true;
   }
   return false;
}


void BoundaryBoxLayout::define( const DisjointBoxLayout&  a_grids,
                                const MultiBlockCoordSys& a_coord_sys,
                                const Box&                a_domain_box,
                                const int&                a_dir,
                                const Side::LoHiSide&     a_side,
                                const IntVect&            a_nghosts )
{
   m_dir = a_dir;
   m_side = a_side;
   m_block = a_coord_sys.whichBlock(a_domain_box);

   std::map<Box,LayoutIterator> box_map;
   m_has_boxes = buildBoxLayout( m_bdry_grids,
                                 box_map,
                                 a_grids,
                                 a_coord_sys,
                                 a_domain_box,
                                 m_dir,
                                 m_side,
                                 a_nghosts );

   if (m_has_boxes) {
      m_box_map.define( m_bdry_grids );
      m_index_map.define( m_bdry_grids );
      for (DataIterator dit( m_bdry_grids.dataIterator() ); dit.ok(); ++dit) {
         LayoutIterator lit( box_map[m_bdry_grids[dit]] );
         m_box_map[dit] = a_grids[lit];
         m_index_map[dit] = DataIndex( lit() );
      }
   }
   m_is_defined = true;
}

#include "NamespaceFooter.H"

