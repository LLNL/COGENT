//#ifndef _BOUNDARY_BOX_DATA_
//#define _BOUNDARY_BOX_DATA_

#include "Box.H"
#include "IntVect.H"
#include "DataIndex.H"
#include "DataIterator.H"
#include "DisjointBoxLayout.H"
#include "LayoutData.H"
#include "LoHiSide.H"
#include "MultiBlockCoordSys.H"

#include "NamespaceHeader.H"

class BoundaryBoxLayout
{
   public:
      BoundaryBoxLayout();

      ~BoundaryBoxLayout() {;}

      BoundaryBoxLayout( const DisjointBoxLayout&  a_grids,
                         const MultiBlockCoordSys& a_coord_sys,
                         const Box&                a_domain_box,
                         const int&                a_dir,
                         const Side::LoHiSide&     a_side,
                         const IntVect&            a_nghosts = IntVect::Zero );

      void define( const DisjointBoxLayout&  a_grids,
                   const MultiBlockCoordSys& a_coord_sys,
                   const Box&                a_domain_box,
                   const int&                a_dir,
                   const Side::LoHiSide&     a_side,
                   const IntVect&            a_nghosts = IntVect::Zero );

      const DisjointBoxLayout& disjointBoxLayout() const;

      const Box& interiorBox( const DataIterator& a_dit ) const;

      const DataIndex& dataIndex( const DataIterator& a_dit ) const;

      const int& dir() const;

      const Side::LoHiSide& side() const;

      const int& block() const;
   
      bool hasBoxes() const;

   private:

      int                   m_dir;
      Side::LoHiSide        m_side;
      int                   m_block;
      bool                  m_has_boxes;
      bool                  m_is_defined;
      DisjointBoxLayout     m_bdry_grids;
      LayoutData<Box>       m_box_map;
      LayoutData<DataIndex> m_index_map;
};


typedef RefCountedPtr<BoundaryBoxLayout> BoundaryBoxLayoutPtr;

typedef Vector<BoundaryBoxLayoutPtr> BoundaryBoxLayoutPtrVect;


inline
BoundaryBoxLayout::BoundaryBoxLayout()
   : m_dir(-1),
     m_side(Side::Invalid),
     m_block(-1),
     m_has_boxes(false),
     m_is_defined(false)
{
}


inline
BoundaryBoxLayout::BoundaryBoxLayout( const DisjointBoxLayout&  a_grids,
                                      const MultiBlockCoordSys& a_coord_sys,
                                      const Box&                a_domain_box,
                                      const int&                a_dir,
                                      const Side::LoHiSide&     a_side,
                                      const IntVect&            a_nghosts )
   : m_dir(-1),
     m_side(Side::Invalid),
     m_block(-1),
     m_has_boxes(false),
     m_is_defined(false)
{
   define( a_grids, a_coord_sys, a_domain_box, a_dir, a_side, a_nghosts );
}


inline
const DisjointBoxLayout& BoundaryBoxLayout::disjointBoxLayout() const
{
   CH_assert(m_is_defined);
   CH_assert(m_has_boxes);
   return m_bdry_grids;
}


inline
const Box& BoundaryBoxLayout::interiorBox( const DataIterator& a_dit ) const
{
   CH_assert(m_is_defined);
   CH_assert(m_has_boxes);
   return m_box_map[a_dit];
}


inline
const DataIndex& BoundaryBoxLayout::dataIndex( const DataIterator& a_dit ) const
{
   CH_assert(m_is_defined);
   CH_assert(m_has_boxes);
   return m_index_map[a_dit];
}

inline
const int& BoundaryBoxLayout::dir() const
{
   CH_assert(m_is_defined);
   return m_dir;
}

inline
const Side::LoHiSide& BoundaryBoxLayout::side() const
{
   CH_assert(m_is_defined);
   return m_side;
}

inline
const int& BoundaryBoxLayout::block() const
{
   CH_assert(m_is_defined);
   return m_block;
}

inline
bool BoundaryBoxLayout::hasBoxes() const
{
   return m_has_boxes;
}

#include "NamespaceFooter.H"

//#endif
