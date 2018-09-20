#include "DataArray.H"
#include "FourthOrderUtil.H"


#include "NamespaceHeader.H"


void
DataArray::setData( const FArrayBox& a_data,
                    const bool&      a_cell_averages )
{
   m_data.define(a_data.box(), a_data.nComp());
   m_data.copy(a_data);
   m_cell_averages = a_cell_averages;
}   


void
DataArray::setData( const DataArray& a_data )
{
   m_data.copy(a_data.m_data);
   m_cell_averages = a_data.m_cell_averages;
}   


void
DataArray::assign( LevelData<FArrayBox>&      a_data,
                   const MultiBlockLevelGeom& a_geometry,
                   const Real&                a_time,
                   const bool&                a_cell_averages ) const
{
   CH_assert(a_data.nComp() <= m_data.nComp());
   CH_assert(a_cell_averages == m_cell_averages);

   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      a_data[dit].copy(m_data);
   }
}
   

void
DataArray::assign( LevelData<FluxBox>&        a_data,
                   const MultiBlockLevelGeom& a_geometry,
                   const Real&                a_time,
                   const bool&                a_cell_averages ) const
{
   MayDay::Error("DataArray::assign() not implemented");
}
   

void
DataArray::assign( FArrayBox&                 a_data,
                   const MultiBlockLevelGeom& a_geometry,
                   const Box&                 a_box,
                   const Real&                a_time,
                   const bool&                a_cell_averages ) const
{
   CH_assert(a_data.nComp() <= m_data.nComp());
   CH_assert(a_cell_averages == m_cell_averages);

   if (a_cell_averages) {
      MayDay::Error("DataArray::assign() is only used for cell centered data");
   }
   else {
      a_data.copy(m_data);  // Copies on the box overlap
   }
}


void
DataArray::assign( FluxBox&                   a_data,
                   const MultiBlockLevelGeom& a_geometry,
                   const Box&                 a_box,
                   const Real&                a_time,
                   const bool&                a_cell_averages ) const
{
   MayDay::Error("DataArray::assign() not implemented");
}


void
DataArray::assign( LevelData<FArrayBox>&      a_data,
                   const MultiBlockLevelGeom& a_geometry,
                   const BoundaryBoxLayout&   a_bdry_layout,
                   const Real&                a_time ) const
{
   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   // NB: This is a cheat - there's one too many cells at the (dir,side) face
   // of the boundary box, but it doesn't matter because one-sided difference
   // will be used at that face to construct the cell average.  We do want the
   // extra cell in all other directions.
   LevelData<FArrayBox> data_tmp( grids, a_data.nComp(), IntVect::Unit );
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      data_tmp[dit].copy(m_data);
   }
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      Box domain_box( data_tmp[dit].box() );
      domain_box.growDir( a_bdry_layout.dir(), a_bdry_layout.side(), -1 );
      ProblemDomain domain( domain_box );
      fourthOrderAverageCell( data_tmp[dit], domain, grids[dit] );
   }
   data_tmp.copyTo( a_data );
   a_data.exchange();
}


#include "NamespaceFooter.H"
