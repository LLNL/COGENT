#include "DataArray.H"

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


#include "NamespaceFooter.H"
