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
DataArray::updateData( const FArrayBox& a_data )
{
   // Copy on overlap
   m_data.copy(a_data);
}   


void DataArray::setPointwise(FArrayBox&                 a_data,
                             const MultiBlockLevelGeom& a_geometry,
                             const FArrayBox&           a_real_coords,
                             const FArrayBox&           a_normalized_flux,
                             const int                  a_block_number) const
{
   a_data.copy(m_data);
}

void DataArray::printParameters() const
{
   if (procID()==0) {
      std::cout << "DataArray grid function parameters:" << std::endl;
   }
}

#include "NamespaceFooter.H"
