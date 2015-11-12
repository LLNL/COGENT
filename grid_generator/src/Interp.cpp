#include "Interp.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"



Interp::Interp( const FArrayBox& a_node_coords,
                const FArrayBox& a_data )
  : m_x(NULL),
    m_y(NULL)
{
   m_data.define(a_data.box(),a_data.nComp());
   m_data.copy(a_data);

   const Box& box = m_data.box();
   CH_assert(box == a_node_coords.box());

   //   m_x = new double[dim(0)];
   m_x = new double[m_data.box().size(0)];
   Box box_x = adjCellLo(box, 1, -1);
   int i = 0;
   BoxIterator bit_x(box_x);
   for (bit_x.begin();bit_x.ok();++bit_x) {
      IntVect iv = bit_x();
      m_x[i] = a_node_coords(iv, 0);
      i++;
   }

   //   m_y = new double[dim(1)];
   m_y = new double[m_data.box().size(1)];
   Box box_y = adjCellLo(box, 0, -1);
   int j = 0;
   BoxIterator bit_y(box_y);
   for (bit_y.begin();bit_y.ok();++bit_y) {
      IntVect iv = bit_y();
      m_y[j] = a_node_coords(iv, 1);
      j++;
   }
}



Interp::~Interp()
{
  if (m_x) delete [] m_x;
  if (m_y) delete [] m_y;
}



#include "NamespaceFooter.H"


