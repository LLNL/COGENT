#include "Interp3D.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"


Interp3D::Interp3D( const FArrayBox&  a_nodes,
                    const FArrayBox&  a_data )
  : m_x(NULL),
    m_y(NULL),
    m_z(NULL)
{
   m_data.define(a_data.box(),a_data.nComp());
   m_data.copy(a_data);

   const Box& box = m_data.box();
   CH_assert(box == a_nodes.box());

   m_x = new double[dim(0)];
   Box box_x = adjCellLo(adjCellLo(box, 1, -1), 2, -1);
   int i = 0;
   for (BoxIterator bit(box_x); bit.ok(); ++bit) {
      m_x[i] = a_nodes(bit(), 0);
      i++;
   }

   m_y = new double[dim(1)];
   Box box_y = adjCellLo(adjCellLo(box, 0, -1), 2, -1);
   int j = 0;
   for (BoxIterator bit(box_y); bit.ok(); ++bit) {
      m_y[j] = a_nodes(bit(), 1);
      j++;
   }

   m_z = new double[dim(2)];
   Box box_z = adjCellLo(adjCellLo(box, 0, -1), 1, -1);
   int k = 0;
   for (BoxIterator bit(box_z); bit.ok(); ++bit) {
      m_z[k] = a_nodes(bit(), 2);
      k++;
   }
}


Interp3D::~Interp3D()
{
  if (m_x) delete [] m_x;
  if (m_y) delete [] m_y;
  if (m_z) delete [] m_z;
}


#include "NamespaceFooter.H"


