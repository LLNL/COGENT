#include "SparseCoupling.H"

#include "BaseNamespaceHeader.H"

#include "NamespaceVar.H"

template < > int linearSize( const CH_XDIR::SparseCoupling& a_input )
{
   return a_input.linearSize();
}

template < > void linearIn( CH_XDIR::SparseCoupling& a_outputT, const void* const a_inBuf )
{
   a_outputT.linearIn( a_inBuf );
}

template < > void linearOut( void* const a_outBuf, const CH_XDIR::SparseCoupling& a_inputT )
{
   a_inputT.linearOut( a_outBuf );
}

#include "BaseNamespaceFooter.H"
#include "NamespaceHeader.H"

int SparseCoupling::linearSize() const
{
   return ( linearListSize( m_index ) + linearListSize( m_weight ) );
}

void SparseCoupling::linearIn( const void* const a_inBuf )
{
   const char* const buf = (char*)a_inBuf;
   linearListIn( m_index, buf );
   linearListIn( m_weight, buf + linearListSize( m_index ) );
}

void SparseCoupling::linearOut( void* const a_outBuf ) const
{
   char* const buf = (char*)a_outBuf;
   linearListOut( buf, m_index );
   linearListOut( buf + linearListSize( m_index ), m_weight );
}


template < >
int BaseFab<SparseCoupling>::preAllocatable()
{
   return 2; // dyanmic allocatable.
}

template < >
size_t BaseFab<SparseCoupling>::size(const Box&      a_R,
                                     const Interval& a_comps) const
{
  int totalsize(0);
  ForAllThisCBNN(SparseCoupling,a_R,a_comps.begin(),a_comps.size())
  {
    totalsize += thisR.linearSize();
  } EndFor;
  return totalsize;
}

template < >
void BaseFab<SparseCoupling>::linearOut(void*           a_buf,
                                        const Box&      a_R,
                                        const Interval& a_comps) const
{
  char* buffer = (char*)a_buf;
  ForAllThisCBNN(SparseCoupling,a_R,a_comps.begin(),a_comps.size())
  {
    thisR.linearOut( buffer );
    buffer += thisR.linearSize();
  } EndFor;
}

template < >
void BaseFab<SparseCoupling>::linearIn(void*           a_buf,
                                       const Box&      a_R,
                                       const Interval& a_comps)
{
  char* buffer = (char*)a_buf;
  ForAllThisBNN(SparseCoupling,a_R,a_comps.begin(),a_comps.size())
  {
    thisR.linearIn( buffer );
    buffer += thisR.linearSize();
  } EndFor;
}


std::ostream& operator<<( std::ostream& a_out, const SparseCoupling& a_spc )
{
   a_out << "[";
   for (int i(0); i<a_spc.m_index.size(); ++i) {
      a_out << "{" << a_spc.m_index[i] << "," << a_spc.m_weight[i] << "}";
      if (i<a_spc.m_index.size()-1) {
         a_out << ", ";
      }
   }
   a_out << "]";
   return a_out;
}

#include "NamespaceFooter.H"
