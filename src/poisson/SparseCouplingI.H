#ifndef _SPARSE_COUPLING_I_H_
#define _SPARSE_COUPLING_I_H_

inline
SparseCoupling::SparseCoupling() {}


inline
SparseCoupling::~SparseCoupling() {}


inline
void SparseCoupling::add( const IntVect& a_index, const Real& a_weight )
{
   const int location( find( a_index ) );
   if (location == m_index.size())
   {
      m_index.push_back( a_index );
      m_weight.push_back( a_weight );
   }
   else
   {
      m_weight[location] += a_weight;
   }
}


inline
int SparseCoupling::size() const
{
   return m_index.size();
}


inline
const IntVect& SparseCoupling::operator[]( const int a_i ) const
{
   CH_assert(a_i>=0);
   CH_assert(a_i<m_index.size());
   return m_index[a_i];
}


inline
const Real& SparseCoupling::weight( const IntVect& a_index ) const
{
   const int location( find( a_index ) );
   CH_assert(location != m_index.size());
   return m_weight[location];
}


inline
const Real& SparseCoupling::weight( const int& a_i ) const
{
   CH_assert(a_i>=0);
   CH_assert(a_i<m_index.size());
   return m_weight[a_i];
}


inline
SparseCouplingIterator SparseCoupling::iterator() const
{
   return SparseCouplingIterator( *this );
}


inline
SparseCoupling::SparseCoupling( const SparseCoupling& a_stencil )
{
   m_index = a_stencil.m_index;
   m_weight = a_stencil.m_weight;
}


inline
const SparseCoupling& SparseCoupling::operator=( const SparseCoupling& a_stencil )
{
   if (*this!=a_stencil) {
      m_index = a_stencil.m_index;
      m_weight = a_stencil.m_weight;
   }
   return (*this);
}


inline
bool SparseCoupling::operator!=( const SparseCoupling& a_stencil ) const
{
   if (m_index.size()!=a_stencil.m_index.size()) {
      return true;
   }
   
   for (int i(0); i<m_index.size(); i++) {
      if (m_index[i]!=a_stencil.m_index[i]) {
         return true;
      }
      if (m_weight[i]!=a_stencil.m_weight[i]) {
         return true;
      }
   }
   return false;
}


inline
bool SparseCoupling::operator==( const SparseCoupling& a_stencil ) const
{
   return !(*this!=a_stencil);
}

   
inline
   void SparseCoupling::print( std::ostream& a_out,
                             const std::string& a_prefix) const
{
   for (int i(0); i<m_index.size(); ++i) {
      a_out << a_prefix
            << m_index[i] << ":  " << m_weight[i]
            << std::endl;
   }
}


inline
int SparseCoupling::find( const IntVect& a_index ) const
{
   int location(0);
   for (; location<m_index.size(); location++)
   {
      if (m_index[location]==a_index) {
         break;
      }
   }
   return location;
}


///////////////////////////////////////////////////////////////

inline
SparseCouplingIterator::SparseCouplingIterator()
   : m_current(-1), m_lo(1), m_hi(0)
{
}


inline
SparseCouplingIterator::SparseCouplingIterator( const SparseCoupling& a_stencil )
{
   define( a_stencil );
}


inline
SparseCouplingIterator::SparseCouplingIterator(
   const SparseCouplingIterator& a_iterator )
{
   m_current = a_iterator.m_current;
   m_lo =   a_iterator.m_lo;
   m_hi =   a_iterator.m_hi;
}


inline
void SparseCouplingIterator::define( const SparseCoupling& a_stencil )
{
   m_lo = 0;
   m_hi = a_stencil.size() - 1;
   m_current = m_lo;
}


inline
void SparseCouplingIterator::begin()
{
   if (m_lo<=m_hi)
      m_current = m_lo;
}


inline
void SparseCouplingIterator::reset()
{
   begin();
}


inline
void SparseCouplingIterator::operator++()
{
   next();
}


inline
void SparseCouplingIterator::next()
{
   m_current++;
}


inline
const int& SparseCouplingIterator::operator()() const
{
   CH_assert(m_current<=m_hi);
   CH_assert(m_current>=m_lo);
   return m_current;
}


inline
bool SparseCouplingIterator::ok()
{
   return (m_current<=m_hi);
}


#endif
