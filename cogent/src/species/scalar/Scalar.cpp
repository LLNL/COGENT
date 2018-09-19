#include "Scalar.H"

#include "NamespaceHeader.H"

Scalar::Scalar( const std::string& a_name,
                const int          a_num_scalars )
   : m_name(a_name)
{
   m_data.resize(a_num_scalars);
}

const Scalar& Scalar::operator=( const Scalar& a_rhs )
{
   if (&a_rhs != this)
   {
      m_name = a_rhs.m_name;
      m_data.resize(a_rhs.m_data.size());
      for (int i=0; i<a_rhs.m_data.size(); ++i) {
         m_data[i] = a_rhs.m_data[i];
      }
   }
   return *this;
}


void Scalar::copy( const Scalar& a_rhs )
{
   if (&a_rhs != this)
   {
      m_name = a_rhs.m_name;
      m_data.resize(a_rhs.m_data.size());
      for (int i=0; i<a_rhs.m_data.size(); ++i) {
         m_data[i] = a_rhs.m_data[i];
      }
   }
}


void Scalar::zeroData()
{
   for (int i=0; i<m_data.size(); ++i) {
      m_data[i] = 0.;
   }
}


void Scalar::addData( const Scalar& a_rhs,
                      const Real    a_factor )
{
   CH_assert(m_data.size()==a_rhs.m_data.size());
   for (int i=0; i<m_data.size(); ++i) {
      m_data[i] += a_factor * a_rhs.m_data[i];
   }
}


bool Scalar::conformsTo( const Scalar& a_rhs,
                         const bool a_include_ghost_cells ) const
{
   return m_data.size() == a_rhs.m_data.size();
}


RefCountedPtr<Scalar>
Scalar::clone( const bool a_copy_soln_data ) const
{
   RefCountedPtr<Scalar> result
      = RefCountedPtr<Scalar>(
        new Scalar( m_name, m_data.size() ) );

   if (a_copy_soln_data) {
      result->copy(*this);
   }

   return result;
}


#include "NamespaceFooter.H"

