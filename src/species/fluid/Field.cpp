#include "Field.H"

#include "NamespaceHeader.H"

Field::Field( const std::string& a_name, const MagGeom& a_geometry )
   : m_geometry(a_geometry),
     m_name(a_name)
{
}

const Field& Field::operator=( const Field& a_rhs )
{
   if (&a_rhs != this)
   {
      CH_assert( (m_data.disjointBoxLayout()).compatible( a_rhs.m_data.disjointBoxLayout() ) );
      m_name = a_rhs.m_name;
      m_data.define( a_rhs.m_data );
   }
   return *this;
}


void Field::copy( const Field& a_rhs )
{
   if (&a_rhs != this)
   {
      CH_assert( (m_data.disjointBoxLayout()).compatible( a_rhs.m_data.disjointBoxLayout() ) );
      m_name = a_rhs.m_name;
      for (DataIterator dit( m_data.dataIterator() ); dit.ok(); ++dit) {
         m_data[dit].copy( a_rhs.m_data[dit] );
      }
   }
}


void Field::zeroData()
{
   DataIterator dit( m_data.dataIterator() );
   for (dit.begin(); dit.ok(); ++dit) {
      m_data[dit].setVal( 0.0 );
   }
}


void Field::addData( const Field& a_rhs,
                            const Real a_factor )
{
   CH_assert( (m_data.disjointBoxLayout()).compatible( a_rhs.m_data.disjointBoxLayout() ) );
   for (DataIterator dit( m_data.dataIterator() ); dit.ok(); ++dit) {
      m_data[dit].plus( a_rhs.m_data[dit], a_factor );
   }
}


bool Field::conformsTo( const Field& a_rhs,
                               const bool a_include_ghost_cells ) const
{
   const LevelData<FArrayBox>& thisData( m_data );
   const LevelData<FArrayBox>& rhsData( a_rhs.m_data );
   
   const DisjointBoxLayout& thisBoxes( thisData.disjointBoxLayout() );
   const DisjointBoxLayout& rhsBoxes( rhsData.disjointBoxLayout() );
   
   bool status( true );
   status &= thisBoxes.compatible( rhsBoxes );
   status &= ( thisData.nComp() == rhsData.nComp() );
   
   if ( a_include_ghost_cells ) {
      status &= ( thisData.ghostVect() == rhsData.ghostVect() );
   }
   
   return status;
}


RefCountedPtr<Field>
Field::clone( const IntVect a_ghost_vect, const bool a_copy_soln_data ) const
{
   RefCountedPtr<Field> result
      = RefCountedPtr<Field>(
         new Field( m_name, m_geometry ) );

   result->m_data.define( m_data.disjointBoxLayout(),
                          m_data.nComp(),
                          a_ghost_vect );

   if (a_copy_soln_data) {
      LevelData<FArrayBox>& result_state( result->m_data );
      for (DataIterator dit( result_state.dataIterator() ); dit.ok(); ++dit) {
         result_state[dit].copy( m_data[dit] );
      }
   }

   return result;
}


#include "NamespaceFooter.H"

