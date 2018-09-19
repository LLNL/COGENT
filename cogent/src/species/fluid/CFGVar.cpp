#include "CFGVar.H"

#include "NamespaceHeader.H"


CFGVar::CFGVar( const string&       a_pp_prefix,
                const std::string&  a_name,
                const MagGeom&      a_geometry )
   : m_pp_prefix(a_pp_prefix),
     m_name(a_name),
     m_geometry(a_geometry)
{
}


const CFGVar& CFGVar::operator=( const CFGVar& a_rhs )
{
   if (&a_rhs != this)
   {
      m_pp_prefix = a_rhs.m_pp_prefix;
      m_name = a_rhs.m_name;
      if ( a_rhs.m_cell_data.isDefined() ) {
         CH_assert( (m_cell_data.disjointBoxLayout()).compatible( a_rhs.m_cell_data.disjointBoxLayout() ) );
         m_cell_data.define( a_rhs.m_cell_data );
      }
      if ( a_rhs.m_face_data.isDefined() ) {
         CH_assert( (m_face_data.disjointBoxLayout()).compatible( a_rhs.m_face_data.disjointBoxLayout() ) );
         m_face_data.define( a_rhs.m_face_data );
      }
   }
   return *this;
}


void CFGVar::copy( const CFGVar& a_rhs )
{
   if (&a_rhs != this)
   {
      m_pp_prefix = a_rhs.m_pp_prefix;
      m_name = a_rhs.m_name;

      if ( a_rhs.m_cell_data.isDefined() ) {
         CH_assert( (m_cell_data.disjointBoxLayout()).compatible( a_rhs.m_cell_data.disjointBoxLayout() ) );
         for (DataIterator dit( m_cell_data.dataIterator() ); dit.ok(); ++dit) {
            m_cell_data[dit].copy( a_rhs.m_cell_data[dit] );
         }
      }
      if ( a_rhs.m_face_data.isDefined() ) {
         CH_assert( (m_face_data.disjointBoxLayout()).compatible( a_rhs.m_face_data.disjointBoxLayout() ) );
         for (DataIterator dit( m_face_data.dataIterator() ); dit.ok(); ++dit) {
            m_face_data[dit].copy( a_rhs.m_face_data[dit] );
         }
      }
   }
}


void CFGVar::zeroData()
{
   if ( m_cell_data.isDefined() ) {
      for (DataIterator dit(m_cell_data.dataIterator()); dit.ok(); ++dit) {
         m_cell_data[dit].setVal( 0.0 );
      }
   }
   if ( m_face_data.isDefined() ) {
      for (DataIterator dit(m_face_data.dataIterator()); dit.ok(); ++dit) {
         m_face_data[dit].setVal( 0.0 );
      }
   }
}


void CFGVar::addData( const CFGVar&  a_rhs,
                      const Real     a_factor )
{
   if ( m_cell_data.isDefined() ) {
      CH_assert( (m_cell_data.disjointBoxLayout()).compatible( a_rhs.m_cell_data.disjointBoxLayout() ) );
      for (DataIterator dit( m_cell_data.dataIterator() ); dit.ok(); ++dit) {
         m_cell_data[dit].plus( a_rhs.m_cell_data[dit], a_factor );
      }
   }
   if ( m_face_data.isDefined() ) {
      CH_assert( (m_face_data.disjointBoxLayout()).compatible( a_rhs.m_face_data.disjointBoxLayout() ) );
      for (DataIterator dit( m_face_data.dataIterator() ); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            m_face_data[dit][dir].plus( a_rhs.m_face_data[dit][dir], a_factor );
         }
      }
   }
}


int CFGVar::size( bool a_count_ghosts )
{
   int size(0);

   if ( m_cell_data.isDefined() ) {
      const DisjointBoxLayout& dbl( m_cell_data.disjointBoxLayout() );
      for (DataIterator dit( dbl.dataIterator() ); dit.ok(); ++dit) {
         const Box& box( a_count_ghosts ? m_cell_data[dit].box() : dbl[dit] );
         size += box.numPts() * m_cell_data.nComp();
      }
   }
   if ( m_face_data.isDefined() ) {
      const DisjointBoxLayout& dbl( m_face_data.disjointBoxLayout() );
      for (DataIterator dit( dbl.dataIterator() ); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            const Box& box( a_count_ghosts ? m_face_data[dit][dir].box() : surroundingNodes(dbl[dit],dir) );
            size += box.numPts() * m_face_data.nComp();
         }
      }
   }

   return size;
}


bool CFGVar::conformsTo( const CFGVar&  a_rhs,
                         const bool     a_include_ghost_cells ) const
{
   bool status( true );

   if ( m_cell_data.isDefined() ) {
      const LevelData<FArrayBox>& thisData( m_cell_data );
      const LevelData<FArrayBox>& rhsData( a_rhs.m_cell_data );
   
      const DisjointBoxLayout& thisBoxes( thisData.disjointBoxLayout() );
      const DisjointBoxLayout& rhsBoxes( rhsData.disjointBoxLayout() );
   
      status &= thisBoxes.compatible( rhsBoxes );
      status &= ( thisData.nComp() == rhsData.nComp() );
   
      if ( a_include_ghost_cells ) {
         status &= ( thisData.ghostVect() == rhsData.ghostVect() );
      }
   }
   if ( m_face_data.isDefined() ) {
      const LevelData<FluxBox>& thisData( m_face_data );
      const LevelData<FluxBox>& rhsData( a_rhs.m_face_data );
   
      const DisjointBoxLayout& thisBoxes( thisData.disjointBoxLayout() );
      const DisjointBoxLayout& rhsBoxes( rhsData.disjointBoxLayout() );
   
      status &= thisBoxes.compatible( rhsBoxes );
      status &= ( thisData.nComp() == rhsData.nComp() );
   
      if ( a_include_ghost_cells ) {
         status &= ( thisData.ghostVect() == rhsData.ghostVect() );
      }
   }
   
   return status;
}


RefCountedPtr<CFGVar>
CFGVar::clone( const IntVect&  a_ghost_vect,
               const bool      a_copy_data ) const
{
   CFGVar* field_ptr = new CFGVar( m_pp_prefix, m_name, m_geometry );

   if ( m_cell_data.isDefined() ) {
      field_ptr->m_cell_data.define( m_cell_data.disjointBoxLayout(),
                                     m_cell_data.nComp(),
                                     a_ghost_vect );

      if (a_copy_data) {
         LevelData<FArrayBox>& result_state( field_ptr->m_cell_data );
         for (DataIterator dit( result_state.dataIterator() ); dit.ok(); ++dit) {
            result_state[dit].copy( m_cell_data[dit] );
         }
      }
   }
   if ( m_face_data.isDefined() ) {
      field_ptr->m_face_data.define( m_face_data.disjointBoxLayout(),
                                     m_face_data.nComp(),
                                     a_ghost_vect );

      if (a_copy_data) {
         LevelData<FluxBox>& result_state( field_ptr->m_face_data );
         for (DataIterator dit( result_state.dataIterator() ); dit.ok(); ++dit) {
            result_state[dit].copy( m_face_data[dit] );
         }
      }
   }

   return RefCountedPtr<CFGVar>(field_ptr);
}


#include "NamespaceFooter.H"

