#include "CFGVars.H"

#include "NamespaceHeader.H"


CFGVars::CFGVars( const string&       a_pp_prefix,
                  const std::string&  a_name,
                  const MagGeom&      a_geometry )
   : m_pp_prefix(a_pp_prefix),
     m_name(a_name),
     m_geometry(a_geometry)
{
}


const CFGVars& CFGVars::operator=( const CFGVars& a_rhs )
{
   if (&a_rhs != this)
   {
      m_pp_prefix = a_rhs.m_pp_prefix;
      m_name = a_rhs.m_name;

      m_cell_data.resize(0);
      for (int i=0; i<a_rhs.num_cell_vars(); ++i) {
         const LevelData<FArrayBox>& data = a_rhs.cell_var(i);
         addCellVar(cell_var_name(i), data.nComp(), data.ghostVect());
         m_cell_data[i]->define(data);
      }

      m_face_data.resize(0);
      for (int i=0; i<a_rhs.num_face_vars(); ++i) {
         const LevelData<FluxBox>& data = a_rhs.face_var(i);
         addFaceVar(face_var_name(i), data.nComp(), data.ghostVect());
         m_face_data[i]->define(data);
      }
   }
   return *this;
}


void CFGVars::copy( const CFGVars& a_rhs )
{
   if ( &a_rhs != this ) {
      m_pp_prefix = a_rhs.m_pp_prefix;
      m_name = a_rhs.m_name;

      int num_cell_data = a_rhs.num_cell_vars();
      CH_assert(num_cell_data == num_cell_vars());
      for (int i=0; i<num_cell_data; ++i) {
         LevelData<FArrayBox>& this_cell_data = cell_var(i);
         const LevelData<FArrayBox>& that_cell_data = a_rhs.cell_var(i);
         CH_assert( (this_cell_data.disjointBoxLayout()).compatible( that_cell_data.disjointBoxLayout() ) );
         for (DataIterator dit( this_cell_data.dataIterator() ); dit.ok(); ++dit) {
            this_cell_data[dit].copy(that_cell_data[dit] );
         }
      }

      int num_face_data = a_rhs.num_face_vars();
      CH_assert(num_face_data == num_face_vars());
      for (int i=0; i<num_face_data; ++i) {
         LevelData<FluxBox>& this_face_data = face_var(i);
         const LevelData<FluxBox>& that_face_data = a_rhs.face_var(i);
         CH_assert( (this_face_data.disjointBoxLayout()).compatible( that_face_data.disjointBoxLayout() ) );
         for (DataIterator dit( this_face_data.dataIterator() ); dit.ok(); ++dit) {
            this_face_data[dit].copy(that_face_data[dit] );
         }
      }
   }
}


void CFGVars::zeroData()
{
   for (int i=0; i<num_cell_vars(); ++i) {
      LevelData<FArrayBox>& this_cell_data = cell_var(i);
      for (DataIterator dit( this_cell_data.dataIterator() ); dit.ok(); ++dit) {
         this_cell_data[dit].setVal(0.);
      }
   }

   for (int i=0; i<num_face_vars(); ++i) {
      LevelData<FluxBox>& this_face_data = face_var(i);
      for (DataIterator dit( this_face_data.dataIterator() ); dit.ok(); ++dit) {
         this_face_data[dit].setVal(0.);
      }
   }
}


void CFGVars::addData( const CFGVars&  a_rhs,
                       const Real      a_factor )
{
   int num_cell_data = a_rhs.num_cell_vars();
   CH_assert(num_cell_data == num_cell_vars());
   for (int i=0; i<num_cell_data; ++i) {
      LevelData<FArrayBox>& this_cell_data = cell_var(i);
      const LevelData<FArrayBox>& that_cell_data = a_rhs.cell_var(i);
      CH_assert( (this_cell_data.disjointBoxLayout()).compatible( that_cell_data.disjointBoxLayout() ) );
      for (DataIterator dit( this_cell_data.dataIterator() ); dit.ok(); ++dit) {
         this_cell_data[dit].plus(that_cell_data[dit], a_factor);
      }
   }

   int num_face_data = a_rhs.num_face_vars();
   CH_assert(num_face_data == num_face_vars());
   for (int i=0; i<num_face_data; ++i) {
      LevelData<FluxBox>& this_face_data = face_var(i);
      const LevelData<FluxBox>& that_face_data = a_rhs.face_var(i);
      CH_assert( (this_face_data.disjointBoxLayout()).compatible( that_face_data.disjointBoxLayout() ) );
      for (DataIterator dit( this_face_data.dataIterator() ); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            this_face_data[dit][dir].plus(that_face_data[dit][dir], a_factor);
         }
      }
   }
}


int CFGVars::size( bool a_count_ghosts )
{
   int size(0);

   for (int i=0; i<num_cell_vars(); ++i) {
      LevelData<FArrayBox>& this_cell_data = cell_var(i);
      const DisjointBoxLayout& dbl( this_cell_data.disjointBoxLayout() );
      for (DataIterator dit( this_cell_data.dataIterator() ); dit.ok(); ++dit) {
         const Box& box( a_count_ghosts ? this_cell_data[dit].box() : dbl[dit] );
         size += box.numPts() * this_cell_data.nComp();
      }
   }

   for (int i=0; i<num_face_vars(); ++i) {
      LevelData<FluxBox>& this_face_data = face_var(i);
      const DisjointBoxLayout& dbl( this_face_data.disjointBoxLayout() );
      for (DataIterator dit( this_face_data.dataIterator() ); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            const Box& box( a_count_ghosts ? this_face_data[dit][dir].box() : surroundingNodes(dbl[dit],dir) );
            size += box.numPts() * this_face_data.nComp();
         }
      }
   }

   return size;
}


bool CFGVars::conformsTo( const CFGVars&  a_rhs,
                          const bool      a_include_ghost_cells ) const
{
   bool status( true );

   for (int i=0; i<num_cell_vars(); ++i) {
      const LevelData<FArrayBox>& this_cell_data = cell_var(i);
      const LevelData<FArrayBox>& that_cell_data = a_rhs.cell_var(i);

      const DisjointBoxLayout& thisBoxes( this_cell_data.disjointBoxLayout() );
      const DisjointBoxLayout& rhsBoxes( that_cell_data.disjointBoxLayout() );
   
      status &= thisBoxes.compatible( rhsBoxes );
      status &= ( this_cell_data.nComp() == that_cell_data.nComp() );
   
      if ( a_include_ghost_cells ) {
         status &= ( this_cell_data.ghostVect() == that_cell_data.ghostVect() );
      }

   }

   for (int i=0; i<num_face_vars(); ++i) {
      const LevelData<FluxBox>& this_face_data = face_var(i);
      const LevelData<FluxBox>& that_face_data = a_rhs.face_var(i);

      const DisjointBoxLayout& thisBoxes( this_face_data.disjointBoxLayout() );
      const DisjointBoxLayout& rhsBoxes( that_face_data.disjointBoxLayout() );
   
      status &= thisBoxes.compatible( rhsBoxes );
      status &= ( this_face_data.nComp() == that_face_data.nComp() );
   
      if ( a_include_ghost_cells ) {
         status &= ( this_face_data.ghostVect() == that_face_data.ghostVect() );
      }
   }

   return status;
}


RefCountedPtr<CFGVars>
CFGVars::clone( const IntVect&  a_ghost_vect,
                const bool      a_copy_data ) const
{
   CFGVars* field_ptr = new CFGVars( m_pp_prefix, m_name, m_geometry );

   for (int i=0; i<num_cell_vars(); ++i) {
      const LevelData<FArrayBox>& this_cell_data = cell_var(i);

      field_ptr->addCellVar(m_cell_data_var_name[i], this_cell_data.nComp(), a_ghost_vect);

      LevelData<FArrayBox>& new_cell_data = field_ptr->cell_var(i);

      if (a_copy_data) {
         for (DataIterator dit( new_cell_data.dataIterator() ); dit.ok(); ++dit) {
            new_cell_data[dit].copy(this_cell_data[dit]);
         }
      }
   }

   for (int i=0; i<num_face_vars(); ++i) {
      const LevelData<FluxBox>& this_face_data = face_var(i);
      
      field_ptr->addFaceVar(m_face_data_var_name[i], this_face_data.nComp(), a_ghost_vect);

      LevelData<FluxBox>& new_face_data = field_ptr->face_var(i);

      if (a_copy_data) {
         for (DataIterator dit( new_face_data.dataIterator() ); dit.ok(); ++dit) {
            new_face_data[dit].copy(this_face_data[dit]);
         }
      }
   }

   return RefCountedPtr<CFGVars>(field_ptr);
}


#include "NamespaceFooter.H"

