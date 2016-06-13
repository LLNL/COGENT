#ifndef _FIELD_H_
#define _FIELD_H_

#include <string>

#include "REAL.H"
#include "Vector.H"
#include "RefCountedPtr.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "MagGeom.H"

#include "NamespaceHeader.H"

/// Field interface.
/**
 * Abstract variable representation.
 */
class Field
{
   public:
      /// Constructor.
      /**
       */
      Field( const std::string& name,
             const MagGeom&     geometry );
   
      /// Copy Constructor.
      /**
       * @param[in] foo Field from which to initialize.
        */
      Field( const Field& foo );

      /// Destructor.
      /**
       */
      ~Field() {;}

      /// Returns species name.
      /**
       */
      virtual const std::string& name() const { return m_name; }

      /// Zeros out the internal data.
      /**
       */
      virtual void zeroData();

      /// Add a multiple of another Field data to this Field.
      /**
       * @param[in] rhs Field to add.
       * @param[in] factor multiple of input Field to add.
       */
      virtual void addData( const Field& rhs,
                            const Real          factor = 1.0 );

      /// Conformity verification
      /**
       * Checks that the rhs and this Field have a similar
       * storage representation and returns true if so, false if not.
       *
       * @param[in] rhs Field with which to compare.
       * @param[in] include_ghost_cells include ghost cells in comparison.
       */
      virtual bool conformsTo( const Field& rhs,
                               const bool include_ghost_cells = true ) const;

      /// Factory method
      /**
       * Dynamically creates a new Field with conforming layout and
       * copies internal parameters.
       *
       * @param[in] copy_soln_data Copies internal LevelDatas; otherwise,
       *                           initial LevelDatas are unititialized.
       * @param[in] include_ghost_cells Create internal LevelDatas with
       *                                ghost cells in comparison.
       */
      virtual
      RefCountedPtr<Field> clone( const IntVect ghost_vect = IntVect::Zero,
                                  const bool include_ghost_cells = true ) const;

      /// Copy
      /**
       * @param[in] rhs Field from which to copy.
       */
      virtual void copy( const Field& rhs );

      /// Assignment operator.
      /**
       * @param[in] rhs Field from which to copy.
       */
      const Field& operator=( const Field& rhs );

      /// Provides a const reference to the configuration space geometry object.
      /**
       */
      const MagGeom& configurationSpaceGeometry() const { return m_geometry; }

      virtual LevelData<FArrayBox>& data() {return m_data;}

      const virtual LevelData<FArrayBox>& data() const {return m_data;}

   inline
   virtual int size( bool a_count_ghosts = false )
   {
      int size(0);
      const DisjointBoxLayout& dbl( m_data.disjointBoxLayout() );
      for (DataIterator dit( dbl.dataIterator() ); dit.ok(); ++dit) {
         const Box& box( a_count_ghosts ? m_data[dit].box() : dbl[dit] );
         size += box.numPts() * m_data.nComp();
      }
      return size;
   }

   private:

      string m_name;
      const MagGeom& m_geometry;
      LevelData<FArrayBox> m_data;

};

typedef RefCountedPtr<Field> FieldPtr;
typedef Vector<FieldPtr> FieldPtrVect;

#include "NamespaceFooter.H"

#endif
