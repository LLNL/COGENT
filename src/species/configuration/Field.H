#ifndef _FIELD_H_
#define _FIELD_H_

#include "Vector.H"
#include "RefCountedPtr.H"
#include "LevelData.H"
#include "FArrayBox.H"

#include "NamespaceHeader.H"

/// Field interface.
/**
 * Abstract variable representation.
 */
class Field
{
   public:

   /// Returns species name.
   /**
    */
   virtual const std::string& name() const = 0;

   virtual LevelData<FArrayBox>& data() = 0;

   const virtual LevelData<FArrayBox>& data() const = 0;

   virtual RefCountedPtr<Field> clone( const IntVect ghostVect = IntVect::Zero,
                                       const bool include_ghost_cells = true ) const = 0;

//   virtual void copy( const Field& rhs ) = 0;

   virtual void zeroData() = 0;

//   virtual void addData( const Field& rhs, const Real factor = 1.0 ) = 0;

//   virtual bool conformsTo( const Field& rhs,
//                            const bool include_ghost_cells = true ) const;


};

typedef RefCountedPtr<Field> FieldPtr;
typedef Vector<FieldPtr> FieldPtrVect;

#include "NamespaceFooter.H"

#endif
