#ifndef _SCALAR_H_
#define _SCALAR_H_

#include <string>

#include "REAL.H"
#include "Vector.H"
#include "RefCountedPtr.H"

#include "NamespaceHeader.H"

/// Scalar interface.
/**
 * Abstract variable representation.
 */
class Scalar
{
   public:
      /// Constructor.
      /**
       */
      Scalar( const std::string& name,
              const int          num_scalars );
   
      /// Copy Constructor.
      /**
       * @param[in] foo Scalar from which to initialize.
        */
      Scalar( const Scalar& foo );

      /// Destructor.
      /**
       */
      virtual ~Scalar() {;}

      /// Returns species name.
      /**
       */
      virtual const std::string& name() const { return m_name; }

      /// Zeros out the internal data.
      /**
       */
      virtual void zeroData();

      /// Add a multiple of another Scalar data to this Scalar.
      /**
       * @param[in] rhs Scalar to add.
       * @param[in] factor multiple of input Scalar to add.
       */
      virtual void addData( const Scalar& rhs,
                            const Real          factor = 1.0 );

      /// Conformity verification
      /**
       * Checks that the rhs and this Scalar have a similar
       * storage representation and returns true if so, false if not.
       *
       * @param[in] rhs Scalar with which to compare.
       * @param[in] include_ghost_cells include ghost cells in comparison.
       */
      virtual bool conformsTo( const Scalar& rhs,
                               const bool include_ghost_cells = true ) const;

      /// Factory method
      /**
       * Dynamically creates a new Scalar with conforming layout and
       * copies internal parameters.
       *
       * @param[in] copy_soln_data Copies internal LevelDatas; otherwise,
       *                           initial LevelDatas are unititialized.
       * @param[in] include_ghost_cells Create internal LevelDatas with
       *                                ghost cells in comparison.
       */
      virtual
      RefCountedPtr<Scalar> clone( const bool copy_soln_data = true ) const;

      /// Copy
      /**
       * @param[in] rhs Scalar from which to copy.
       */
      virtual void copy( const Scalar& rhs );

      /// Assignment operator.
      /**
       * @param[in] rhs Scalar from which to copy.
       */
      const Scalar& operator=( const Scalar& rhs );

      virtual Vector<Real>& data() {return m_data;}

      const virtual Vector<Real>& data() const {return m_data;}

      inline
      virtual int size() const
        {
           return m_data.size();
        }

   private:

      std::string m_name;
      Vector<Real> m_data;

};

typedef RefCountedPtr<Scalar> ScalarPtr;
typedef Vector<ScalarPtr> ScalarPtrVect;

#include "NamespaceFooter.H"

#endif
