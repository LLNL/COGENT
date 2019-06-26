#ifndef _CONSTANT_H_
#define _CONSTANT_H_

#include "GridFunction.H"

#include "NamespaceHeader.H"

/**
 * Generic constant function class.
 * Assigns LevelData of any dimension to a constant.
 *
 * Required input keys and data types:
 * \b value
 * real constant value to set
 *
 * The following represents a sample input entry for this function choice:
 *
 * \verbatim
 * type = "constant"
 * value = 1.0
*/
class Constant : public GridFunction
{
   public:

      /// Constructor using a value.
      /**
       * Constructs the object using a passed-in parameter.
       *
       * @param[in] value Real value to set.
       * @param[in] verbosity Flag to turn on additional levels of output
       */
      Constant( const Real& value, const int& verbosity )
         : GridFunction(verbosity), m_value(value) {;}

      /// Constructor using a database.
      /**
       * Constructs the object by getting parameter from a ParmParse database.
       *
       * @param[in] pp ParmParse database that contains the parameters for
       *               the new GridFunction object.
       * @param[in] verbosity Flag to turn on additional levels of output
       */
      Constant( ParmParse& pp, const int& verbosity )
         : GridFunction(verbosity), m_value(0.0)
      {
         parseParameters( pp );
      }

      /// Destructor.
      /**
       */
      virtual ~Constant() {;}

      /// Print object parameters.
      /**
       */
      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Constant grid function parameters:" << std::endl;
            std::cout << "  value: "  << m_value  << std::endl << std::endl;
         }
      }

   private:

      // prohibit copying
      Constant( const Constant& );
      Constant& operator=( const Constant& );

      /// Parse the input database for parameters.
      /**
       */
      void parseParameters( ParmParse& pp )
      {
         pp.get( "value", m_value );
         if (m_verbosity) {
            printParameters();
         }
      }

      /// Set pointwise values at cell centers.
      /**
       */
      virtual void setPointwise(FArrayBox&                  dfn,
                                const MultiBlockLevelGeom&  geometry,
                                const FArrayBox&            real_coords,
                                const FArrayBox&            normalized_flux,
                                const int                   block_number ) const
      {
         dfn.setVal(m_value);
      }

      /// Check geometry validity.
      /**
       */
      virtual void checkGeometryValidity( const MultiBlockLevelGeom& geometry ) const {;}
   
      Real m_value;
};

#include "NamespaceFooter.H"

#endif

