#ifndef _CONSTANT_H_
#define _CONSTANT_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "DataIterator.H"
#include "LevelData.H"
#include "MultiBlockLevelGeom.H"

#include <iostream>

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
         : m_value(value), m_verbosity(verbosity) {;}

      /// Constructor using a database.
      /**
       * Constructs the object by getting parameter from a ParmParse database.
       *
       * @param[in] pp ParmParse database that contains the parameters for
       *               the new GridFunction object.
       * @param[in] verbosity Flag to turn on additional levels of output
       */
      Constant( ParmParse& pp, const int& verbosity )
         : m_value(0.0), m_verbosity(verbosity)
      {
         parseParameters( pp );
      }

      /// Destructor.
      /**
       */
      virtual ~Constant() {;}

      /// Set the grid function.
      /**
       * Sets the grid function on the LevelData for the given geometry
       * at the given time.  It is assumed that derived classes will check
       * that the specification is suitable for the given geometry.  Cell
       * average values are set unless the optional parameter is false, in
       * which case, pointwise values are set.
       *
       * @param[out] data LevelData to set.
       * @param[in] geometry Corresponding level geometry.
       * @param[in] time Time value at which to evaluate the function.
       * @param[in] as_cell_name Boolean switch for cell average values or
       *                         point values.
       */
      virtual void assign( LevelData<FArrayBox>& data,
                           const MultiBlockLevelGeom& geometry,
                           const Real& time,
                           const bool& cell_averages = true ) const
      {
         const DisjointBoxLayout& grid( data.getBoxes() );
         for (DataIterator dit( grid.dataIterator() ); dit.ok(); ++dit) {
            data[dit].setVal( m_value );
         }
      }

      /// Set the grid function.
      /**
       * Sets the grid function on the box for the given geometry
       * at the given time.  It is assumed that derived classes will check
       * that the specification is suitable for the given geometry.  Cell
       * average values are set unless the optional parameter is false, in
       * which case, pointwise values are set.
       *
       * @param[out] data LevelData to set.
       * @param[in] geometry Corresponding level geometry.
       * @param[in] box Box from which to get coordinate system.
       * @param[in] time Time value at which to evaluate the function.
       * @param[in] as_cell_name Boolean switch for cell average values or
       *                         point values.
       */
       virtual void assign( FArrayBox& data,
                            const MultiBlockLevelGeom& geometry,
                            const Box& box,
                            const Real& time,
                            const bool& cell_averages = true ) const
      {
         data.setVal( m_value );
      }

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

      Real m_value;
      int m_verbosity;
};

#include "NamespaceFooter.H"

#endif

