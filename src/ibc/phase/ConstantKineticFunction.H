#include "BoundaryBoxLayout.H.multidim"
#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "KineticSpecies.H"
#include "KineticFunction.H"

#include <iostream>
#include <string>

#include "NamespaceHeader.H"

/**
 * Constant distibution function initial condition class.
 * Sets distribution function of a KineticSpecies to a constant.
 *
 * Required input keys and data types:
 * \b value
 * real constant value to set
 *
 * The following represents a sample input entry for this Function choice:
 *
 * \verbatim
 * type = "constant"
 * value = 1.0
*/
class ConstantKineticFunction : public KineticFunction
{
   public:

      /// Constructor using a value.
      /**
       * Constructs the object using a passed-in parameter.
       *
       * @param[in] value Real value to set.
       * @param[in] verbosity Flag to turn on additional levels of output
       */
     ConstantKineticFunction( const Real& value, const int& verbosity )
         : m_value(value), m_verbosity(verbosity) {;}

      /// Constructor using a database.
      /**
       * Constructs the object by getting parameter from a ParmParse database.
       *
       * @param[in] pp ParmParse database that contains the parameters for
       *               the new InitialCondition object.
       * @param[in] verbosity Flag to turn on additional levels of output
       */
     ConstantKineticFunction( ParmParse& pp, const int& verbosity )
         : m_value(0.0), m_verbosity(verbosity)
      {
         parseParameters( pp );
      }

      /// Destructor.
      /**
       */
      virtual ~ConstantKineticFunction() {;}

      /// Set the function on the interior.
      /**
       * Sets the cell-averaged function value for the kinetic species
       * at the given time.  It is assumed that derived classes will check
       * that the specification is suitable for the given geometry embedded in
       * the KineticSpecies object.
       *
       * NB: The cell-averaged values include the J and BStarParallel factors,
       * i.e., the internal LevelData contains the average of fJB^*_{\parallel}
       *
       * @param[inout] species Kinetic species to intiialize.
       * @param[in] time Time value at which to evaluate the function.
       */
      virtual void assign( KineticSpecies& species, const Real& time ) const;

      /// Set the function on the boundary boxes.
      /**
       * Sets the cell-averaged function value for the kinetic species
       * at the given time on the boxes specified in the BoundaryBoxLayout.
       * It is assumed that derived classes will check that the specification is
       * suitable for the given geometry embedded in the KineticSpecies object.
       *
       * NB: The cell-averaged values include the J and BStarParallel factors,
       * i.e., the internal LevelData contains the average of fJB^*_{\parallel}
       *
       * @param[inout] species Kinetic species to intiialize.
       * @param[in] bdry_layout BoundaryBoxLayout defining boxes to fill.
       * @param[in] time Time value at which to evaluate the function.
       */
      virtual void assign( KineticSpecies& species,
                           const BoundaryBoxLayout& bdry_layout,
                           const Real& time ) const;

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const {
         if (procID()==0) {
            std::cout << "Constant kinetic function parameters:" << std::endl;
            std::cout << "  value: "  << m_value  << std::endl << std::endl;
         }
      }

   private:

      // prohibit copying
      ConstantKineticFunction( const ConstantKineticFunction& );
      ConstantKineticFunction& operator=( const ConstantKineticFunction& );

      /// Parse the input database for parameters.
      /**
       */
      void parseParameters( ParmParse& a_pp )
      {
         a_pp.get( "value", m_value );
         if (m_verbosity) {
            printParameters();
         }
      }

      /// Convert cell center pointwise values to cell averaged values.
      /**
       */
      void convertToCellAverage( const MultiBlockCoordSys&  a_coord_sys,
                                 LevelData<FArrayBox>&      a_dfn ) const;

      Real m_value;
      int m_verbosity;
};

#include "NamespaceFooter.H"

