#ifndef _TENSOR_PRODUCT_SINE_KINETIC_FUNCTION_H_
#define _TENSOR_PRODUCT_SINE_KINETIC_FUNCTION_H_

#include "REAL.H"
#include "BoundaryBoxLayout.H.multidim"
#include "KineticSpecies.H"
#include "KineticFunction.H"
#include "LevelData.H"
#include "PhaseGeom.H"
#include "ParmParse.H"
#include "RealVect.H"

#include <string>

#include "NamespaceHeader.H"

/**
 * Localized distribution function initial condition class.
 * Initializes Gaussian functions of all space variables in the form
 *
 *   f = A \prod_{i=1}^D \sin(  2\pi k_i x_i )
 *
 * where A is the peak scale and the k is the wavevector.
 *
 * Implemented for both annular and single null geometries.
 *
 * Required input keys and data types:
 * \b amplitude
 * real scale of the peak of the Gaussian function
 *
 * \b wavevector
 * real vector of wave numbers
 *
 * The following represents a sample input entry for this Function choice.
 *
 * \verbatim
 * type = "localized"
 * amplitude = 1.0
 * wavevector = 1.0 1.0 1.0 1.0
 *
*/
class TensorProductSineKineticFunction : public KineticFunction
{
   public:

      /// Constructor using a database.
      /**
       * Constructs the object by getting parameter from a ParmParse database.
       *
       * @param[in] pp ParmParse database that contains the parameters for
       *               the new InitialCondition object.
       * @param[in] verbosity Flag to turn on additional levels of output
       */
      TensorProductSineKineticFunction( ParmParse& pp, const int& verbosity );

      /// Destructor.
      /**
       */
      virtual ~TensorProductSineKineticFunction() {;}

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
      inline void printParameters() const;

   private:

      // prohibit copying
      TensorProductSineKineticFunction( const TensorProductSineKineticFunction& );
      TensorProductSineKineticFunction& operator=( const TensorProductSineKineticFunction& );

      /// Set pointwise values at cell centers.
      /**
       */
      inline void setPointValues( FArrayBox&                a_dfn,
                                  const PhaseBlockCoordSys& a_coord_sys,
                                  const Real& a_time ) const;

      /// Parse the input database for parameters.
      /**
       */
      inline void parseParameters( ParmParse& a_pp );

      /// Check that the geometry is annular or single null.
      /**
       */
      inline void checkGeometryValidity( const PhaseGeom& a_geometry ) const;

      /// Convert cell center pointwise values to cell averaged values.
      /**
       */
      void convertToCellAverage( const MultiBlockCoordSys&  a_coord_sys,
                                 LevelData<FArrayBox>&      a_dfn ) const;

      int m_verbosity;

      Real m_amplitude;
      RealVect m_wavevector;

};

#include "NamespaceFooter.H"

#endif
