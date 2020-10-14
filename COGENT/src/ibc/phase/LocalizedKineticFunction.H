#ifndef _LOCALIZED_KINETIC_FUNCITON_H_
#define _LOCALIZED_KINETIC_FUNCITON_H_

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
 *   f = A \exp( \suum_{i=1}^D ((x_i-c_i)/w_i)^2 ) + floor
 *
 * where A is the peak scale, floor is the background floor, 
 * the vector c is the center of the pulse, and the vector w
 * give the width scales for each direction.
 *
 * Implemented for both annular and single null geometries.
 *
 * Required input keys and data types:
 * \b amplitude
 * real scale of the peak of the Gaussian function
 *
 * \b location
 * real coordinates of peak of Gaussian in physical coordinates (R,Z)
 *
 * \b width
 * real vector of 1/e widths of Gaussian in each direction
 *
 * \b floor
 * constant real floor
 *
 * The following represents a sample input entry for this Function choice.
 *
 * \verbatim
 * type = "localized"
 * amplitude = 1.0
 * location = 0.5 1.0 1.0 1.0
 * width = 0.05 0.1 0.1 0.1
 * floor = 1.0
 *
*/
class LocalizedKineticFunction : public KineticFunction
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
      LocalizedKineticFunction( ParmParse& pp, const int& verbosity );

      /// Destructor.
      /**
       */
      virtual ~LocalizedKineticFunction() {;}

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
      LocalizedKineticFunction( const LocalizedKineticFunction& );
      LocalizedKineticFunction& operator=( const LocalizedKineticFunction& );

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

      int m_verbosity;

      Real m_amplitude;
      RealVect m_location;
      RealVect m_width;
      Real m_floor;

	  bool m_useSG;
};

#include "NamespaceFooter.H"

#endif
