#ifndef _KINETIC_FUNCTION_H_
#define _KINETIC_FUNCTION_H_

#include "BoundaryBoxLayout.H.multidim"
#include "REAL.H"
#include "KineticSpecies.H"
#include "MultiBlockLevelGeom.H"
#include "PhaseGeom.H"

#include <string>

#include "NamespaceHeader.H"

/**
 * Generic initial condition abstract interface for kinetic species.
*/
class KineticFunction
{
   public:

      /// Destructor.
      /**
       */
      virtual ~KineticFunction() {;}

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
      virtual void assign( KineticSpecies& species, const Real& time ) const = 0;

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
                           const Real& time ) const = 0;

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const = 0;
};

#include "NamespaceFooter.H"

#endif
