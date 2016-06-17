#ifndef _ARBITRARY_KINETIC_FUNCITON_H_
#define _ARBITRARY_KINETIC_FUNCITON_H_

#include "REAL.H"
#include "BoundaryBoxLayout.H.multidim"
#include "KineticSpecies.H"
#include "KineticFunction.H"
#include "LevelData.H"
#include "PhaseGeom.H"
#include "ParmParse.H"
#include "RealVect.H"

#include <string>

#include "ParsingCore.H"

#include "NamespaceHeader.H"

/**
 * Arbitrary distribution function initial condition class.
 * Initializes user defined functions of all space variables in the form
 *
 * The following represents a sample input entry for this Function choice.
 *
 * \verbatim
 * type = "arbitrary"
 * function = "1.0"
 *
*/
class ArbitraryKineticFunction : public KineticFunction
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
      ArbitraryKineticFunction( ParmParse& pp, const int& verbosity );

      /// Destructor.
      /**
       */
      virtual ~ArbitraryKineticFunction() {delete m_pscore;}

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

      ParsingCore *m_pscore;

   private:

      // prohibit copying
      ArbitraryKineticFunction( const ArbitraryKineticFunction& );
      ArbitraryKineticFunction& operator=( const ArbitraryKineticFunction& );

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

      std::string m_function;

};

#include "NamespaceFooter.H"

#endif
