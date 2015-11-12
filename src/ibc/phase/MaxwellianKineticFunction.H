#ifndef _MAXWELLIAN_KINETIC_FUNCTION_H_
#define _MAXWELLIAN_KINETIC_FUNCTION_H_

#include "BoundaryBoxLayout.H.multidim"
#include "REAL.H"
#include "RefCountedPtr.H"
#include "KineticSpecies.H"
#include "KineticFunction.H"
#include "LevelData.H"
#include "PhaseGeom.H"
#include "ParmParse.H"
#include "RealVect.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "GridFunction.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include <string>

#include "NamespaceHeader.H"

/**
 * Maxwellian distribution function initial condition class.
 * Initializes Maxwellian distribution function given density and
 * temperature fields.
 *
 * Required input keys and data types:
 * \b vparallel
 * real offset in the parallel velocity
 *
 * \b density
 * sub-database prefix of initial condition specification of density field
 *
 * \b temperature
 * sub-database prefix of initial condition specification of temperature field
 *
 * The following represents a sample input entry for this Function choice.
 *
 * \verbatim
 * type = "maxwellian"
 * vparallel.type = "type0"
 * vparallel.* = ...
 * temperature.type = "type1"
 * temperature.* = ...
 * density.type = "type2"
 * density.* = ...
 *
*/
class MaxwellianKineticFunction : public KineticFunction
{
   public:

      /// Constructor using a database.
      /**
       * Constructs the object by getting parameter from a ParmParse database.
       *
       * @param[in] pp ParmParse database that contains the parameters for
       *               the new GridFunction object.
       * @param[in] verbosity Flag to turn on additional levels of output
       */
     MaxwellianKineticFunction( ParmParse& pp, const int& verbosity );

      /// Destructor.
      /**
       */
     virtual ~MaxwellianKineticFunction() {;}

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
      MaxwellianKineticFunction( const MaxwellianKineticFunction& );
      MaxwellianKineticFunction& operator=( const MaxwellianKineticFunction& );

      /// Set pointwise values at cell centers.
      /**
       */
      inline void setPointValues( FArrayBox&                a_dfn,
                                  const Box&                a_box,
                                  const PhaseBlockCoordSys& a_coord_sys,
                                  const FArrayBox&          a_density,
                                  const FArrayBox&          a_temperature,
                                  const FArrayBox&          a_vparallel,
                                  const FArrayBox&          a_B,
                                  const Real&               a_mass ) const;

      /// Initialize density or temperature field on a LevelData.
      /**
       */
      inline void initializeField( LevelData<FArrayBox>& a_field,
                                   const CFG::GridFunction& a_ic,
                                   const PhaseGeom& a_phase_geometry,
                                   const Real& a_time ) const;

      /// Initialize density or temperature field on a box.
      /**
       */
      inline void initializeField( FArrayBox& a_field,
                                   const CFG::GridFunction& a_ic,
                                   const PhaseGeom& a_geometry,
                                   const Box& a_interior_box,
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

      RefCountedPtr<CFG::GridFunction> m_ic_vparallel;
      RefCountedPtr<CFG::GridFunction> m_ic_temperature;
      RefCountedPtr<CFG::GridFunction> m_ic_density;
};

#include "NamespaceFooter.H"

#endif
