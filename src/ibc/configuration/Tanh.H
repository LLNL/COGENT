#ifndef _TANH_H_
#define _TANH_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "GridFunction.H"
#include "LevelData.H"
#include "ParmParse.H"
#include "MagBlockCoordSys.H"
#include "MultiBlockLevelGeom.H"
#include "RealVect.H"

#include "NamespaceHeader.H"

/**
 * Tanh function class.
 * Assigns hyperbolic tangent in flux coordnates across an annular
 * geometry:
 *
 *   f = a \tanh( (\psi-\psi_m) / w ) + b
 *
 * where
 *
 *   a = (f_i-f_o) / (\tanh(\phi_i-\phi_m)-\tanh(\phi_o-\phi_m)),
 *   b = f_o - a \tanh(\phi_o-\phi_m),
 *
 * f_i and f_o are the values at the inner (\psi_i) and outer (\psi_o)
 * flux surfaces, and \psi_m is the flux coordinate at which the inflection
 * occurs.
 *
 * This function has only been implemented for annular geometries.
 *
 * Required input keys and data types:
 * \b inner_radial_value
 * real value at inner flux surface
 *
 * \b outer_radial_value
 * real value at outer flux surface
 *
 * \b radial_midpoint
 * real value of the inflection flux surface in normalized flux coordinates
 *
 * \b radial_width
 * real width scale in normalized flux coordinates
 *
 * The following represents a sample input entry for this function choice.
 *
 * \verbatim
 * type = "tanh"
 * inner_radial_value = 0.0
 * outer_radial_value = 1.0
 * radial_midpoint = 0.5
 * radial_width = 0.25
 *
*/
class Tanh : public GridFunction
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
      Tanh( ParmParse& a_pp, const int& a_verbosity );

       /// Destructor.
      /**
       */
      virtual ~Tanh() {;}

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
      virtual void assign( LevelData<FArrayBox>& a_dfn,
                           const MultiBlockLevelGeom& a_geometry,
                           const Real& a_time,
                           const bool& a_cell_averages = true ) const;

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
      virtual void assign( FArrayBox& a_dfn,
                           const MultiBlockLevelGeom& a_geometry,
                           const Box& a_box,
                           const Real& a_time,
                           const bool& cell_averages = true ) const;

      /// Set the function on the boundary boxes.
      /**
       * Sets the cell-averaged function value for the grid function
       * at the given time on the boxes specified in the BoundaryBoxLayout.
       *
       * NB: The cell-averaged values include the J and BStarParallel factors,
       * i.e., the internal LevelData contains the average of fJB^*_{\parallel}
       *
       * @param[inout] data  Grid function to set.
       * @param[in] bdry_layout BoundaryBoxLayout defining boxes to fill.
       * @param[in] time Time value at which to evaluate the function.
       */
      virtual void assign( LevelData<FArrayBox>& data,
                           const MultiBlockLevelGeom& geometry,
                           const BoundaryBoxLayout& bdry_layout,
                           const Real& time ) const;

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const;

   private:

      // prohibit copying
      Tanh( const Tanh& );
      Tanh& operator=( const Tanh& );

      /// Set cell-average values.
      /**
       */
      inline void setCellAverages( FArrayBox&              a_dfn,
                                   const MagBlockCoordSys& a_coord_sys ) const;

      /// Set pointwise values at cell centers.
      /**
       */
      inline void setPointwise( FArrayBox&              a_dfn,
                                const MagBlockCoordSys& a_coord_sys ) const;

      /// Parse the input database for parameters.
      /**
       */
      inline void parseParameters( ParmParse& a_pp );

      /// Check that the geometry is annular.
      /**
       */
      inline void checkGeometryValidity( const MultiBlockLevelGeom& a_geometry ) const;

      int m_verbosity;

      Real m_inner_radial_value;
      Real m_outer_radial_value;
      Real m_radial_midpoint;
      Real m_radial_width;
};

#include "NamespaceFooter.H"

#endif
