#ifndef _GRID_FUNCTION_H_
#define _GRID_FUNCTION_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "MultiBlockLevelGeom.H"
#include "BoundaryBoxLayout.H.multidim"

#include "NamespaceHeader.H"

/**
 * Generic grid function abstract interface.
*/
class GridFunction
{
   public:

      /// Destructor.
      /**
       */
      virtual ~GridFunction() {;}

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
                           const bool& as_cell_averages = true ) const = 0;

      virtual void assign( LevelData<FluxBox>& data,
                           const MultiBlockLevelGeom& geometry,
                           const Real& time,
                           const bool& as_cell_averages = true ) const = 0;

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
      virtual void assign( FArrayBox& dfn,
                           const MultiBlockLevelGeom& geometry,
                           const Box& box,
                           const Real& time,
                           const bool& cell_averages = true ) const = 0;

      virtual void assign( FluxBox& dfn,
                           const MultiBlockLevelGeom& geometry,
                           const Box& box,
                           const Real& time,
                           const bool& cell_averages = true ) const = 0;

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
                           const Real& time ) const = 0;

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const = 0;
};

#include "NamespaceFooter.H"

#endif
