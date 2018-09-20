#ifndef _LOCALIZED_H_
#define _LOCALIZED_H_

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
 * Localized function class.
 * Assigns Gaussian functions of all space variables in the form
 *
 *   f = A \exp( \sum_{i=1}^D ((x_i-c_i)/w_i)^2 ) + floor
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
 * The following represents a sample input entry for this function choice.
 *
 * \verbatim
 * type = "localized"
 * amplitude = 1.0
 * location = 0.5 1.0
 * width = 0.025 0.1
 * floor = 1.0
 *
*/
class Localized : public GridFunction
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
      Localized( ParmParse& pp, const int& verbosity );

      /// Destructor.
      /**
       */
      virtual ~Localized() {;}

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
                           const bool& cell_averages = true ) const;

      virtual void assign( LevelData<FluxBox>& data,
                           const MultiBlockLevelGeom& geometry,
                           const Real& time,
                           const bool& cell_averages = true ) const;

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
                           const bool& cell_averages = true ) const;

      virtual void assign( FluxBox& data,
                           const MultiBlockLevelGeom& geometry,
                           const Box& box,
                           const Real& time,
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
      Localized( const Localized& );
      Localized& operator=( const Localized& );

      /// Set cell-average values.
      /**
       */
      inline void setCellAverages( FArrayBox&              data,
                                   const MagBlockCoordSys& coord_sys ) const;

      /// Set pointwise values at cell centers.
      /**
       */
      inline void setPointwise( FArrayBox&              data,
                                const MagBlockCoordSys& coord_sys ) const;

      /// Parse the input database for parameters.
      /**
       */
      inline void parseParameters( ParmParse& pp );

      /// Check that the geometry is annular or single null.
      /**
       */
      inline void checkGeometryValidity( const MultiBlockLevelGeom& geometry ) const;

      int m_verbosity;

      Real        m_amplitude;
      RealVect    m_location;
      RealVect    m_width;
      Real        m_floor;
};

#include "NamespaceFooter.H"

#endif
