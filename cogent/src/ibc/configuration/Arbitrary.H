#ifndef _ARBITRARY_H_
#define _ARBITRARY_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "GridFunction.H"
#include "LevelData.H"
#include "ParmParse.H"
#include "MagBlockCoordSys.H"
#include "MultiBlockLevelGeom.H"
#include "RealVect.H"

#include "ParsingCore.H"

#include "NamespaceHeader.H"

/**
 * Arbitrary function class.
 * Assigns arbitrary function in configuration space:
 * 
 * coordinate_type="mapped" uses the mapped r and theta coordinates as "x" and "y"  (default)
 * coordinate_type="physical" uses the physical R and Z coordinates as "x" and "y"
 * coordinate_type="flux" uses the normalized flux as "x" and mapped theta coordinate as "y"
 * coordinate_type="outer_midplane" uses the shifter outer radial coordinate (R-Rsep) as "x" and mapped theta coordinate as "y"
 *
 * The following represents a sample input entry for this function choice.
 *
 * \verbatim
 * type = "arbitrary"
 * coordinate_type="flux"
 * function="sin(x)+sin(y)"
 *
*/
class Arbitrary : public GridFunction
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
      Arbitrary( ParmParse& pp, const int& verbosity );

      /// Destructor.
      /**
       */
      virtual ~Arbitrary() {delete m_pscore;}

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

      ParsingCore *m_pscore;
      ParsingCore *m_pscore2;

   private:

      // prohibit copying
      Arbitrary( const Arbitrary& );
      Arbitrary& operator=( const Arbitrary& );

      /// Set cell-average values.
      /**
       */
      inline void setCellAverages( FArrayBox&              a_dfn,
                                   const MagBlockCoordSys& a_coord_sys,
				   const int               a_block_number ) const;

      /// Set pointwise values at cell centers.
      /**
       */
      inline void setPointwise( FArrayBox&              a_dfn,
                                const MagBlockCoordSys& a_coord_sys,
				const int               a_block_number ) const;

      /// Parse the input database for parameters.
      /**
       */
      inline void parseParameters( ParmParse& a_pp );

      /// Check that the geometry is annular.
      /**
       */
      inline void checkGeometryValidity( const MultiBlockLevelGeom& a_geometry ) const;

      int m_verbosity;
      std::string m_function;
      std::string m_function2;
      std::string m_coord_type;

};

#include "NamespaceFooter.H"

#endif
