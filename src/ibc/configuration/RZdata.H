#ifndef _RZDATA_H_
#define _RZDATA_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "GridFunction.H"
#include "LevelData.H"
#include "ParmParse.H"
#include "MagBlockCoordSys.H"
#include "MultiBlockLevelGeom.H"
#include "RealVect.H"

#include <cmath>
#include "mba.hpp"

#include <iostream>
#include <fstream>

#include "NamespaceHeader.H"

/**
 * RZdata function class that interpolates scattered data onto the COGENT grid for 
 * an arbitrary choise of geometry.
 *
 * Subtype RBF:
 * Uses the radial basis function interpolation (see  rbf_interp_2d, r8lib,
 * test_interp_2d and 
 * http://people.sc.fsu.edu/~jburkardt/cpp_src/rbf_interp_2d/rbf_interp_2d.html for details).
 * NB: works slow, but provides smooth interpolation function (thus, reasonble behavior outise
 * the interpolated data range. Use for a small set of points, e.g., ~ 100.
 *
 * Subtype MBA:
 * Uses the multilevel B-spline approximation (see mba.hpp and the reference wherein for control
 * parameters). NB: works very fast, but provides only C2 interpolation. Behaves badly outside
 * the interpolated data range. Use for a large set of data points.
 *
 * The data file should be written as (R, Z, Data), and should not contain empty lines
 * (the code assinges the number of the data points to the number of the lines in the data
 * file). 
 *
 *
 * Required input keys and data types:
 * \b data_file
 * name of the data file 
 *
 * The following represents a sample input entry for this function choice.
 *
 * \verbatim
 * type = "RZdata"
 * subtype = "MBA"
 * data_file = "data.txt"
 *
*/
class RZdata : public GridFunction
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
      RZdata( ParmParse& pp, const int& verbosity );

      /// Destructor.
      /**
       */
      ~RZdata();

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
      //RZdata( const RZdata& );
      //RZdata& operator=( const RZdata& );

      /// Set cell-average values.
      /**
       */
      inline void setCellAverages( FArrayBox&              data,
                                   const MultiBlockCoordSys& coords,
                                   const MagBlockCoordSys& block_coord_sys,
                                   const int               block_number ) const;

      /// Set pointwise values at cell centers.
      /**
       */
      inline void setPointwise( FArrayBox&              data,
                                const MagBlockCoordSys& block_coord_sys) const;

      /// Get geoemtry parameters
      void getGeomParameters(const MultiBlockCoordSys& a_coord_sys) const;

      /// Create RBF inerpolation
      void createInterpolationRBF();
    
      /// Create RBF inerpolation
      void createInterpolationMBA();

      /// Parse the input database for parameters.
      /**
       */
      inline void parseParameters( ParmParse& pp );

      int m_verbosity;
      std::string m_data_file;
      std::string m_subtype;

      int m_data_npoints;
      double m_r0;
      double* m_data_grid;
      double* m_weights;
    
      mba::cloud<2> *m_MBA;
      int m_init_lattice_MBA;
    
};

#include "NamespaceFooter.H"

#endif
