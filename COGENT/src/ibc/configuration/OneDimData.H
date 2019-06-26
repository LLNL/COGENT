#ifndef _ONEDIMDATA_H_
#define _ONEDIMDATA_H_

#include "GridFunction.H"

#include <cmath>
#include "mba.hpp"

#include <iostream>
#include <fstream>

#include "NamespaceHeader.H"

/**
 * OneDimData function class that interpolates scattered data onto the COGENT grid for 
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
 * The data file should be written as (coord, data), and should not contain empty lines
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
 * type = "OneDimData"
 * subtype = "MBA"
 * coordinate_type="flux"
 * data_file = "data.txt"
 *
*/
class OneDimData : public GridFunction
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
      OneDimData( ParmParse& pp, const int& verbosity );

      /// Destructor.
      /**
       */
      ~OneDimData();

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const;

   private:

      // prohibit copying
      //OneDimData( const OneDimData& );
      //OneDimData& operator=( const OneDimData& );

      /// Set pointwise values at cell centers.
      /**
       */
      virtual void setPointwise(FArrayBox&                  dfn,
                                const MultiBlockLevelGeom&  geometry,
                                const FArrayBox&            real_coords,
                                const FArrayBox&            normalized_flux,
                                const int                   block_number ) const;
   
      /// Parse the input database for parameters.
      /**
       */
      inline void parseParameters( ParmParse& pp );
   
      /// Check geometry validity.
      /**
       */
      virtual void checkGeometryValidity( const MultiBlockLevelGeom& geometry ) const;

      /// Create RBF inerpolation
      void createInterpolationRBF();
    
      /// Create RBF inerpolation
      void createInterpolationMBA();

      std::string m_data_file;
      std::string m_subtype;
      std::string m_coord_type;

      int m_data_npoints;
      double m_r0;
      double* m_coords;
      double* m_weights;
    
      mba::cloud<1> *m_MBA;
      int m_init_lattice_MBA;
    
};

#include "NamespaceFooter.H"

#endif
