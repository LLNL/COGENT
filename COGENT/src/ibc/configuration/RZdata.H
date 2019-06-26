#ifndef _RZDATA_H_
#define _RZDATA_H_

#include "GridFunction.H"

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

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const;

   private:

      // prohibit copying
      //RZdata( const RZdata& );
      //RZdata& operator=( const RZdata& );

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
      virtual void checkGeometryValidity( const MultiBlockLevelGeom& geometry ) const {;}

      /// Create RBF inerpolation
      void createInterpolationRBF();
    
      /// Create RBF inerpolation
      void createInterpolationMBA();

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
