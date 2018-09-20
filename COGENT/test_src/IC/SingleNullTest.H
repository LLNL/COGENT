#ifndef  _SINGLE_NULL_TEST_H_
#define  _SINGLE_NULL_TEST_H_

#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "MagCoordSys.H"
#include "MultiBlockLevelGeom.H"
#include "ParmParse.H"
#include "ProblemDomain.H"

#include "InitialCondition.H.multidim"

#include <vector>
#include <string>

#include "NamespaceHeader.H"


class SingleNullTest
{
   public:

      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       */
      SingleNullTest( ParmParse& pp );

      /// Destructor.
      /**
       */
      ~SingleNullTest();

      /// Initialize the data
      /**
       * Applies each initial condition to an associatred LevelData
       */
      void initialize();

      /// Write plot file.
      /**
       * Write plotting data to HDF5 files with a common prefix.
       *
       *  @param[in] prefix the string used to prefix the visit files
       *                    and directory.
       */
      void writePlotFile( const std::string& prefix );

   private:

      // To add a new IC to test, add a new enum name before "NUM_IC"
      enum {INVALID=-1,
            ZERO,
            CONSTANT,
            LOCALIZED,
            NUM_IC};

      /// Apply initial conditions.
      /**
       * Private method that applies initial conditions.
       */
      void setInitialConditions();

      /// Create configuration space geometry and coordinates.
      /**
       * Private method that creates configuration space geometry and
       * coordinates.
       */
      void createConfigurationSpace();

      /// Create grid.
      /**
       * Private method that creates the DisjointBoxLayout and distributes it
       * across processors.
       */
      void createDisjointBoxLayout( DisjointBoxLayout& grids );

      /// Parse parameters.
      /**
       * Private method to obtain control parameters from
       * ParmParse database
       *
       * \param[in] pp input database, properly delimited with a prefix.
       */
      void parseParameters( ParmParse& pp );

      MultiBlockLevelGeom* m_mag_geom;
      MagCoordSys*         m_mag_geom_coords;
      ProblemDomain        m_domain;

      std::vector<LevelData<FArrayBox>*> m_solution;

      std::vector<InitialCondition*> m_ic;

      std::vector<std::string> m_ic_name;
};

#include "NamespaceFooter.H"

#endif
