#ifndef  _SINGLE_NULL_SPECIES_TEST_H_
#define  _SINGLE_NULL_SPECIES_TEST_H_

#define CFG_NAMESPACE D2
#define PDIM 4
#define VEL_DIM 2
#define VEL_NAMESPACE D2

#include "Dimensions.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "IntVect.H"
#include "KineticSpecies.H"
#include "KineticSpeciesIC.H"
#include "LevelData.H"
#include "PhaseGeom.H"
#include "ParmParse.H"
#include "ProblemDomain.H"
#include "SingleNullPhaseCoordSys.H"

#include <vector>
#include <string>

#undef CH_SPACEDIM
#define CH_SPACEDIM VEL_DIM
#include "DisjointBoxLayout.H"
#include "InitialCondition.H.multidim"
#include "IntVect.H"
#include "VelCoordSys.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "DisjointBoxLayout.H"
#include "InitialCondition.H.multidim"
#include "IntVect.H"
#include "MagCoordSys.H"
#include "MultiBlockLevelGeom.H"
#include "ProblemDomain.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"
namespace CFG = CFG_NAMESPACE;
namespace VEL = VEL_NAMESPACE;

class SingleNullSpeciesTest
{
   public:

      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       */
      SingleNullSpeciesTest( ParmParse& pp );

      /// Destructor.
      /**
       */
      ~SingleNullSpeciesTest();

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
            MAXWELLIAN,
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
      void createConfigurationSpaceDisjointBoxLayout( CFG::DisjointBoxLayout& grids );

      /// Create configuration space geometry and coordinates.
      /**
       * Private method that creates configuration space geometry and
       * coordinates.
       */
      void createVelocitySpace();

      /// Returns the problem domain for velocity space.
      /**
       * Private method that computes the problem domain for velocity space.
       */
      VEL::ProblemDomain getVelocityDomain() const;

      /// Create grid.
      /**
       * Private method that creates the DisjointBoxLayout and distributes it
       * across processors.
       */
      void createVelocitySpaceDisjointBoxLayout( VEL::DisjointBoxLayout& grids,
                                                 const VEL::ProblemDomain& domain );

      /// Create configuration space geometry and coordinates.
      /**
       * Private method that creates configuration space geometry and
       * coordinates.
       */
      void createPhaseSpace( ParmParse& a_pp );

      /// Create grid.
      /**
       * Private method that creates the DisjointBoxLayout and distributes it
       * across processors.
       */
      void createPhaseSpaceDisjointBoxLayout( DisjointBoxLayout& grids,
                                              const Vector<ProblemDomain>& a_domains );

      /// Parse parameters.
      /**
       * Private method to obtain control parameters from
       * ParmParse database
       *
       * \param[in] pp input database, properly delimited with a prefix.
       */
      void parseParameters( ParmParse& pp );

      CFG::MultiBlockLevelGeom* m_mag_geom;
      CFG::MagCoordSys*         m_mag_geom_coords;
      CFG::ProblemDomain        m_domain;

      VEL::VelCoordSys*    m_velocity_coords;

      SingleNullPhaseCoordSys* m_phase_coords;
      PhaseGeom*               m_phase_geom;

      std::vector<KineticSpecies*>   m_solution;

      std::vector<KineticSpeciesIC*> m_ic;

      std::vector<std::string> m_ic_name;

      std::vector<int> m_num_velocity_cells;
      std::vector<int> m_velocity_decomposition;
      std::vector<int> m_fixed_plotindices;
};

#include "NamespaceFooter.H"

#endif
