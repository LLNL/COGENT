#ifndef _FLUID_SPECIES_H_
#define _FLUID_SPECIES_H_

#include <string>

#include "REAL.H"
#include "Dimensions.H"
#include "Vector.H"
#include "IntVect.H"
#include "RefCountedPtr.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "MagGeom.H"
#include "Variable.H"

#include "NamespaceHeader.H"

/// Fluid species interface.
/**
 * Concrete reification of abstract SpeciesModel interface for kinetic
 * species.
 */
class FluidSpecies : public Variable
{
   public:

      /// Constructor.
      /**
       */
      FluidSpecies( const std::string& name,
                    const Real         mass,
                    const Real         charge,
                    const std::string& model,
                    const MagGeom&     geometry );

      /// Copy Constructor.
      /**
       * @param[in] foo FluidSpecies from which to initialize.
        */
      FluidSpecies( const FluidSpecies& foo );

      /// Destructor.
      /**
       */
      ~FluidSpecies() {;}

      /// Returns species name.
      /**
       */
      virtual const std::string& name() const { return m_name; }

      /// Returns species mass.
      /**
       */
      virtual Real mass() const { return m_mass; }

      /// Returns species charge.
      /**
       */
      virtual Real charge() const { return m_charge; }

      /// Returns species number density.
      /**
       */
      virtual void numberDensity( LevelData<FArrayBox>& rho ) const;

      /// Returns species mass density.
      /**
       */
      virtual void massDensity( LevelData<FArrayBox>& rho ) const;

      /// Returns species charge density.
      /**
       */
      virtual void chargeDensity( LevelData<FArrayBox>& rho ) const;

      /// Returns species pressure.
      /**
       */
      virtual void velocity( LevelData<FArrayBox>& velocity ) const;

      /// Returns species pressure.
      /**
       */
      virtual void pressure( LevelData<FArrayBox>& pressure ) const;

      /// Returns species pressure.
      /**
       */
      virtual void temperature( LevelData<FArrayBox>& pressure ) const;


      /// Returns true if the species is name.
      /**
       */
      virtual bool isSpecies( const string& name ) const;

      /// Zeros out the internal data.
      /**
       */
      virtual void zeroData();

      /// Add a multiple of another Species Model data to this SpeciesModel.
      /**
       * @param[in] rhs SpeciesModel to add.
       * @param[in] factor multiple of input SpeciesModel to add.
       */
      virtual void addData( const FluidSpecies& rhs,
                            const Real          factor = 1.0 );

      /// Conformity verification
      /**
       * Checks that the rhs and this SpeciesModel have a similar
       * storage representation and returns true if so, false if not.
       *
       * @param[in] rhs SpeciesModel with which to compare.
       * @param[in] include_ghost_cells include ghost cells in comparison.
       */
      virtual bool conformsTo( const FluidSpecies& rhs,
                               const bool include_ghost_cells = true ) const;

      /// Factory method
      /**
       * Dynamically creates a new FluidSpecies with conforming layout and
       * copies internal parameters.
       *
       * @param[in] copy_soln_data Copies internal LevelDatas; otherwise,
       *                           initial LevelDatas are unititialized.
       * @param[in] include_ghost_cells Create internal LevelDatas with
       *                                ghost cells in comparison.
       */
      virtual
      RefCountedPtr<FluidSpecies> clone( const IntVect ghostVect = IntVect::Zero,
                                         const bool include_ghost_cells = true ) const;

      /// Copy
      /**
       * @param[in] rhs FluidSpecies from which to copy.
       */
      virtual void copy( const FluidSpecies& rhs );

      /// Assignment operator.
      /**
       * @param[in] rhs FluidSpecies from which to copy.
       */
      const FluidSpecies& operator=( const FluidSpecies& rhs );

      /// Provides the full fluid state.
      /**
       */
      LevelData<FArrayBox>& fluidState() { return m_state; }

      /// Provides a const reference to the full fluid state.
      /**
       */
      const LevelData<FArrayBox>& fluidState() const { return m_state; }

      /// Provides a const reference to the configuration space geometry object.
      /**
       */
      const MagGeom& congifurationSpaceGeometry() const { return m_geometry; }

      virtual LevelData<FArrayBox>& data() {return m_state;}

      const virtual LevelData<FArrayBox>& data() const {return m_state;}

   private:

      enum {INVALID=-1, MASS_DENSITY, MOMENTUM_DENSITY, TOTAL_ENERGY_DENSITY };

      enum {ISOTHERMAL, EULER};

      const MagGeom& m_geometry;
      string m_name;
      string m_model;
      Real m_mass;
      Real m_charge;
      LevelData<FArrayBox> m_state;
};

#include "NamespaceFooter.H"

#endif
