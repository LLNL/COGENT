#ifndef _FLUID_SPECIES_H_
#define _FLUID_SPECIES_H_

#include "CFGVar.H"

#include "NamespaceHeader.H"


/// Fluid species interface.
/**
 * Fluid species interface
 */
class FluidSpecies 
   : public CFGVar
{
   public:

      /// Constructor.
      /**
       */
      FluidSpecies( const string&      pp_prefix,
                    const std::string& name,
                    const MagGeom&     geometry,
                    const IntVect&     ghost_vect );

      /// Copy Constructor.
      /**
       * @param[in] foo FluidSpecies from which to initialize.
        */
      FluidSpecies( const FluidSpecies& foo );

      /// Destructor.
      /**
       */
      ~FluidSpecies() {;}

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

      /// Returns species velocity.
      /**
       */
      virtual void velocity( LevelData<FArrayBox>& velocity ) const;

      /// Returns species velocity.
      /**
       */
   //      virtual LevelData<FluxBox>& velocity() {return face_data();}
      virtual LevelData<FluxBox>& velocity() {return m_velocity;}

   //      virtual const LevelData<FluxBox>& velocity() const {return face_data();}
      virtual const LevelData<FluxBox>& velocity() const {return m_velocity;}

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

      virtual
      RefCountedPtr<CFGVar> convertToPhysical( const IntVect& a_ghost_vect ) const;

      virtual void convertFromPhysical();

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
      RefCountedPtr<CFGVar> clone( const IntVect& ghost_vect = IntVect::Zero,
                                   const bool copy_soln_data = true ) const;

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

      Real maxValue() const;
      Real minValue() const;

   private:

      LevelData<FluxBox> m_velocity;
   
      enum {INVALID=-1, MASS_DENSITY, MOMENTUM_DENSITY, TOTAL_ENERGY_DENSITY };

      enum {ISOTHERMAL, EULER};

      string m_op_type;
      Real m_mass;
      Real m_charge;
};

#include "NamespaceFooter.H"

#endif
