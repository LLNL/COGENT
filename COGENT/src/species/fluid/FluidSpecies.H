#ifndef _FLUID_SPECIES_H_
#define _FLUID_SPECIES_H_

#include "CFGVars.H"
#include "NodeFArrayBox.H"
#include "EdgeDataBox.H"

#include "NamespaceHeader.H"


/// Fluid species interface.
/**
 * Fluid species interface
 */
class FluidSpecies 
   : public CFGVars
{
   public:

      /// Constructor.
      /**
       */
      FluidSpecies( const string&       pp_prefix,
                    const std::string&  name,
                    const MagGeom&      geometry,
                    const IntVect&      ghost_vect );

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
      
      /// Returns species velocity in virtual direction
      /**
       */
      virtual void velocity_virtual( LevelData<FArrayBox>& velocity ) const;

      /// Returns species mean energy.
      /**
       */
      virtual void meanEnergy( LevelData<FArrayBox>& meanEnergy ) const;

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
      virtual LevelData<FluxBox>& velocity() {return m_velocity;}

      virtual const LevelData<FluxBox>& velocity() const {return m_velocity;}

      /// Returns species pressure.
      /**
       */
      virtual void pressure( LevelData<FArrayBox>& pressure ) const;
      
      /// Returns species magnetic pressure.
      /**
       */
      virtual void Bpressure( LevelData<FArrayBox>& Bpressure ) const;

      /// Returns species temperature.
      /**
       */
      virtual void temperature( LevelData<FArrayBox>& temperature ) const;


      void interpFaceVarToCell( LevelData<FArrayBox>&  a_cell_var,
                          const string&                a_face_var_name ) const;
      
      void interpEdgeVarToCell( LevelData<FArrayBox>&  a_cell_var,
                          const string&                a_edge_var_name ) const;

      /// Returns true if the species is name.
      /**
       */
      virtual bool isSpecies( const string& name ) const;

      virtual void convertToPhysical();
      
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
      RefCountedPtr<CFGVars> clone( const IntVect& ghost_vect = IntVect::Zero,
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

   //      int cell_component(const string& variable_name) const;

      Real maxValue() const;
      Real minValue() const;

      int m_evolve_momentumDensity=0, m_evolve_momentumDensity_virtual=0; 
      int m_evolve_energyDensity=0;
      int m_evolve_magneticField=0, m_evolve_magneticField_virtual=0;
      int m_evolve_electricField=0, m_evolve_electricField_virtual=0;
      int m_evolve_currentDensity=0, m_evolve_currentDensity_virtual=0;
      
      Real m_gamma = 5.0/3.0;

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
