#ifndef _KINETIC_SPECIES_H_
#define _KINETIC_SPECIES_H_

#include "Dimensions.H"
#include "PhaseGeom.H"
//#include "ICInterface.H.multidim"
//#include "BCInterface.H.multidim"
#include "Vector.H"
#include "RefCountedPtr.H"
#define CH_SPACEDIM PDIM

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FluxSurface.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"
namespace CFG = CFG_NAMESPACE;
namespace VEL = VEL_NAMESPACE;

class Kernel;
class MomentOp;

/// Kinetic species interface.
/**
 * Concrete reification of abstract SpeciesModel interface for kinetic
 * species.
 */
class KineticSpecies
{
   public:

      /// Constructor.
      /**
       */
      KineticSpecies(
         const string&    name,
         const Real       mass,
         const Real       charge,
         const PhaseGeom& geometry );

      /// Copy Constructor.
      /**
       * @param[in] foo KineticSpecies from which to initialize.
        */
      KineticSpecies( const KineticSpecies& foo );

      /// Destructor.
      /**
       */
      ~KineticSpecies() {;}

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
      virtual void numberDensity( CFG::LevelData<CFG::FArrayBox>& rho ) const;

      /// Returns species mass density.
      /**
       */
      virtual void massDensity( CFG::LevelData<CFG::FArrayBox>& rho ) const;

      /// Returns species charge density.
      /**
       */
      virtual void chargeDensity( CFG::LevelData<CFG::FArrayBox>& rho ) const;

      /// Returns species momentum density.
      /**
       */
      virtual void momentumDensity( CFG::LevelData<CFG::FArrayBox>& momentum ) const;

      /// Returns species parallel momentum (n * Vpar).
      /**
       */
      virtual void ParallelMomentum( CFG::LevelData<CFG::FArrayBox>& Parallel_Vel ) const;

      /// Returns species poloidal momentum (n * Vpol).
      /**
       */
      virtual void PoloidalMomentum( CFG::LevelData<CFG::FArrayBox>& Poloidal_Vel,
                                     const LevelData<FluxBox>& field,
                                     const double larmor  ) const; 

      /// Returns species radial particle flux.
      /**
       */
      virtual void ParticleFlux( CFG::LevelData<CFG::FArrayBox>& Particle_Flux,
                                 const LevelData<FluxBox>& field  ) const;

      /// Returns species radial heat flux.
      /**
       */
      virtual void HeatFlux( CFG::LevelData<CFG::FArrayBox>& Heat_Flux,
                                   const LevelData<FluxBox>& field,
                                   const LevelData<FArrayBox>& phi ) const;

      /// Returns species parallel heat flux.
      /**
       */
      virtual void parallelHeatFluxMoment( CFG::LevelData<CFG::FArrayBox>& parallelHeatFlux,
                                           CFG::LevelData<CFG::FArrayBox>& vparshift ) const;

      /// Returns species pressure.
      /**
       */
      virtual void pressureMoment( CFG::LevelData<CFG::FArrayBox>& pressure,
                                   CFG::LevelData<CFG::FArrayBox>& vparshift ) const;

      /// Returns species fourth velocity space moment.
      /**
       */
      virtual void fourthMoment( CFG::LevelData<CFG::FArrayBox>& fourthMoment ) const;

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
      virtual void addData( const KineticSpecies& rhs,
                            const Real          factor = 1.0 );

      /// Conformity verification
      /**
       * Checks that the rhs and this SpeciesModel have a similar
       * storage representation and returns true if so, false if not.
       *
       * @param[in] rhs SpeciesModel with which to compare.
       * @param[in] include_ghost_cells include ghost cells in comparison.
       */
      virtual bool conformsTo( const KineticSpecies& rhs,
                               const bool include_ghost_cells = true ) const;

      /// Factory method
      /**
       * Dynamically creates a new KineticSpecies with conforming layout and
       * copies internal parameters.
       *
       * @param[in] copy_soln_data Copies internal LevelDatas; otherwise,
       *                           initial LevelDatas are unititialized.
       * @param[in] include_ghost_cells Create internal LevelDatas with
       *                                ghost cells in comparison.
       */
      virtual
      RefCountedPtr<KineticSpecies> clone( const IntVect ghostVect = IntVect::Zero,
                                           const bool include_ghost_cells = true ) const;

      /// Copy
      /**
       * @param[in] rhs KineticSpecies from which to copy.
       */
      virtual void copy( const KineticSpecies& rhs );

      /// Assignment operator.
      /**
       * @param[in] rhs KineticSpecies from which to copy.
       */
      const KineticSpecies& operator=( const KineticSpecies& rhs );

      /// Provides the distribution function.
      /**
       */
      LevelData<FArrayBox>& distributionFunction() { return m_dist_func; }

      /// Provides a const reference to the distribution function.
      /**
       */
      const LevelData<FArrayBox>& distributionFunction() const { return m_dist_func; }

      /// Provides a const reference to the phase space geometry object.
      /**
       */
      const PhaseGeom& phaseSpaceGeometry() const { return m_geometry; }

      Real maxValue() const;
      Real minValue() const;

      void computeVelocity( LevelData<FluxBox>& a_velocity,
                            const LevelData<FluxBox>& a_E_field ) const;

      void computeMappedVelocity( LevelData<FluxBox>& a_velocity,
                                  const LevelData<FluxBox>& a_E_field ) const;

   private:

      const PhaseGeom& m_geometry;
      string m_name;
      Real m_mass;
      Real m_charge;
      LevelData<FArrayBox> m_dist_func;
      MomentOp& m_moment_op;
};

typedef Vector<RefCountedPtr<KineticSpecies> > KineticSpeciesPtrVect;

typedef RefCountedPtr<KineticSpecies> KineticSpeciesPtr;

#include "NamespaceFooter.H"

#endif
