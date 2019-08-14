#ifndef _KINETIC_SPECIES_H_
#define _KINETIC_SPECIES_H_

#include "Dimensions.H"
#include "PhaseGeom.H"
//#include "ICInterface.H.multidim"
//#include "BCInterface.H.multidim"
#include "Vector.H"
#include "RefCountedPtr.H"

#include "GyroaverageOperator.H"

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
         const string&                   name,
         const Real                      mass,
         const Real                      charge,
         const RefCountedPtr<PhaseGeom>& geometry,
         const bool                      is_gyrokinetic = false );

      /// Copy Constructor.
      /**
       * @param[in] foo KineticSpecies from which to initialize.
        */
      KineticSpecies( const KineticSpecies& foo );

      /// Destructor.
      /**
       */
      virtual ~KineticSpecies() {;}

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

      /// Returns gyroaveraged species charge density.
      /**
       */
      virtual void gyroaveragedChargeDensity( CFG::LevelData<CFG::FArrayBox>& rho ) const;

      /// Returns species momentum density.
      /**
       */
      virtual void momentumDensity( CFG::LevelData<CFG::FArrayBox>& momentum,
				    const LevelData<FluxBox>& field,
				    const double larmor ) const;

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

      /// Returns species perpendicular current density.
      /**
       */
      virtual void perpCurrentDensity( CFG::LevelData<CFG::FArrayBox>& Perp_Current_Density,
                                       const LevelData<FluxBox>& field  ) const;
   
      /// Returns species pressure.
      /**
       */
      virtual void pressureMoment( CFG::LevelData<CFG::FArrayBox>& pressure,
                                   CFG::LevelData<CFG::FArrayBox>& vparshift ) const;

      /// Returns species energy.
      /**
       */
      virtual void energyMoment( CFG::LevelData<CFG::FArrayBox>& energy ) const;

      /// Returns species energy for an arbitrary FArrayBox
      /**
       */
      virtual void energyMoment( CFG::LevelData<CFG::FArrayBox>& energy,
                                 const LevelData<FArrayBox>& dfn) const;

      /// Returns species fourth velocity space moment.
      /**
       */
      virtual void fourthMoment( CFG::LevelData<CFG::FArrayBox>& fourthMoment ) const;

      /// Returns maxwellian fit with FS-averaged density and temperature.
      /**
       */
      virtual void computeFSavgMaxwellian( LevelData<FArrayBox>&  a_F0 ) const;
   
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

      /// Provides a const reference to a DisjointBoxLayout constructed by extending
      /// the distribution function DisjointBoxLayout grown in the configuration
      /// space directions at block boundaries by an amount given by the 
      /// distribution function ghostVect.  This is used to construct a LevelData
      /// for use in moment operations that include ghost cells and use a ReductionCopier
      /// (which, per Chombo policy, does not copy from ghost cells).
      /**
       */
      const DisjointBoxLayout& getGhostDBL() const;

      /// Provides a const reference to the phase space geometry object.
      /**
       */
      const PhaseGeom& phaseSpaceGeometry() const { return *m_geometry; }

      const RefCountedPtr<PhaseGeom>& phaseSpaceGeometryPtr() const { return m_geometry; }

      Real maxValue() const;
      Real minValue() const;

      void computeVelocity( LevelData<FluxBox>& velocity,
                            const LevelData<FluxBox>& E_field,
                            const int  velocity_option,
                            const Real& time,
                            const bool  apply_axisymmetric_correction = true  ) const;

      void computeMappedVelocityNormals( LevelData<FluxBox>& velocity,
                                         const CFG::LevelData<CFG::FArrayBox>& E_field_cell,
                                         const CFG::LevelData<CFG::FArrayBox>& phi_node,
                                         const bool fourth_order_Efield,
                                         const int  velocity_option ) const;

      void computeMappedVelocity( LevelData<FluxBox>& velocity,
                                  const LevelData<FluxBox>& E_field,
                                  const Real& time) const;

      inline
      virtual int size( bool a_count_ghosts = false )
      {
         int size(0);
         const DisjointBoxLayout& dbl( m_dist_func.disjointBoxLayout() );
         for (DataIterator dit( dbl.dataIterator() ); dit.ok(); ++dit) {
            const Box& box( a_count_ghosts ? m_dist_func[dit].box() : dbl[dit] );
            size += box.numPts() * m_dist_func.nComp();
         }
         return size;
      }

      inline void gyroaverageOp(const GyroaverageOperator* const a_op)
      {
        if (!m_is_gyrokinetic) {
          MayDay::Error("KineticSpecies::gyroaverageOp() called for non-gyrokinetic species!");
        }
        m_gyroavg_op = a_op;
        (*m_geometry).gyroaverageOp(a_op);
      }

      inline const GyroaverageOperator* const gyroaverageOp() const
      {
        if (!m_is_gyrokinetic) {
          MayDay::Error("KineticSpecies::gyroaverageOp() called for non-gyrokinetic species!");
        }
        return m_gyroavg_op;
      }

      inline void gyroaverage(LevelData<FArrayBox>& a_phi_bar,
                              const CFG::LevelData<CFG::FArrayBox>& a_phi) const
      {
        if (!m_is_gyrokinetic) {
          MayDay::Error("KineticSpecies::gyroaverage() called for non-gyrokinetic species!");
        }
        if (m_gyroavg_op == NULL) {
          MayDay::Error("Error in KineticSpecies::gyroaverage() - m_gyroavg_op is NULL.");
        }
        m_gyroavg_op->applyOp(a_phi_bar, a_phi);
      }

      inline void gyroaveragedEField( LevelData<FArrayBox>&                 a_E,
                                      const CFG::LevelData<CFG::FArrayBox>& a_phi,
                                      const int                             a_order ) const
      {
        if (!m_is_gyrokinetic) {
          MayDay::Error("KineticSpecies::gyroaveragedEfield() called for non-gyrokinetic species!");
        }
        if (m_gyroavg_op == NULL) {
          MayDay::Error("Error in KineticSpecies::gyroaveragedEfield() - m_gyroavg_op is NULL.");
        }
        m_gyroavg_op->gyroaveragedEField(a_E, a_phi, a_order);
      }

      inline void gyroaveragedEField( LevelData<FluxBox>&                   a_E,
                                      const CFG::LevelData<CFG::FArrayBox>& a_phi,
                                      const int                             a_order ) const
      {
        if (!m_is_gyrokinetic) {
          MayDay::Error("KineticSpecies::gyroaveragedEfield() called for non-gyrokinetic species!");
        }
        if (m_gyroavg_op == NULL) {
          MayDay::Error("Error in KineticSpecies::gyroaveragedEfield() - m_gyroavg_op is NULL.");
        }
        m_gyroavg_op->gyroaveragedEField(a_E, a_phi, a_order);
      }

      inline bool isGyrokinetic() const
      {
        return m_is_gyrokinetic;
      }

   private:

      const RefCountedPtr<PhaseGeom> m_geometry;
      string m_name;
      Real m_mass;
      Real m_charge;
   
      LevelData<FArrayBox> m_dist_func;
      MomentOp& m_moment_op;
      mutable DisjointBoxLayout m_ghost_dbl;
      mutable LevelData<FluxBox> m_velocity;
      mutable Real m_time;
      mutable int m_velocity_option;
   
      GyroaverageOperator const *m_gyroavg_op;
      bool m_is_gyrokinetic;
};

typedef Vector<RefCountedPtr<KineticSpecies> > KineticSpeciesPtrVect;

typedef RefCountedPtr<KineticSpecies> KineticSpeciesPtr;

#include "NamespaceFooter.H"

#endif
