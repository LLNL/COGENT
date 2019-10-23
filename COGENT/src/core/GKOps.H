#ifndef  _GKOPS_H_
#define  _GKOPS_H_

#include <vector>
#include <map>
#include <string>

#include "GKVlasov.H"
#include "GKCollisions.H"
#include "GKNeutrals.H"
#include "GKTransport.H"
#include "GKUnit.H"
#include "GKSystemIC.H"
#include "GKSystemBC.H"
#include "CH_HDF5.H"
#include "GKVector.H"
#include "GKScalarOp.H"
#include "GyroaverageOperatorFactory.H"
#include "Checkpoints.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "GKFluidOp.H"
#include "EFieldAmpere.H"
#include "VorticityOp.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

typedef enum {unsplit, imex_exp, imex_imp} opType;

/**
 * Gyrokinetic System class.
 *
 * A concrete implementation of the ODEInterface<Vector<SpeciesModel>>
 * interface, this class owns the simulation data and the Vlasov and Poisson
 * operator classes and orchestrates the time advancement.
 *
 * GKSystem implements parsing of some of the input data for gksystem.
 * Presently the list
 * for gksystem includes "verbosity", controlling the amount of diagnostic output,
 * "num_cells", the number of cells in the radial, poloidal, vparallel, and mu directions,
 * "is_periodic", an integer array of zeros and ones for periodicity in the coordinate
 * directions (same order). some parameters to control specification of what hdf5
 * files to write for plotting (see below), and "magneic_geometry_mapping", the name of the
 * magnetic geometry.   Specific parameters of the geometry are parsed in the geometry
 * implementation.  The plotting parameters are hdf_density (default true),
 * df_vpartheta (default false), and hdf_frtheta (default false), to control creation
 * of hdf5 files for density, f versus vparallel and theta,  and f versus r and theta, as well
 * as m_fixed_plotindices, an array of length 5 denoting the radial, poloidal, toroidal,
 * vparallel and mu indices
 *
 * Sample input:
 * \verbatim
 * gksystem.verbosity   = 1
 * gksystem.num_cells   = 32 128 4 4
 * gksystem.is_periodic =  0  1  0  0
 * gksystem.magnetic_geometry_mapping = "miller"
 * gksystem.hdf_density = true
 * gkysytem.hdf_vpartheta = false
 * gkysytem.hdf_frtheta = true
 * gksystem.fixed_plot_indices = 3 1 0 2 2
 *
 */

class GKOps
{
public:

   /// Constructor with initialization.
   /**
    *  Pure virtual function that will be resposible for setting
    *  initial conditions on the soln data.
    *
    *  @param[in] pp the input database.
    */

   void initializeElectricField( const GKState& state, const int cur_step, const double cur_time );
 
   void initializeFluidSpeciesPhysical( const CFG::FluidSpeciesPtrVect&  a_species_comp );

   inline
   void initializeState( GKState& a_state, const int a_cur_step, const Real a_time )
   {
      if ( a_cur_step > 0 && m_old_vorticity_model ) {
         // If restarting the old vorticity model, then the potential has already been set
         m_initial_conditions->initializeSystem( a_state, a_time );
      }
      else {
         m_initial_conditions->initializeSystem( a_state, m_phi, a_time );
      }
      m_fluidOp->initialize( a_state.dataFluid(), a_time );
      enforceQuasiNeutrality( a_state.dataKinetic(), m_phi );
   }

  inline
   void initializePotential( const Real a_time )
   {
      m_initial_conditions->initializePotentialField( m_phi, a_time );
   }
   
   /// Destructor.
   /**
    */
   virtual ~GKOps();

   /// Returns the Vlasov time scale.
   /**
    * Returns the Vlasov time scale.
    */
   Real dtScaleVlasov( const GKState& state, const int step_number )
   {
      return m_time_scale_vlasov;
   }

   /// Returns the collisions time scale.
   /**
    * Returns the collisions time scale.
    */
   Real dtScaleCollisions( const GKState& a_state, const int step_number )
   {
      return m_time_scale_collisions;
   }

   /// Returns the transport time scale.
   /**
    * Returns the transport time scale.
    */
   Real dtScaleTransport( const GKState& a_state,  const int step_number )
   {
      return m_time_scale_transport;
   }

   /// Returns the time scale for the fluid components.
   /**
    * Returns the time scale for the fluid components.
    */
   Real dtScaleFluids( const GKState& a_state,  const int step_number )
   {
      return m_time_scale_fluids;
   }
   
   /// Returns the time scale for the scalar components.
   /**
    * Returns the time scale for the scalar components.
    */
   Real dtScaleScalars( const GKState& a_state,  const int step_number )
   {
      return m_time_scale_scalars;
   }
   
   /// Returns the neutrals time scale.
   /**
    * Returns the neutrals time scale.
    */
   Real dtScaleNeutrals( const GKState& a_state, const int step_number )
   {
      return m_time_scale_neutrals;
   }

   /// Compute a stable time step for explicit time integration
   /**
    * Computes and returns a stable time step estimate.
    */
   Real stableDtExpl( const GKState& a_state, const int step_number );

   /// Compute a stable time step for IMEX time integration
   /**
    * Computes and returns a stable time step estimate.
    */
   Real stableDtImEx( const GKState& a_state, const int step_number );

   /// Evaluates the RHS of the ODE.
   /**
    * Concrete implementation of pure virtual base class member function.
    * Given the current time and solution, evaluates the RHS of the ODE.
    *
    * @param[out] rhs data holder for rhs (without ghost cells).
    * @param[in] soln current solution with which to evaluate RHS.
    * @param[in] time current time.
    * @param[in] stage  stage number
    */
   
   GKOps() :
      m_vlasov(NULL),
      m_collisions(NULL),
      m_neutrals(NULL),
      m_transport(NULL),
      m_fluidOp(NULL),
      m_scalarOp(NULL),
      m_boltzmann_electron(NULL),
      m_units(NULL),
      m_initial_conditions(NULL),
      m_boundary_conditions(NULL),
      m_dt(0.0),
      m_dt_scale(1.0),
      m_is_defined(false),
      m_ghost_vect(4*IntVect::Unit),
      m_kinetic_ghost_vect(4*IntVect::Unit),
      m_cfg_nghosts(2*CFG::IntVect::Unit),
      m_fixed_efield(false),
      m_zero_efield(false), // JRA
      m_transport_model_on(false),     
      m_enforce_quasineutrality(false),
      m_step_const_coef(false),
      m_history(false),
      m_hist_freq(1),
      m_last_hist(0),
      m_hist_count(0),
      m_expand_incr(100),
      m_hist_fieldname("field"),
      m_dt_vlasov_expl(DBL_MAX),
      m_dt_collisions_expl(DBL_MAX),
      m_dt_transport_expl(DBL_MAX),
      m_dt_fluids_expl(DBL_MAX),
      m_dt_scalars_expl(DBL_MAX),
      m_dt_neutrals_expl(DBL_MAX),
      m_dt_vlasov_imex(DBL_MAX),
      m_dt_collisions_imex(DBL_MAX),
      m_dt_transport_imex(DBL_MAX),
      m_dt_fluids_imex(DBL_MAX),
      m_dt_scalars_imex(DBL_MAX),
      m_dt_neutrals_imex(DBL_MAX),
      m_time_scale_vlasov(-1),
      m_time_scale_collisions(-1),
      m_time_scale_transport(-1),
      m_time_scale_fluids(-1),
      m_time_scale_scalars(-1),
      m_time_scale_neutrals(-1),
      m_count_vlasov(0),
      m_count_collision(0),
      m_count_transport(0),
      m_count_fluids(0),
      m_count_scalars(0),
      m_count_neutrals(0),
      m_verbosity(0) {;}   
   
   void define( const GKState& a_state, const Real a_dt );
   
   bool isDefined() const { return m_is_defined; }
   
   void resetDt( const Real a_dt )
   {
      CH_assert( a_dt>0.0 );
      CH_assert( isDefined() );
      m_dt = a_dt;
   }

   void preTimeStep  (const int,const Real,const GKState&,const GKState&);
   void postTimeStep (const int,const Real, const Real,const GKState&);

   void postTimeStage(const int,const Real,const GKState&, const int);
   inline void postTimeStage( const int        a_step,
                              const Real       a_time,
                              const GKVector&  a_vec,
                              const int        a_stage )
   {
      m_Y.copyFrom(a_vec.data());
      postTimeStage(a_step,a_time,m_Y,a_stage);
   }
   
   void setElectricField( const GKState&                   state_comp,
                          const int                        step,
                          CFG::LevelData<CFG::FArrayBox>&  phi,
                          CFG::EField&                     E_field) const;

   inline bool trivialSolutionOp() const { return m_trivial_solution_op; }

   void preOpEval( const GKState&, const Real, const Checkpoint&);
   inline void preOpEval(const GKVector&    a_state,
                         const Real         a_time,
                         const Checkpoint&  a_chkpt )
   {
      m_Y.copyFrom(a_state.data());
      preOpEval(m_Y, a_time, a_chkpt);
   }

   void preSolutionOpEval( const GKState&, const Real, const Checkpoint&);
   inline void preSolutionOpEval(const GKVector&    a_state,
                                 const Real         a_time,
                                 const Checkpoint&  a_chkpt )
   {
      m_Y.copyFrom(a_state.data());
      preSolutionOpEval(m_Y, a_time, a_chkpt);
   }

   void explicitOp( GKRHSData&, const Real, const GKState&);
   inline void explicitOp(GKVector&       a_rhs, 
                          const Real      a_time, 
                          const GKVector& a_state_comp)
   {
      m_Y.copyFrom(a_state_comp.data());
      explicitOp(m_rhs,a_time,m_Y);
      m_rhs.copyTo(a_rhs.data());
   }

   void explicitOpImEx( GKRHSData&, const Real, const GKState&);
   inline void explicitOpImEx(GKVector&        a_rhs,
                              const Real       a_time,
                              const GKVector&  a_state_comp )
   {
      m_Y.copyFrom(a_state_comp.data());
      explicitOpImEx(m_rhs,a_time,m_Y);
      m_rhs.copyTo(a_rhs.data());
   }

   void implicitOpImEx( GKRHSData&, const Real, const GKState&);
   inline void implicitOpImEx(GKVector&        a_rhs,
                              const Real       a_time,
                              const GKVector&  a_state_comp )
   {
      m_Y.copyFrom(a_state_comp.data());
      implicitOpImEx(m_rhs, a_time, m_Y);
      m_rhs.copyTo(a_rhs.data());
   }

   void explicitPC( GKRHSData&, const Real, const GKState&);
   inline void explicitPC(GKVector&       a_rhs,
                          const Real      a_time,
                          const GKVector& a_state )
   {
      m_Y.copyFrom(a_state.data());
      explicitPC(m_rhs,a_time,m_Y);
      m_rhs.copyTo(a_rhs.data());
   }

   void implicitPCImEx( GKRHSData&, const Real, const GKState&);
   inline void implicitPCImEx(GKVector&        a_rhs,
                              const Real       a_time,
                              const GKVector&  a_state_comp )
   {
      m_Y.copyFrom(a_state_comp.data());
      implicitPCImEx(m_rhs,a_time,m_Y);
      m_rhs.copyTo(a_rhs.data());
   }

   void solutionOp( GKRHSData&, const Real, const GKState&);
   inline void solutionOp(GKVector&       a_rhs,
                          const Real      a_time,
                          const GKVector& a_state )
   {
      m_Y.copyFrom(a_state.data());
      solutionOp(m_rhs,a_time,m_Y);
      m_rhs.copyTo(a_rhs.data());
   }


   void solutionPC( GKRHSData&, const Real, const GKState&);
   inline void solutionPC(GKVector&       a_rhs,
                          const Real      a_time,
                          const GKVector& a_state )
   {
      m_Y.copyFrom(a_state.data());
      solutionPC(m_rhs,a_time,m_Y);
      m_rhs.copyTo(a_rhs.data());
   }

   void defineMultiPhysicsPC( std::vector<Preconditioner<GKVector,GKOps>*>&,
                              std::vector<DOFList>&,
                              const GKVector&,
                              const std::string&,
                              const std::string&,
                              bool);

   void updateMultiPhysicsPC( std::vector<Preconditioner<GKVector,GKOps>*>&,
                              const GKVector&,
                              const Real,
                              const Real,
                              const bool );

   bool preallocatePCMatrixImEx (void*,const GKVector&);
   void assemblePCImEx( const bool, void*, const GKVector&, const Real, const Real );

   bool preallocateSolPCMatrix (void*,const GKVector&);
   void assembleSolPC( const bool, void*, const GKVector&, const Real );
   
   void solveSolutionPC( GKVector&, const GKVector&, int a_idx=-1);
   void solvePCImEx( GKVector&, const GKVector&, int a_idx=-1);

   void divideJ( const KineticSpeciesPtrVect& soln_mapped,
                 KineticSpeciesPtrVect&       soln_physical ) const;

   void convertToPhysical( const GKState& soln_mapped, GKState& soln_physical );

   inline
   bool usingBoltzmannElectrons() const
   {
      return (m_boltzmann_electron!=NULL);
   };

   inline
   bool fixedEField() const { return m_fixed_efield; };

   inline bool usingAmpereLaw() const { return m_ampere_law; }

   double getLoRadialField(const GKState& state) const;
   void setLoRadialField(GKState& state, const double Er_lo) const;
   double getHiRadialField(const GKState& state) const;
   void setHiRadialField(GKState& state, const double Er_hi) const;

   const CFG::LevelData<CFG::FluxBox>& getErAverageFace() const
   {
      CFG::EFieldAmpere& efield = static_cast<CFG::EFieldAmpere&>(*m_E_field);
      return efield.getErAverageFace();
   }

   const CFG::LevelData<CFG::FArrayBox>& getErAverageCell() const
   {
      CFG::EFieldAmpere& efield = static_cast<CFG::EFieldAmpere&>(*m_E_field);
      return efield.getErAverageCell();
   }

   const CFG::LevelData<CFG::FluxBox>& getETildeFace() const
   {
      CFG::EFieldAmpere& efield = static_cast<CFG::EFieldAmpere&>(*m_E_field);
      return efield.getETildeFace();
   }

   const CFG::LevelData<CFG::FArrayBox>& getETildeCell() const
   {
      CFG::EFieldAmpere& efield = static_cast<CFG::EFieldAmpere&>(*m_E_field);
      return efield.getETildeCell();
   }

   const CFG::LevelData<CFG::FArrayBox>& getPhi() const
   {
      return m_phi;
   }

   void setErAverage(const LevelData<FluxBox>& Er_face_injected);
   void setErAverage(const LevelData<FArrayBox>& Er_cell_injected);
   void setETilde( const LevelData<FluxBox>& E_tilde_face_injected );
   void setETilde( const LevelData<FArrayBox>& E_tilde_cell_injected);
   void setPhi( const LevelData<FArrayBox>& phi_injected);
    
   void plotDistributionFunction( const std::string&    filename,
                                  const KineticSpecies& soln_species,
                                  const double&         time ) const;

   void plotFluid( const std::string&       filename,
                   const CFG::FluidSpecies& fluid_species,
                   const string&            var_name,
                   const double&            time ) const;

   void plotDeltaF( const std::string&    filename,
                    const KineticSpecies& soln_species,
                    const double&         time ) const;

   void plotDistributionFunctionAtMu( const std::string&    filename,
                                      const KineticSpecies& soln_species,
                                      const int             mu,
                                      const double&         time ) const;

   void plotBStarParallel( const std::string&    filename,
                           const KineticSpecies& soln_species,
                           const double&         time ) const;

   void plotVParallelTheta( const std::string&    filename,
                            const KineticSpecies& soln_species,
                            const int             radial_index,
                            const int             toroidal_index,
                            const int             mu_index,
                            const double&         time ) const;

   void plotRTheta( const std::string&    filename,
                    const KineticSpecies& soln_species,
                    const int             vpar_index,
                    const int             mu_index,
                    const double&         time ) const;

   void plotVParallelMu( const std::string&    filename,
                         const KineticSpecies& soln_species,
                         const int             radial_index,
                         const int             poloidal_index,
                         const int             toroidal_index,
                         const double&         time ) const;

   void plotChargeDensity( const std::string&    filename,
                           const KineticSpecies& soln_species,
                           const double&         time ) const;

   void plotChargeDensity( const std::string&     filename,
                           const KineticSpeciesPtrVect& species,
                           const double&          time ) const;

   void plotMomentum( const std::string&     filename,
		      const KineticSpecies&  soln_species,
		      const double&          time ) const;

   void plotParallelMomentum( const std::string&    filename,
                              const KineticSpecies& soln_species,
                              const double&         time ) const;

   void plotPoloidalMomentum( const std::string&    filename,
                              const KineticSpecies& kinetic_species,
                              const double&         time ) const;
   
   void plotPressure( const std::string&    filename,
                      const KineticSpecies& soln_species,
                      const double&         time ) const;
   
   void plotParallelHeatFlux( const std::string&    filename,
                              const KineticSpecies& soln_species,
                              const double&         time ) const;
   

   void plotTemperature( const std::string&    filename,
                         const KineticSpecies& soln_species,
                         const double&         time ) const;

   void plotFourthMoment( const std::string&    filename,
                          const KineticSpecies& soln_species,
                          const double&         time ) const;

   void plotParticleFlux( const std::string&    filename,
                          const KineticSpecies& soln_species,
                          const double&         time ) const;

   void plotHeatFlux( const std::string&    filename,
                      const KineticSpecies& soln_species,
                      const double&         time ) const;

   void plotPotential( const std::string& filename,
                       const bool&        non_zonal_comp,
                       const double&      time ) const;

   void plotEField( const std::string& filename,
                    const double&      time ) const;

   void plotAmpereErIncrement( const std::string&                filename,
                               const CFG::FluidSpeciesPtrVect&   soln_fluid,
                               const ScalarPtrVect&              soln_scalar,
			       const double&                     time ) const;

   void setupFieldHistories( ParmParse& a_ppsim );

   void writeFieldHistory( int cur_step, double cur_time, bool startup_flag );
  
      /// Write checkpoint file.
      /**
       * Write checkpoint data to an output HDF5 file.
       *
       * @param[in] handle a reference to an HDF5 plot file.
       */
   void writeCheckpointFile( HDF5Handle& handle ) const;
   
      /// Read checkpoint file.
      /**
       * Read checkpoint data from an output HDF5 file and reinitialize.
       *
       * @param[in] handle a reference to an HDF5 plot file.
       */
   void readCheckpointFile( HDF5Handle& handle, const int& cur_step  );

   bool stateIncludesPotential( const GKState& a_state ) const;

   //   inline bool isLinear() { return(m_collisions->isLinear()); }
   inline bool isLinear() { return false; }

   void printFunctionCounts()
   {
    if (!procID()) {
      cout << "  Function counts:-\n";
      cout << "    Vlasov     : " << m_count_vlasov     << "\n";
      cout << "    Collision  : " << m_count_collision  << "\n";
      cout << "    Transport  : " << m_count_transport  << "\n";
      cout << "    Fluids     : " << m_count_fluids     << "\n";
      cout << "    Scalars    : " << m_count_scalars    << "\n";
      cout << "    Neutrals   : " << m_count_neutrals   << "\n";
    }
   }

   inline void gyroaverage(LevelData<FArrayBox>&                  a_phi_bar,
                           const CFG::LevelData<CFG::FArrayBox>&  a_phi,
                           const std::string&                     a_species_name )
   {
     std::map<std::string, GyroaverageOperator*>::iterator it;
     it = m_gyroavg_ops.find(a_species_name);
     if ((it == m_gyroavg_ops.end()) || (it->second == NULL)) {
       MayDay::Error("Error in GKOps::gyroaverage() - cannot find gyroaveraging operator specified species name.");
     }
     m_gyroavg_ops[a_species_name]->applyOp(a_phi_bar, a_phi);
   }

   inline const std::map<std::string,GyroaverageOperator*>& getGyroaverageOps() const
   {
     return m_gyroavg_ops;
   }
   inline std::map<std::string,GyroaverageOperator*>& getGyroaverageOps()
   {
     return m_gyroavg_ops;
   }

   /// Return kinetic species ghost vect size require for methods used.
   /**
    * Private method to return kinetic species ghost vect.
    */
   void getKineticSpeciesGhostVect(IntVect& ghost_vect);

   inline bool stepConstCoef() const { return m_step_const_coef; }
   
private:
   
   /// Parse parameters.
   /**
    * Private method to obtain control parameters from "gksystem" section
    * of ParmParse database.
    *
    * \param[in] ppgksys gksystem input database.
    */
   void parseParameters( ParmParse& ppgksys );

   void createPhysicalSpeciesVector( KineticSpeciesPtrVect&        out,
                                     const KineticSpeciesPtrVect&  in,
                                     const double                  time );
   
   void updatePhysicalSpeciesVector( const CFG::FluidSpeciesPtrVect&  in,
                                     const double                     time );
   
   void applyVlasovOperator( GKRHSData&                    rhs,
                             const KineticSpeciesPtrVect&  kinetic_phys,
                             const opType                  op_type,
                             const Real&                   time );

   void applyKineticSpeciesSolutionOperator( KineticSpeciesPtrVect&       rhs,
                                             const KineticSpeciesPtrVect& soln,
                                             const Real&                  time );

   void applyFluidSpeciesSolutionOperator( CFG::FluidSpeciesPtrVect&        rhs,
                                           const KineticSpeciesPtrVect&     kinetic_species_comp,
                                           const KineticSpeciesPtrVect&     kinetic_species_phys,
                                           const CFG::FluidSpeciesPtrVect&  fluid_species_comp,
                                           const CFG::FluidSpeciesPtrVect&  fluid_species_phys,
                                           const ScalarPtrVect&             scalars,
                                           const Real&                      time );

   void applyScalarSolutionOperator( ScalarPtrVect&             a_rhs,
                                     const ScalarPtrVect&       a_soln,
                                     const Real&                a_time );

   void applyCollisionOperator( KineticSpeciesPtrVect&       rhs,
                                const KineticSpeciesPtrVect& soln,
                                const opType                 op_type,
                                const Real&                  time );

   void applyNeutralsOperator(   KineticSpeciesPtrVect&      rhs,
                                const KineticSpeciesPtrVect& soln,
                                const Real&                  time );
    

   void applyTransportOperator( KineticSpeciesPtrVect&       rhs,
                                const KineticSpeciesPtrVect& soln,
                                const Real&                  time );

   void applyFluidOperator( CFG::FluidSpeciesPtrVect&            rhs,
                            const KineticSpeciesPtrVect&         kinetic_species_phys,
                            const CFG::FluidSpeciesPtrVect&      fluid_species,
                            const ScalarPtrVect&                 scalars,
                            const CFG::EField&                   E_field,
                            const opType                         op_type,
                            const Real&                          time);
   
   void applyScalarOperator( GKRHSData&                           rhs,
                             const KineticSpeciesPtrVect&         kinetic_species,
                             const CFG::FluidSpeciesPtrVect&      fluid_species,
                             const ScalarPtrVect&                 scalars,
                             const CFG::EField&                   E_field,
                             const opType                         op_type,
                             const bool                           recompute_kinetic_terms,
                             const Real&                          time);

   /// Compute the signed charge densities
   /**
    * Compute the signed charge densities
    *
    * @param[out] pos_charge_density    Ion charge density
    * @param[out] neg_charge_density    Electron charge density
    * @param[in]  species               Species vector
    */

   void computeSignedChargeDensities( CFG::LevelData<CFG::FArrayBox>& pos_charge_density,
                                      CFG::LevelData<CFG::FArrayBox>& neg_charge_density,
                                      const KineticSpeciesPtrVect&    species ) const;

   void enforceQuasiNeutrality( KineticSpeciesPtrVect&          species,
                                CFG::LevelData<CFG::FArrayBox>& potential ) const;
   
   inline
   void computeQuasiNeutralElectronDensity(
      CFG::LevelData<CFG::FArrayBox>&       a_quasineutral_density,
      CFG::LevelData<CFG::FArrayBox>&       a_potential,
      const CFG::EllipticOpBC&              a_bc, 
      const CFG::LevelData<CFG::FArrayBox>& a_ion_density) const;
   
   inline
   void renormalizeElectronDistributionFunction(
      KineticSpecies&                 a_electron_species,
      CFG::LevelData<CFG::FArrayBox>& a_initial_density,
      CFG::LevelData<CFG::FArrayBox>& a_quasineutral_density ) const;

   GKVlasov*                      m_vlasov;
   GKCollisions*                  m_collisions;
   GKNeutrals*                    m_neutrals;
   GKTransport*                   m_transport;
   CFG::GKFluidOp*                m_fluidOp;
   GKScalarOp*                    m_scalarOp;
   CFG::VorticityOp*              m_VorticityOp;

   CFG::EField*                   m_E_field;
   CFG::EField*                   m_E_field_ImOpImEx;
   CFG::LevelData<CFG::FArrayBox> m_phi;

   CFG::BoltzmannElectron*        m_boltzmann_electron;

   KineticSpeciesPtrVect          m_kinetic_species_phys;
   CFG::FluidSpeciesPtrVect       m_fluid_species_phys;

   const GKUnits*                 m_units;
   const GKSystemIC*              m_initial_conditions;
   const GKSystemBC*              m_boundary_conditions;
   RefCountedPtr<PhaseGeom>       m_phase_geometry;

   Real m_dt;
   Real m_dt_scale;
   bool m_is_defined;
   
   IntVect                        m_ghost_vect;
   IntVect                        m_kinetic_ghost_vect;
   const CFG::IntVect             m_cfg_nghosts;
      
   bool                           m_consistent_potential_bcs;
   bool                           m_ampere_post_step_update;
   bool                           m_ampere_law;
   bool                           m_ampere_cold_electrons;
   bool                           m_fixed_efield;
   bool                           m_zero_efield; // JRA
   bool                           m_transport_model_on;
   bool                           m_neutrals_model_on;
   bool                           m_enforce_quasineutrality;
   bool                           m_Esol_extrapolation;
   bool                           m_dealignment_corrections;
   bool                           m_old_vorticity_model;
   
   /* The following 4 flags are relevant only for IMEX time integration */
   bool                           m_transport_imex_implicit;
   bool                           m_neutrals_imex_implicit;

   GKState                        m_Y;
   GKRHSData                      m_rhs;

   bool m_state_includes_potential;
   bool m_trivial_solution_op;
   bool m_step_const_coef;

   bool m_history;   // whether to write out histories
   int m_hist_freq;  // how often to write out field history
   int m_last_hist;  // moribund, not used
   int m_hist_count; // how many watchpoints
   int m_expand_incr; // expansion increment for vectors in FieldHist
   CFG::IntVect m_hist_indices; // for indexing watchpoints, now a temporary
   std::string m_hist_fieldname; // type of field, e.g. "potential", now a temporary

   /* stable time step sizes and time scales */

   Real m_dt_vlasov_expl;
   Real m_dt_collisions_expl;
   Real m_dt_transport_expl;
   Real m_dt_fluids_expl;
   Real m_dt_scalars_expl;
   Real m_dt_neutrals_expl;

   Real m_dt_vlasov_imex;
   Real m_dt_collisions_imex;
   Real m_dt_transport_imex;
   Real m_dt_fluids_imex;
   Real m_dt_scalars_imex;
   Real m_dt_neutrals_imex;

   Real m_time_scale_vlasov;
   Real m_time_scale_collisions;
   Real m_time_scale_transport;
   Real m_time_scale_fluids;
   Real m_time_scale_scalars;
   Real m_time_scale_neutrals;

   // stores data pertaining to a watchpoint history
   typedef struct {
      int hist_index; // numerical index of this record, zero-based
      CFG::IntVect grid_indices; // location of watchpoint in grid
      Vector<Real> *fieldvals; // field value, e.g., potential
      Vector<Real> *timevals; // time corresponding to above
      Vector<int> *timestep; // timestep corresponding to above
      int valsize; // length of above three fields
      int cur_entry; //index into fieldvars and timevals
      string fieldname; // e.g., "potential"
   } FieldHist;
   
   Vector<FieldHist> m_fieldHistLists; // for saving all histories internally
   
   /* Function counters */
   int m_count_vlasov;
   int m_count_collision;
   int m_count_transport;
   int m_count_fluids;
   int m_count_scalars;
   int m_count_neutrals;

   int m_verbosity;

   /* gyroaverage operators */
   std::map<std::string,GyroaverageOperator*> m_gyroavg_ops;
};

#include "NamespaceFooter.H"

#endif
