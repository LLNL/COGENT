#ifndef  _GKSYSTEM_H_
#define  _GKSYSTEM_H_

#include "AppCtxt.H"
#include "GKOps.H"
#include "evalOpMatrix.H"
#include "PositivityPostProcessor.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "MagGeom.H"
#include "FloorPostProcessor.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

typedef Vector<RefCountedPtr<LevelData<FArrayBox> > > VelocityPtrVect;
typedef RefCountedPtr<LevelData<FArrayBox> > VelocityPtr;

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

class GKSystem : public AppCtxt
{
   public:

      /// Constructor with initialization.
      /**
       *  Pure virtual function that will be resposible for setting
       *  initial conditions on the soln data.
       *
       *  @param[in] pp                     the input database.
       *  @param[in] a_useExternalTI   Flag to indicate if using native (false, default) 
       *                                    or external (true) time integrators. 
       */
      GKSystem( ParmParse& pp );

      /// Destructor.
      /**
       */
      virtual ~GKSystem();

      virtual void initialize( const int cur_step, const double cur_time );
   
      /// Compute a stable time step.
      /**
       * Computes and returns a stable time step estimate.
       */
      Real stableDt( const int step_number, 
                     const TimeIntegratorType );

      /// Write plot file.
      /**
       * Write plotting data to HDF5 files with a common prefix.
       *
       */
      void writePlotFile( const char *prefix, const int cur_step, const double& cur_time );

      /// Augment a history file with a value at a spatial index
      /**
       *Augment a history file with a value at a spatial index
       */

      void writeFieldHistory(int cur_step, double cur_time, bool startup_flag);

      /// Write checkpoint file.
      /**
       * Write checkpoint data to an output HDF5 file.
       *
       * @param[in] handle a reference to an HDF5 plot file.
       */
      void writeCheckpointFile( HDF5Handle& handle, const int cur_step,
                                const double cur_time, const double cur_dt );

      /// Read checkpoint file.
      /**
       * Read checkpoint data from an output HDF5 file and reinitialize.
       *
       * @param[in] handle a reference to an HDF5 plot file.
       */
      void readCheckpointFile( HDF5Handle& handle, int& cur_step,
                               double& cur_time, double& cur_dt );

      void printDiagnostics();

      void preTimeStep  ( const ODEVector&, int, Real );
   
      void postTimeStep ( ODEVector&, int, Real, Real );

      inline void postTimeStage( const int     a_step,
                                 const Real    a_time,
                                 const Real    a_dt,
                                 ODEVector&    a_vec,
                                 const int     a_stage )
      {
         copyStateFromArray( a_vec );
         postTimeStage(a_step, a_dt, a_time, a_stage);
         copyStateToArray( a_vec );
      }

      void postTimeStage( const int, 
                          const Real, 
                          const Real, 
                          const int );

      int  getVectorSize()  { return m_state_comp.getVectorSize(); }

      inline void copyStateToArray (ODEVector& a_Y) const
      { 
        m_state_comp.copyTo(a_Y.data());   
      }

      inline void copyStateFromArray (const ODEVector& a_Y) 
      { 
        m_state_comp.copyFrom(a_Y.data()); 
      }

      inline void copyRHSToArray (ODEVector& a_Y) const
      { 
        m_rhs.copyTo(a_Y.data());   
      }
      
      inline void copyRHSFromArray (const ODEVector& a_Y) 
      { 
        m_rhs.copyFrom(a_Y.data()); 
      }

      inline void preOpEval(const ODEVector&   a_Y,
                            const Real         a_time,
                            const Checkpoint&  a_chkpt )
      {
        copyStateFromArray( a_Y );
        m_gk_ops->preOpEval(m_state_comp, a_time, a_chkpt);
      }

      inline void preSolutionOpEval(const ODEVector&   a_Y,
                                    const Real         a_time,
                                    const Checkpoint&  a_chkpt )
      {
         copyStateFromArray( a_Y );
         m_gk_ops->preSolutionOpEval(m_state_comp, a_time, a_chkpt);
      }

      inline void explicitOp(ODEVector&        a_F, 
                             const Real        a_time, 
                             const ODEVector&  a_Y)
      {
         copyStateFromArray( a_Y );
         m_gk_ops->explicitOp( m_rhs, a_time, m_state_comp );
         copyRHSToArray( a_F );
      }

      inline void explicitOpImEx(ODEVector&         a_F,
                                 const Real         a_time,
                                 const ODEVector&   a_Y )
      {
         copyStateFromArray( a_Y );
         m_gk_ops->explicitOpImEx( m_rhs, a_time, m_state_comp );
         copyRHSToArray( a_F );
      }

      inline void implicitOpImEx(ODEVector&         a_F,
                                 const Real         a_time,
                                 const ODEVector&   a_Y )
      {
         copyStateFromArray( a_Y );
         m_gk_ops->implicitOpImEx( m_rhs, a_time, m_state_comp );
         copyRHSToArray( a_F );
      }

      inline void solutionOp(ODEVector&        a_F,
                             const Real        a_time,
                             const ODEVector&  a_Y )
      {
         copyStateFromArray( a_Y );
         m_gk_ops->solutionOp( m_rhs, a_time, m_state_comp );
         copyRHSToArray( a_F );
      }

      inline void explicitPC(ODEVector&        a_F,
                             const Real        a_time,
                             const ODEVector&  a_Y )
      {
         copyStateFromArray( a_Y );
         m_gk_ops->explicitPC( m_rhs, a_time, m_state_comp );
         copyRHSToArray( a_F );
      }

      inline void implicitPCImEx(ODEVector&         a_F,
                                 const Real         a_time,
                                 const ODEVector&   a_Y )
      {
         copyStateFromArray( a_Y );
         m_gk_ops->implicitPCImEx( m_rhs, a_time, m_state_comp );
         copyRHSToArray( a_F );
      }

      inline void solutionPC(ODEVector&        a_F,
                             const Real        a_time,
                             const ODEVector&  a_Y )
      {
         copyStateFromArray( a_Y );
         m_gk_ops->solutionPC( m_rhs, a_time, m_state_comp );
         copyRHSToArray( a_F );
      }

      inline
      void defineMultiPhysicsPC(std::vector<Preconditioner<ODEVector,AppCtxt>*>&  a_pc,
                                std::vector<DOFList>&                             a_dof_list,
                                const ODEVector&                                  a_Y,
                                const std::string&                                a_out_string,
                                const std::string&                                a_opt_string,
                                bool                                              a_im)
      {
         copyStateFromArray( a_Y );
         m_gk_ops->defineMultiPhysicsPC(  a_pc,
                                          a_dof_list,
                                          a_Y,
                                          m_state_comp,
                                          (void*) this,
                                          a_out_string,
                                          a_opt_string,
                                          a_im );
      }

      void updateMultiPhysicsPC(std::vector<Preconditioner<ODEVector,AppCtxt>*>&  a_pc,
                                const ODEVector&                                  a_Y,
                                const Real                                        a_time,
                                const int                                         a_step, 
                                const int                                         a_stage,
                                const Real                                        a_shift,
                                const bool                                        a_im)
      {
        copyStateFromArray( a_Y );
        m_gk_ops->updateMultiPhysicsPC( a_pc,
                                        m_state_comp,
                                        a_time,
                                        a_step,
                                        a_stage,
                                        a_shift,
                                        a_im );
      }

      void solveSolutionPC( ODEVector&        a_Z, 
                            const ODEVector&  a_R, 
                            int               a_i=-1)
      {
        copyStateFromArray( a_R );
        copyRHSFromArray( a_Z );
        m_gk_ops->solveSolutionPC(  m_rhs,
                                    m_state_comp,
                                    a_i );
        copyRHSToArray( a_Z );
      }
   
      void solveFluidOpPCImEx(  ODEVector&        a_Z, 
                                const ODEVector&  a_R, 
                                int               a_i=-1)
      {
        copyStateFromArray( a_R );
        copyRHSFromArray( a_Z );
        m_gk_ops->solveFluidOpPCImEx( m_rhs,
                                      m_state_comp,
                                      a_i );
        copyRHSToArray( a_Z );
      }
   
      void solveVlasovOpPCImEx(Preconditioner<ODEVector,AppCtxt>* const a_pc,
                               ODEVector&                             a_Z, 
                               const ODEVector&                       a_R, 
                               int                                    a_i=-1 )
      {
        copyStateFromArray( a_R );
        copyRHSFromArray( a_Z );
        m_gk_ops->solveVlasovOpPCImEx(  a_pc,
                                        m_rhs,
                                        m_state_comp,
                                        a_i );
        copyRHSToArray( a_Z );
      }

      inline bool isLinear()  { return m_gk_ops->isLinear(); }
      inline void printFunctionCounts() { m_gk_ops->printFunctionCounts(); }

      inline
      void defineODEVector( ODEVector& a_vec) const
      {
        a_vec.define(m_state_comp);
      }

      inline
      const std::string& tiClass() const { return m_ti_class; }

      inline
      const std::string& tiMethod() const { return m_ti_method; }

      inline
      void useScales( const bool a_use_scales = true )
      {
        m_use_scales = a_use_scales;
      }
   
      inline
      void setScaleTolerance(const Real a_tol)
      {
        m_scale_tol = a_tol;
      }
   
      inline
      void computeUScale( const ODEVector&  a_vec )
      {
        static bool first_call = true;
        if (m_use_scales) {
          m_state_comp.copyFrom(a_vec.data());
          m_state_comp.setVarScaling(m_scale_u.data(), m_scale_tol);
        } else {
          if (first_call) m_scale_u.ones();
        }
        first_call = false;
      }
   
      inline
      void computeRHSOpScale( const ODEVector&  a_vec )
      {
        static bool first_call = true;
        if (m_use_scales) {
          m_rhs.copyFrom(a_vec.data());
          m_rhs.setVarScaling(m_scale_rhsop.data(), m_scale_tol);
        } else {
          if (first_call) m_scale_rhsop.ones();
        }
        first_call = false;
      }
   
      inline
      void computeLHSOpScale( const ODEVector&  a_vec )
      {
        static bool first_call = true;
        if (m_use_scales) {
          m_rhs.copyFrom(a_vec.data());
          m_rhs.setVarScaling(m_scale_lhsop.data(), m_scale_tol);
        } else {
          if (first_call) m_scale_lhsop.ones();
        }
        first_call = false;
      }
   
      inline void scaleU(ODEVector& a_vec) 
      { 
        a_vec.elementwiseMultiply(m_scale_u); 
      }
      
      inline void scaleRHSOp(ODEVector& a_vec) 
      { 
        a_vec.elementwiseMultiply(m_scale_rhsop); 
      }
      
      inline void scaleLHSOp(ODEVector& a_vec) 
      { 
        a_vec.elementwiseMultiply(m_scale_lhsop); 
      }
   
      inline void deScaleU(ODEVector& a_vec) 
      { 
        a_vec.elementwiseDivide(m_scale_u); 
      }
   
      inline void deScaleRHSOp(ODEVector& a_vec) 
      { 
        a_vec.elementwiseDivide(m_scale_rhsop); 
      }
      
      inline void deScaleLHSOp(ODEVector& a_vec) 
      { 
        a_vec.elementwiseDivide(m_scale_lhsop); 
      }
   
      inline const ODEVector& uScale() { return m_scale_u; }
      inline const ODEVector& rhsOpScale() { return m_scale_rhsop; }
      inline const ODEVector& lhsOpScale() { return m_scale_lhsop; }
   
      inline Real scaledVecNormU( const ODEVector& a_vec) 
      { 
        return a_vec.computeNorm(m_scale_u); 
      }
      
      inline Real scaledVecNormRHSOp( const ODEVector& a_vec) 
      { 
        return a_vec.computeNorm(m_scale_rhsop); 
      }
      
      inline Real scaledVecNormLHSOp( const ODEVector& a_vec) 
      { 
        return a_vec.computeNorm(m_scale_lhsop); 
      }

      inline bool trivialSolutionOp() const 
      { 
        return m_gk_ops->trivialSolutionOp(); 
      }

   private:

      void setupFieldHistories();

      void getConfigurationSpaceDisjointBoxLayout( CFG::DisjointBoxLayout& );

      void getVelocitySpaceDisjointBoxLayout( const VEL::ProblemDomain&,
                                              VEL::DisjointBoxLayout& );

      void createState();
   
      void createFluidSpecies( CFG::FluidSpeciesPtrVect& fluid_species );

      void createScalars( ScalarPtrVect& scalars );

      void createKineticSpecies( KineticSpeciesPtrVect& kinetic_species );

      inline
      void createGlobalDOF()
      {
        m_global_dof.define( m_state_comp.getVectorSize(),
                             m_state_comp.dataKinetic(),
                             m_state_comp.dataFluid(),
                             m_state_comp.dataScalar());
      }

      VEL::ProblemDomain getVelocityDomain() const;
   
      void createConfigurationSpace();

      void createVelocitySpace();

      void createPhaseSpace( ParmParse& ppgksys );

      void enforcePositivity( KineticSpeciesPtrVect& a_soln );

      void computeError(const double                a_time,
                        const int                   a_step,
                        const LevelData<FArrayBox>& a_solution,
                        double&                     a_l1_error,
                        double&                     a_l2_error,
                        double&                     a_max_error);

      /// Print parameters.
      /**
       * Private method to print parameters.
       */
      void printParameters() const;

      /// Parse parameters.
      /**
       * Private method to obtain control parameters from "gksystem" section
       * of ParmParse database.
       *
       * \param[in] ppgksys gksystem input database.
       */
      void parseParameters( ParmParse& ppgksys );

      void giveSpeciesTheirGyroaverageOps();

      bool m_using_electrons; // true for either boltzmann or kinetic
      bool m_enforce_stage_positivity;
      bool m_enforce_step_positivity;
      bool m_enforce_step_floor;
   
      int               m_max_grid_size;
      std::vector<int>  m_num_cells;
      std::vector<int>  m_num_velocity_cells;
      std::vector<bool> m_is_periodic;
      int               m_kinetic_ghosts;
      int               m_fluid_ghosts;

      std::vector<int> m_configuration_decomposition;
      std::vector<int> m_velocity_decomposition;
      std::vector<int> m_phase_decomposition;

      Vector<ProblemDomain> m_domains;
   
      std::string m_mag_geom_type;

      GKUnits* m_units;

      RefCountedPtr<CFG::MagGeom>     m_mag_geom;
      RefCountedPtr<CFG::MagCoordSys> m_mag_geom_coords;
      RefCountedPtr<VEL::VelCoordSys> m_velocity_coords;
   
      VEL::RealVect            m_dv;

      RefCountedPtr<PhaseGeom> m_phase_geom;
      RefCountedPtr<PhaseCoordSys> m_phase_coords;

      RefCountedPtr<PhaseGrid>  m_phase_grid;
      VelocityPtrVect          m_velocity;

      GKOps* m_gk_ops;

      GKState m_state_comp;
      GKState m_state_phys;
      GKRHSData m_rhs;

      GlobalDOF m_global_dof;

     // Parameters to control plotting
     /**
      * Private variables to control what hdf5 files to generate
      */
      bool m_hdf_potential;
      bool m_hdf_potential_non_zonal;
      bool m_hdf_efield;
      bool m_hdf_ExBdata;
      bool m_hdf_density;
      bool m_hdf_momentum;
      bool m_hdf_ParallelMomentum;
      bool m_hdf_PoloidalMomentum;
      bool m_hdf_ParallelVelocity;

      bool m_hdf_energyDensity;
      bool m_hdf_kineticEnergyDensity;
      bool m_hdf_perpEnergyDensity;
      bool m_hdf_parallelEnergyDensity;
      bool m_hdf_pressure;
      bool m_hdf_perpPressure;
      bool m_hdf_parallelPressure;
      bool m_hdf_gradPoverN;
      bool m_hdf_temperature;
      bool m_hdf_perpTemperature;
      bool m_hdf_parallelTemperature;

      bool m_hdf_parallelHeatFlux;
      bool m_hdf_totalParallelHeatFlux;
      bool m_hdf_fourthMoment;
      bool m_hdf_ParticleFlux;
      bool m_hdf_HeatFlux;
      
      bool m_hdf_vlasov_divergence;
      bool m_hdf_AmpereErIncrement;
      bool m_hdf_total_density;
   
      bool m_hdf_dfn;
      bool m_hdf_deltaF;
      bool m_hdf_dfn_at_mu;
      bool m_hdf_fluids;
   
      bool m_hdf_vpartheta;
      bool m_hdf_bfvpartheta;
      bool m_hdf_rtheta;
      bool m_hdf_vparmu;

      std::vector<int> m_fixed_plotindices;

      PositivityPostProcessor m_positivity_post_processor;
      CFG::FloorPostProcessor m_floor_post_processor;

      bool m_old_vorticity_model;

      int m_verbosity;
      
      std::string m_ti_class;
      std::string m_ti_method;

      /* Scales */
      bool m_use_scales;
      Real m_scale_tol;
      ODEVector m_scale_u, m_scale_rhsop, m_scale_lhsop;
};

#include "NamespaceFooter.H"

#endif
