#ifndef  _GKSYSTEM_H_
#define  _GKSYSTEM_H_

#include "AppCtxt.H"
#include "GKOps.H"
#include "evalOpMatrix.H"
#include "PositivityPostProcessor.H"

#include "GKSystemParameters.H"
#include "GKDiagnostics.H"

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
      GKSystem( int sys_id = -1 );

      GKSystem( ParmParse&, int sys_id = -1 );

      GKSystem( const GKSystemParameters&, int sys_id = -1 );

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

      inline int sysID() const { return m_sys_id; }
      inline void sysID( const int a_id ) { m_sys_id = a_id; }

      /// Write plot file.
      /**
       * Write plotting data to HDF5 files with a common prefix.
       *
       */
      void writePlotFile( const char *prefix, 
                          const int cur_step, 
                          const double& cur_time );

      /// Augment a history file with a value at a spatial index
      /**
       *Augment a history file with a value at a spatial index
       */

      void writeFieldHistory( int cur_step, 
                              double cur_time, 
                              bool startup_flag);

      /// Write checkpoint file.
      /**
       * Write checkpoint data to an output HDF5 file.
       *
       * @param[in] handle a reference to an HDF5 plot file.
       */
      void writeCheckpointFile( const std::string&, 
                                const int cur_step,
                                const double cur_time, 
                                const double cur_dt );

      /// Read checkpoint file.
      /**
       * Read checkpoint data from an output HDF5 file and reinitialize.
       *
       * @param[in] handle a reference to an HDF5 plot file.
       */
      void readCheckpointFile( const std::string&, 
                               int& cur_step,
                               double& cur_time, 
                               double& cur_dt );

      void getDiagnostics(std::vector<Real>&, std::vector<Real>& );
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
      void defineMultiPhysicsPC(std::vector<
                                  Preconditioner< ODEVector,
                                                  AppCtxt >* >&  a_pc,
                                std::vector<DOFList>& a_dof_list,
                                const ODEVector&      a_Y,
                                const std::string&    a_out_string,
                                const std::string&    a_opt_string,
                                bool                  a_im,
                                const int             a_sys_id = -1)
      {
         copyStateFromArray( a_Y );
         m_gk_ops->defineMultiPhysicsPC(  a_pc,
                                          a_dof_list,
                                          a_Y,
                                          m_state_comp,
                                          (void*) this,
                                          a_out_string,
                                          a_opt_string,
                                          a_im,
                                          a_sys_id );
      }

      inline
      void defineMultiPhysicsPC(std::vector<
                                  Preconditioner< ODEVector,
                                                  AppCtxt >* >&  a_pc,
                                std::vector<DOFList>& a_dof_list,
                                const ODEVector&      a_Y,
                                void*                 a_ops,
                                const std::string&    a_out_string,
                                const std::string&    a_opt_string,
                                bool                  a_im,
                                const int             a_sys_id = -1)
      {
         copyStateFromArray( a_Y );
         m_gk_ops->defineMultiPhysicsPC(  a_pc,
                                          a_dof_list,
                                          a_Y,
                                          m_state_comp,
                                          a_ops,
                                          a_out_string,
                                          a_opt_string,
                                          a_im,
                                          a_sys_id );
      }

      void updateMultiPhysicsPC(std::vector<
                                  Preconditioner< ODEVector,
                                                  AppCtxt >* >&  a_pc,
                                const ODEVector&  a_Y,
                                const Real        a_time,
                                const int         a_step, 
                                const int         a_stage,
                                const Real        a_shift,
                                const bool        a_im)
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
                            const int         a_sys_id,
                            int               a_i=-1)
      {
        CH_assert(a_sys_id == m_sys_id);
        copyStateFromArray( a_R );
        copyRHSFromArray( a_Z );
        m_gk_ops->solveSolutionPC(  m_rhs,
                                    m_state_comp,
                                    a_i );
        copyRHSToArray( a_Z );
      }
   
      void solveFluidOpPCImEx(  ODEVector&          a_Z, 
                                const ODEVector&    a_R, 
                                const std::string&  a_name,
                                const int           a_sys_id,
                                int                 a_i=-1)
      {
        CH_assert(a_sys_id == m_sys_id);
        copyStateFromArray( a_R );
        copyRHSFromArray( a_Z );
        m_gk_ops->solveFluidOpPCImEx( m_rhs,
                                      m_state_comp,
                                      a_name,
                                      a_i );
        copyRHSToArray( a_Z );
      }
   
      void solveVlasovOpPCImEx(Preconditioner<ODEVector,AppCtxt>* const a_pc,
                               ODEVector&                             a_Z, 
                               const ODEVector&                       a_R, 
                               const int                              a_sid,
                               int                                    a_i=-1 )
      {
        CH_assert(a_sid == m_sys_id);
        copyStateFromArray( a_R );
        copyRHSFromArray( a_Z );
        m_gk_ops->solveVlasovOpPCImEx(  a_pc,
                                        m_rhs,
                                        m_state_comp,
                                        a_i );
        copyRHSToArray( a_Z );
      }

      inline bool isLinear()  { return m_gk_ops->isLinear(); }

      inline void printFunctionCounts() 
      { 
        if (!procID() && m_sys_id >= 0) {
          printf("  System %d", m_sys_id);
        }
        m_gk_ops->printFunctionCounts(); 
      }

      inline
      void defineODEVector( ODEVector& a_vec)
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
        if (m_use_scales) {
          m_state_comp.copyFrom(a_vec.data());
          m_state_comp.setVarScaling(m_scale_u.data(), m_scale_tol);
        } else {
          m_scale_u.ones();
        }
      }
   
      inline
      void computeRHSOpScale( const ODEVector&  a_vec )
      {
        if (m_use_scales) {
          m_rhs.copyFrom(a_vec.data());
          m_rhs.setVarScaling(m_scale_rhsop.data(), m_scale_tol);
        } else {
          m_scale_rhsop.ones();
        }
      }
   
      inline
      void computeLHSOpScale( const ODEVector&  a_vec )
      {
        if (m_use_scales) {
          m_rhs.copyFrom(a_vec.data());
          m_rhs.setVarScaling(m_scale_lhsop.data(), m_scale_tol);
        } else {
          m_scale_lhsop.ones();
        }
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

      inline const GKState& stateComp() const { return m_state_comp; }
      inline GKState& stateComp() { return m_state_comp; }
      inline const GKState& statePhys() const { return m_state_phys; }
      inline GKState& statePhys() { return m_state_phys; }

      inline void copyDOFs( ODEVector&        a_Y,
                            const ODEVector&  a_X,
                            const DOFList&    a_dof_list,
                            const int         a_sys_id ) const
      {
        CH_assert(a_sys_id == m_sys_id);
        a_Y.copy(a_X, a_dof_list);
      }

      inline void getVec( const ODEVector&  a_X,
                          const int,
                          ODEVector&        a_Y ) const
      {
        CH_assert(!a_Y.isDefined());
        a_Y.define( a_X.getVectorSize(),
                    a_X.data(),
                    a_X.getState() );
      }

      virtual void defineShell( const GKSystemParameters& );
      virtual void initializeShell( const int, const double );

      inline 
      void getTimeScales( std::vector<std::pair<std::string,Real> >& a_T) const
      {
        a_T.clear();
        for (int i = 0; i < m_dt_vlasov.size(); i++) {
          std::string name_str = "Vlasov (" + m_dt_vlasov[i].first + ")"; 
          a_T.push_back( std::pair<std::string,Real>( name_str, m_dt_vlasov[i].second ) );
        }
        a_T.push_back( std::pair<std::string,Real>( "Collisions", 
                                                    m_dt_collisions) );
        a_T.push_back( std::pair<std::string,Real>( "Transport", 
                                                    m_dt_transport) );
        a_T.push_back( std::pair<std::string,Real>( "Neutrals", 
                                                    m_dt_neutrals) );
      }

      inline
      const GKOps* getOps() const { return m_gk_ops; }

      inline
      const GKDiagnostics* getDiagnostics() const { return m_diagnostics; }

      inline
      void ARKFirstStageExplicit( const bool a_flag ) { m_gk_ops->ARKFirstStageExplicit(a_flag); }

   private:

      void define( const GKSystemParameters& );


      void setupFieldHistories();

      void getConfigurationSpaceDisjointBoxLayout( CFG::DisjointBoxLayout& );

      void getVelocitySpaceDisjointBoxLayout( const VEL::ProblemDomain&,
                                              VEL::DisjointBoxLayout& );

      void createState( GKState&,
                        const GKSystemParameters&,
                        const int,
                        const int );
   
      void createFluidSpecies( CFG::FluidSpeciesPtrVect& fluid_species,
                               const GKSystemParameters& );

      void createScalars( ScalarPtrVect& scalars,
                          const GKSystemParameters& );

      void createKineticSpecies(  KineticSpeciesPtrVect& kinetic_species,
                                  const GKSystemParameters& );

      inline
      void createGlobalDOF()
      {
        m_global_dof.define( m_state_comp.getVectorSize(),
                             m_state_comp.dataKinetic(),
                             m_state_comp.dataFluid(),
                             m_state_comp.dataScalar());
      }

      VEL::ProblemDomain getVelocityDomain() const;
   
      void createConfigurationSpace( const GKSystemParameters& );

      void createVelocitySpace();

      void createPhaseSpace( const GKSystemParameters& );

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

      void setParameters( const GKSystemParameters& );

      void giveSpeciesTheirGyroaverageOps();

      inline void prependSysID(std::string& a_str)
      {
        if (m_sys_id >= 0) {
          a_str = m_sys_id_str + "." + a_str;
        }
      }

      inline void appendSysID(std::string& a_str)
      {
        if (m_sys_id >= 0) {
          a_str = a_str + "." + m_sys_id_str;
        }
      }

      int m_sys_id;
      std::string m_sys_id_str;

      bool m_using_electrons; // true for either boltzmann or kinetic
      bool m_enforce_stage_positivity;
      bool m_enforce_step_positivity;
      bool m_enforce_positivity_pointwise;
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
      GKDiagnostics* m_diagnostics;

      GKState m_state_comp;
      GKState m_state_phys;
      GKRHSData m_rhs;

      GlobalDOF m_global_dof;

      std::vector<std::string> m_diagnostics_dfn_plots;
      std::vector<std::string> m_diagnostics_cfg_field_vars;
      std::vector<std::string> m_diagnostics_cfg_vars_kin_spec;
      std::vector<std::string> m_diagnostics_cfg_vars_total_kin_spec;
      bool m_diagnostics_fluids;

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

      /* Time scales */
      std::vector<std::pair<std::string,Real> > m_dt_vlasov;
      Real m_dt_collisions, m_dt_transport, m_dt_neutrals;
};

#include "NamespaceFooter.H"

#endif
