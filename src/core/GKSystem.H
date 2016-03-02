#ifndef  _GKSYSTEM_H_
#define  _GKSYSTEM_H_

#undef BLOB_ERROR_TEST

#if CFG_DIM==2
#define CFG_NAMESPACE D2
#define PDIM 4
#endif
#if CFG_DIM==3
#define CFG_NAMESPACE D3
#define PDIM 5
#endif
#define VEL_DIM 2
#define VEL_NAMESPACE D2


#include "KineticSpecies.H"
#include "GKVlasov.H"
#include "PhaseGrid.H"
#include "RefCountedPtr.H"

#include "GKCollisions.H"
#include "GKNeutrals.H"
#include "GKTransport.H"

#include "GKUnit.H"
#include "PhaseCoordSys.H"
#include "GKSystemIC.H"
#include "GKSystemBC.H"
#include "GKOps.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM VEL_DIM
#include "VelCoordSys.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
//#include "BoltzmannElectron.H"
#include "MagCoordSys.H"
#include "MagGeom.H"
#include "MultiBlockLevelGeom.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "ParmParse.H"
#include "Vector.H"
#include "CH_HDF5.H"
#include "RefCountedPtr.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"
#include "PositivityPostProcessor.H"

#include "GKTimeIntegration.H"
#include "TimeIntegrator.H"
#include "TiRK.H"
#include "TiARK.H"

#include "NamespaceHeader.H"
namespace CFG = CFG_NAMESPACE;
namespace VEL = VEL_NAMESPACE;

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
 * gksystem.magnetic_geometry_mapping = "Miller"
 * gksystem.hdf_density = true
 * gkysytem.hdf_vpartheta = false
 * gkysytem.hdf_frtheta = true
 * gksystem.fixed_plot_indices = 3 1 0 2 2
 *
 */

class GKSystem 
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
      GKSystem( ParmParse& pp, bool a_useExternalTI = false);

      /// Destructor.
      /**
       */
      virtual ~GKSystem();

      virtual void initialize(const int cur_step);

      /// Define Vector<SpeciesModel> data holder for data.
      /**
       * Concrete implementation of pure virtual base class member function.
       * Defines the Vector<SpeciesModel> data holder for data without ghost
       * cells.
       *
       * @param[out] rhs data holder for rhs (without ghost cells).
       * @param[in] prototype prototype data holder.
       */
      virtual void defineRHSData( KineticSpeciesPtrVect&       rhs,
                                  const KineticSpeciesPtrVect& prototype );

      /// Defines the Vector<SpeciesModel> data holder for data with ghost cells.
      /**
       * Concrete implementation of pure virtual base class member function.
       * Defines the Vector<SpeciesModel> data holder for data with ghost cells.
       *
       * @param[out] soln data holder for solution with ghost cells.
       * @param[in] prototype prototype data holder.
       */
      virtual void defineSolnData( KineticSpeciesPtrVect&       soln,
                                   const KineticSpeciesPtrVect& prototype );

      /// Checks for conformity between soln and protoSoln.
      /**
       * Concrete implementation of pure virtual base class member function.
       * Checks for conformity between soln and protoSoln.
       *
       * @param[out] soln data holder.
       * @param[in] protoSoln prototype soln data holder.
       */
      virtual bool validSolnData( const KineticSpeciesPtrVect& soln,
                                  const KineticSpeciesPtrVect& protoSoln );

      /// Checks for conformity between rhs and protoSoln.
      /**
       * Concrete implementation of pure virtual base class member function.
       * Checks for conformity between soln and protoSoln ignoring ghost cells.
       *
       * @param[out] rhs data holder.
       * @param[in] protoSoln prototype soln data holder.
       */
      virtual bool validRHSData( const KineticSpeciesPtrVect& rhs,
                                 const KineticSpeciesPtrVect& protoSoln );

      /// Sets the Vector<SpeciesModel> data to zero.
      /**
       * Concrete implementation of pure virtual base class member function.
       * Sets the Vector<SpeciesModel> data to zero.
       *
       * @param[out] soln data holder.
       */
      virtual void zeroSolnData( KineticSpeciesPtrVect& soln );

      /// Adds (scale * increment) to the Vector<SpeciesModel> data holder.
      /**
       * Concrete implementation of pure virtual base class member function.
       * Adds (scale * increment) to the Vector<SpeciesModel> data holder.
       *
       * @param[out] soln solution data to which is added.
       * @param[in] increment amount to add.
       * @param[in] scale multiplier by which to scale increment.
       */
      virtual void addSolnData( KineticSpeciesPtrVect&       soln,
                                const KineticSpeciesPtrVect& increment,
                                const Real                   scale );

      /// Copies all data from srcSoln to dstSoln.
      /**
       * Concrete implementation of pure virtual base class member function.
       * Copies all data from srcSoln to dstSoln.
       *
       * @param[out] dstSoln destination data holder.
       * @param[in] srcSoln source data holder.
       */
      virtual void copySolnData( KineticSpeciesPtrVect&       dstSoln,
                                 const KineticSpeciesPtrVect& srcSoln );

      /// Compute a stable time step.
      /**
       * Computes and returns a stable time step estimate.
       */
      Real stableDt(const int step_number);

      /// Advances the solution.
      /**
       * Advance the solution by a single time step.  Returns the current
       * simulation time.
       */
      void advance( Real* cur_time, Real* dt, int* step_number );

      /// Callback to allow for additional operations after stage advance.
      /**
       * Callback to allow for additional operations after stage advance.
       *
       * @param[in,out] soln data holder.
       */
      virtual void postStageAdvance( KineticSpeciesPtrVect& soln );

      /// Write plot file.
      /**
       * Write plotting data to HDF5 files with a common prefix.
       *
       */
      void writePlotFile(const char *prefix, const int cur_step, const double& cur_time );

      // Read history parameters from input file (see writeFieldHistory)

      void setupFieldHistories();

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

      void computeEField(const int step_number, KineticSpeciesPtrVect& soln,
                         CFG::LevelData<CFG::FluxBox>& E_field,
                         CFG::LevelData<CFG::FArrayBox>& E_field_cell);

#ifdef BLOB_ERROR_TEST
      void computeError(const double                 time,
                        const int                    step,
                        const KineticSpeciesPtrVect& species,
                        double&                      l1_error,
                        double&                      l2_error,
                        double&                      max_error);
#endif

      void getConfigurationSpaceDisjointBoxLayout( CFG::DisjointBoxLayout& grids );

      void getVelocitySpaceDisjointBoxLayout( const VEL::ProblemDomain& domain,
                                              VEL::DisjointBoxLayout& grids );

      void printDiagnostics();

      void preTimeStep  (int,Real);
      void postTimeStep (int,Real,Real);
      inline void postTimeStage(const int a_step,const Real a_time)
        { m_gk_ops->postTimeStage(a_step,a_time,m_state); }

      int  getVectorSize()  { return m_state.getVectorSize(); }
      int  getNComponents() { return m_state.getNComponents();}

      inline void constructPCMatIMEX(void* a_Pmat,Real a_time)
        { m_gk_ops->constructPCMatIMEX(a_Pmat,a_time,m_state); }

      inline void copyStateToArray   (Real *Y) { m_state.copyTo(Y);   }
      inline void copyStateFromArray (Real *Y) { m_state.copyFrom(Y); }
      inline void copyRHSToArray     (Real *Y) { m_rhs.copyTo(Y);   }
      inline void copyRHSFromArray   (Real *Y) { m_rhs.copyFrom(Y); }
      inline void addStateFromArray  (Real *Y, Real a_scale = 1.0) { m_state.addFrom(Y,a_scale); }

      inline void computeRHSFunctionExp (Real t)            { m_gk_ops->explicitOp     (m_rhs,t,m_state); }
      inline void computeRHSFunctionImEx(Real t)            { m_gk_ops->explicitOpImEx (m_rhs,t,m_state); }
      inline void computeIFunctionImEx  (Real t,int flag=1) { m_gk_ops->implicitOpImEx (m_rhs,t,m_state,1,flag); }

      inline bool isLinear()  { return m_gk_ops->isLinear(); }

   private:

      void createSpecies();

      VEL::ProblemDomain getVelocityDomain() const;

      void createConfigurationSpace();

      void createVelocitySpace();

      void createPhaseSpace( ParmParse& ppgksys );

      void enforcePositivity( KineticSpeciesPtrVect& a_soln );

      double sumDfn( const LevelData<FArrayBox>& dfn );

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

      bool m_using_electrons; // true for either boltzmann or kinetic
      bool m_enforce_stage_positivity;
      bool m_enforce_step_positivity;

      int               m_max_grid_size;
      std::vector<int>  m_num_cells;
      std::vector<int>  m_num_velocity_cells;
      std::vector<bool> m_is_periodic;
      IntVect m_ghostVect;

      std::vector<int> m_configuration_decomposition;
      std::vector<int> m_velocity_decomposition;
      std::vector<int> m_phase_decomposition;

      std::string m_mag_geom_type;

      std::string m_ti_class;
      std::string m_ti_method;
      short int   m_ti_type;

      GKUnits*                  m_units;
      CFG::MagGeom*             m_mag_geom;
      CFG::MagCoordSys*         m_mag_geom_coords;
      VEL::VelCoordSys*         m_velocity_coords;

      RefCountedPtr<PhaseGeom>  m_phase_geom;
      PhaseCoordSys* m_phase_coords;

      PhaseGrid*                     m_phase_grid;
      KineticSpeciesPtrVect          m_species;
      KineticSpeciesPtrVect          m_kinetic_species;
      VelocityPtrVect                m_velocity;

      GKOps* m_gk_ops;
      GKSystemIC* m_initial_conditions;
      GKSystemBC* m_boundary_conditions;

      TimeIntegrator<GKState, GKRHSData, GKOps> *m_integrator;
      GKState m_state;
      GKRHSData m_rhs;

     // Parameters to control plotting
     /**
      * Private variables to control what hdf5 files to generate
      */
      bool m_hdf_potential;
      bool m_hdf_efield;
      bool m_hdf_density;
      bool m_hdf_vpartheta;
      bool m_hdf_rtheta;
      bool m_hdf_pressure;
      bool m_hdf_temperature;
      bool m_hdf_fourthMoment;
      bool m_hdf_ParallelMomentum;
      bool m_hdf_PoloidalMomentum;
      bool m_hdf_ParticleFlux;
      bool m_hdf_HeatFlux;
      bool m_hdf_vparmu;
      bool m_hdf_vlasov_divergence;
      bool m_hdf_AmpereErIncrement;
      bool m_hdf_total_density;
      bool m_hdf_dfn;
      bool m_hdf_deltaF;
      bool m_hdf_dfn_at_mu;
      std::vector<int> m_fixed_plotindices;

      PositivityPostProcessor m_positivity_post_processor;

      int m_verbosity;
      bool m_useNativeTimeIntegrator;
};

#include "NamespaceFooter.H"

#endif
