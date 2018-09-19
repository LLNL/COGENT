#ifndef  _XPSYSTEM_H_
#define  _XPSYSTEM_H_

#include "RK4Integrator.H"
#include "ODEInterface.H"
#include "CFGBlock.H"
#include "MagCoordSys.H"
#include "MultiBlockLevelGeom.H"
#include "MultiBlockLevelExchangeAverage.H"
#include "Vector.H"
#include "CH_HDF5.H"
#include "LoadBalance.H"
#include "LocSpatialIBC.H"
#include "Advect.H"
#include "PositivityPostProcessor.H"

#include "NamespaceHeader.H"


class XPSystem : public ODEInterface<LevelData<FArrayBox> >
{
   public:

      /// Constructor with initialization.
      /**
       *  Pure virtual function that will be resposible for setting
       *  initial conditions on the soln data.
       *
       *  @param[in] pp the input database.
       */
      XPSystem( ParmParse& pp );

      /// Destructor.
      /**
       */
      virtual ~XPSystem();

      void createBlocks( Vector<CFGBlock>& block_vec ) const;

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
      virtual void defineRHSData( LevelData<FArrayBox>&       rhs,
                                  const LevelData<FArrayBox>& prototype );

      /// Defines the Vector<SpeciesModel> data holder for data with ghost cells.
      /**
       * Concrete implementation of pure virtual base class member function.
       * Defines the Vector<SpeciesModel> data holder for data with ghost cells.
       *
       * @param[out] soln data holder for solution with ghost cells.
       * @param[in] prototype prototype data holder.
       */
      virtual void defineSolnData( LevelData<FArrayBox>&       soln,
                                   const LevelData<FArrayBox>& prototype );

      /// Checks for conformity between soln and protoSoln.
      /**
       * Concrete implementation of pure virtual base class member function.
       * Checks for conformity between soln and protoSoln.
       *
       * @param[out] soln data holder.
       * @param[in] protoSoln prototype soln data holder.
       */
      virtual bool validSolnData( const LevelData<FArrayBox>& soln,
                                  const LevelData<FArrayBox>& protoSoln );

      /// Checks for conformity between rhs and protoSoln.
      /**
       * Concrete implementation of pure virtual base class member function.
       * Checks for conformity between soln and protoSoln ignoring ghost cells.
       *
       * @param[out] rhs data holder.
       * @param[in] protoSoln prototype soln data holder.
       */
      virtual bool validRHSData( const LevelData<FArrayBox>& rhs,
                                 const LevelData<FArrayBox>& protoSoln );

      /// Sets the Vector<SpeciesModel> data to zero.
      /**
       * Concrete implementation of pure virtual base class member function.
       * Sets the Vector<SpeciesModel> data to zero.
       *
       * @param[out] soln data holder.
       */
      virtual void zeroSolnData( LevelData<FArrayBox>& soln );

      /// Adds (scale * increment) to the Vector<SpeciesModel> data holder.
      /**
       * Concrete implementation of pure virtual base class member function.
       * Adds (scale * increment) to the Vector<SpeciesModel> data holder.
       *
       * @param[out] soln solution data to which is added.
       * @param[in] increment amount to add.
       * @param[in] scale multiplier by which to scale increment.
       */
      virtual void addSolnData( LevelData<FArrayBox>&       soln,
                                const LevelData<FArrayBox>& increment,
                                const Real                   scale );

      /// Copies all data from srcSoln to dstSoln.
      /**
       * Concrete implementation of pure virtual base class member function.
       * Copies all data from srcSoln to dstSoln.
       *
       * @param[out] dstSoln destination data holder.
       * @param[in] srcSoln source data holder.
       */
      virtual void copySolnData( LevelData<FArrayBox>&       dstSoln,
                                 const LevelData<FArrayBox>& srcSoln );

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
      virtual void evalRHS( LevelData<FArrayBox>&       rhs,
                            const LevelData<FArrayBox>& soln,
                            const int             step_number,
                            const Real            time,
                            const int             stage);

      void getConditionNumbers(const MultiBlockLevelExchangeAverage* a_mblexPtr,
                               const MultiBlockLevelGeom&            a_geom ) const;

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
      Real advance( const Real cur_time, const Real dt, const int step_number );

      /// Callback to allow for additional operations after stage advance.
      /**
       * Callback to allow for additional operations after stage advance.
       *
       * @param[in,out] soln data holder.
       */
      virtual void postStageAdvance( LevelData<FArrayBox>& soln );

      /// Write plot file.
      /**
       * Write plotting data to HDF5 files with a common prefix.
       *
       */
      void writePlotFile(const char *prefix, const int cur_step );

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
      void readCheckpointFile( HDF5Handle& handle, int* cur_step, double* cur_time, double* cur_dt );

      void writeFieldHistory(int cur_step, double cur_time, bool startup_flag);

      void printDiagnostics();

   private:

      void setInitialConditions( LevelData<FArrayBox>& soln );

      void setBoundaryData( LevelData<FArrayBox>& soln );

      void multJ( const LevelData<FArrayBox>& soln_physical,
                  LevelData<FArrayBox>&       soln_mapped );

      void divideJ( const LevelData<FArrayBox>& soln_mapped,
                    LevelData<FArrayBox>&       soln_physical );

      void createConfigurationSpace(const Vector<CFGBlock>& block_vec,
                                    Vector<ProblemDomain>& cs_domains);

      void enforcePositivity( LevelData<FArrayBox>& a_soln );

      int m_verbosity;

      MultiBlockLevelGeom* m_mag_geom;
      MagCoordSys*         m_mag_geom_coords;
      ProblemDomain        m_domain;

      LevelData<FArrayBox> m_solution;
      LevelData<FArrayBox> m_solution_mapped;
      LevelData<FArrayBox> m_solution_mapped_old;
      LevelData<FArrayBox> m_initial_condition;

      RK4Integrator< LevelData<FArrayBox> >* m_integrator;

      Advect* m_operator;

      LocSpatialIBC* m_ibc;

      // Parameters to control plotting
      /**
       * Private variables to control what hdf5 files to generate
       */
      bool m_hdf_density;
      bool m_hdf_vpartheta;
      bool m_hdf_frtheta;
      bool m_hdf_pressure;
      bool m_hdf_vparmu;
      std:: vector<int> m_fixed_plotindices;

      /// Print parameters.
      /**
       * Private method to print parameters.
       */
      void printParameters(const Vector<CFGBlock>& block_vec);

      /// Parse parameters.
      /**
       * Private method to obtain control parameters from "gksystem" section
       * of ParmParse database, and also generic (not-species or potential-related)
       * parts of "IBC" section of ParmParse database
       *
       * \param[in] ppgksys gksystem input database.
       *\param[in] ppibc IBC input database.
       */
      void parseParameters( ParmParse& ppgksys , Vector<CFGBlock>& block_vec );

      bool m_enforce_stage_positivity;
      bool m_enforce_step_positivity;

      IntVect m_ghostVect;

      PositivityPostProcessor m_positivity_post_processor;

      ParmParse m_ppgksys;
};

#include "NamespaceFooter.H"

#endif
