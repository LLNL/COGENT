#ifndef  _SNSYSTEM_H_
#define  _SNSYSTEM_H_

#include "MagCoordSys.H"
#include "MultiBlockLevelGeom.H"
#include "MultiBlockLevelExchangeAverage.H"
#include "Vector.H"
#include "CH_HDF5.H"
#include "LoadBalance.H"
#include "PositivityPostProcessor.H"
#include "FluidSpecies.H"
#include "SNOps.H"
#include "SNSystemIC.H"
#include "SNSystemBC.H"
#include "SNVector.H"

#include "SNTimeIntegration.H"
#include "TimeIntegrator.H"
#include "TiRK.H"

#include "NamespaceHeader.H"


class SNSystem 
{
   public:

      /// Constructor with initialization.
      /**
       *  Pure virtual function that will be resposible for setting
       *  initial conditions on the soln data.
       *
       *  @param[in] pp the input database.
       */
      SNSystem( ParmParse& pp, bool a_useExternalTI = false );

      /// Destructor.
      /**
       */
      virtual ~SNSystem();

      virtual void initialize( const int cur_step );

      /// Compute a stable time step.
      /**
       * Computes and returns a stable time step estimate.
       */
      Real stableDt( const int step_number );

      /// Advances the solution.
      /**
       * Advance the solution by a single time step.  Returns the current
       * simulation time.
       */
      void advance( Real& cur_time, Real& dt, int& step_number );

      /// Callback to allow for additional operations after stage advance.
      /**
       * Callback to allow for additional operations after stage advance.
       *
       * @param[in,out] soln data holder.
       */
      virtual void postStageAdvance( FluidSpeciesPtrVect& soln );

      /// Write plot file.
      /**
       * Write plotting data to HDF5 files with a common prefix.
       *
       */
      void writePlotFile( const char *prefix, const int cur_step, const double& cur_time );

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
      void readCheckpointFile( HDF5Handle& handle, int cur_step, double cur_time, double cur_dt );

      void printDiagnostics();

      void writeFieldHistory(int&, Real&, bool&) {};
   
      void writeHistFile(bool) {};

      void preTimeStep( int, Real );
   
      void postTimeStep( int, Real, Real );

      inline
      void postTimeStage( const int a_step, const Real a_time, const int a_stage )
        {
           m_sn_ops->postTimeStage( a_step, a_time, m_state_comp, a_stage );
        }

      inline
      int getVectorSize() { return 0; }

      inline
      void getFunctionCounts( int& a, int& b, int& c, int& d, int& e, int& f )
      {}
   
      void printFunctionCounts() {};

      void printTimeIntegratorCounts() {};

   private:

      void createDisjointBoxLayout( DisjointBoxLayout& grids );

      void createState();
   
      void setInitialConditions( LevelData<FArrayBox>& soln );

      void setBoundaryData( LevelData<FArrayBox>& soln );

      void multJ( const LevelData<FArrayBox>& soln_physical,
                  LevelData<FArrayBox>&       soln_mapped );

      void divideJ( const LevelData<FArrayBox>& soln_mapped,
                    LevelData<FArrayBox>&       soln_physical,
                    bool                        restrict_to_valid = false);

      void createConfigurationSpace();

      void enforcePositivity( FluidSpeciesPtrVect& a_soln );

      int m_verbosity;

      RefCountedPtr<MagGeom>       m_mag_geom;
      MagCoordSys*         m_mag_geom_coords;

      SNOps*      m_sn_ops;
      SNSystemIC* m_initial_conditions;
      SNSystemBC* m_boundary_conditions;

      TimeIntegrator<SNVector, SNOps> *m_integrator;
      SNState m_state_comp;
      SNState m_state_phys;
      SNRHSData m_rhs;
      SNVector m_serialized_vector;

      /// Print parameters.
      /**
       * Private method to print parameters.
       */
      void printParameters();

      bool m_enforce_stage_positivity;
      bool m_enforce_step_positivity;
      bool m_use_native_time_integrator;
      std::string m_ti_class;
      std::string m_ti_method;
   
      /// Parse parameters.
      /**
       * Private method to obtain control parameters from "gksystem" section
       * of ParmParse database, and also generic (not-species or potential-related)
       * parts of "IBC" section of ParmParse database
       *
       * \param[in] ppgksys gksystem input database.
       *\param[in] ppibc IBC input database.
       */
      void parseParameters( ParmParse& ppgksys );

      IntVect m_ghostVect;

      PositivityPostProcessor m_positivity_post_processor;

      ParmParse m_ppsnsys;
};

#include "NamespaceFooter.H"
#endif
