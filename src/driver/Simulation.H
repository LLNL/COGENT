#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <string>
#include <vector>

#include "RefCountedPtr.H"
#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "REAL.H"
#ifdef CH_USE_TIMER
#include "CH_Timer.H"
#endif // CH_USE_TIMER


#ifdef USE_ARRAYVIEW
#include "UsingNamespace.H"
extern "C" {
#include <fpu_control.h>
}
/* IM: Invalid operation mask
 * DM: Denormalized operand mask
 * ZM: Zero-divide mask
 * OM: Overflow mask
 * UM: Underflow mask
 * PM: Precision (inexact result) mask
  ---(pm is kinda stupid)
*/
static void __attribute__ ((constructor)) trapfpe(void)
{
   fpu_control_t cw =
      _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM | _FPU_MASK_UM);
   _FPU_SETCW(cw);
}
#endif // USE_ARRAYVIEW

#include "NamespaceHeader.H"

/**
 * Simulation driver class.
 *
 * Responsible for generic initialization, advancement, and finalization of
 * simulation.  This class specifically controls the monitoring of the number
 * of steps and of the simulation time in order to execute plotfile and
 * checkpoint file creation as well as simulation termination.
 *
 * Required input keys and data types:
 *
 *    -\b max_step
 *      positiveinteger value maximum number of simulation steps
 *
 *    -\b max_time
 *      positive real value final time at which simulation terminates
 *
 *    -\b fixed_dt
 *      positive real value, if set, for the time step to be used
 *
 *    -\b max_dt_grow
 *      positive real value of the multiple of the current time step to use for
 *      the next time step.  Must be grater than 1.  Multually exclusive with
 *      fixed_dt
 *
 *    -\b initial_dt_fraction
 *      positive real value fraction of the initial stable time step to use.
 *      Must be less than 1.  Multually exclusive with fixed_dt
 *
 *    -\b checkpoint_interval
 *      integer value specifying the number of steps between checkpoint dumps
 *
 *    -\b checkpoint_prefix
 *      string used as prefix for checkpoint file names ["chk"]
 *
 *    -\b plot_interval
 *      integer value specifying the number of steps between plot dumps
 *
 *    -\b plot_prefix
 *      string used as prefix for plot file names ["plt"]
 *
 *    -\b verbosity
 *      integer flag specifying the verbosity of logging output; zero turns
 *      of output, and increasing values produce more detailed logs
 *
 * The following represents a sample input entry:
 *
 * \verbatim
 *  simulation.max_step = 1000
 *  simulation.max_time = 1.0
 *  simulation.max_dt_grow = 1.1
 *  simultaion.initial_dt_fraction = 0.1
 *  simulation.checkpoint_interval = 50
 *  simulation.checkpoint_prefix = "chk"
 *  simulation.plot_interval = 10
 *  simulation.plot_prefix = "plt"
 *  simulation.verbosity = 1
 * \endverbatim
*/
template <class SYSTEM>
class Simulation
{
   public:

      /// Constructor with definition.
      /**
       * \param[in] parm_parse input database from which to intialize.
       */
      Simulation( ParmParse& pp );

      /// Destructor.
      /**
       */
      ~Simulation();

      /// Returns true if simulation is incomplete, false otherwise.
      /**
       */
      bool notDone();

      /// Advances the simulation a single timestep.
      /**
       * Advances the simulation a single timestep and writes
       * restart/diagnostics as necessary.
       */
      void advance();

      /// Finalizes the simulation.
      /**
       * Cleans up and closes open resources at completion of simulation.
       */
      void finalize();

   private:

      /// Initialize timers.
      /**
       * Private method to initialize timers.
       */
      void initializeTimers();

      /// Parse parameters.
      /**
       * Private method to obtain control parameters from "simulation" section
       * of ParmParse database.
       *
       * \param[in] ppsim simulation input database.
       */
      void parseParameters( ParmParse& ppsim );

      /// Print parameters.
      /**
       * Private method to print parameters.
       */
      void printParameters();

      /// Load from restart file.
      /**
       * Private method to initialize from restart file specified in input
       * database.
       *
       * \param[in] ppsim simulation input database.
       */
      void loadRestartFile( ParmParse& ppsim );

      /// Write data to plot file.
      /**
       * Private method to write plot file.
       */
      void writePlotFile();

      /// Write data to history file.
      /**
       * Private method to write history file.
       */
      void writeHistFile(bool startup_flag);

      /// Write data to checkpoint file.
      /**
       * Private method to write checkpoint file.
       */
      void writeCheckpointFile();

      /// Set fixed time step.
      /**
       * Private method to set fixed time step.
       *
       * \param[in] dt_stable maximum stable time step.
       */
      void setFixedTimeStep( const Real& dt_stable );

      /// Select time step.
      /**
       * Private method to set time step based on multiple controls.
       *
       * \param[in] dt_stable maximum stable time step.
       */
      void selectTimeStep();

      int m_verbosity;

      int m_cur_step;
      int m_max_step;

      Real m_cur_time;
      Real m_max_time;

      Real m_cur_dt;
      Real m_fixed_dt;
      Real m_max_dt_grow;
      Real m_init_dt_frac;
      Real m_cfl;
      bool m_adapt_dt;
      static const Real s_DT_EPS;

      int         m_checkpoint_interval;
      int         m_last_checkpoint;
      std::string m_checkpoint_prefix;

      int         m_plot_interval;
      int         m_last_plot;
      std::string m_plot_prefix;


      SYSTEM* m_system;

#ifdef CH_USE_TIMER
      Timer* m_all_timer;
      Timer* m_setup_timer;
      Timer* m_solve_timer;
      Timer* m_shutdown_timer;
#endif
};

#include "NamespaceFooter.H"

#endif // _SIMULATION_H_
