#ifndef _COGENTTimeIntegration_H_
#define _COGENTTimeIntegration_H_

#include <string>
#include <vector>
#include <time.h>
#include <limits>

#include "RefCountedPtr.H"
#include "Simulation.H"

#ifdef CH_USE_TIMER
#include "CH_Timer.H"
#endif // CH_USE_TIMER

#include "evalOpMatrix.H"
#include "TimeIntegratorLibrary.H"

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
 * COGENTTimeIntegration driver class.
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
 *    -\b plot_time_interval
 *      Real value specifying the time interval between plot dumps
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
template <class Vector, class System>
class COGENTTimeIntegration : public Simulation<Vector,System>
{
   public:

      /// Constructor with definition.
      /**
       * \param[in] parm_parse input database from which to intialize.
       */
      COGENTTimeIntegration<Vector,System>();

      /// Destructor.
      /**
       */
      ~COGENTTimeIntegration<Vector,System>();

      void initialize( ParmParse&, System* const );

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

      /// Runs the simulation by advancing until done
      inline  void solve()
      {
        while ( notDone() ) advance();
      }

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

      /// Set fixed time step
      /**
       * Private method to set the fixed time step
       *
       * \param[in] a_dt_stable stable time step size
       */
      void setFixedTimeStep(const Real&);

      /// Pre time step operations.
      /**
       * Private method called before a time step is started.
       */
      inline void preTimeStep();

      /// Post time step operations.
      /**
       * Private method called after time step has been completed.
       */
      inline void postTimeStep();

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
      Real        m_plot_time_interval;
      Real        m_plot_time;
      int         m_last_plot;
      std::string m_plot_prefix;
      bool        m_write_odevec;

      System* m_system;
      TimeIntegrator<Vector, System>* m_integrator;

      Vector m_U, m_U_old;
      Real m_step_norm_rel, m_step_norm_abs;

#ifdef CH_USE_TIMER
      Timer* m_all_timer;
      Timer* m_setup_timer;
      Timer* m_solve_timer;
      Timer* m_shutdown_timer;
#endif

      /* Walltime measurements */
      clock_t m_main_start, m_solve_start, m_main_end, m_solve_end;

      bool m_compute_op_matrices;
      Real m_op_matrices_tolerance;

      EvalOpMatrix<Vector,System> m_op_matrices;
};

template <class Vector, class System> const Real
COGENTTimeIntegration<Vector,System>::s_DT_EPS = std::numeric_limits<Real>::epsilon();

template <class Vector, class System> inline void
COGENTTimeIntegration<Vector,System>::printParameters()
{
   pout() << "maximum step = " << m_max_step << endl;
   pout() << "maximum time = " << m_max_time << endl;
   pout() << "checkpoint interval = " << m_checkpoint_interval << endl;
   pout() << "plot interval = " << m_plot_interval << endl;
   if(m_plot_time_interval>0.0) pout() << "plot time interval = " << m_plot_time_interval << endl;
}


template <class Vector, class System> inline void
COGENTTimeIntegration<Vector,System>::loadRestartFile( ParmParse& a_ppsim )
{
   std::string restartFile;
   a_ppsim.query( "restart_file", restartFile );
   m_system->readCheckpointFile(  restartFile, 
                                  m_cur_step, 
                                  m_cur_time, 
                                  m_cur_dt );
}

static int plottedonce = 0;

template <class Vector, class System> inline void
COGENTTimeIntegration<Vector,System>::writePlotFile()
{
   if (m_verbosity >= 3) {
      pout() << "COGENTTimeIntegration<Vector,System>::writePlotFile" << endl;
   }

   if (m_write_odevec) {

     char iter_str[10];
     sprintf(iter_str, "%05d", m_cur_step);
     std::string fname = "odevec_" + std::string(iter_str);

     m_U.writeToFile( fname );

   }

#ifdef CH_USE_HDF5
   ParmParse ppsim( "simulation" );
   parseParameters( ppsim );

   char iter_str[100];

   if ( ppsim.contains( "restart_file" )&& plottedonce == 0 ) {
      sprintf( iter_str, "restart.%s", m_plot_prefix.c_str());
      plottedonce = 1;
   } else {
      sprintf( iter_str, "%s", m_plot_prefix.c_str());
   }

   m_system->writePlotFile( iter_str, m_cur_step, m_cur_time);
#else
   MayDay::Error( "plotting only defined with hdf5" );
#endif
}

template <class Vector, class System> inline void
COGENTTimeIntegration<Vector,System>::writeHistFile(bool startup_flag)
{
  // the startup_flag is used to force a write on start-up
   if (m_verbosity>=3) {
      pout() << "COGENTTimeIntegration<Vector,System>::writeHistoryFile" << endl;
   }
   m_system->writeFieldHistory(m_cur_step, m_cur_time, startup_flag);
}

template <class Vector, class System> inline void
COGENTTimeIntegration<Vector,System>::writeCheckpointFile()
{
   m_system->writeCheckpointFile( m_checkpoint_prefix, 
                                  m_cur_step, 
                                  m_cur_time, 
                                  m_cur_dt );
}

template <class Vector, class System> inline void
COGENTTimeIntegration<Vector,System>::initializeTimers()
{
#ifdef CH_USE_TIMER
   m_all_timer      = new Timer("All",0);
   m_setup_timer    = new Timer("Setup",*m_all_timer);
   m_solve_timer    = new Timer("Solve",*m_all_timer,1);
   m_shutdown_timer = new Timer("Shutdown",*n_all_timer);
   Timer::TimerInit( 0 );
#endif
}


template <class Vector, class System>
COGENTTimeIntegration<Vector,System>::COGENTTimeIntegration()
   :   m_verbosity(0),
       m_cur_step(0),
       m_max_step(0),
       m_cur_time(0.0),
       m_max_time(0.0),
       m_cur_dt(-1.0),
       m_fixed_dt(-1.0),
       m_max_dt_grow(1.1),
       m_init_dt_frac(0.1),
       m_cfl(1.0),
       m_adapt_dt(true),
       m_checkpoint_interval(0),
       m_last_checkpoint(0),
       m_checkpoint_prefix( "chk" ),
       m_plot_interval(0),
       m_plot_time_interval(0.0),
       m_plot_time(0.0),
       m_last_plot(0),
       m_plot_prefix( "plt" ),
       m_write_odevec(false),
       m_system( NULL ),
       m_compute_op_matrices(false),
       m_op_matrices_tolerance(1e-6)
#ifdef CH_USE_TIMER
,
       m_all_timer( NULL ),
       m_setup_timer( NULL ),
       m_solve_timer( NULL ),
       m_shutdown_timer( NULL )
#endif
{
}

template <class Vector, class System>
void COGENTTimeIntegration<Vector,System>::initialize(  ParmParse&    a_pp,
                                                        System* const a_sys )
{
   m_system = a_sys;
#ifdef CH_MPI
   MPI_Barrier(MPI_COMM_WORLD);
#endif
   m_main_start = clock();
#ifdef CH_USE_TIMER
   initializeTimers();
   m_all_timer->start();
   m_setup_timer->start();
#endif

#ifdef USE_ARRAYVIEW
   trapfpe();
#endif

   ParmParse ppsim( "simulation" );
   parseParameters( ppsim );

   if ( ppsim.contains( "restart_file" ) ) {
      loadRestartFile( ppsim );
   }
   m_system->initialize(m_cur_step, m_cur_time);
   m_system->printDiagnostics();

   m_system->defineODEVector(m_U);
   m_U.printLoadBalanceInfo();

   m_system->copyStateToArray( m_U );
   m_system->defineODEVector(m_U_old);

   if (m_system->tiClass() == _TI_RK_) {
      m_integrator = new TiRK<Vector, AppCtxt>;
   } else if (m_system->tiClass() == _TI_ARK_) {
      m_integrator = new TiARK<Vector, AppCtxt>;
   } else {
      MayDay::Error("Unrecognized input for ti_class.");
   }
   m_integrator->define(  a_pp, 
                          m_system->tiMethod(), 
                          m_U, 
                          m_system );
   m_integrator->initialize( m_cur_time, m_U );

   if ( m_plot_interval>=0 || m_plot_time_interval>=0.0) {
      writePlotFile();
      m_last_plot = m_cur_step;
      if ( m_plot_time_interval>=0.0 ) m_plot_time = m_plot_time_interval;
   }

   writeHistFile(true); // "true" forces a write on start-up

   if (m_compute_op_matrices) {
     m_op_matrices.define(m_system, m_U, m_op_matrices_tolerance);
   }

#ifdef CH_USE_TIMER
   setup_timer->stop();
#endif

#ifdef CH_USE_TIMER
   solve_timer->start();
#endif
#ifdef CH_MPI
   MPI_Barrier(MPI_COMM_WORLD);
#endif
   m_solve_start = clock();
}


template <class Vector, class System>
COGENTTimeIntegration<Vector,System>::~COGENTTimeIntegration<Vector,System>()
{
   delete m_integrator;
#ifdef CH_USE_TIMER
   delete m_all_timer;
   delete m_setup_timer;
   delete m_solve_timer;
   delete m_shutdown_timer;
#endif
}


template <class Vector, class System>
bool COGENTTimeIntegration<Vector,System>::notDone()
{
   return ( (m_cur_step < m_max_step) && (m_cur_time < m_max_time) );
}


template <class Vector, class System>
void COGENTTimeIntegration<Vector,System>::advance()
{
   CH_TIMERS("COGENTTimeIntegration<Vector,System>::advance");
   CH_TIMER("print_diagnostics",t_print_diagnostcs);

   preTimeStep();
   if (m_verbosity >= 1) {
      if (procID()==0) {
        printf("\nStep %4d, dt = %1.6e\n", m_cur_step+1, m_cur_dt);
      }
   }
   m_integrator->setTimeStepSize( m_cur_dt );

   m_U_old = m_U;
   Real ref_norm = m_U_old.computeNorm();

   m_integrator->advance( m_cur_time, m_U );

   m_U_old -= m_U;
   m_step_norm_abs = m_U_old.computeNorm();
   if (ref_norm > 1e-15) {
     m_step_norm_rel = m_step_norm_abs / ref_norm;
   } else {
     m_step_norm_rel = m_step_norm_abs;
   }

   m_integrator->getCurrentTime( m_cur_time );
   m_integrator->getTimeStep( m_cur_step );
   m_integrator->getTimeStepSize( m_cur_dt );

   postTimeStep();

#ifdef CH_MPI
   MPI_Barrier(MPI_COMM_WORLD);
#endif
   clock_t m_now = clock();
   double walltime, walltime_g;
   walltime = ((double) (m_now - m_solve_start)) / CLOCKS_PER_SEC;
#ifdef CH_MPI
   MPI_Allreduce(&walltime,&walltime_g,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#else
   walltime_g = walltime;
#endif

   CH_START(t_print_diagnostcs);
   if (m_verbosity >= 1) {
      m_system->printDiagnostics();
      if(procID()==0) {
        printf( "  Step norm: %1.4e (rel), %1.4e (abs)\n", 
                m_step_norm_rel, m_step_norm_abs);
        printf("\n");
        printf( "Step %4d completed, simulation time is %f, ",
                m_cur_step, m_cur_time );
        printf("solver wall time is %f seconds\n", walltime_g);
        printf("----\n");
      }
   }
   CH_STOP(t_print_diagnostcs);

   if ( m_plot_time_interval>0.0 ) {
      if ( m_cur_time>m_plot_time-m_cur_dt/1000.0 ) {
         writePlotFile();
         m_last_plot = m_cur_step;
         m_plot_time = m_plot_time + m_plot_time_interval;
      } 
   } 
   else {
      if ( (m_cur_step % m_plot_interval)==0 ) {
         writePlotFile();
         m_last_plot = m_cur_step;
      }
   }

   writeHistFile(false);

   if ( (m_cur_step % m_checkpoint_interval)==0 ) {
      if (!procID()) {
         printf("----\n");
      }
      m_system->printFunctionCounts();
      if (!procID()) {
         printf("----\n");
      }
      m_integrator->printCounts();
      if (!procID()) {
         printf("----\n");
      }
      writeCheckpointFile();
      m_last_checkpoint = m_cur_step;
   }
}


template <class Vector, class System>
void COGENTTimeIntegration<Vector,System>::finalize()
{
#ifdef CH_MPI
   MPI_Barrier(MPI_COMM_WORLD);
#endif
  m_solve_end = clock();
#ifdef CH_USE_TIMER
   solve_timer->stop();
#endif
#ifdef CH_USE_TIMER
   shutdown_timer->start() ;
#endif

   m_system->copyStateFromArray( m_U );
   
   if (     (     m_plot_interval >= 0 ||  m_plot_time_interval >= 0.0  ) 
        &&  ( m_last_plot!=m_cur_step )   ) {
      writePlotFile();
   }

   if ( m_last_checkpoint!=m_cur_step ) {
      if (!procID()) {
         cout << "----\n";
      }
      m_system->printFunctionCounts();
      if (!procID()) {
        cout << "----\n";
      }
      m_integrator->printCounts();
      if (!procID()) {
        cout << "----\n";
      }
      writeCheckpointFile();
   }

#ifdef CH_MPI
   MPI_Barrier(MPI_COMM_WORLD);
#endif
   m_main_end = clock();

   double main_walltime = ((double) (m_main_end - m_main_start)) / CLOCKS_PER_SEC;
   double solve_walltime = ((double) (m_solve_end - m_solve_start)) / CLOCKS_PER_SEC;
   double main_walltime_g, solve_walltime_g;
#ifdef CH_MPI
   MPI_Allreduce(&main_walltime,&main_walltime_g,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
   MPI_Allreduce(&solve_walltime,&solve_walltime_g,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#else
   main_walltime_g = main_walltime;
   solve_walltime_g = solve_walltime;
#endif
   if (!procID()) {
     cout << "Solve wall time (in seconds): " << solve_walltime_g << "\n";
     cout << "Total wall time (in seconds): " << main_walltime_g << "\n";
   }

#ifdef CH_USE_TIMER
   m_shutdown_timer->stop() ;
   m_all_timer->stop() ;
   Timer::TimerSummary();
#endif
}


template <class Vector, class System>
void COGENTTimeIntegration<Vector,System>::parseParameters( ParmParse& a_ppsim )
{
   // This determines the amount of diagnositic output generated
   a_ppsim.query( "verbosity", m_verbosity );
   CH_assert( m_verbosity >= 0 );

   // Stop after this number of steps
   a_ppsim.query( "max_step", m_max_step );
   CH_assert( m_max_step >= 0 );

   // Stop when the simulation time get here
   a_ppsim.query( "max_time", m_max_time );
   CH_assert( m_max_time >= 0.0 );

   // If set, use as the fixed time step size
   if ( a_ppsim.query( "fixed_dt", m_fixed_dt ) ) {
      CH_assert( m_fixed_dt>0.0 );
      m_adapt_dt = false;
   }

   // Multiply by which to increase dt each step
   if ( a_ppsim.query( "max_dt_grow", m_max_dt_grow ) ) {
      CH_assert( m_max_dt_grow>1.0 );
      if (!m_adapt_dt) {
         MayDay::Error( "fixed_dt and max_dt_grow are mutually exclusive!" );
      }
   }

   // If set, reduce initial stable time step to this fraction
   if ( a_ppsim.query( "initial_dt_fraction", m_init_dt_frac ) ) {
      CH_assert( m_init_dt_frac>0.0 && m_init_dt_frac<=1.0 );
      if (!m_adapt_dt) {
         MayDay::Error( "fixed_dt and initial_dt_frac are mutually exclusive!" );
      }
   }

   // set cfl number for the case of dynamic timestep selection
   if ( a_ppsim.query( "cfl_number", m_cfl ) ) {
      CH_assert( m_cfl>0.0 && m_cfl<=2.0 );
      if (!m_adapt_dt) {
         MayDay::Error( "fixed_dt and cfl are mutually exclusive!" );
      }
   }

   // Set up checkpointing
   a_ppsim.query( "checkpoint_interval", m_checkpoint_interval );
   a_ppsim.query( "checkpoint_prefix", m_checkpoint_prefix );

   // Set up plot file writing
   a_ppsim.query( "plot_interval", m_plot_interval );
   a_ppsim.query( "plot_time_interval", m_plot_time_interval );
   if( m_plot_time_interval>0.0 ) m_plot_interval = 0;
   a_ppsim.query( "plot_prefix", m_plot_prefix );
   a_ppsim.query( "write_ode_vec", m_write_odevec );


   // compute linearized operator matrices?
   // Caution: a *very* slow process
   a_ppsim.query("compute_op_matrices",m_compute_op_matrices);
   a_ppsim.query("op_matrices_tolerance",m_op_matrices_tolerance);

   if (m_verbosity) {
      printParameters();
   }
}


template <class Vector, class System>
inline void COGENTTimeIntegration<Vector,System>::setFixedTimeStep( const Real& )
{
   m_cur_dt = m_fixed_dt; 
}


template <class Vector, class System>
void COGENTTimeIntegration<Vector,System>::postTimeStep()
{
   m_system->postTimeStep(  m_U, 
                            m_cur_step, 
                            m_cur_dt, 
                            m_cur_time );
}

template <class Vector, class System>
void COGENTTimeIntegration<Vector,System>::preTimeStep()
{
   m_integrator->setCurrentTime( m_cur_time );
   m_integrator->setTimeStep( m_cur_step );

   m_system->preTimeStep( m_U, m_cur_step, m_cur_time );

   Real dt_stable = m_system->stableDt( m_cur_step, m_integrator->type() );
   dt_stable *= m_cfl;
   CH_assert( dt_stable > 1.0e-16 );

   if ( m_cur_time>0.0 ) { 
      // not initial time step
      if ( m_adapt_dt ) { 
         // adjustable time step
         m_cur_dt = std::min( dt_stable, m_max_dt_grow * m_cur_dt );
      } else {                 
         // fixed time step
         setFixedTimeStep( dt_stable );
      }
   } else { 
      // initial time step
      if ( m_adapt_dt ) { 
         // adjustable time step
         m_cur_dt = m_init_dt_frac * dt_stable;
      } else {                 
         // fixed time step
         setFixedTimeStep( dt_stable );
      }
   }
   
   // If less than a time step from the final time, adjust time step
   // to end just over the final time.
   Real timeRemaining = m_max_time - m_cur_time;
   if ( m_cur_dt > timeRemaining ) {
      m_cur_dt = timeRemaining + m_max_time * s_DT_EPS;
   }

}

#include "NamespaceFooter.H"
#endif // _COGENTTimeIntegration_H_
