#include "Simulation.H"
#include <limits>

#include "NamespaceHeader.H"

template <class SYSTEM> const Real
Simulation<SYSTEM>::s_DT_EPS = std::numeric_limits<Real>::epsilon();

template <class SYSTEM> inline void
Simulation<SYSTEM>::printParameters()
{
   pout() << "maximum step = " << m_max_step << endl;
   pout() << "maximum time = " << m_max_time << endl;
   pout() << "checkpoint interval = " << m_checkpoint_interval << endl;
   pout() << "plot interval = " << m_plot_interval << endl;
}


template <class SYSTEM> inline void
Simulation<SYSTEM>::loadRestartFile( ParmParse& a_ppsim )
{
   std::string restartFile;
   a_ppsim.query( "restart_file", restartFile );
#ifdef CH_USE_HDF5
   HDF5Handle handle( restartFile, HDF5Handle::OPEN_RDONLY );
   m_system->readCheckpointFile( handle, m_cur_step, m_cur_time, m_cur_dt );
   handle.close();
#else
   MayDay::Error("restart only defined with hdf5");
#endif
}

static int plottedonce = 0;

template <class SYSTEM> inline void
Simulation<SYSTEM>::writePlotFile()
{
   if (m_verbosity >= 3) {
      pout() << "Simulation<SYSTEM>::writePlotFile" << endl;
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

   //if (m_verbosity >= 2) {
   //   pout() << "plot file name = " << iter_str << endl;
   //}

   // HDF5Handle handle( iter_str, HDF5Handle::CREATE );
   // Instead of dummy file with prefix, use prefix for all plot files.
   m_system->writePlotFile( iter_str, m_cur_step, m_cur_time);
   // handle.close();
#else
   MayDay::Error( "plotting only defined with hdf5" );
#endif
}

template <class SYSTEM> inline void
Simulation<SYSTEM>::writeHistFile(bool startup_flag)
{
  // the startup_flag is used to force a write on start-up
   if (m_verbosity>=3) {
      pout() << "Simulation<SYSTEM>::writeHistoryFile" << endl;
   }
   m_system->writeFieldHistory(m_cur_step, m_cur_time, startup_flag);
}

template <class SYSTEM> inline void
Simulation<SYSTEM>::writeCheckpointFile()
{
   if (m_verbosity>=3) {
      pout() << "Simulation<SYSTEM>::writeCheckpointFile" << endl;
   }

#ifdef CH_USE_HDF5
   char iter_str[100];

   sprintf( iter_str, "%s%04d.%dd.hdf5",
            m_checkpoint_prefix.c_str(), m_cur_step, SpaceDim );

   if (m_verbosity>=2) {
      pout() << "checkpoint file name = " << iter_str << endl;
   }

   HDF5Handle handle( iter_str, HDF5Handle::CREATE );
   m_system->writeCheckpointFile( handle, m_cur_step, m_cur_time, m_cur_dt );
   handle.close();
#else
   MayDay::Error( "restart only defined with hdf5" );
#endif
}

template <class SYSTEM> inline void
Simulation<SYSTEM>::initializeTimers()
{
#ifdef CH_USE_TIMER
   m_all_timer      = new Timer("All",0);
   m_setup_timer    = new Timer("Setup",*m_all_timer);
   m_solve_timer    = new Timer("Solve",*m_all_timer,1);
   m_shutdown_timer = new Timer("Shutdown",*n_all_timer);
   Timer::TimerInit( 0 );
#endif
}


template <class SYSTEM>
Simulation<SYSTEM>::Simulation( ParmParse& a_pp )
   :   m_verbosity(0),
       m_cur_step(0),
       m_max_step(0),
       m_cur_time(0.0),
       m_max_time(0.0),
       m_cur_dt(-1.0),
       m_fixed_dt(-1.0),
       m_fixed_dt_subiteration(false),
       m_sub_dt(-1.0),
       m_max_dt_grow(1.1),
       m_init_dt_frac(0.1),
       m_cfl(1.0),
       m_adapt_dt(true),
       m_checkpoint_interval(0),
       m_last_checkpoint(0),
       m_checkpoint_prefix( "chk" ),
       m_plot_interval(0),
       m_last_plot(0),
       m_plot_prefix( "plt" ),
       m_system( NULL )
#ifdef CH_USE_TIMER
,
       m_all_timer( NULL ),
       m_setup_timer( NULL ),
       m_solve_timer( NULL ),
       m_shutdown_timer( NULL )
#endif
{
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

   m_system = new SYSTEM( a_pp );
   if ( ppsim.contains( "restart_file" ) ) {
      loadRestartFile( ppsim );
   }
   m_system->initialize(m_cur_step);
   m_system->printDiagnostics();

   if ( m_plot_interval>=0 ) {
      writePlotFile();
      m_last_plot = m_cur_step;
   }

   writeHistFile(true); // "true" forces a write on start-up


#ifdef CH_USE_TIMER
   setup_timer->stop();
#endif

#ifdef CH_USE_TIMER
   solve_timer->start();
#endif
}


template <class SYSTEM>
Simulation<SYSTEM>::~Simulation<SYSTEM>()
{
   delete m_system;
#ifdef CH_USE_TIMER
   delete m_all_timer;
   delete m_setup_timer;
   delete m_solve_timer;
   delete m_shutdown_timer;
#endif
}


template <class SYSTEM>
bool Simulation<SYSTEM>::notDone()
{
   return ( (m_cur_step<m_max_step) && (m_cur_time<m_max_time) );
}


template <class SYSTEM>
void Simulation<SYSTEM>::advance()
{
   if (m_verbosity >= 1) {
      pout() << endl << "Step " << m_cur_step+1 << endl;
      if(procID()==0) cout << endl << "Step " << m_cur_step+1 << endl;
   }
   selectTimeStep();

   m_cur_step++;
   if (m_fixed_dt_subiteration) {
      for (int sub_iter=1; sub_iter<=m_subiterations; sub_iter++) {
         if (!procID()) {
            cout << "  --\n";
            cout << "  Subiteration:" << sub_iter << " of " << m_subiterations << endl;
         }
         m_cur_time = m_system->advance( m_cur_time, m_sub_dt, m_cur_step );
      }
   } else {
      m_cur_time = m_system->advance( m_cur_time, m_cur_dt, m_cur_step );
   }

   postTimeStep();

   if (m_verbosity >= 1) {
      m_system->printDiagnostics();
      pout() << "Step " << m_cur_step << " completed, simulation time is "
             << m_cur_time << endl;
      if(procID()==0) cout << "Step " << m_cur_step << " completed, simulation time is "
                           << m_cur_time << endl << "----" << endl;
   }

   if ( (m_cur_step % m_plot_interval)==0 ) {
      writePlotFile();
      m_last_plot = m_cur_step;
   }

   writeHistFile(false);

   if ( (m_cur_step % m_checkpoint_interval)==0 ) {
      writeCheckpointFile();
      m_last_checkpoint = m_cur_step;
   }
}


template <class SYSTEM>
void Simulation<SYSTEM>::finalize()
{
#ifdef CH_USE_TIMER
   solve_timer->stop();
#endif
#ifdef CH_USE_TIMER
   shutdown_timer->start() ;
#endif

   if ( (m_plot_interval >= 0) && (m_last_plot!=m_cur_step) ) {
      writePlotFile();
   }

   if ( m_last_checkpoint!=m_cur_step ) {
      writeCheckpointFile();
   }

#ifdef CH_USE_TIMER
   m_shutdown_timer->stop() ;
   m_all_timer->stop() ;
   Timer::TimerSummary();
#endif
}


template <class SYSTEM>
void Simulation<SYSTEM>::parseParameters( ParmParse& a_ppsim )
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
      a_ppsim.query( "fixed_dt_subiteration", m_fixed_dt_subiteration ); //this member has false value by default in constructor 
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

   /* possibly redundant check */
   if (m_adapt_dt) m_fixed_dt_subiteration = false;

   // Set up checkpointing
   a_ppsim.query( "checkpoint_interval", m_checkpoint_interval );
   a_ppsim.query( "checkpoint_prefix", m_checkpoint_prefix );

   // Set up plot file writing
   a_ppsim.query( "plot_interval", m_plot_interval );
   a_ppsim.query( "plot_prefix", m_plot_prefix );

   // History parameter parsing moved to GKSystem.cpp

   if (m_verbosity) {
      printParameters();
   }
}


template <class SYSTEM>
inline void Simulation<SYSTEM>::setFixedTimeStep( const Real& a_dt_stable )
{
   if (m_fixed_dt_subiteration) {
      /* subiteration to satisfy stable dt */
      m_subiterations = ceil(m_fixed_dt/a_dt_stable);
      m_sub_dt = m_fixed_dt/m_subiterations;
      if (m_fixed_dt > a_dt_stable) {
         if (!procID()) {
            cout << "  --\n";
            cout << "  Specified time step is higher than stable time step. Using subiterations.\n";
            cout << "  " << m_subiterations << " subiterations will be made with sub_dt = ";
            cout << m_sub_dt << endl;
            cout << "  You may disable this by setting simulation.fixed_dt_subiteration = false.\n";
            cout << "  --\n";
         }
      }
      m_cur_dt = m_sub_dt; 
   } else {
      if (m_fixed_dt > a_dt_stable) {
         if (!procID()) {
            cout << "  --\n";
            cout << "  Warning: fixed time step may be higher than the stable time step.\n";
            cout << "  Stable time step = " << a_dt_stable << ".\n";
            cout << "  You may reduce the specified dt or enable subiteration by setting\n";
            cout << "  simulation.fixed_dt_subiteration = true (default:false).\n";
            cout << "  --\n";
         }
      }
      m_cur_dt = m_fixed_dt; 
   }
}


template <class SYSTEM>
void Simulation<SYSTEM>::selectTimeStep()
{
   Real dtStable = m_system->stableDt(m_cur_step) * m_cfl;
   CH_assert(dtStable > 1.0e-16);

   if ( m_cur_time>0.0 ) { 
      // not initial time step
      if ( m_adapt_dt ) { 
         // adjustable time step
         m_cur_dt = std::min( dtStable, m_max_dt_grow * m_cur_dt );
      } else {                 
         // fixed time step
         setFixedTimeStep(dtStable);
      }
   } else { 
      // initial time step
      if ( m_adapt_dt ) { 
         // adjustable time step
         m_cur_dt = m_init_dt_frac * dtStable;
      } else {                 
         // fixed time step
         setFixedTimeStep(dtStable);
      }
   }

   // If less than a time step from the final time, adjust time step
   // to end just over the final time.
   Real timeRemaining = m_max_time - m_cur_time;
   if ( m_cur_dt > timeRemaining ) {
      m_cur_dt = timeRemaining + m_max_time * s_DT_EPS;
   }
}


template <class SYSTEM>
void Simulation<SYSTEM>::postTimeStep()
{
   m_system->printTimeStep(m_cur_step,m_cur_dt);
}

#include "NamespaceFooter.H"
