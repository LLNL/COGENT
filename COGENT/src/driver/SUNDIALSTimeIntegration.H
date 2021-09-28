#ifdef with_sundials

#ifndef _SUNDIALSTIMEINTEGRATOR_H_
#define _SUNDIALSTIMEINTEGRATOR_H_

#include <string>
#include <time.h>
#include "Simulation.H"

#include "ODEVector.H"
#include "AppCtxt.H"
#include "MultiPhysicsPreconditioner.H"
#include "TiButcherTables.H"

/* include SUNDIALS header files */
#include "arkode/arkode_arkstep.h"
#include "sunnonlinsol/sunnonlinsol_newton.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "nvector/nvector_parallel.h"

#include "NamespaceHeader.H"

template <class AppVec, class System>
class SundialsTimeIntegrator : public Simulation<AppVec,System>
{
  public: 
    SundialsTimeIntegrator<AppVec,System>();
    ~SundialsTimeIntegrator<AppVec,System>();

    void initialize(ParmParse&, System* const);
    void finalize();
    void solve();

    inline
    void cpAppVecToSundials(  N_Vector&     a_Y, 
                              const AppVec& a_U )
    {
      realtype *Yarr = N_VGetArrayPointer(a_Y);
      a_U.copyTo( (Real* const) Yarr );
    }

    inline
    void cpAppVecFromSundials(  const N_Vector& a_Y, 
                                AppVec&         a_U )
    {
      const realtype *Yarr = N_VGetArrayPointer(a_Y);
      a_U.copyFrom( (const Real* const) Yarr );
    }

    inline System*    getSystem()                 { return m_system; }
    inline int        getVerbosity()              { return m_verbosity; }
    inline int        getPlotInterval()           { return m_plot_interval; }
    inline int        getChkptInterval()          { return m_checkpoint_interval; }
    inline Real       getTime()                   { return m_cur_time; }
    inline Real       getShift()                  { return m_shift; }
    inline bool       isLinear()                  { return m_is_linear; }
    inline void*      getTIObj()                  { return m_ti_obj; }

    inline void    setLastPlot  (int a_int) { m_last_plot = a_int; }      
    inline void    setLastChkpt (int a_int) { m_last_checkpoint = a_int; }      

    inline void setTime           (Real a_t)    { m_cur_time = a_t; }
    inline void setShift          (Real a_a)    { m_shift = a_a; }
    inline void setTimeStep       (Real a_dt)   { m_dt = a_dt; }
    inline void setCFL            (Real a_cfl)  { m_cfl = a_cfl; }

    inline Real getCFL()  { return m_cfl; }

    void writeHistFile(bool);
    void writePlotFile();
    void writeCheckpointFile();

    void preTimeStep(Real);
    void postTimeStep(Real);

    AppVec m_U, m_Uold, m_F;
    Preconditioner<AppVec,System>* m_preCond;

  private:

    void parseParametersSimulation( ParmParse& );
    void parseParametersSUNDIALS( ParmParse& );
    void printTIStats();

    System*   m_system;
    long int  m_cur_step;
    long int  m_max_step;
    Real      m_cur_time;
    Real      m_max_time;
    Real      m_dt;
    Real      m_cfl;
    Real      m_shift;
    int       m_verbosity;
    int       m_plottedonce;

    int  m_checkpoint_interval;
    int  m_last_checkpoint;

    int  m_plot_interval;
    int  m_last_plot;

    int  m_ark_order;
    Real m_ark_rtol;
    Real m_ark_atol;

    std::string m_checkpoint_prefix;
    std::string m_plot_prefix;

    int  m_local_size,
         m_global_size;

    bool  m_is_linear;

    N_Vector m_Y;
    void* m_ti_obj;
    SUNNonlinearSolver m_nonlinsolver;
    SUNLinearSolver m_linsolver;

    bool m_linsolver_verbose;
    bool m_nonlinsolver_verbose;
    int m_linsolver_maxits;
    int m_nonlinsolver_maxits;

    bool m_adapt_step;

    Real m_step_norm_abs, m_step_norm_rel;

    /* Walltime measurements */
    clock_t m_main_start, m_solve_start, m_main_end, m_solve_end;
};

static int RHSFunction( realtype  a_t,
                        N_Vector  a_Y,
                        N_Vector  a_F, 
                        void*     a_ctxt )
{
  SundialsTimeIntegrator<ODEVector,AppCtxt> *context 
      = (SundialsTimeIntegrator<ODEVector,AppCtxt>*) a_ctxt;
  AppCtxt* system = context->getSystem();

  context->cpAppVecFromSundials( a_Y, context->m_U );
  system->preOpEval( context->m_U, (Real)a_t, chkpt_naive );
  system->explicitOp( context->m_F, (Real)a_t, context->m_U );
  context->cpAppVecToSundials( a_F, context->m_F);

  return 0;
}

static int RHSFunctionIMEX( realtype  a_t,
                            N_Vector  a_Y,
                            N_Vector  a_F, 
                            void*     a_ctxt )
{
  SundialsTimeIntegrator<ODEVector,AppCtxt> *context 
      = (SundialsTimeIntegrator<ODEVector,AppCtxt>*) a_ctxt;
  AppCtxt* system = context->getSystem();

  context->cpAppVecFromSundials( a_Y, context->m_U );
  system->preOpEval( context->m_U, (Real)a_t, chkpt_naive );
  system->explicitOpImEx( context->m_F, (Real)a_t, context->m_U );
  context->cpAppVecToSundials( a_F, context->m_F);

  return 0;
}

static int IFunctionIMEX( realtype  a_t,
                          N_Vector  a_Y,
                          N_Vector  a_F, 
                          void*     a_ctxt )
{
  SundialsTimeIntegrator<ODEVector,AppCtxt> *context 
      = (SundialsTimeIntegrator<ODEVector,AppCtxt>*) a_ctxt;
  AppCtxt* system = context->getSystem();

  context->cpAppVecFromSundials( a_Y, context->m_U );
  system->preOpEval( context->m_U, (Real)a_t, chkpt_naive );
  system->implicitOpImEx( context->m_F, (Real)a_t, context->m_U );
  context->cpAppVecToSundials( a_F, context->m_F);

  return 0;
}

static int PCUpdateIMEX(  realtype      a_t,
                          N_Vector      a_Y,
                          N_Vector      a_F,
                          booleantype   a_jok,
                          booleantype*  a_jcurPtr,
                          realtype      a_gamma,
                          void*         a_ctxt )
{
  SundialsTimeIntegrator<ODEVector,AppCtxt> *context 
      = (SundialsTimeIntegrator<ODEVector,AppCtxt>*) a_ctxt;
  void* ti_obj = context->getTIObj();

  long int cur_step, stage = 0;
  ARKStepGetNumSteps( ti_obj, &cur_step );

  Real cur_dt;
  ARKStepGetCurrentStep( ti_obj, &cur_dt );
  Real shift = 1.0 / (a_gamma * cur_dt);
  context->setShift( shift );

  context->cpAppVecFromSundials( a_Y, context->m_U );
  context->m_preCond->update( context->m_U,
                              (Real) a_t,
                              cur_step+1,
                              stage,
                              shift );

  *a_jcurPtr = true;
  return 0;
}

static int PCFunctionIMEX(  realtype  a_t,
                            N_Vector  a_Y0,
                            N_Vector  a_F0,
                            N_Vector  a_X,
                            N_Vector  a_Y,
                            realtype  a_gamma,
                            realtype  a_delta,
                            int       a_lr,
                            void*     a_ctxt )
{
  SundialsTimeIntegrator<ODEVector,AppCtxt>* context
   = (SundialsTimeIntegrator<ODEVector,AppCtxt>*) a_ctxt;

  context->cpAppVecFromSundials( a_X, context->m_U );
  context->m_F.zero();
  context->m_preCond->apply(context->m_F, context->m_U);
  context->m_F.scale( context->getShift() );
  context->cpAppVecToSundials( a_Y, context->m_F );

  return 0;
}

template <class AppVec, class System>
void SundialsTimeIntegrator<AppVec,System>::preTimeStep(Real a_time)
{
  m_Uold = m_U;
  m_system->preTimeStep(m_U, m_cur_step, a_time);

  if (m_cfl > 0) {

    Real dt_stable = 0.0;
    if ( m_system->tiClass() == _TI_ARK_ ) {
      dt_stable = m_system->stableDt( m_cur_step, ti_imex );
    } else if ( m_system->tiClass() == _TI_RK_ ) {
      dt_stable = m_system->stableDt( m_cur_step, ti_explicit );
    } else {
      MayDay::Error("SUNDIALS interface: time integrator type not supported!");
    }
    CH_assert( dt_stable > 1.0e-16 );
    m_dt = dt_stable*m_cfl;

  }

  if (!m_adapt_step) {
    ARKStepSetFixedStep( m_ti_obj, m_dt ); 
  }

  return;
}

template <class AppVec, class System>
void SundialsTimeIntegrator<AppVec,System>::postTimeStep(Real a_time)
{
  m_cur_time = a_time;
  ARKStepGetNumSteps( m_ti_obj, &m_cur_step );
  ARKStepGetCurrentStep( m_ti_obj, &m_dt );
  m_system->postTimeStep(m_U, a_time, m_dt, a_time);

  Real ref_norm = m_Uold.computeNorm();
  m_Uold -= m_U;
  m_step_norm_abs = m_Uold.computeNorm();
  if (ref_norm > 1e-15) {
    m_step_norm_rel = m_step_norm_abs / ref_norm;
  } else {
    m_step_norm_rel = m_step_norm_abs;
  }

  if ( ((m_cur_step) % m_plot_interval)==0 ) {
    writePlotFile();
    setLastPlot(m_cur_step);
  }

  writeHistFile(false);

  if ( ((m_cur_step) % getChkptInterval())==0 ) {
    writeCheckpointFile();
    setLastChkpt(m_cur_step);
  }

  return;
}

static int postTimeStageSD( realtype  a_t, 
                            N_Vector  a_Y,
                            void*     a_ctxt )
{
  SundialsTimeIntegrator<ODEVector,AppCtxt>* context
   = (SundialsTimeIntegrator<ODEVector,AppCtxt>*) a_ctxt;
  AppCtxt* system = context->getSystem();
  void* ti_obj = context->getTIObj();

  long int cur_step;
  ARKStepGetNumSteps( ti_obj, &cur_step );
  Real dt;
  ARKStepGetCurrentStep( ti_obj, &dt );

  static int step_number  = 0;
  static int stage_number = -1;

  if (cur_step == step_number+1) {
    step_number++;
    stage_number = -1;
  } else if (cur_step == step_number) {
    stage_number++;
  } else {
    if (!procID()) {
      std::cout << "Error in postTimeStage() "
                << "in src/driver/SUNDIALSTimeIntegration.H: cur_step=" 
                << cur_step 
                << ", step_number=" 
                << step_number << "."
                << " Something has caused stage-counting to go haywire.\n";
    }
  }

  context->cpAppVecFromSundials( a_Y, context->m_U );
  system->postTimeStage(  cur_step,
                          a_t,
                          dt,
                          context->m_U,
                          stage_number );
  context->cpAppVecToSundials( a_Y, context->m_U );

  return 0;
}

static int postTimeStepSD(realtype  a_t, 
                          N_Vector  a_Y,
                          void*     a_ctxt )
{
  SundialsTimeIntegrator<ODEVector,AppCtxt>* context
   = (SundialsTimeIntegrator<ODEVector,AppCtxt>*) a_ctxt;

  context->cpAppVecFromSundials( a_Y, context->m_U );
  context->postTimeStep(a_t);
  context->preTimeStep(a_t);

  return 0;
}

template <class AppVec, class System>
void SundialsTimeIntegrator<AppVec,System>::writeHistFile(bool flag)
{
  m_system->writeFieldHistory(m_cur_step, m_cur_time,flag);
}

template <class AppVec, class System> 
void SundialsTimeIntegrator<AppVec,System>::writeCheckpointFile()
{
   m_system->writeCheckpointFile( m_checkpoint_prefix, 
                                  m_cur_step, 
                                  m_cur_time, 
                                  m_dt );
}

template <class AppVec, class System> 
void SundialsTimeIntegrator<AppVec,System>::writePlotFile()
{
#ifdef CH_USE_HDF5
   ParmParse ppsim( "simulation" );
   parseParametersSimulation( ppsim );
   ParmParse ppsun( "sundials" );
   parseParametersSUNDIALS( ppsun );

   char iter_str[100];
   if ( ppsim.contains( "restart_file" )&& m_plottedonce == 0 ) {
      sprintf( iter_str, "restart.%s", m_plot_prefix.c_str());
      m_plottedonce = 1;
   } else {
      sprintf( iter_str, "%s", m_plot_prefix.c_str());
   }
   m_system->writePlotFile( iter_str, m_cur_step, m_cur_time);
#else
   MayDay::Error( "plotting only defined with hdf5" );
#endif
}

template <class AppVec, class System>
SundialsTimeIntegrator<AppVec,System>::SundialsTimeIntegrator( )
  : m_preCond(NULL),
    m_system(NULL),
    m_cur_step(0),
    m_max_step(0),
    m_cur_time(0.0),
    m_max_time(0.0),
    m_dt(0.0),
    m_cfl(-1.0),
    m_shift(0.0),
    m_verbosity(0),
    m_plottedonce(0),
    m_checkpoint_interval(0),
    m_last_checkpoint(0),
    m_plot_interval(0),
    m_last_plot(0),
    m_ark_order(-1),
    m_ark_rtol(1e-6),
    m_ark_atol(1e-12),
    m_checkpoint_prefix( "chk" ),
    m_plot_prefix( "plt" ),
    m_local_size(0),
    m_global_size(0),
    m_is_linear("true"),
    m_Y(NULL),
    m_ti_obj(NULL),
    m_nonlinsolver(NULL),
    m_linsolver(NULL),
    m_linsolver_verbose(false),
    m_nonlinsolver_verbose(false),
    m_linsolver_maxits(50),
    m_nonlinsolver_maxits(10),
    m_adapt_step(false)
{
}

template <class AppVec, class System>
void SundialsTimeIntegrator<AppVec,System>::initialize( ParmParse&    a_pp,
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


  m_main_start = clock();
  ParmParse ppsim( "simulation" );
  parseParametersSimulation( ppsim );

  if ( ppsim.contains( "restart_file" ) ) {
     //loadRestartFile( ppsim );
  }
  m_system->initialize(m_cur_step, m_cur_time);
  m_system->printDiagnostics();

  m_system->defineODEVector( m_U );
  m_U.printLoadBalanceInfo();

  m_system->copyStateToArray( m_U );
  m_F.define( m_U );

  if ( m_plot_interval >= 0 ) {
    writePlotFile();
    m_last_plot = m_cur_step;
  }
  writeHistFile(true);

  m_is_linear = m_system->isLinear();
  
  /* create and initialize the SUNDIALS solution vector */
  m_local_size = m_U.getVectorSize();
  m_global_size = 0;
#if CH_MPI
  MPI_Allreduce(  &m_local_size,
                  &m_global_size,
                  1,
                  MPI_INT,
                  MPI_SUM,
                  MPI_COMM_WORLD );
#else
  m_global_size = m_local_size;
#endif
  m_Y = N_VNew_Parallel(  MPI_COMM_WORLD, 
                          m_local_size, 
                          m_global_size );

  /* copy initial solution to SUNDIALS vector */
  cpAppVecToSundials( m_Y, m_U );

  /* create and initialize the SUNDIALS time integration object */
  if (m_system->tiClass() == _TI_RK_) {
    m_ti_obj = ARKStepCreate( RHSFunction,
                              NULL,
                              m_cur_time,
                              m_Y );
    if (m_ark_order > 0) {
      ARKStepSetOrder(  m_ti_obj, m_ark_order );
      if (!procID()) {
        std::cout << "SUNDIALS: using SUNDIALS method of order "
                  << m_ark_order << "\n";
      }
    } else {
      Real  *A = NULL, 
            *b = NULL, 
            *c = NULL, 
            *bemb = NULL;
      int nstages, order, order_emb;
      TiButcherTables::getButcherTableRK( m_system->tiMethod(),
                                          nstages,
                                          order,
                                          order_emb,
                                          &A, &b, &c, &bemb );
      ARKodeButcherTable butcher_tab = NULL;
      butcher_tab = ARKodeButcherTable_Create(  nstages,
                                                order,
                                                order_emb,
                                                c, A, b, bemb );
      ARKStepSetTables( m_ti_obj, 
                        order,
                        order_emb,
                        NULL,
                        butcher_tab );

      if (!procID()) {
        std::cout << "SUNDIALS: using provided method "
                  << m_system->tiClass()
                  << "-"
                  << m_system->tiMethod() << "\n";
      }
    }

  } else if (m_system->tiClass() == _TI_ARK_) {
    /* linear solver */
    m_linsolver = SUNLinSol_SPGMR(  m_Y,
                                    PREC_RIGHT,
                                    m_linsolver_maxits );
    SUNLinSolSetPrintLevel_SPGMR( m_linsolver, 
                                  (m_linsolver_verbose ? 1 : 0) );
    /* nonlinear solver */
    m_nonlinsolver = SUNNonlinSol_Newton( m_Y );
    SUNNonlinSolSetPrintLevel_Newton( m_nonlinsolver, 
                                      (m_nonlinsolver_verbose ? 1 : 0));
    /* ARK time integrator */
    m_ti_obj = ARKStepCreate( RHSFunctionIMEX,
                              IFunctionIMEX,
                              m_cur_time,
                              m_Y );

    ARKStepSetNonlinearSolver( m_ti_obj, m_nonlinsolver );
    ARKStepSetMaxNonlinIters( m_ti_obj, m_nonlinsolver_maxits );
    ARKStepSetLinearSolver( m_ti_obj, m_linsolver, NULL);
    ARKStepSetPreconditioner( m_ti_obj, PCUpdateIMEX, PCFunctionIMEX );
    ARKStepSetLSetupFrequency( m_ti_obj, -1 );
    if (m_is_linear) ARKStepSetLinear( m_ti_obj, 1);

    if (m_ark_order > 0) {
      ARKStepSetOrder(  m_ti_obj, m_ark_order );
      if (!procID()) {
        std::cout << "SUNDIALS: using SUNDIALS method of order "
                  << m_ark_order << "\n";
      }
    } else {
      Real  *Ai = NULL, 
            *Ae = NULL, 
            *bi = NULL, 
            *be = NULL, 
            *c = NULL, 
            *bembi = NULL,
            *bembe = NULL,
            *binterpi = NULL,
            *binterpe = NULL;
      int nstages, order, order_emb, pinterp;
      TiButcherTables::getButcherTableARK(  m_system->tiMethod(),
                                            nstages,
                                            order,
                                            order_emb,
                                            &Ai, &Ae,
                                            &bi, &be,
                                            &c, 
                                            &bembi, &bembe,
                                            pinterp,
                                            &binterpi, &binterpe );
      ARKodeButcherTable butcher_tab_i = NULL;
      ARKodeButcherTable butcher_tab_e = NULL;
      butcher_tab_i = ARKodeButcherTable_Create(  nstages,
                                                  order,
                                                  order_emb,
                                                  c, Ai, bi, bembi );
      butcher_tab_e = ARKodeButcherTable_Create(  nstages,
                                                  order,
                                                  order_emb,
                                                  c, Ae, be, bembe );
      ARKStepSetTables( m_ti_obj, 
                        order,
                        order_emb,
                        butcher_tab_i,
                        butcher_tab_e );

      if (!procID()) {
        std::cout << "SUNDIALS: using provided method "
                  << m_system->tiClass()
                  << "-"
                  << m_system->tiMethod() << "\n";
      }
    }

    m_preCond = new MultiPhysicsPreconditioner<AppVec,System>;
    m_preCond->define(  m_U,
                        (void*) m_system,
                        std::string("ARK.PC"),
                        std::string("ark_pc"),
                        true );
    m_preCond->setIsLinear(m_is_linear);
  } else {
    if (!procID()) {
      printf("--------\n");
      printf(" SUNDIALS Interface: ");
      printf("Could not set time integration class through ");
      printf("COGENT input file. (%s)\n", m_system->tiClass().c_str());
      printf("--------\n");
    }
    MayDay::Error("see error message!");
  }

  ARKStepSetUserData( m_ti_obj, (void*) this );
  ARKStepSStolerances(  m_ti_obj, m_ark_rtol, m_ark_atol );

  ARKStepSetMaxNumSteps( m_ti_obj, m_max_step);
  ARKStepSetStopTime( m_ti_obj, m_max_time );

  ARKStepSetPostprocessStageFn( m_ti_obj, postTimeStageSD );
  ARKStepSetPostprocessStepFn( m_ti_obj, postTimeStepSD );

  preTimeStep(m_cur_time);

  if (m_adapt_step) {
    int adapt_method = 0;
    if (m_system->tiClass() == _TI_RK_) adapt_method = 3;
    else if (m_system->tiClass() == _TI_ARK_) adapt_method = 5;
    ARKStepSetAdaptivityMethod( m_ti_obj, 
                                adapt_method,
                                SUNTRUE,
                                SUNFALSE, 
                                NULL );
  }

  m_solve_start = clock();
}

template <class AppVec, class System>
void SundialsTimeIntegrator<AppVec,System>::solve()
{
#ifdef CH_MPI
   MPI_Barrier(MPI_COMM_WORLD);
#endif
   m_solve_start = clock();

  while ( (m_cur_step<m_max_step) && (m_cur_time<m_max_time) ) {

    if (!procID()) {
      std::cout << endl << "Step " << m_cur_step+1 << endl;
    }

    ARKStepEvolve( m_ti_obj, 
                   m_cur_time+m_dt, 
                   m_Y,
                   &m_cur_time,
                   ARK_ONE_STEP );

#ifdef CH_MPI
   MPI_Barrier(MPI_COMM_WORLD);
#endif
   clock_t now = clock();
   double walltime, walltime_g;
   walltime = ((double) (now - m_solve_start)) / CLOCKS_PER_SEC;
#ifdef CH_MPI
   MPI_Allreduce( &walltime,
                  &walltime_g,
                  1,
                  MPI_DOUBLE,
                  MPI_MAX,
                  MPI_COMM_WORLD);
#else
   walltime_g = walltime;
#endif

    m_system->printDiagnostics();
    if(procID()==0) {
      printf( "  Step norm: %1.4e (rel), %1.4e (abs)\n", 
              m_step_norm_rel, m_step_norm_abs);
      printf("\n");
      printf( "Step %4d completed, simulation time is %f, ",
              (int)m_cur_step, m_cur_time );
      printf("solver wall time is %f seconds\n", walltime_g);
      printf("----\n");
    }
  }

#ifdef CH_MPI
   MPI_Barrier(MPI_COMM_WORLD);
#endif
  m_solve_end = clock();

  return;
}

template <class AppVec, class System>
void SundialsTimeIntegrator<AppVec,System>::parseParametersSimulation( ParmParse& a_pp )
{
   // This determines the amount of diagnositic output generated
   a_pp.query( "verbosity", m_verbosity );
   CH_assert( m_verbosity >= 0 );

   // Stop after this number of steps
   int tmp_max_step;
   a_pp.query( "max_step", tmp_max_step );
   CH_assert( m_max_step >= 0 );
   m_max_step = tmp_max_step;

   // Stop when the simulation time get here
   a_pp.query( "max_time", m_max_time );
   CH_assert( m_max_time >= 0.0 );

   //get time step size or CFL
   if (a_pp.contains("fixed_dt")) {
     a_pp.get( "fixed_dt", m_dt );
     CH_assert( m_dt>0.0 );
     m_cfl = -1;
   } else {
     m_cfl = 1.0;
     a_pp.query( "cfl_number", m_cfl );
     CH_assert( m_cfl>0.0 );
     m_dt = -1;
   }

   // Set up checkpointing
   a_pp.query( "checkpoint_interval", m_checkpoint_interval );
   a_pp.query( "checkpoint_prefix", m_checkpoint_prefix );

   // Set up plot file writing
   a_pp.query( "plot_interval", m_plot_interval );
   a_pp.query( "plot_prefix", m_plot_prefix );
}

template <class AppVec, class System>
void SundialsTimeIntegrator<AppVec,System>::parseParametersSUNDIALS( ParmParse& a_pp )
{
   a_pp.query( "ti_order", m_ark_order );
   a_pp.query( "ti_atol", m_ark_atol );
   a_pp.query( "ti_rtol", m_ark_rtol );
   a_pp.query( "linsolver_maxits", m_linsolver_maxits );
   a_pp.query( "linsolver_verbose", m_linsolver_verbose );
   a_pp.query( "nonlinsolver_maxits", m_nonlinsolver_maxits );
   a_pp.query( "nonlinsolver_verbose", m_nonlinsolver_verbose );
   a_pp.query( "adapt_step", m_adapt_step );
}

template <class AppVec, class System>
void SundialsTimeIntegrator<AppVec,System>::finalize()
{
  m_solve_end = clock();
  cpAppVecFromSundials( m_Y, m_U );

  if ( (m_plot_interval >= 0) && (m_last_plot!=m_cur_step) ) {
     writePlotFile();
  }

  if ( m_last_checkpoint!=m_cur_step ) {
     writeCheckpointFile();
  }

  if (!procID()) {
    std::cout << "SUNDIALS time integration completed: ";
    std::cout << "final time=" << m_cur_time << ", ";
    std::cout << "steps=" << m_cur_step << ".\n";
  }

  if (!procID()) {
    std::cout << "----\n";
  }
  m_system->printFunctionCounts();

  printTIStats();

#ifdef CH_MPI
   MPI_Barrier(MPI_COMM_WORLD);
#endif
   m_main_end = clock();
  double main_walltime = ((double) (m_main_end - m_main_start)) / CLOCKS_PER_SEC;
  double solve_walltime = ((double) (m_solve_end - m_solve_start)) / CLOCKS_PER_SEC;
  double main_walltime_g, solve_walltime_g;
#ifdef CH_MPI
  MPI_Allreduce(  &main_walltime,
                  &main_walltime_g,
                  1,
                  MPI_DOUBLE,
                  MPI_MAX,
                  MPI_COMM_WORLD );
  MPI_Allreduce(  &solve_walltime,
                  &solve_walltime_g,
                  1,
                  MPI_DOUBLE,
                  MPI_MAX,
                  MPI_COMM_WORLD );
#else
  main_walltime_g = main_walltime;
  solve_walltime_g = solve_walltime;
#endif
  if (!procID()) {
    std::cout << "Solve wall time (in seconds): " << solve_walltime_g << "\n";
    std::cout << "Total wall time (in seconds): " << main_walltime_g << "\n";
  }
}

template <class AppVec, class System>
SundialsTimeIntegrator<AppVec,System>::~SundialsTimeIntegrator()
{
  ARKStepFree(&m_ti_obj);
  SUNLinSolFree(m_linsolver);
  SUNNonlinSolFree(m_nonlinsolver);
  N_VDestroy(m_Y);
  if (m_preCond) delete m_preCond;
}

template <class AppVec, class System>
void SundialsTimeIntegrator<AppVec,System>::printTIStats()
{
  // Get integrator and solver stats
  long int  nst = 0, 
            nst_a = 0, 
            netf = 0, 
            nfe = 0, 
            nfi = 0, 
            nni = 0, 
            ncfn = 0, 
            nli = 0, 
            nlcf = 0, 
            nsetups = 0, 
            nfi_ls = 0, 
            nJv = 0;

  ARKStepGetNumSteps(m_ti_obj, &nst);
  ARKStepGetNumStepAttempts(m_ti_obj, &nst_a);
  ARKStepGetNumErrTestFails(m_ti_obj, &netf);
  ARKStepGetNumRhsEvals(m_ti_obj, &nfe, &nfi);
  if (!procID()) {
    std::cout << fixed;
    std::cout << setprecision(6);
    std::cout << "  Steps                    = " << nst     << endl;
    std::cout << "  Step attempts            = " << nst_a   << endl;
    std::cout << "  Error test fails         = " << netf    << endl;
    std::cout << "  Explicit RHS evals       = " << nfe     << endl;
    std::cout << endl;
  }
  
  if (m_system->tiClass() == _TI_ARK_) {
    ARKStepGetNumNonlinSolvIters(m_ti_obj, &nni);
    ARKStepGetNumNonlinSolvConvFails(m_ti_obj, &ncfn);
    ARKStepGetNumLinIters(m_ti_obj, &nli);
    ARKStepGetNumLinConvFails(m_ti_obj, &nlcf);
    ARKStepGetNumLinSolvSetups(m_ti_obj, &nsetups);
    ARKStepGetNumLinRhsEvals(m_ti_obj, &nfi_ls);
    ARKStepGetNumJtimesEvals(m_ti_obj, &nJv);

    if (!procID()) {
      std::cout << fixed;
      std::cout << setprecision(6);
      std::cout << "  Implicit RHS evals       = " << nfi     << endl;
      std::cout << "  Nonlinear solver iters   = " << nni     << endl;
      std::cout << "  Nonlinear solver fails   = " << ncfn    << endl;
      std::cout << "  Linear solver iters      = " << nli     << endl;
      std::cout << "  Linear solver fails      = " << nlcf    << endl;
      std::cout << "  Linear solver setups     = " << nsetups << endl;
      std::cout << "  Linear solver RHS evals  = " << nfi_ls  << endl;
      std::cout << "  Jacobian-vector products = " << nJv     << endl;
      std::cout << endl;
    }
    // Compute average nls iters per step attempt and ls iters per nls iter
    realtype avgnli = (realtype) nni / (realtype) nst_a;
    realtype avgli  = (realtype) nli / (realtype) nni;
    if (!procID()) {
      std::cout << "  Avg nonlinear solver iters per step attempt = " 
                << avgnli 
                << endl;
      std::cout << "  Avg linear solver iters per NLS iter        = " 
                << avgli  
                << endl;
      std::cout << endl;
    }
    
    // Get preconditioner stats
    long int npe, nps;
    ARKStepGetNumPrecEvals(m_ti_obj, &npe);
    ARKStepGetNumPrecSolves(m_ti_obj, &nps);
    if (!procID()) {
      std::cout << "  Preconditioner setups = " << npe << endl;
      std::cout << "  Preconditioner solves = " << nps << endl;
      std::cout << endl;
    }
  }

  return;
}

#include "NamespaceFooter.H"

#endif /*_SUNDIALSTIMEINTEGRATOR_H_ */
#endif /* with_sundials */
