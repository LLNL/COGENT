#ifdef with_petsc

#ifndef _PETSCTIMEINTEGRATOR_H_
#define _PETSCTIMEINTEGRATOR_H_

#include <string>
#include <time.h>
#include "Simulation.H"

#include "ODEVector.H"
#include "AppCtxt.H"
#include "MultiPhysicsPreconditioner.H"

/* include PETSc header files */
#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscsnes.h>
#include <petscts.h>

#include "NamespaceHeader.H"

template <class AppVec, class System>
class PetscTimeIntegrator : public Simulation<AppVec,System>
{
  public: 
    PetscTimeIntegrator<AppVec,System>();
    ~PetscTimeIntegrator<AppVec,System>();

    void initialize(ParmParse&, System* const);
    void finalize();
    void solve() { TSSolve(m_ts, NULL); }

    void  cpAppVecToPetsc    (Vec&, const AppVec&);
    void  cpAppVecFromPetsc  (const Vec&, AppVec&);

    inline System*    getSystem()                 { return m_system; }
    inline int        getVerbosity()              { return m_verbosity; }
    inline int        getPlotInterval()           { return m_plot_interval; }
    inline int        getChkptInterval()          { return m_checkpoint_interval; }
    inline Real       getShift()                  { return m_shift; }
    inline Real       getTime()                   { return m_cur_time; }
    inline Real       getOldTime()                { return m_old_time; }
    inline Real       getTimeStep()               { return m_dt; }
    inline Real       getEpsilonJFNK()            { return m_epsilon_jfnk; }
    inline bool       isLinear()                  { return m_is_linear; }

    inline void    setLastPlot  (int a_int) { m_last_plot = a_int; }      
    inline void    setLastChkpt (int a_int) { m_last_checkpoint = a_int; }      

    inline void setTime           (Real a_t)    { m_cur_time = a_t; }
    inline void setOldTime        (Real a_t)    { m_old_time = a_t; }
    inline void setTimeStep       (Real a_dt)   { m_dt = a_dt; }
    inline void setCFL            (Real a_cfl)  { m_cfl = a_cfl; }
    inline void setTimeStepNumber (int a_iter)  { m_cur_step = a_iter; }
    inline void setShift          (Real a_a)    { m_shift = a_a; }

    inline Real getCFL()  { return m_cfl; }

    void writeHistFile(bool);
    void writePlotFile();
    void writeCheckpointFile();

    AppVec m_U, m_F, m_U0, m_F0;
    Preconditioner<AppVec,System>* m_preCond;

  private:

    void parseParametersSimulation( ParmParse& );

    System* m_system;
    int     m_cur_step;
    int     m_newton_its;
    int     m_gmres_its;
    int     m_max_step;
    Real    m_cur_time, m_old_time;
    Real    m_max_time;
    Real    m_dt;
    Real    m_cfl;
    int     m_verbosity;
    int     m_plottedonce;

    int     m_checkpoint_interval;
    int     m_last_checkpoint;

    int     m_plot_interval;
    int     m_last_plot;

    std::string m_checkpoint_prefix;
    std::string m_plot_prefix;

    int   m_local_size;
    Real  m_shift,
          m_epsilon_jfnk;

    bool  m_is_linear,
          m_is_Jac_allocated;

    Vec           m_Y;
    TS            m_ts;
    Mat           m_Jac, m_Pre;

    /* Walltime measurements */
    clock_t m_main_start, m_solve_start, m_main_end, m_solve_end;
};

static PetscErrorCode RHSFunction(  TS        ts, 
                                    PetscReal t, 
                                    Vec       Y, 
                                    Vec       F, 
                                    void*     ctxt  )
{
  PetscTimeIntegrator<ODEVector,AppCtxt> *context 
      = (PetscTimeIntegrator<ODEVector,AppCtxt>*) ctxt;
  AppCtxt* system = context->getSystem();

  context->cpAppVecFromPetsc( Y, context->m_U );
  system->preOpEval( context->m_U, (Real)t, chkpt_naive );
  system->explicitOp( context->m_F, (Real)t, context->m_U );
  context->cpAppVecToPetsc( F, context->m_F);

  return 0;
}

static PetscErrorCode RHSFunctionIMEX(  TS        ts, 
                                        PetscReal t, 
                                        Vec       Y, 
                                        Vec       F, 
                                        void*     ctxt )
{
  PetscTimeIntegrator<ODEVector,AppCtxt> *context 
      = (PetscTimeIntegrator<ODEVector,AppCtxt>*) ctxt;
  AppCtxt* system = context->getSystem();

  context->cpAppVecFromPetsc( Y, context->m_U );
  system->preOpEval( context->m_U, (Real)t, chkpt_naive );
  system->explicitOpImEx( context->m_F, (Real)t, context->m_U );
  context->cpAppVecToPetsc( F, context->m_F);

  return 0;
}

static PetscErrorCode IFunctionIMEX(  TS        ts, 
                                      PetscReal t, 
                                      Vec       Y, 
                                      Vec       Ydot, 
                                      Vec       F, 
                                      void*     ctxt  )
{
  PetscTimeIntegrator<ODEVector,AppCtxt> *context 
      = (PetscTimeIntegrator<ODEVector,AppCtxt>*) ctxt;
  AppCtxt* system = context->getSystem();

  context->cpAppVecFromPetsc( Y, context->m_U );
  system->preOpEval( context->m_U, (Real)t, chkpt_naive );
  system->implicitOpImEx( context->m_F, (Real)t, context->m_U );
  context->cpAppVecToPetsc( F, context->m_F);

  context->m_U0 = context->m_U;
  context->m_F0 = context->m_F;
  context->m_preCond->setBaseSolution(context->m_U);
  context->m_preCond->setBaseRHS(context->m_F, context->m_U);

  VecAYPX(F,-1.0,Ydot);

  return 0;
}

static PetscErrorCode IJacobianIMEX(  TS        ts, 
                                      PetscReal t, 
                                      Vec       Y, 
                                      Vec       Ydot, 
                                      PetscReal a,
                                      Mat       A, 
                                      Mat       B, 
                                      void*     ctxt  )
{
  PetscTimeIntegrator<ODEVector,AppCtxt> *context 
      = (PetscTimeIntegrator<ODEVector,AppCtxt>*) ctxt;

  int cur_step, stage = 0;
  TSGetStepNumber( ts, &cur_step );

  context->setShift((Real)a);
  context->setTime((Real)t);

  context->cpAppVecFromPetsc( Y, context->m_U );
  context->m_preCond->update( context->m_U,
                              (Real) t,
                              cur_step+1,
                              stage,
                              (Real) a );

  return 0;
}

static PetscErrorCode JacobianFunctionIMEX( Mat A, 
                                            Vec Y, 
                                            Vec F )
{
  PetscTimeIntegrator<ODEVector,AppCtxt> *context;
  MatShellGetContext(A,&context);
  AppCtxt* system = context->getSystem();

  PetscReal normY;
  VecNorm( Y, NORM_2, &normY);
  Real normY0 = context->m_U0.computeNorm();

  if (normY < 1e-15) {

    VecZeroEntries(F);
    VecAXPBY( F, context->getShift(), 0, Y );

  } else {

    Real epsilon;
    if (context->isLinear()) {
      epsilon = 1.0;
    } else {
      Real epsJFNK = context->getEpsilonJFNK();
      epsilon = epsJFNK * sqrt(1 + normY0) / normY;
    }

    context->cpAppVecFromPetsc( Y, context->m_U );
    context->m_U *= epsilon;
    context->m_U += context->m_U0;

    system->preOpEval(  context->m_U, 
                        context->getTime(), 
                        chkpt_naive );
    system->implicitOpImEx( context->m_F, 
                            context->getTime(), 
                            context->m_U );
    context->m_F -= context->m_F0;

    context->cpAppVecToPetsc( F, context->m_F );
    VecAXPBY( F, 
              (PetscReal) context->getShift(), 
              (-1.0/epsilon), 
              Y );
  }

  return 0;
}

static PetscErrorCode PCFunctionIMEX( PC  pc,
                                      Vec X,
                                      Vec Y )
{
  void* ctxt;
  PCShellGetContext(pc, &ctxt);
  PetscTimeIntegrator<ODEVector,AppCtxt>* context
   = (PetscTimeIntegrator<ODEVector,AppCtxt>*) ctxt;

  context->cpAppVecFromPetsc( X, context->m_U );
  context->m_F.zero();
  context->m_preCond->apply(context->m_F, context->m_U);
  context->cpAppVecToPetsc( Y, context->m_F );

  return 0;
}

static PetscErrorCode PreTimeStep(TS ts)
{
  PetscTimeIntegrator<ODEVector,AppCtxt> *context;
  TSGetApplicationContext( ts, &context );
  AppCtxt* system = context->getSystem();

  int cur_step;
  Real cur_time, dt;

  TSGetTime( ts, &cur_time );
  TSGetTimeStep( ts, &dt );
  TSGetStepNumber( ts, &cur_step );

  TSType ts_type; 
  TSGetType( ts, &ts_type );

  if (context->getVerbosity() >= 1) {
    if (procID()==0) cout << endl << "Step " << cur_step+1 << endl;
  }

  Vec Y; 
  TSGetSolution( ts, &Y );

  context->cpAppVecFromPetsc( Y, context->m_U );
  system->preTimeStep(context->m_U, cur_step, cur_time);

  TSAdapt ts_adapt;
  TSGetAdapt(ts, &ts_adapt);
  TSAdaptType ts_adapt_type;
  TSAdaptGetType(ts_adapt, &ts_adapt_type);

  if (!strcmp(ts_adapt_type, TSADAPTNONE)) {

    if (context->getCFL() > 0) {
  
      Real dt_stable = 0.0;
      if (!strcmp(ts_type, TSARKIMEX)) {
        dt_stable = system->stableDt( cur_step, ti_imex );
      } else if (!strcmp(ts_type, TSRK)) {
        dt_stable = system->stableDt( cur_step, ti_explicit );
      } else {
        MayDay::Error("PETSc interface: ts_type not supported!");
      }
      CH_assert( dt_stable > 1.0e-16 );
  
      TSSetTimeStep( ts, dt_stable*context->getCFL() );
  
    } else if (dt <= 0.0) {
  
      MayDay::Error("PETSc interface: non-positive dt!");
  
    }

  }

  context->setOldTime( cur_time );

  return 0;
}

static PetscErrorCode PostTimeStep(TS ts)
{
  PetscTimeIntegrator<ODEVector,AppCtxt> *context;
  TSGetApplicationContext( ts, &context );
  AppCtxt* system = context->getSystem();

  int cur_step, 
      verbosity = context->getVerbosity(),
      plot_interval = context->getPlotInterval();
  Real cur_time, dt;

  TSGetTime (ts, &cur_time );
  dt = cur_time - context->getOldTime();
  //TSGetTimeStep   (ts, &dt );
  TSGetStepNumber (ts, &cur_step );

  context->setTime( cur_time );
  context->setTimeStepNumber( cur_step );
  context->setTimeStep( dt );

  Vec Y; 
  TSGetSolution( ts, &Y );

  context->cpAppVecFromPetsc( Y, context->m_U );
  system->postTimeStep( context->m_U, 
                        cur_step, 
                        dt, 
                        cur_time );

  if (verbosity >= 1) {
    system->printDiagnostics();
    if(procID()==0) {
      cout  << "Step " 
            << cur_step 
            << " completed, simulation time is "
            << cur_time 
            << endl 
            << "----" 
            << endl;
    }
  }

  if ( (cur_step % plot_interval)==0 ) {
    context->writePlotFile();
    context->setLastPlot(cur_step);
  }

  context->writeHistFile(false);

  if ( (cur_step % context->getChkptInterval())==0 ) {
    context->writeCheckpointFile();
    context->setLastChkpt(cur_step);
  }

  return 0;
}

static PetscErrorCode PostTimeStage(  TS        ts, 
                                      PetscReal stagetime, 
                                      PetscInt  stageindex, 
                                      Vec* Y )
{
  PetscTimeIntegrator<ODEVector,AppCtxt> *context;
  TSGetApplicationContext( ts, &context );
  AppCtxt* system = context->getSystem();

  int cur_step;
  Real cur_time, dt;

  TSGetTime       ( ts, &cur_time );
  TSGetTimeStep   ( ts, &dt );
  TSGetStepNumber ( ts, &cur_step );

  static int step_number  = 0;
  static int stage_number = -1;

  if (cur_step == step_number+1) {
    step_number++;
    stage_number = -1;
  } else if (cur_step == step_number) {
    stage_number++;
  } else {
    if (!procID()) {
      std::cout << "Error in PostTimeStage() in src/driver/PETScTimeIntegration.H: cur_step=" 
                << cur_step 
                << ", step_number=" 
                << step_number << "."
                << " Something has caused stage-counting to go haywire.\n";
    }
  }

  context->cpAppVecFromPetsc( Y[stageindex], context->m_U );

  system->postTimeStage(  cur_step,
                          stagetime,
                          dt,
                          context->m_U,
                          stage_number );

  return 0;
}

template <class AppVec, class System>
inline void PetscTimeIntegrator<AppVec,System>::writeHistFile(bool flag)
{
  m_system->writeFieldHistory(m_cur_step, m_cur_time,flag);
}

template <class AppVec, class System> 
inline void PetscTimeIntegrator<AppVec,System>::writeCheckpointFile()
{
   m_system->writeCheckpointFile( m_checkpoint_prefix, 
                                  m_cur_step, 
                                  m_cur_time, 
                                  m_dt );
}

template <class AppVec, class System> 
inline void PetscTimeIntegrator<AppVec,System>::writePlotFile()
{
#ifdef CH_USE_HDF5
   ParmParse ppsim( "simulation" );
   parseParametersSimulation( ppsim );

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
PetscTimeIntegrator<AppVec,System>::PetscTimeIntegrator( )
  : m_preCond(NULL),
    m_system(NULL),
    m_cur_step(0),
    m_max_step(0),
    m_cur_time(0.0),
    m_old_time(0.0),
    m_max_time(0.0),
    m_dt(0.0),
    m_cfl(0.0),
    m_verbosity(0),
    m_plottedonce(0),
    m_checkpoint_interval(0),
    m_last_checkpoint(0),
    m_plot_interval(0),
    m_last_plot(0),
    m_checkpoint_prefix( "chk" ),
    m_plot_prefix( "plt" ),
    m_local_size(0),
    m_shift(1.0),
    m_epsilon_jfnk(1e-6),
    m_is_linear("true"),
    m_is_Jac_allocated(false)
{
}

template <class AppVec, class System>
void PetscTimeIntegrator<AppVec,System>::initialize(ParmParse&    a_pp,
                                                    System* const a_sys ) 
{
  m_system = a_sys;

  m_main_start = clock();
  ParmParse ppsim( "simulation" );
  parseParametersSimulation( ppsim );

  if ( ppsim.contains( "restart_file" ) ) {
     //loadRestartFile( ppsim );
  }
  m_system->initialize(m_cur_step, m_cur_time);
  m_old_time = m_cur_time;
  m_system->printDiagnostics();

  m_system->defineODEVector( m_U );
  m_U.printLoadBalanceInfo();

  m_system->copyStateToArray( m_U );
  m_F.define( m_U );
  m_U0.define( m_U );
  m_F0.define( m_U );

  if ( m_plot_interval >= 0 ) {
    writePlotFile();
    m_last_plot = m_cur_step;
  }
  writeHistFile(true);

  m_is_linear = m_system->isLinear();
  
  /* create and initialize the PETSc solution vector */
  m_local_size = m_system->getVectorSize();
  VecCreate( MPI_COMM_WORLD, &m_Y );
  VecSetSizes( m_Y, m_local_size, PETSC_DECIDE );
  VecSetUp( m_Y );

  /* copy initial solution to PETSc vector */
  cpAppVecToPetsc( m_Y, m_U );

  /* create and initialize the PETSc time integration object */
  TSCreate( MPI_COMM_WORLD, &m_ts );
  TSSetMaxSteps( m_ts, m_max_step );
  TSSetTimeStep( m_ts, m_dt );
  TSSetMaxTime( m_ts, m_max_time );
  TSSetTime( m_ts, 0.0 );
  TSSetExactFinalTime( m_ts, TS_EXACTFINALTIME_MATCHSTEP );

  if (m_system->tiClass() == _TI_RK_) {
    TSSetType( m_ts, TSRK );
    if (m_system->tiMethod() == "1fe") {
      TSRKSetType( m_ts, TSRK1FE );
    } else  if (m_system->tiMethod() == "2") {
      TSRKSetType( m_ts, TSRK2A );
    } else  if (m_system->tiMethod() == "2a") {
      TSRKSetType( m_ts, TSRK2A );
    } else  if (m_system->tiMethod() == "3") {
      TSRKSetType( m_ts, TSRK3 );
    } else  if (m_system->tiMethod() == "4") {
      TSRKSetType( m_ts, TSRK4 );
    }
  } else if (m_system->tiClass() == _TI_ARK_) {
    TSSetType( m_ts, TSARKIMEX );
    if (m_system->tiMethod() == "1bee") {
      TSARKIMEXSetType( m_ts, TSARKIMEX1BEE );
    } else if (m_system->tiMethod() == "2a") {
      TSARKIMEXSetType( m_ts, TSARKIMEXA2 );
    } else if (m_system->tiMethod() == "2c") {
      TSARKIMEXSetType( m_ts, TSARKIMEX2C );
    } else if (m_system->tiMethod() == "3") {
      TSARKIMEXSetType( m_ts, TSARKIMEX3 );
    } else if (m_system->tiMethod() == "4") {
      TSARKIMEXSetType( m_ts, TSARKIMEX4 );
    } else if (m_system->tiMethod() == "5") {
      TSARKIMEXSetType( m_ts, TSARKIMEX5 );
    }
  } else {
    if (!procID()) {
      printf("--------\n");
      printf(" PETSc Interface: ");
      printf("Could not set time integration method through ");
      printf("COGENT input file.\n");
      printf("Please set it through PETSc inputs (-ts_type, etc.)\n");
      printf("Default: RK4\n");
      printf("--------\n");
    }
    TSSetType( m_ts, TSRK );
    TSRKSetType( m_ts, TSRK4 );
  }

  TSAdapt adapt; 
  TSGetAdapt( m_ts,&adapt );
  TSAdaptSetType( adapt, TSADAPTNONE );

  TSSetSolution( m_ts, m_Y );
  TSSetApplicationContext( m_ts, this );
  if (m_is_linear) {
    TSSetProblemType( m_ts, TS_LINEAR );
  }

  TSSetFromOptions( m_ts );

  PetscOptionsGetReal(NULL,NULL,"-jfnk_epsilon",&m_epsilon_jfnk,NULL);
  
  TSType ts_type; 
  TSGetType( m_ts, &ts_type );

  if (!strcmp(ts_type,TSARKIMEX)) {

    TSSetRHSFunction( m_ts, PETSC_NULL, RHSFunctionIMEX, this);
    TSSetIFunction  ( m_ts, PETSC_NULL, IFunctionIMEX  , this);

    MatCreateShell( MPI_COMM_WORLD,
                    m_local_size,
                    m_local_size,
                    PETSC_DETERMINE,
                    PETSC_DETERMINE,
                    this,
                    &m_Jac );

    MatShellSetOperation( m_Jac,
                          MATOP_MULT,
                          (void (*)(void))JacobianFunctionIMEX );

    MatSetUp(m_Jac);
    m_is_Jac_allocated = true;

    SNES  snes; 
    TSGetSNES( m_ts, &snes );
    if (m_is_linear) {
      SNESSetType( snes, SNESKSPONLY );
    }

    /* Set the IJacobian function for TS */
    TSSetIJacobian( m_ts, m_Jac, m_Jac, IJacobianIMEX, this );

    /* Set PC (preconditioner) to PCSHELL */
    KSP ksp;  
    SNESGetKSP( snes, &ksp );
    PC pc;   
    KSPGetPC( ksp, &pc );

    PCSetType( pc, PCSHELL );
    PCShellSetApply(pc, PCFunctionIMEX);
    PCShellSetContext(pc, this);

    m_preCond = new MultiPhysicsPreconditioner<AppVec,System>;
    m_preCond->define(  m_U,
                        (void*) m_system,
                        std::string("ARK.PC"),
                        std::string("ark_pc"),
                        true );
    m_preCond->setJFNKEps(m_epsilon_jfnk);
    m_preCond->setIsLinear(m_is_linear);

  } else {

    TSSetRHSFunction( m_ts, PETSC_NULL, RHSFunction, this );

  }

  /* set pre- and post-step functions, and post-stage function */
  TSSetPreStep( m_ts, PreTimeStep );
  TSSetPostStep( m_ts, PostTimeStep );
  TSSetPostStage( m_ts, PostTimeStage );

  TSSetUp( m_ts );

  m_solve_start = clock();
}

template <class AppVec, class System>
void PetscTimeIntegrator<AppVec,System>::cpAppVecToPetsc( Vec&          a_Y, 
                                                          const AppVec& a_U )
{
  PetscErrorCode  ierr;
  PetscScalar     *Yarr;

  ierr = VecGetArray(a_Y,&Yarr);
  a_U.copyTo( (Real* const) Yarr );
  ierr = VecRestoreArray(a_Y,&Yarr);
}

template <class AppVec, class System>
void PetscTimeIntegrator<AppVec,System>::cpAppVecFromPetsc( const Vec&  a_Y, 
                                                            AppVec&     a_U )
{
  PetscErrorCode    ierr;
  const PetscScalar *Yarr;

  ierr = VecGetArrayRead(a_Y,&Yarr);
  a_U.copyFrom( (const Real* const) Yarr );
  ierr = VecRestoreArrayRead(a_Y,&Yarr);
}


template <class AppVec, class System>
void PetscTimeIntegrator<AppVec,System>::parseParametersSimulation( ParmParse& a_ppsim )
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

   //get time step size or CFL
   if (a_ppsim.contains("fixed_dt")) {
     a_ppsim.get( "fixed_dt", m_dt );
     CH_assert( m_dt>0.0 );
     m_cfl = -1;
   } else {
     m_cfl = 1.0;
     a_ppsim.query( "cfl_number", m_cfl );
     CH_assert( m_cfl>0.0 );
     m_dt = -1;
   }

   // Set up checkpointing
   a_ppsim.query( "checkpoint_interval", m_checkpoint_interval );
   a_ppsim.query( "checkpoint_prefix", m_checkpoint_prefix );

   // Set up plot file writing
   a_ppsim.query( "plot_interval", m_plot_interval );
   a_ppsim.query( "plot_prefix", m_plot_prefix );
}

template <class AppVec, class System>
void PetscTimeIntegrator<AppVec,System>::finalize()
{
  m_solve_end = clock();
  /* Get the final solution and copy from PETSc vector to system */
  TSGetSolution( m_ts, &m_Y );
  cpAppVecFromPetsc( m_Y, m_U );

  TSGetTime( m_ts, &m_cur_time );
  TSGetTimeStep( m_ts, &m_dt );
  TSGetStepNumber( m_ts, &m_cur_step );

  TSGetSNESIterations( m_ts, &m_newton_its );
  TSGetKSPIterations( m_ts, &m_gmres_its );

  if ( (m_plot_interval >= 0) && (m_last_plot!=m_cur_step) ) {
     writePlotFile();
  }

  if ( m_last_checkpoint!=m_cur_step ) {
     writeCheckpointFile();
  }

  m_main_end = clock();

  if (!procID()) {
    cout << "PETSc time integration completed: ";
    cout << "final time=" << m_cur_time << ", ";
    cout << "steps=" << m_cur_step << ".\n";
  }

  if (!procID()) {
    cout << "----\n";
  }
  m_system->printFunctionCounts();
  if (!procID()) {
    cout << "----\n";
    cout << "  Time integrator counts:-\n";
    cout << "    Time steps          : " << m_cur_step << "\n";
    cout << "    Nonlinear iterations: " << m_newton_its << "\n";
    cout << "    Linear iterations   : " << m_gmres_its << "\n";
    cout << "----\n";
  }

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
}

template <class AppVec, class System>
PetscTimeIntegrator<AppVec,System>::~PetscTimeIntegrator()
{
  if (m_is_Jac_allocated) MatDestroy(&m_Jac);
  TSDestroy (&m_ts);
  VecDestroy(&m_Y);
  if (m_preCond) delete m_preCond;
}

#include "NamespaceFooter.H"

#endif /*_PETSCTIMEINTEGRATOR_H_ */
#endif /* with_petsc */
