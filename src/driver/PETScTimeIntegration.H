#ifdef with_petsc

#ifndef _PETSCTIMEINTEGRATOR_H_
#define _PETSCTIMEINTEGRATOR_H_

#include <string>
#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "REAL.H"

/* include PETSc header files */
#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscsnes.h>
#include <petscts.h>

#include "NamespaceHeader.H"

template <class System>
class PetscTimeIntegrator
{
  public: 
    PetscTimeIntegrator<System>(ParmParse&,bool);
    ~PetscTimeIntegrator<System>();
    void finalize();

  private:

    void writeHistFile(bool);
    void writePlotFile();
    void writeCheckpointFile();
    void parseParametersSimulation( ParmParse& );
    void parseParametersPetsc     ( ParmParse& );

    System      *m_system;
    int         m_cur_step;
    int         m_max_step;
    Real        m_cur_time;
    Real        m_max_time;
    Real        m_dt;
    int         m_verbosity;
    int         m_plottedonce;

    int         m_checkpoint_interval;
    int         m_last_checkpoint;
    std::string m_checkpoint_prefix;

    int         m_plot_interval;
    int         m_last_plot;
    std::string m_plot_prefix;

    int         m_local_size;

    std::string m_ts_type, m_ts_subtype;
    
    Vec   m_Y;
    TS    m_ts;

    void  copyStateToPetsc    (Vec& a_Y, System* a_System);
    void  copyStateFromPetsc  (Vec& a_Y, System* a_System);

};

template <class System>
inline void PetscTimeIntegrator<System>::writeHistFile(bool flag)
{
  m_system->writeFieldHistory(m_cur_step, m_cur_time,flag);
}

template <class System> 
inline void PetscTimeIntegrator<System>::writeCheckpointFile()
{
#ifdef CH_USE_HDF5
   char iter_str[100];

   sprintf( iter_str, "%s%04d.%dd.hdf5",
            m_checkpoint_prefix.c_str(), m_cur_step, SpaceDim );

   HDF5Handle handle( iter_str, HDF5Handle::CREATE );
   m_system->writeCheckpointFile( handle, m_cur_step, m_cur_time, m_dt );
   handle.close();
#else
   MayDay::Error( "restart only defined with hdf5" );
#endif
}

template <class System> 
inline void PetscTimeIntegrator<System>::writePlotFile()
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

template <class System>
PetscTimeIntegrator<System>::PetscTimeIntegrator( ParmParse &a_pp, bool usePetsc)
  : m_cur_step(0),
    m_max_step(0),
    m_cur_time(0.0),
    m_max_time(0.0),
    m_dt(-1.0),
    m_verbosity(0),
    m_checkpoint_interval(0),
    m_last_checkpoint(0),
    m_checkpoint_prefix( "chk" ),
    m_plot_interval(0),
    m_last_plot(0),
    m_plot_prefix( "plt" ),
    m_plottedonce(0),
    m_system(NULL),
    m_local_size(0),
    m_ts_type("rk"),
    m_ts_subtype("4")
{
  ParmParse ppsim( "simulation" );
  parseParametersSimulation( ppsim );
  ParmParse pppet( "petsc" );
  parseParametersPetsc( pppet );

  m_system = new System(a_pp, usePetsc);
  m_system->initialize(m_cur_step);
  m_system->printDiagnostics();
  if ( m_plot_interval >= 0 ) {
    writePlotFile();
    m_last_plot = m_cur_step;
  }
  writeHistFile(true);
  
  if (!procID()) {
    cout << "Using PETSc time integrator: " << m_ts_type;
    cout << " (" << m_ts_subtype << ")\n";
  }

  /* create and initialize the PETSc solution vector */
  m_local_size = m_system->getVectorSize();
  VecCreate(MPI_COMM_WORLD,&m_Y);
  VecSetSizes(m_Y,m_local_size,PETSC_DECIDE);
  VecSetUp(m_Y);

  /* copy initial solution to PETSc vector */
  copyStateToPetsc(m_Y,m_system);

  /* create and initialize the PETSc time integration object */
  TSCreate(MPI_COMM_WORLD,&m_ts);
  TSSetDuration(m_ts,m_max_step,m_max_time);
  TSSetInitialTimeStep(m_ts,0.0,m_dt);
  TSSetExactFinalTime(m_ts,TS_EXACTFINALTIME_MATCHSTEP);
  TSSetSolution(m_ts,m_Y);
  TSSetApplicationContext(m_ts,m_system);

  TSSetType(m_ts,m_ts_type.c_str());
  if (!strcmp(m_ts_type.c_str(),TSRK)) {
    TSRKSetType(m_ts,m_ts_subtype.c_str());
    //TSSetRHSFunction(m_ts,PETSC_NULL,RHSFunction,NULL);
  } else {
    if (!procID()) {
      cout << "ts_type (" << m_ts_type << ") is invalid or ";
      cout << "not yet supported.\n";
    }
  }
  TSSetUp(m_ts);
}

template <class System>
void PetscTimeIntegrator<System>::copyStateToPetsc( Vec& a_Y, System* a_System )
{
  PetscErrorCode  ierr;
  PetscScalar     *Yarr;

  ierr = VecGetArray(a_Y,&Yarr);
  a_System->copyStateToArray((Real*)Yarr);
  ierr = VecRestoreArray(a_Y,&Yarr);
}

template <class System>
void PetscTimeIntegrator<System>::copyStateFromPetsc( Vec& a_Y, System* a_System )
{
  PetscErrorCode  ierr;
  const PetscScalar     *Yarr;

  ierr = VecGetArrayRead(a_Y,&Yarr);
  a_System->copyStateFromArray((Real*)Yarr);
  ierr = VecRestoreArrayRead(a_Y,&Yarr);
}

template <class System>
void PetscTimeIntegrator<System>::parseParametersSimulation( ParmParse& a_ppsim )
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

   a_ppsim.query( "fixed_dt", m_dt );
   CH_assert( m_dt>0.0 );

   // Set up checkpointing
   a_ppsim.query( "checkpoint_interval", m_checkpoint_interval );
   a_ppsim.query( "checkpoint_prefix", m_checkpoint_prefix );

   // Set up plot file writing
   a_ppsim.query( "plot_interval", m_plot_interval );
   a_ppsim.query( "plot_prefix", m_plot_prefix );
}

template <class System>
void PetscTimeIntegrator<System>::parseParametersPetsc( ParmParse& a_pppet )
{
   // This determines the amount of diagnositic output generated
   a_pppet.query( "ts_type", m_ts_type );
   a_pppet.query( "ts_subtype", m_ts_subtype );
}

template <class System>
void PetscTimeIntegrator<System>::finalize()
{
  /* copy final solution from PETSc vector to system */
  copyStateFromPetsc(m_Y,m_system);

  TSGetTime(m_ts,&m_cur_time);
  TSGetTimeStep(m_ts,&m_dt);
  TSGetTimeStepNumber(m_ts,&m_cur_step);

  if ( (m_plot_interval >= 0) && (m_last_plot!=m_cur_step) ) {
     writePlotFile();
  }

  if ( m_last_checkpoint!=m_cur_step ) {
     writeCheckpointFile();
  }

  if (!procID()) {
    cout << "PETSc time integration completed: ";
    cout << "final time=" << m_cur_time << ", ";
    cout << "steps=" << m_cur_step <<".\n";
  }
}

template <class System>
PetscTimeIntegrator<System>::~PetscTimeIntegrator()
{
  TSDestroy(&m_ts);
  VecDestroy(&m_Y);
}

#include "NamespaceFooter.H"

#endif /*_PETSCTIMEINTEGRATOR_H_ */
#endif /* with_petsc */
