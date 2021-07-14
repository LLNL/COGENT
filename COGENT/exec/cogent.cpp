#include <string.h>
#include "ParmParse.H"

#define CH_SPACEDIM CFG_DIM+2

#include "ODEVector.H"
#include "AppCtxtLibrary.H"
#include "SimulationLibrary.H"

#include "parstream.H"
#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#include "UsingNamespace.H"

inline int checkCommandLineArgs( int a_argc, char* a_argv[] )
{
   // Check for an input file
   if (a_argc<=1) {
      pout() << "Usage:  cogent...ex <inputfile>" << endl;
      pout() << "No input file specified" << endl;
      return -1;
   }
   return 0;
}

#ifdef with_petsc
static const char help[] = "COGENT";
#endif

typedef enum {native, petsc, sundials} TimeIntegrationImplementation;

TimeIntegrationImplementation parseParams( ParmParse &a_pp)
{
  /* default: use native time integrator */
  TimeIntegrationImplementation ti_impl = native;

  std::string ti_impl_name = "native";
  a_pp.query("ti_implementation", ti_impl_name);

  if      (ti_impl_name == "native")    ti_impl = native;
  else if (ti_impl_name == "petsc")     ti_impl = petsc;
  else if (ti_impl_name == "sundials")  ti_impl = sundials;
  else {
    MayDay::Error("Invalid value for ti_implementation");
  }

  return ti_impl;
}

int main( int a_argc, char* a_argv[] )
{
#ifdef CH_MPI
   // Start MPI
   MPI_Init( &a_argc, &a_argv );
   setChomboMPIErrorHandler();
#endif

#ifdef with_petsc
  PetscInitialize(&a_argc,&a_argv,(char*)0,help);
#endif

   int status = checkCommandLineArgs( a_argc, a_argv );

   if (status==0) {

      ParmParse pp( a_argc-2, a_argv+2, NULL, a_argv[1] );
      auto ti_impl = parseParams(pp);

      AppCtxt* system = new GKSystem(pp);

      Simulation<ODEVector,AppCtxt>* simulation(NULL); 
      if (ti_impl == petsc) {
#ifdef with_petsc
        simulation = new PetscTimeIntegrator<ODEVector,AppCtxt>;
#else
        MayDay::Error("Not compiled with PETSc");
#endif
      } else if (ti_impl == sundials) {
#ifdef with_sundials
        simulation = new SundialsTimeIntegrator<ODEVector,AppCtxt>;
#else
        MayDay::Error("Not compiled with SUNDIALS");
#endif
      } else {
        simulation = new COGENTTimeIntegration<ODEVector,AppCtxt>;
      }

      simulation->initialize( pp, system );
      simulation->solve();
      simulation->finalize();

      delete simulation;
      delete system;
   }

#ifdef with_petsc
  PetscFinalize();
#endif

#ifdef CH_MPI
   MPI_Finalize();
#endif

   return status;
}
