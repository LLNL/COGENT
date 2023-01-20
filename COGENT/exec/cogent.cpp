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

typedef enum {  single_domain, 
                multidomain, 
                sparse_grids } SimulationType;

typedef enum {  native, 
                petsc, 
                sundials  } TimeIntegrationImplementation;

TimeIntegrationImplementation parseParamsTI( ParmParse& a_pp)
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

SimulationType parseParamsSimType( ParmParse& a_pp)
{
  /* default: single system */
  SimulationType sim_type = single_domain;

  std::string sim_type_name = "single_domain";
  a_pp.query("simulation_type", sim_type_name);

  if      (sim_type_name == "single_domain")  sim_type = single_domain;
  else if (sim_type_name == "multidomain")    sim_type = multidomain;
  else if (sim_type_name == "sparse_grids")   sim_type = sparse_grids;
  else {
    MayDay::Error("Invalid value for simulation_type");
  }

  return sim_type;
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
      auto ti_impl = parseParamsTI(pp);
      auto sim_type = parseParamsSimType(pp);

      AppCtxt* system(NULL);
      if (sim_type == single_domain) {
        system = new GKSystem(pp);
      } else if (sim_type == multidomain) {
        system = new MultidomainSimulation(pp);
      } else if (sim_type == sparse_grids) {
        system = new SparseGridsSimulation(pp);
      } else {
        MayDay::Error("Invalid value for sim_type");
      }

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
