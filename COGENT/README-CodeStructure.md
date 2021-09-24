# Cogent code structure

The `main` entry-point function is in
[exec/cogent.cpp](exec/cogent.cpp#52). It parses the command-line
arguments using Chombo's
[ParmParse](../Chombo/lib/src/BaseTools/ParmParse.H), creates a
[GKSystem](src/core/GKSystem.H) to represent the simulation state, and
a [Simulation](src/driver/Simulation.H) implementation to perform the
time integration. Currently the `Simulation` implementations available
are [PetscTimeIntegrator](src/driver/PETScTimeIntegration.H),
[SundialsTimeIntegrator](src/driver/SUNDIALSTimeIntegration.H) or
[COGENTTimeIntegration](src/driver/COGENTTimeIntegration.H).  The
`Simulation` is initialized by passing a pointer to the `GKSystem` to
`Simulation::initialize()`, and the simulation is then run by calling
`Simulation::solve()`.

## Time integrators

The time integrators expect a `System` which implements the
[AppCtxt](src/driver/AppCtxt.H) interface.

The [COGENTTimeIntegration](src/driver/COGENTTimeIntegration.H) class
is a wrapper around an integrator [TiRK](src/time/TiRK.H) or
[TiARK](src/time/TiARK.H). 


