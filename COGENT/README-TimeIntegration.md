# Cogent time integration

COGENT can be used with the following external time integration libraries:
- SUNDIALS
- PETSc
By default, it will use native time integrators.

## Compiling with SUNDIALS

See SUNDIALS (https://github.com/LLNL/sundials) documentation and repo 
for instructions on how to dowload and compile it. A quick set of 
instructions is available here:
https://debog.github.io/codes/sundials.html

- Define environment variable SUNDIALS_DIR pointing the installation
location of SUNDIALS (i.e., where the directories "include" and "lib"
or "lib64" are located).
- Compile COGENT as usual - it will use the SUNDIALS_DIR variable
to find SUNDIALS and compile with the SUNDIALS interface.

### Running with SUNDIALS Time Integration

In the input file, add the following line:

    ti_implementation = sundials

To use SUNDIALS' own time integrators, add the line:

    sundials.ti_order = <n>

where <n> is the desired order. Otherwise, it will
use the COGENT method specified through
gksystem.ti_class and gksystem.ti_method
(default: RK4)

For other available options (setting tolerances,
verbosity, etc.), see 
SundialsTimeIntegration::parseParametersSUNDIALS()

## Compiling with PETSc

See PETSc (https://gitlab.com/petsc/petsc) documentation
and repo for instructions on how to download and compile
it. A quick set of instructions is available here:
https://debog.github.io/codes/petsc.html

- PETSc's compilation requires setting the environment
variables PETSC_DIR and PETSC_ARCH. Make sure they are 
set to meaningful values.
- Compile COGENT as usual - it will use the PETSC_DIR 
and PETSC_ARCH variables to find PETSc and compile with
the PETSc interface.

### Running with PETSc Time Integration

In the input file, add the following line:

    ti_implementation = petsc

The PETSc interface will use the inputs related
to time integration from the input file to set
the corresponding PETSc parameters.

However, if already familiar with PETSc, options
can be provided through .petscrc or the command line
(eg., -ts_type, -snes_monitor, -ksp_monitor, etc)
and these will overwrite input file values.
