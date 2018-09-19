COGENT is a tokamak edge plasma gyro-kinetic simulation code. 
This repository contains the COGENT code (cogent/) as well as
Chombo (Chombo/), the mapped, multiblock PDE library from 
Berkeley Lab on which COGENT is built.

- The copy of Chombo in this repository is only for easily
  installing COGENT. It is developed and maintained by 
  LBL-ANAG, and the original copy and documentation can be
  downloaded from: 
  https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations

Install Instructions
--------------------

It is assumed the C++ and FORTRAN compilers and their MPI wrappers are already
available in the environment.

- Step 1: Make sure the following environment variables are defined and point
  to the location of the serial and parallel HDF5 libraries compiled with
  the **same** compilers that you intend to compile COGENT with.
  - **HDF5_DIR_SERIAL** pointing to the serial HDF5 library
  - **HDF5_DIR_PARALLEL** pointing to the parallel HDF5 library

    In bash, this can be done by:
  
        export HDF5_DIR_SERIAL=/path/to/hdf5-serial
        export HDF5_DIR_PARALLEL=/path/to/hdf5-parallel

  - LC/Quartz: running the script *scripts/lc_quartz_hdf5.sh* will do
    this. 

        ./scripts/lc_quartz_hdf5.sh
    
    If you want to do it manually, the following commands should work
    when using GNU compilers v4.9.3 and MVAPICH2 v2.2, for example.

        export HDF5_DIR_SERIAL=/usr/tce/packages/hdf5/hdf5-serial-1.8.18-gcc-4.9.3
        export HDF5_DIR_PARALLEL=/usr/tce/packages/hdf5/hdf5-parallel-1.8.18-gcc-4.9.3-mvapich2-2.2

After this step, either run the install script (*install.sh*) 

    ./install.sh

**or** follow the steps below.

- Step 2: Copy the file "Make.defs.local" from the directory Chombo_Makefile to
  the directory Chombo/lib/mk/

      cp Chombo_Makefile/Make.defs.local Chombo/lib/mk/

- Step 3: Compile hypre - go to cogent/hypre-2.9.0b/ and run the following commands

      ./doconfig-lc-opt
      ./doinstall

- Step 4: Compile COGENT - go to cogent/exec/ and run

      make -j all OPT=TRUE DEBUG=FALSE

Quick Install on Specific Platforms
-----------------------------------

- LLNL-LC-Quartz:

      ./scripts/lc_quartz_hdf5.sh && ./install.sh

