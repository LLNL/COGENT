COGENT is a continuum (Eulerian) plasma simulation code.  It is primarily focused on tokamak edge plasma geometries, but includes options for, and is extensible to, other configurations. This repository contains the COGENT code (COGENT/) as well as Chombo (Chombo/), the adaptive mesh refinement application framework from Lawrence Berkeley National Laboratory upon which COGENT is built.

- The copy of Chombo in this repository is only for easily
  installing COGENT. It is developed and maintained by 
  LBL-ANAG, and the original copy and documentation can be
  downloaded from: 
  https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations

Quick Install on Specific Platforms
-----------------------------------

- On the following machines, just run the provided install script:
  - LLNL-LC-Quartz 
  - NERSC-Cori
  - NERSC-Edison

        ./install.sh

One can look inside this script to figure out the steps necessary for
compiling COGENT. Or read on ...

Install Instructions
--------------------

- It is assumed the C++ and FORTRAN compilers and their MPI wrappers are already
  available in the environment. The following environment variables must be 
  defined and set to the correct compilers:
  - CC (C compiler)
  - CXX (C++ compiler)
  - F77 (FORTRAN 77 compiler)
  - FC (FORTRAN 90 compiler)
  - MPICC (MPI C wrapper)
  - MPICXX (MPI C++ wrapper)
  - MPIF77 (MPI F77 wrapper)
  - MPIF90 (MPI F90 wrapper)

- Step 1: Set environment variables pointing to the HDF5 installation that is compiled
  with **the same compilers**  that COGENT will be compiled with.
  
    - On LC-Quartz, the following environment variables are needed
        - **HDF5_DIR_SERIAL** pointing to the serial HDF5 library
        - **HDF5_DIR_PARALLEL** pointing to the parallel HDF5 library

      In bash, this can be done by:
    
          export HDF5_DIR_SERIAL=/path/to/hdf5-serial
          export HDF5_DIR_PARALLEL=/path/to/hdf5-parallel
  
      For example, when using GNU compilers v4.9.3 and MVAPICH2 v2.2:
  
          export HDF5_DIR_SERIAL=/usr/tce/packages/hdf5/hdf5-serial-1.8.18-gcc-4.9.3
          export HDF5_DIR_PARALLEL=/usr/tce/packages/hdf5/hdf5-parallel-1.8.18-gcc-4.9.3-mvapich2-2.2

      should work.

    - On NERSC machines, the following environment variable is needed:
        - **HDF5_DIR** pointing to the parallel HDF5 library

      For example:
  
          export HDF5_DIR=/opt/cray/pe/hdf5-parallel/1.10.2.0/gnu/8.2/


      should work.

- Step 2: Copy the file appropriate "Make.defs.<machine_name>" from the directory 
  Chombo_Makefile to the directory Chombo/lib/mk/ and rename as "Make.defs.local"

      cp Chombo_Makefile/Make.defs.<machine_name> Chombo/lib/mk/Make.defs.local

  For example, on LC-Quartz,

      cp Chombo_Makefile/Make.defs.quartz Chombo/lib/mk/Make.defs.local

  and on NERSC machines (Cori/Edison),

      cp Chombo_Makefile/Make.defs.nersc Chombo/lib/mk/Make.defs.local

- Step 3: Compile hypre - go to COGENT/hypre-2.9.0b/ and run the following commands

      ./doconfig-opt
      ./doinstall

- Step 4: Compile COGENT - go to COGENT/exec/ and run

      make -j all OPT=TRUE DEBUG=FALSE

