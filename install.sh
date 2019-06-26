#!/bin/bash

########################
# set up the environment
########################

# Load GNU compilers
echo "---------------------------------------"
if [[ $HOSTNAME == "quartz"* ]]; then
  echo "Loading GNU compilers for Quartz (LLNL)"
  module load gcc
  module load mvapich2/2.2
  export CC=$(which gcc)
  export CXX=$(which g++)
  export CPP=$(which cpp)
  export F77=$(which gfortran)
  export FC=$(which gfortran)
  export MPICC=$(which mpicc)
  export MPICXX=$(which mpicxx)
  export MPIF77=$(which mpif77)
  export MPIF90=$(which mpif90)
  export MACHINE="quartz"
else
  if [[ $HOSTNAME == "cori"* || $HOSTNAME == "edison"* ]]; then
    echo "Loading GNU compilers for Cori/Edison (NERSC)"
    module swap PrgEnv-intel PrgEnv-gnu
    export CC=$(which cc)
    export CXX=$(which CC)
    export CPP=$(which cpp)
    export F77=$(which ftn)
    export FC=$(which ftn)
    export MPICC=$(which cc)
    export MPICXX=$(which CC)
    export MPIF77=$(which ftn)
    export MPIF90=$(which ftn)
    export MACHINE="nersc"
  fi
fi
echo "CC is $CC"
echo "CXX is $CXX"
echo "CPP is $CPP"
echo "F77 is $F77"
echo "FC is $FC"
echo "MPICXX is $MPICXX"
echo "---------------------------------------"

# Load HDF5 libraries and set 
# HDF5-related environment vars

H5DIFF_SUFFIX=/bin/h5diff
if [[ $HOSTNAME == "quartz"* ]]; then
  echo "Loading HDF5-serial for Quartz (LLNL)"
  module load hdf5-serial/1.8.18
  H5DIFF_PATH=$(which h5diff)
  echo "h5diff is $H5DIFF_PATH"
  HDF5_DIR=${H5DIFF_PATH%$H5DIFF_SUFFIX}
  export HDF5_DIR_SERIAL="${HDF5_DIR}"
  echo "Serial HDF5: $HDF5_DIR_SERIAL"
  echo "Loading HDF5-parallel for Quartz (LLNL)"
  module load hdf5-parallel/1.8.18
  H5DIFF_PATH=$(which h5diff)
  echo "h5diff is $H5DIFF_PATH"
  HDF5_DIR=${H5DIFF_PATH%$H5DIFF_SUFFIX}
  export HDF5_DIR_PARALLEL="${HDF5_DIR}"
  echo "Parallel HDF5: $HDF5_DIR_PARALLEL"
else
  if [[ $HOSTNAME == "cori"* || $HOSTNAME == "edison"* ]]; then
    echo "Loading HDF5 for Cori/Edison (NERSC)"
    export HDF5_DIR="/opt/cray/pe/hdf5-parallel/1.10.0.3/GNU/5.1"
    echo "HDF5: $HDF5_DIR"
  fi
fi
echo "---------------------------------------"

################
# Compile COGENT
################

CH_MAKEFILE_DIR=Chombo_Makefile
CH_MAKEFILE=Make.defs.$MACHINE
CH_LIB_DIR=Chombo/lib/mk

COGENT_DIR=COGENT
HYPRE_DIR=hypre-2.9.0b
COGENT_EXEC_DIR=exec

ROOT_DIR=$PWD

#copy Chombo makefile
cp $CH_MAKEFILE_DIR/$CH_MAKEFILE $CH_LIB_DIR/Make.defs.local

# build Hypre
cd $COGENT_DIR/$HYPRE_DIR/
./doconfig-opt
./doinstall
cd $ROOT_DIR

# build COGENT
cd $COGENT_DIR/$COGENT_EXEC_DIR
make -j all OPT=TRUE DEBUG=FALSE
cd $ROOT_DIR
