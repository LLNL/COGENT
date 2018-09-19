#!/bin/bash

# For use on LC/Quartz only!
# Once the correct compilers are loaded, this
# script will load the corresponding HDF5 libraries
# and set the HDF5-related environment variables
# needed for compiling COGENT.

H5DIFF_SUFFIX=/bin/h5diff
module load hdf5-serial
H5DIFF_PATH=$(which h5diff)
HDF5_DIR=${H5DIFF_PATH%$H5DIFF_SUFFIX}
export HDF5_DIR_SERIAL="${HDF5_DIR}"
module load hdf5-parallel
H5DIFF_PATH=$(which h5diff)
HDF5_DIR=${H5DIFF_PATH%$H5DIFF_SUFFIX}
export HDF5_DIR_PARALLEL="${HDF5_DIR}"
