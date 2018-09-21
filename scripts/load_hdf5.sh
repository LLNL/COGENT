#!/bin/bash

# For use on LC/Quartz only!
# Once the correct compilers are loaded, this
# script will load the corresponding HDF5 libraries
# and set the HDF5-related environment variables
# needed for compiling COGENT.

H5DIFF_SUFFIX=/bin/h5diff

echo "Loading serial HDF5"
if [[ $HOSTNAME == "quartz"* ]]; then
  module load hdf5-serial
fi
H5DIFF_PATH=$(which h5diff)
echo "h5diff is $H5DIFF_PATH"
HDF5_DIR=${H5DIFF_PATH%$H5DIFF_SUFFIX}
export HDF5_DIR_SERIAL="${HDF5_DIR}"
echo "Serial HDF5: $HDF5_DIR_SERIAL"

echo "Loading parallel HDF5"
if [[ $HOSTNAME == "quartz"* ]]; then
  module load hdf5-parallel
fi
H5DIFF_PATH=$(which h5diff)
echo "h5diff is $H5DIFF_PATH"
HDF5_DIR=${H5DIFF_PATH%$H5DIFF_SUFFIX}
export HDF5_DIR_PARALLEL="${HDF5_DIR}"
echo "Parallel HDF5: $HDF5_DIR_PARALLEL"
