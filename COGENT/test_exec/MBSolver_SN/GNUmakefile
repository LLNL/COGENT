# -*- Mode: Makefile; -*-

# the location of Chombo lib dir
CHOMBO_HOME = ../../../Chombo/lib

ebase = mbsolver_sn

MINDIM = 1
MAXDIM = 4

export MAXDIM 

# this is the local GNUmakefile which contains this example's multidim
# build info -- libraries and source directory information
MULTIDIM_MAKEFILE = GNUmakefile.multidim

all:  all-multidim

noLibs: all-multidim-nolibs

# this file contains the basic rules used to build multidim codes (using the 
# GNUmakefile.multidim in this directory), including the shell script
# which orchestrates the make process 
include $(CHOMBO_HOME)/mk/Make.multidim.basic
