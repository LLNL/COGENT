# -*- Mode: Makefile; -*-

# the location of Chombo lib dir
CHOMBO_HOME = $(shell pwd)/../../../Chombo/lib

ebase = test_Annulus

MINDIM = 1
MAXDIM = 2

# this is the local GNUmakefile which contains this example's multidim
# build info -- libraries and source directory information
MULTIDIM_MAKEFILE = GNUmakefile.multidim

all: move-script all-multidim

noLibs: all-multidim-nolibs

#all: all-multidim

move-script: 
	cd $(CHOMBO_HOME)/../example/fourthOrderMappedGrids; \
	./moveToBoxTools.sh;

# this file contains the basic rules used to build multidim codes (using the 
# GNUmakefile.multidim in this directory), including the shell script
# which orchestrates the make process 
include $(CHOMBO_HOME)/mk/Make.multidim.basic
