#----------------------------------------------------------------------------
#          Makefile for the EFITXS package
#
#----------------------------------------------------------------------------

ifndef $(CPUNAME)
  CPUNAME = $(shell uname -n)
endif

ifneq ($(findstring gre,$(CPUNAME)),)  #grendel.llnl.gov
  PDB_PATH = /usr/local/pact/pact04_05_11/
  IDL_DIR = /usr/local/rsi/idl
  CC = cc
  FF = pgf90
  LINKER = pgf90
endif

ifneq ($(findstring sma,$(CPUNAME)),)  #smaug.llnl.gov
  #use /afs/localcell/usr/rsi/idl70/bin/idl
  PDB_PATH = /afs/localcell/usr/pact/@sys/pact07_07_18/
  IDL_DIR = /afs/localcell/usr/rsi/idl70
  CC = cc
  FF = ifort ##f95
  LINKER = ifort ##f95
endif

ifneq ($(findstring hop,$(CPUNAME)),)  # hopper.nersc.gov
  IDL_DIR =  /usr/common/usg/idl/idl80/
  CC = cc
  FF = gfortran
  LINKER = gfortran
endif

ifneq ($(findstring cori,$(CPUNAME)),)  # cori.nersc.gov
  IDL_DIR =  /usr/common/software/idl/idl83/
  CC = cc
  FF = gfortran
  LINKER = gfortran
endif

ifneq ($(findstring edis,$(CPUNAME)),)  # edison.nersc.gov
  ##-do first: module load gcc
  IDL_DIR =  /usr/common/usg/idl/idl82/
  CC = cc
  FF = gfortran
  LINKER = gfortran
endif

ifneq ($(findstring cab,$(CPUNAME)),)  # cab.llnl.gov
  IDL_DIR =  /usr/local/tools/idl-85/idl85/
  CC = cc
  FF = gfortran
  LINKER = gfortran
endif

ifneq ($(findstring sunfire,$(CPUNAME)),)  # sunfire.pppl.gov
  IDL_DIR =  /usr/pppl/itt/idl/idl82
  CC = cc
  #FF = pathf95
  #LINKER = pathf95
  FF = gfortran
  LINKER = gfortran
endif

CFILES = io_gridue_wrap.c
FFILES = io_gridue.f90

COBJECTS = io_gridue_wrap.o
FOBJECTS = io_gridue.o
EXECFIL = io_gridue.so

X_LIBS =
LIBS =

CFLAGS = -I$(IDL_DIR)/external/include -fPIC -D UNDERSCORE
FFLAGS = -I$(IDL_DIR)/external/include -fPIC

LDFLAGS = -shared
X_LD_FLAGS = -I$(IDL_DIR)/external/include
EXPORT_CFLG = 
F_LD_POST= 



OBJECTS = $(COBJECTS) $(FOBJECTS)


$(EXECFIL) : $(OBJECTS)
	$(LINKER)  $(X_LD_FLAGS) $(LDFLAGS) $(EXPORT_FLAGS) -o $(EXECFIL) \
	$(OBJECTS) $(LIBS)
	$(F_LD_POST)



io_gridue_wrap.o: io_gridue_wrap.c
	$(CC) $(CFLAGS) -c io_gridue_wrap.c

io_gridue.o: io_gridue.f90
	$(FF) $(FFLAGS) -c io_gridue.f90


clean: 
	rm -f *.o *.so core.* *.mod *.export
