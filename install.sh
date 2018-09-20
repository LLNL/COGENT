#!/bin/bash

CH_MAKEFILE_DIR=Chombo_Makefile
CH_MAKEFILE=Make.defs.local
CH_LIB_DIR=Chombo/lib/mk

COGENT_DIR=COGENT
HYPRE_DIR=hypre-2.9.0b
COGENT_EXEC_DIR=exec

ROOT_DIR=$PWD

#copy Chombo makefile
cp $CH_MAKEFILE_DIR/$CH_MAKEFILE $CH_LIB_DIR/

# build Hypre
cd $COGENT_DIR/$HYPRE_DIR
./doconfig-lc-opt
./doinstall
cd $ROOT_DIR

# build COGENT
cd $COGENT_DIR/$COGENT_EXEC_DIR
make -j all OPT=TRUE DEBUG=FALSE
cd $ROOT_DIR
