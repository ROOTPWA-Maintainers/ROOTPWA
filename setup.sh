#!/usr/bin/bash

export LIBCONFIG=$HOME/scratch/COMPASS/libconfig-1.4b4
export PWA2000=$HOME/scratch/COMPASS/BNL
export ROOTPWA=$HOME/scratch/COMPASS/rootpwa/trunk
export QTROOT=/afs/e18/sw/e18/lib/qt4/4.5.1/
export QTINCL=$QTROOT/include
export QTLIB=$QTROOT/lib
export QTBIN=$QTROOT/bin

. prepare_root 5.24.
#export ROOTSYS=/opt/sw/ROOT/root_v5.24.00.gcc412/
export LD_LIBRARY_PATH=$PWA2000/lib:$LIBCONFIG/lib:$ROOTPWA/build/lib:$QTLIB:$LD_LIBRARY_PATH
export PATH=$ROOTPWA/build/bin:$PWA2000/bin:$QTBIN:$PATH


