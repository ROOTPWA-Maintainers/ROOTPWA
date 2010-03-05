#!/usr/bin/bash

export LIBCONFIG=$HOME/scratch/COMPASS/libconfig-1.4b4
export ROOTPWA=$HOME/scratch/COMPASS/rootpwa/trunk


. prepare_root 5.24.
#export ROOTSYS=/opt/sw/ROOT/root_v5.24.00.gcc412/
export LD_LIBRARY_PATH=$LIBCONFIG/lib:$ROOTPWA/build/lib:$LD_LIBRARY_PATH
export PATH=$ROOTPWA/build/bin:$PATH

