#!/usr/bin/bash

export PWA2000=$HOME/scratch/COMPASS/BNL
export ROOTPWA=$HOME/scratch/COMPASS/rootpwa/trunk


. prepare_root 5.24.
#export ROOTSYS=/opt/sw/ROOT/root_v5.24.00.gcc412/
export LD_LIBRARY_PATH=$PWA2000/lib:$ROOTPWA/build/lib:$LD_LIBRARY_PATH
export PATH=$ROOTPWA/build/bin:$PWA2000/bin:$PATH


