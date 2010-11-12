#!/usr/bin/bash

export LIBCONFIG=$HOME/scratch/COMPASS/libconfig-1.4b4
export ROOTPWA=$HOME/scratch/COMPASS/rootpwa/branches/gpu-test
export CMAKE=$HOME/scratch/COMPASS/cmake-2.8.1
export BOOST_ROOT=$HOME/scratch/COMPASS/boost_1_42_0 

. prepare_root 5.26
#export ROOTSYS=/opt/sw/ROOT/root_v5.24.00.gcc412/
export LD_LIBRARY_PATH=$LIBCONFIG/lib:$ROOTPWA/build/lib:$LD_LIBRARY_PATH
export PATH=$CMAKE/bin:$ROOTPWA/build/bin:$PATH

