#!/usr/bin/bash

export ROOTSYS=/afs/cern.ch/sw/lcg/external/root/5.20.00/slc4_amd64_gcc34/root/
export LIBCONFIG=$HOME/scratch/COMPASS/libconfig-1.4b4
export ROOTPWA=$HOME/scratch/COMPASS/rootpwa/trunk

export LD_LIBRARY_PATH=$LIBCONFIG/lib:$ROOTPWA/build/lib:$LD_LIBRARY_PATH
export PATH=$ROOTPWA/build/bin:$PATH




