#!/usr/bin/bash

export LIBCONFIG=${HOME}/scratch/COMPASS/libconfig-1.4b4
export BOOST_ROOT=${HOME}/scratch/COMPASS/boost_1_42_0 
export ROOTPWA=${HOME}/scratch/COMPASS/rootpwa/trunk
export CMAKE=${HOME}/scratch/COMPASS/cmake-2.8.1

# setup ROOT
export ROOTSYS=/afs/cern.ch/sw/lcg/external/root/5.20.00/slc4_amd64_gcc34/root/
source ${ROOTSYS}/bin/thisroot.sh

export LD_LIBRARY_PATH=${LIBCONFIG}/lib:${ROOTPWA}/build/lib:${LD_LIBRARY_PATH}
export PATH=${CMAKE}/bin:${ROOTPWA}/build/bin:${PATH}
