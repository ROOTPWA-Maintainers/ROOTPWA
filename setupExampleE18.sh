#!/usr/bin/bash

export LIBCONFIG=${HOME}/scratch/COMPASS/libconfig-1.4b4
export BOOST_ROOT=${HOME}/scratch/COMPASS/boost_1_42_0 
export ROOTPWA=${HOME}/scratch/COMPASS/rootpwa/trunk
export CMAKE=${HOME}/scratch/COMPASS/cmake-2.8.1


. prepare_root 5.26.
export LD_LIBRARY_PATH=${LIBCONFIG}/lib:${ROOTPWA}/build/lib:${LD_LIBRARY_PATH}
export PATH=${CMAKE}/bin:${ROOTPWA}/build/bin:${PATH}
