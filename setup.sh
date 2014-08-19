#!/usr/bin/bash

export LIBCONFIG=${HOME}/scratch/COMPASS/libconfig-1.4b4
export BOOST_ROOT=${HOME}/scratch/COMPASS/boost_1_42_0
#export PWA2000=${HOME}/scratch/COMPASS/BNL
export ROOTPWA=${HOME}/scratch/COMPASS/rootpwa/trunk
export CMAKE=${HOME}/scratch/COMPASS/cmake-2.8.1
#export QTROOT=/afs/e18/sw/e18/lib/qt4/4.5.1/
#export QTINCL=${QTROOT}/include
#export QTLIB=${QTROOT}/lib
#export QTBIN=${QTROOT}/bin


. prepare_root 5.26
#export ROOTSYS=/opt/sw/ROOT/root_v5.24.00.gcc412/
#export LD_LIBRARY_PATH=${PWA2000}/lib:${QTLIB}:${LD_LIBRARY_PATH}
#export PATH=${PWA2000}/bin:${QTBIN}:${PATH}
export LD_LIBRARY_PATH=${LIBCONFIG}/lib:${ROOTPWA}/build/lib:${LD_LIBRARY_PATH}
export PATH=${CMAKE}/bin:${ROOTPWA}/build/bin:${PATH}
