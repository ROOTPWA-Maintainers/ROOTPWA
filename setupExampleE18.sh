#!/usr/bin/bash

export PATH=$PATH:$HOME/scratch/SWIPL/pl-5.6.52/bin

. prepare_root 5.20.

export PWA2000=$HOME/scratch/COMPASS/BNL
export ROOTPWA=$PWD;
export FITS=/afs/e18/compass/analysis/sneubert/PWAFITS

export PATH=$PATH:$ROOTPWA/bin:$PWA2000/bin


export BAT=/afs/e18/compass/analysis/sneubert/bat-0.1
export BATHEADERS=$BAT/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BAT:$ROOTPWA/lib