#!/usr/bin/bash

export PATH=$PATH:$HOME/scratch/SWIPL/pl-5.6.52/bin:$HOME/scratch/COMPASS/PWA2000/bin


. prepare_root 5.20.

export PWA2000=$PWD
export FITS=/afs/e18/compass/analysis/sneubert/PWAFITS

export BAT=/afs/e18/compass/analysis/sneubert/bat-0.1
export BATHEADERS=$BAT/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BAT