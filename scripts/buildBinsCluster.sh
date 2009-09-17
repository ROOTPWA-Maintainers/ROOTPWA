#!/bin/bash

export KEYDIR=$HOME/scratch/COMPASS/PWA2000/key5pi/SET5;
echo "---- STARTING AMPLITUDE CALCULATIONS"; 
echo "---- WAVESET: $KEYDIR";
echo;

test -s $KEYDIR || exit
ls $KEYDIR

export DATADIR=/nfs/data/user/sneubert/5PiData2

for bin in $DATADIR/*;
  do echo "---- Submitting job for bin $bin";
     qsub -q medium -v KEY=$KEYDIR,BIN=$bin buildBinCluster.sh;
done;

echo "SUBMITTED JOBS"
