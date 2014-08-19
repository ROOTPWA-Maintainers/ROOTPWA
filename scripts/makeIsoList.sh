#!/bin/bash


cd /afs/e18/compass/analysis/sneubert/5PiLTData/1000.1120/AMPS

export KEYS=( 1-3????pi-*a2*.amp );

export NG=${#KEYS[@]};

echo NG=$NG;

for ((  i = 0 ;  i < $NG;  i++  )) ; do
    for (( j = i+1 ; j < $NG; j++ )) ; do
	addisokey ${KEYS[$i]} ${KEYS[$j]};
    done;
done
