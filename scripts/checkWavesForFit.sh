#!/bin/bash

# checks if all waves in a waveset are in the integrals
# usage: checkWavesInIntegrals <Wavesetfile> <DataDir>

SET=$1
ROOTDIR=$2

ERRCOUNT=0

while read amp; do
    echo "Checking wave $amp ...";
    for bin in $ROOTDIR/*; do
	if [ ! -s $bin/AMPS/$amp ]; then
	    echo "$amp not found in $bin/AMPS";
	    let ERRCOUNT=$ERRCOUNT+1
	    echo $ERRCOUNT
	fi;
	if ! grep -q $amp $bin/AMPS/norm.int; then
	    echo "$amp not found in $bin/AMPS/norm.int";
	    let ERRCOUNT=$ERRCOUNT+1
	    echo $ERRCOUNT
	fi;
	if ! grep -q $amp $bin/AMPS/accnorm.int; then
	    echo "$amp not found in $bin/AMPS/accnorm.int";
	    let ERRCOUNT=$ERRCOUNT+1
	    echo $ERRCOUNT
	fi;
    done
done < $SET

echo "Number of Errors=$ERRCOUNT"

exit $ERRCOUNT
