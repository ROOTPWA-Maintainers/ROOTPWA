#!/bin/bash

cd $ROOTPWA/src/pwafitTest

./runGamp.sh

./runPwaFit.sh

root -b -q -l ./runCompareTFitBins.C
# function returns true if bins match; invert logic for exit status
if [[ "$?" -eq "1" ]]
then
    exit 0
fi

exit 1
