#!/bin/bash

export FITDIR=$ROOTPWA_DATA_DIR/fits
export FITFILE=$FITDIR/fit.root
export PLOTFILE=${FITFILE/.root/.plots.root}

echo FITDIR=$FITDIR
echo FITFILE=$FITFILE
echo PLOTFILE=$PLOTFILE

cd $FITDIR
root -l -b -q "$ROOTPWA/generators/rootlogon.C" "$ROOTPWA/generators/plotGlobalWeightedEvts_3pin.C+(\"fit.plots.root\", \"weightedMC.globalplots.root\")"

exit
