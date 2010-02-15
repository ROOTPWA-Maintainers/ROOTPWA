#!/bin/bash

export WORKDIR=/afs/e18/compass/analysis/sneubert/
export FITDIR=$WORKDIR/PWAFITS/LOWT/fit11/
export FITFILE=$FITDIR/fit11.acc.root
export PLOTFILE=${FITFILE/.root/.plots.root}
export BOOKY=${FITFILE/.root/.booky.pdf}
export DATADIR=$WORKDIR/5PiLTData3/

echo WORKDIR=$WORKDIR
echo FITDIR=$FITDIR
echo FITFILE=$FITFILE
echo PLOTFILE=$PLOTFILE
echo BOOKY=$BOOKY
echo DATADIR=$DATADIR


cd $DATADIR


for i in *; do
    echo "$i";
    cd $i;
    # convert events to root tree if not already done
    test -s $i.root || cat $i.evt | evt2tree $i.root;
    # run evtweight on accepted events:
    cd ACCAMPS
    WEIGHTEDFILE=${FITFILE/.root/.kineval.$i.root}
    test -s $WEIGHTEDFILE || evtweight -e ../$i.acc.evt -o $WEIGHTEDFILE  -w $FITFILE -i accnorm.int -m $i  
    # produce nice plots
    root -b -q "$ROOTPWA/generators/doPlotWEvts.C(\"../$i.root\",\"$WEIGHTEDFILE\",\"$PLOTFILE\",\"$i\")"

    cd $DATADIR;
done;

echo "CREATING BOOKY $BOOKY ..."
# collect all plots into a booky:
cd $FITDIR
for i in *plots*.ps; do ps2pdf $i; rm -f $i; done;
pdftk *plots*.pdf cat output $BOOKY
for i in *plots*.pdf; do rm -f $i; done;


