#!/bin/bash

export FITDIR=$PWA_DATA_DIR/fits
#export FITDIR=$WORKDIR/PWAFITS/GENETICS/ltRUN21/gen39/set42
export FITFILE=$FITDIR/fit.root
export PLOTFILE=${FITFILE/.root/.plots.root}
export BOOKY=${FITFILE/.root/.booky.pdf}
export DATADIR=$PWA_DATA_DIR

echo FITDIR=$FITDIR
echo FITFILE=$FITFILE
echo PLOTFILE=$PLOTFILE
echo BOOKY=$BOOKY
echo DATADIR=$DATADIR


cd $DATADIR

rm $PLOTFILE;

# here generate sequence of the mass bin folders only ...

for i in `ls | awk '/[0-9].[0-9]/ {print $1}'`; do
#for i in "1300.1340"; do
    echo "MassBin: $i";
    cd $i;
    # convert events to root tree if not already done
    test -s $i.neubert.root || cat $i.evt | evt2tree $i.neubert.root;
    # run evtweight on accepted events:
    #cd ACCAMPS
    cd PSPAMPS
    WEIGHTEDFILE=${FITFILE/.root/.kineval.$i.root}
    test -s $WEIGHTEDFILE || evtweight -e ../$i.ps.evt -o $WEIGHTEDFILE  -w $FITFILE -i norm.int -m $i;
    echo "created weightfile for massbin $i";
    # produce nice plots
    root -b -q "$ROOTPWA/generators/rootlogon.C" "$ROOTPWA/generators/doPlotWEvts.C(\"../$i.neubert.root\",\"$WEIGHTEDFILE\",\"$PLOTFILE\",\"$i\")";

    cd $DATADIR;
done;


#echo "CREATING BOOKY $BOOKY ..."
## collect all plots into a booky:
#cd $FITDIR
#for i in *plots*.ps; do ps2pdf $i; rm -f $i; done;
#pdftk *plots*.pdf cat output $BOOKY
#for i in *plots*.pdf; do rm -f $i; done;
