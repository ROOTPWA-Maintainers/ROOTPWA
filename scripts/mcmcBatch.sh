#!/bin/bash

#export FITDIR=/afs/e18/compass/analysis/sneubert/PWAFITS/fit6
test -z $FITDIR && echo "FITDIR not set!" && exit;
test -z $OUTFILE && echo "OUTFILE not set!" && exit;
test -z $START && echo "STARTFILE not set!" ;
test -z $1 && echo "No bin dir given!" && exit;

echo "FITDIR=$FITDIR"

test -s $FITDIR/$OUTFILE && mv $FITDIR/$OUTFILE $FITDIR/$OUTFILE.previous

for i in $1;
do
  export BIN=`echo $i | gawk -F"/" '{ print $(NF) }' `;
  echo "---- Massbin: $BIN";
  LOWER=`echo $BIN | gawk -F"." '{ print $1 }' `;
  UPPER=`echo $BIN | gawk -F"." '{ print $2 }' `;
  echo "$LOWER...$UPPER"
  cd $i/AMPS;
  time hamMCpwa -w $FITDIR/$WLIST -o $FITDIR/$OUTFILE -q -r 2 -l $LOWER -u $UPPER -S $START -s 0.001 -i 10000

  #cd -;
done ;
