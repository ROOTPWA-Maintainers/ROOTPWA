#!/bin/bash

test -z $FITDIR && echo "FITDIR not set!" && exit;
test -z $OUTFILE && echo "OUTFILE not set!" && exit;
test -z $1 && echo "No bin dir given!" && exit;
test -z $RANK && export RANK=2;

echo "FITDIR=$FITDIR"

test -s $FITDIR/$OUTFILE && mv $FITDIR/$OUTFILE $FITDIR/$OUTFILE.previous

for i in $1/*;
do
  export BIN=`echo $i | gawk -F"/" '{ print $(NF) }' `;
  echo "---- Massbin: $BIN";
  LOWER=`echo $BIN | gawk -F"." '{ print $1 }' `;
  UPPER=`echo $BIN | gawk -F"." '{ print $2 }' `;
  echo "$LOWER...$UPPER"
  cd $i/AMPS;
  time pwaConstraintFit -q -w $FITDIR/$WLIST -c $FITDIR/$CONSTRFILE -o $FITDIR/$OUTFILE -r $RANK -l $LOWER -u $UPPER $RENORM -a $ACC -n $NORM

  #cd -;
done ;

#-S $FITDIR/$OUTFILE
