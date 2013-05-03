#!/bin/bash

# script to check if all waves from a given set have been calculated
# takes as a parameter the data directory and the waveset directory

DATA=$1
SET=$2

# build list
cd $SET
echo $PWD

rm -f $DATA/amplist
ls -1 $DATA/*.amp > $DATA/amplist
ls -1 $DATA/OLD/*.amp >> $DATA/amplist

for j in F2BINS ; do
cd $j;
echo "Checking Set $j"
  FOUND=0;
  for i in *.key ; do
      THISFOUND=0;
      grep -q ${i/key/amp} $DATA/amplist || THISFOUND=1;
      if [ $THISFOUND -gt 0 ]; then
	FOUND=1;
	echo "$i missing from $j"
      fi
  done
  if [ $FOUND -gt 0 ]; then echo "+++ Waveset $j incomplete!"; 
  else echo "Waveset $j fine!"; 
  fi
cd - 
done
cat $DATA/amplist | wc