#!/bin/bash
# requires DATAROOT to be set to directory which contains mass bin directories

export KEYDIR=$HOME/scratch/COMPASS/PWA2000/key5pi/SET7;

echo "---- STARTING AMPLITUDE CALCULATIONS"; 
echo "---- WAVESET: $KEYDIR";
echo;

for bin in $DATAROOT/*; 
do
  echo "---- PROCESSING MASS BIN: $bin";
  echo "---- Starting time: ";
  date;
  export AMPDIR=$bin/AMPS;
  export FILE=`echo $bin | gawk -F"/" '{ print $(NF) }' `;
  export FILE=$bin/$FILE.evt;
  echo "---- input: $FILE";
  ./doamps.sh $FILE;

  # Do monte carlo
  export FILE=${FILE/evt/genbod.evt};
  echo "---- mc-input: $FILE";
  export AMPDIR=$bin/PSPAMPS;
  ./doamps.sh $FILE;

  # Do integration
  cd $AMPDIR
  int *.amp > norm.int;
  cd -;

  echo;
done;


