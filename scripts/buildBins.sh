#!/bin/bash
# requires DATAROOT to be set to directory which contains mass bin directories

#export KEYDIR=$HOME/scratch/COMPASS/PWA2000/key5pi/SET7;
export KEYDIR=$ROOTPWA/keyfiles/keyKpipi/SETWA3;
export DATAROOT=$HOME/disk2/analysis/Kp_Kppipi/pwa/WA03_like;

echo "---- STARTING AMPLITUDE CALCULATIONS"; 
echo "---- WAVESET: $KEYDIR";
echo;

for bin in $DATAROOT/*; 
do
  echo "---- PROCESSING MASS BIN: $bin";
  echo "---- Starting time: ";
  date;
  export AMPDIR=$bin/AMPS;
#  export AMPDIR=$bin/amps;
  export FILE=`echo $bin | gawk -F"/" '{ print $(NF) }' `;
  export FILE=$bin/$FILE.evt;
  echo "---- input: $FILE";
  ./doamps.sh $FILE;

  # Do monte carlo
  export FILE=${FILE/evt/genbod.evt};
  echo "---- mc-input: $FILE";
  export AMPDIR=$bin/PSPAMPS;
#  export AMPDIR=$bin/pspamps;
  ./doamps.sh $FILE;

  # Do integration
  cd $AMPDIR
  int *.amp > norm.int;
  cd -;

  echo;
done;


