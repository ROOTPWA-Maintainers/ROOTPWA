#!/bin/bash

export bindir=$1;
export bin=`echo $bindir | gawk -F"/" '{ print $(NF) }' `;
echo "---- BIN-DIRECTORY: $bindir";
echo "---- PROCESSING MASS BIN: $bin";
echo "---- Starting time: ";
date;
if [ -s $bindir/$bin.evt ]; then
    test -s $bindir/AMPS || mkdir $bindir/AMPS;
    export AMPDIR=$bindir/AMPS;
    export FILE=$bindir/$bin.evt;
    echo "---- input: $FILE";
    ./doamps.sh $FILE;
fi;

  # Do monte carlo

export FILE=$bindir/$bin.genbod.evt;
if [ -s $FILE ]; then
    echo "---- mc-input: $FILE";
    test -s "$bindir/PSPAMPS" || mkdir $bindir/PSPAMPS;
    export AMPDIR=$bindir/PSPAMPS;
    ./doamps.sh $FILE;

  # Do integration
   # cd $AMPDIR
   # int *.amp > norm.int;
   # /bin/cp norm.int $bindir/AMPS;
   # cd -;
fi;

  # Do accepted monte carlo
export FILE=$bindir/$bin.acc.evt;
if [ -s $FILE ]; then
    echo "---- mc-acc-input: $FILE";
    test -s "$bindir/ACCAMPS" || mkdir $bindir/ACCAMPS;
    export AMPDIR=$bindir/ACCAMPS;
    ./doamps.sh $FILE;

  # Do integration
  #  cd $AMPDIR
  #  int *.amp > accnorm.int;
  #  /bin/cp accnorm.int $bindir/AMPS;
  #  cd -;
fi;

echo;