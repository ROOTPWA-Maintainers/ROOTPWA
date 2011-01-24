#!/bin/bash
#///////////////////////////////////////////////////////////////////////////
#//
#//    Copyright 2009 Sebastian Neubert
#//
#//    This file is part of rootpwa
#//
#//    rootpwa is free software: you can redistribute it and/or modify
#//    it under the terms of the GNU General Public License as published by
#//    the Free Software Foundation, either version 3 of the License, or
#//    (at your option) any later version.
#//
#//    rootpwa is distributed in the hope that it will be useful,
#//    but WITHOUT ANY WARRANTY; without even the implied warranty of
#//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#//    GNU General Public License for more details.
#//
#//    You should have received a copy of the GNU General Public License
#//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
#//
#///////////////////////////////////////////////////////////////////////////


# Master Script for performing a genetic waveset run
# This will do NG generations of wavesets. Each generation
# will constitute a population of NW wavesets. 
# The mutants are created with mutator from a MOTHER waveset.
# Available waves are listed in the WPOOL file
# Data is stored into a directory tree rooted in GENEROOT
# All wavesets are fitted in parallel on the cluster.
# t.b.c.

export DATAROOT=/afs/e18/compass/analysis/sneubert/5PiLTData3/
#export GENEROOT=/afs/e18/compass/analysis/sneubert/PWAFITS/GENETICS/ltRUN22Extended
export GENEROOT=/nfs/nas/user/sneubert/PWAFITS/GENETICS/ltRUN30cont

#export MOTHER=$GENEROOT/gen2/results.dat
export MOTHER=$GENEROOT/parents.list
export STARTGEN=50

# in later generations several ancestors are possible
export NANC=50 # number of ancestors per generation (top)
export WPOOL=$GENEROOT/wavepool
export NW=50 # number of wavesets in one generation
export NG=100 # number of generations
export FIX=1   # first N waves to be fixed
export EXCH=1; # number of waves to exchange per mutation
export ADD=0;  # number of waves to add ...
export DROP=0; # number of waves to drop ...
export VARY=5; # range -V..+V by which to vary number of waves per mutation
export PRESSURE=1.7; # selective pressure 1..2
export CROSS=0.85; # crossover probability 0..1

# define additional options to pass to qsub
# for example to exclude single nodes:
# export SGEOPT="-l h='!node59|!node60|'"


# check if we are ready to go
if ! $ROOTPWA/scripts/checkWavesForFit.sh $WPOOL $DATAROOT; then 
    echo "Go and check your amplitudes!";
    exit 1;
fi


# loop over generations
for ((  g = $STARTGEN ;  g < $NG;  g++  )) ; do
echo ;
GENDIR=$GENEROOT/gen$g;
mkdir $GENDIR;
cp $MOTHER $GENDIR/ancestors.dat;
export MOTHER=$GENDIR/ancestors.dat;

echo "---------------- Generation $g -----------------";
echo "------------------------------------------------";
echo "------------------------------------------------";

# generate mutant population and start fits
for((  i = 0 ;  i < $NW;  i++  ))  ; do
    SETDIR=$GENDIR/set$i;
    mkdir $SETDIR;
    WLIST=wavelist;
    MUTATION=0;
    echo "mutator -S$i -E$EXCH -A$ADD -D$DROP -V$VARY -P$WPOOL -F$FIX -p$PRESSURE -c$CROSS $MOTHER > $SETDIR/$WLIST"
    if mutator -S$i -E$EXCH -A$ADD -D$DROP -V$VARY -P$WPOOL -F$FIX -p$PRESSURE -c$CROSS $MOTHER > $SETDIR/$WLIST; then
    # start fits. Pass WLIST and FITDIR to script
    
	echo "Starting pwafit in directory $SETDIR with wavelist $WLIST"
	QMES=(`qsub $SGEOPT -v "WLIST=$WLIST,FITDIR=$SETDIR" -t 2-7:1 fitgene.sge | tr ' ' ' '`)
	JOBNR=${QMES[2]}
	echo "Jobnumber: $JOBNR"
	
    else
	echo "++++++ Mutator failed!"
    fi
    
done; # end loop over populations

# wait for fits to be finished
date
echo -n "Generation $g: Waiting for fits to finish ..."
OUTFILE=`mktemp -t qstat_resXXX` 
while [ 1 ]; do
    qstat &> $OUTFILE;
    # restart pending error jobs
    cat $OUTFILE | grep Eqw | while read a; do
	qmod -cj ${a:0:7} &> /dev/null;
    done
	
    if [ $?=0 ]; then
	if grep -q error $OUTFILE; then 
	    echo "qstat Error!";
	else
	    if ! grep -q fitgene $OUTFILE; then
		echo " no jobs running anymore";
		break;
	    fi;
	fi;
    fi;
    if [ -n "$RES" ]; then
	echo -n ".";
    fi;
    sleep 30;
done;
date

# collect fit results and determine the best fit (use all generations)

selector 28 $NANC $MOTHER $GENDIR/set*/ > $GENDIR/results.dat
export MOTHER=$GENDIR/results.dat



done; # end loop over generations
