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

export GENEROOT=/afs/e18/compass/analysis/sneubert/PWAFITS/GENETICS/ltRUN5
export MOTHER=$GENEROOT/parents.dat 
# in later generations several ancestors are possible
export NANC=40 # number of ancestors per generation
export WPOOL=$GENEROOT/wavepool
export NW=50 # number of wavesets in one generation
export NG=40 # number of generations
export FIX=1   # first N waves to be fixed
export EXCH=2; # number of waves to exchange per mutation
export ADD=0;  # number of waves to add ...
export DROP=1; # number of waves to drop ...

# loop over generations
for ((  g = 0 ;  g < $NG;  g++  )) ; do
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
    #echo "mutator -S$i -E$EXCH -A$ADD -D$DROP -P$WPOOL -F$FIX $MOTHER > $SETDIR/$WLIST"
    mutator -S$i -E$EXCH -A$ADD -D$DROP -P$WPOOL -F$FIX $MOTHER > $SETDIR/$WLIST
    # start fits. Pass WLIST and FITDIR to script
    
	echo "Starting pwafit in directory $SETDIR with wavelist $WLIST"
	qsub -v "WLIST=$WLIST,FITDIR=$SETDIR" -t 2-7:1 fitgene.sge
   
done; # end loop over populations

# wait for fits to be finished

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


# collect fit results and determine the best fit

selector 14 $NANC $GENDIR/set* > $GENDIR/results.dat
export MOTHER=$GENDIR/results.dat



done; # end loop over generations
