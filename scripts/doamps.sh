#!/bin/bash

WORKDIR=$PWD
echo $KEYDIR
cd $KEYDIR

COUNTER=0;

echo "MAXKEYPERJOB = $MAXKEYPERJOB"

for i in *.key ; 
  do echo "Processing keyfile: $i" ;
     date;
     outfile=${i/key/amp};
     # check if amplitude is already there! 
     test -s $AMPDIR/$outfile && echo "File $AMPDIR/$outfile already exists! Skipping!"
     test -s $AMPDIR/$outfile && continue;
     test -s $AMPDIR/OLD/$outfile && echo "File $AMPDIR/OLD/$outfile already exists! Skipping!"
     test -s $AMPDIR/OLD/$outfile && continue;

     test -s $AMPDIR/$outfile.inprogress && echo "File $AMPDIR/$outfile is being processed already! Skipping!"
     test -s $AMPDIR/$outfile.inprogress && continue;
     date > $AMPDIR/$outfile.inprogress
     test -s $AMPDIR/$outfile || cat $1 | gamp -P $PDG $i > $AMPDIR/$outfile ;
     rm -f $AMPDIR/$outfile.inprogress
     let COUNTER=$COUNTER+1;
     if [ $COUNTER -ge $MAXKEYPERJOB ] ; then 
	 echo "Reached max key per job limit! Aborting";
	 break; 
     fi;

     #test -s $AMPDIR/$outfile || echo do it > $AMPDIR/$outfile;
done

# now do adding of amps
#for i in *list.dat ; 
#  do echo "Adding amplitudes as specified in $i" ;
#  addamp $i $AMPDIR/OLD/ ;
#done;

cd $WORKDIR
