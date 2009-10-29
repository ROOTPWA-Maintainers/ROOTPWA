#!/bin/bash

WORKDIR=$PWD
cd $KEYDIR

for i in *.key ; 
  do echo "Processing keyfile: $i" ;
     date;
     outfile=${i/key/amp};
     # check if amplitude is already there! 
     test -s $AMPDIR/$outfile && echo "File $AMPDIR/$outfile already exists! Skipping!"
     test -s $AMPDIR/$outfile && continue;
     test -s $AMPDIR/OLD/$outfile && echo "File $AMPDIR/OLD/$outfile already exists! Skipping!"
     test -s $AMPDIR/OLD/$outfile && continue;
     test -s $AMPDIR/$outfile || cat $1 | gamp -P ../pdgTable.txt $i > $AMPDIR/$outfile ; 
     #test -s $AMPDIR/$outfile || echo do it > $AMPDIR/$outfile;
done

# now do adding of amps
for i in *list.dat ; 
  do echo "Adding amplitudes as specified in $i" ;
  addamp $i $AMPDIR/OLD/ ;
done;

cd $WORKDIR