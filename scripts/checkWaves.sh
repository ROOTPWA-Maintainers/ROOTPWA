#!/bin/bash

# $2 is reference bin
# $1 gives root-directory with bins to compare
# e.g.  checkWaves.sh ~/mydirector 1000.1120/AMPS

cd $1/$2

SELECT=$3

if [ -z $SELECT ]; then SELECT=amp ; fi

ls *.amp > /tmp/refamps
cat /tmp/refamps

cd $1

COUNT=1

for i in $1/*; do
    cd $i/PSPAMPS/
    ls *.amp > /tmp/diffamps
    echo "Checking bin $COUNT: "
    let COUNT=$COUNT+1
    diff /tmp/diffamps /tmp/refamps | grep $SELECT
    
done;

