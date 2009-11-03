#!/bin/bash

# $2 is reference bin
# $1 gives root-directory with bins to compare
# e.g.  checkWaves.sh ~/mydirector 1000.1120/AMPS

cd $1/$2

ls *.amp > /tmp/refamps
cat /tmp/refamps

cd $1
for i in $1/*; do
    cd $i/AMPS/
    ls *.amp > /tmp/diffamps
    echo "Diffing bin $i to bin $2: "
    diff /tmp/diffamps /tmp/refamps
    
done;

