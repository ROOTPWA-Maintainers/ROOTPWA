#!/bin/bash
#PBS -m be
#PBS -e localhost:/dev/null
#PBS -o localhost:/dev/null
#PBS -M sneubert@hamlet.e18.physik.tu-muenchen.de
#PBS -l nodes=1:Opteron:etch

# this script requires KEY to be set to the keydirectory!
# and BIN to be set to the mass bin directory

exec >/nfs/data/user/sneubert/logs/$PBS_JOBID.log
exec 2>/nfs/data/user/sneubert/logs/$PBS_JOBID.err

cd ~/scratch/COMPASS/PWA2000
source setup.sh
cd scripts/

export KEYDIR=$KEY
./buildBin.sh $BIN

cd
