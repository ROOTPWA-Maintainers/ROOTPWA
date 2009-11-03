#!/bin/bash
# start 4 fitting jobs simultaneusly!

#IDS=( 1 2 3 4 )
IDS=( 5 6 7 8 )
#IDS=( 9 10 11 12 )

for i in 0 1 2 3 ; do
export SGE_TASK_ID=${IDS[$i]}
echo "Starting Job $SGE_TASK_ID"
./rootfitBatch.sge > batch$SGE_TASK_ID.log &
done