#!/bin/bash
#
# Script to get workdir from a list of directories
# SGE ArrayJob index of the task we are in
# Input: $1=LIST $2=JobID

LIST=$1
JOBID=$2

MYDIR=`head -n $JOBID $LIST | tail -n 1`
echo $MYDIR