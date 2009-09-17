#!/bin/bash

#
# runs pwafit with various tools of valgrind using test data set
#
# assumes that valgirnd is installed and in path
#

function usage {
    cat << EOF
Usage: $(basename ${0}) [OPTIONS]...

Runs pwafit with various tools of valgrind using test data set

OPTIONS:
  -${CALLGRIND_OPT}                if set, callgrind is run instead of memcheck
  -${HELP_OPT}                display this help and exit
EOF
    echo
}

# option variables
CALLGRIND_OPT="C"
HELP_OPT="h"

VALGRIND_TOOL="memcheck"
# parse command line options
while getopts "${CALLGRIND_OPT}${HELP_OPT}" OPTION
do
    case ${OPTION} in
        ${CALLGRIND_OPT}) VALGRIND_TOOL="callgrind";;
        ${HELP_OPT})      usage && exit 0;;
    esac
done
shift $((${OPTIND} - 1))  # remove used optins and leave remaining arguments in $*

# pwafit parameters
PROGRAM="../../bin/pwafit"
AMP_DIR="./amplitudes"
WAVE_LIST="./keyfiles/wavelist"
OUT_FILE="./result.root"
LOG_FILE="./pwafit_${VALGRIND_TOOL}.log"
BIN_MASS_MIN=2100
BIN_MASS_MAX=2140

# get absolute paths
AMP_DIR=$(readlink -f "${AMP_DIR}")
WAVE_LIST=$(readlink -f "${WAVE_LIST}")
OUT_FILE=$(readlink -f "${OUT_FILE}")
LOG_FILE=$(readlink -f "${LOG_FILE}")
CURRENT_DIR=$(pwd)

# run fit
echo ">>> ${0} started on $(date)"
echo ">>> running valgrind tool ${VALGRIND_TOOL} on ${PROGRAM}"
echo

# prepare command line
case ${VALGRIND_TOOL} in
    "memcheck")  # run memcheck tool 
	#VALGRIND_OPT="--leak-check=full --freelist-vol=500000000 --show-reachable=yes --verbose"
	VALGRIND_OPT="--leak-check=full --freelist-vol=500000000 --verbose"
	;;
    "callgrind")  # run callgrind tool
	VALGRIND_OPT="--tool=callgrind --simulate-cache=yes --verbose"
	;;
    *)
	echo "Unknown valgrind tool ${VALGRIND_TOOL}. Exiting."
	exit 1
	;;
esac
PROGRAM_OPT="-q -w ${WAVE_LIST} -o ${OUT_FILE} -r 2 -l ${BIN_MASS_MIN} -u ${BIN_MASS_MAX} -N"

# do it
cd ${AMP_DIR}
COMMAND="valgrind ${VALGRIND_OPT} ${PROGRAM} ${PROGRAM_OPT} &> ${LOG_FILE}"
echo "${COMMAND}"
time eval ${COMMAND}
cd ${CURRENT_DIR}
echo

echo ">>> ${0} finished on $(date)"
exit 0
