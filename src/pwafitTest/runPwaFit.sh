#!/bin/bash

#
# performs PWA fit
# reads amplitudes in AMP_DIR and writes results to results.root
#
# assumes that pwafit is in path
#

#parameters
AMP_DIR="./amplitudes"
WAVE_LIST="./keyfiles/wavelist"
OUT_FILE="./result.root"
LOG_FILE="./pwafit.log"
BIN_MASS_MIN=2100
BIN_MASS_MAX=2140

# get absolute paths
AMP_DIR=$(readlink -f "${AMP_DIR}")
WAVE_LIST=$(readlink -f "${WAVE_LIST}")
OUT_FILE=$(readlink -f "${OUT_FILE}")
LOG_FILE=$(readlink -f "${LOG_FILE}")
CURRENT_DIR=$(pwd)

if [[ -s ${OUT_FILE} ]]
then
    rm ${OUT_FILE}
fi  

# run fit
echo ">>> ${0} started on $(date)"
echo ">>> fitting amplitude data in ${AMP_DIR} using wave list ${WAVE_LIST}"
cd ${AMP_DIR}
COMMAND="$ROOTPWA/src/pwafit -q -w ${WAVE_LIST} -o ${OUT_FILE} -r 2 -l ${BIN_MASS_MIN} -u ${BIN_MASS_MAX} -N &> ${LOG_FILE}"
echo "${COMMAND}"
time eval ${COMMAND}
cd ${CURRENT_DIR}
echo

echo ">>> ${0} finished on $(date)"
exit 0
