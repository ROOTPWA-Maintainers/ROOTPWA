#!/usr/bin/bash

# setup data directories
export PWA_DATA_DIR=$(readlink --canonicalize $(dirname ${BASH_SOURCE}))
export PWA_KEYS_DIR="${ROOTPWA}/keyfiles/key3pi/SET?_new"
export PWA_WAVE_LIST="${ROOTPWA}/keyfiles/key3pi/wavelist"
export PWA_FITS_DIR="${PWA_DATA_DIR}/fits"
export PWA_LOGS_DIR="${PWA_DATA_DIR}/logs"

# make sure directories exist
if [[ ! -d ${PWA_FITS_DIR} ]]
then
    mkdir --parents --verbose ${PWA_FITS_DIR}
fi
if [[ ! -d ${PWA_LOGS_DIR} ]]
then
    mkdir --parents --verbose ${PWA_LOGS_DIR}
fi

if [[ "${PWA_ENV_SET}" == "TRUE" ]]
then
    # set flag that indicates that environment was set successfully
    export PWA_DATA_ENV_SET="TRUE"
else
    echo "??? warning: ROOTPWA environment was not set. some paths may be invalid."
fi

echo ">>> info: setting PWA data environment: "
echo "        PWA_DATA_DIR  = '${PWA_DATA_DIR}'"
echo "        PWA_KEYS_DIR  = '${PWA_KEYS_DIR}'"
echo "        PWA_WAVE_LIST = '${PWA_WAVE_LIST}'"
echo "        PWA_FITS_DIR  = '${PWA_FITS_DIR}'"
echo "        PWA_LOGS_DIR  = '${PWA_LOGS_DIR}'"
