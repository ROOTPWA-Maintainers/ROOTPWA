#!/usr/bin/bash
##########################################################################
#
#    Copyright 2010
#
#    This file is part of rootpwa
#
#    rootpwa is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    rootpwa is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################
#-------------------------------------------------------------------------
# File and Version Information:
# $Rev::                             $: revision of last commit
# $Author::                          $: author of last commit
# $Date::                            $: date of last commit
#
# Description:
#      sets up environment variables that are used by ROOTPWA scripts
#      to process the data; this script should be copied into the
#      directory where the mass bins are residing
#
#      uses ROOTPWA environment variables
#      PWA_ENV_SET
#      ROOTPWA
#
#      defines ROOTPWA environment variables
#      PWA_DATA_DIR
#      PWA_KEYS_DIR
#      PWA_WAVE_LIST
#      PWA_FITS_DIR
#      PWA_LOGS_DIR
#
#
# Author List:
#      Boris Grube          TUM            (original author)
#
#
#-------------------------------------------------------------------------


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
    echo "??? warning: ROOTPWA environment was not set. some paths may be invalid. source setupThisProg.sh in ROOTPWA directory."
fi

echo ">>> info: setting ROOTPWA data environment: "
echo "        PWA_DATA_DIR  = '${PWA_DATA_DIR}'"
echo "        PWA_KEYS_DIR  = '${PWA_KEYS_DIR}'"
echo "        PWA_WAVE_LIST = '${PWA_WAVE_LIST}'"
echo "        PWA_FITS_DIR  = '${PWA_FITS_DIR}'"
echo "        PWA_LOGS_DIR  = '${PWA_LOGS_DIR}'"
