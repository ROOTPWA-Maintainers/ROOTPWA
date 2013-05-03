#!/bin/bash
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
#
# Description:
#      converts all .evt files in the given mass bin directories into ROOT files
#
#      uses PWA environment variable(s)
#      ROOTPWA_ENV_SET
#      ROOTPWA
#      ROOTPWA_DATA_DIR
#
#
# Author List:
#      Boris Grube          TUM            (original author)
#
#
#-------------------------------------------------------------------------


echo ">>> info: ${0} started on $(date)"
echo ">>> info: converting .evt files in ${ROOTPWA_DATA_DIR}"
echo
if [[ "${ROOTPWA_DATA_ENV_SET}" != "true" ]]
then
    echo "!!! error: ROOTPWA data environment is not setup. cannot convert .evt files. please source the setup script for the data environment first."
    exit 1
fi


#internal parameters
#REGENERATE_EVT_FILES="true"
REGENERATE_EVT_FILES="false"
PDG_TABLE="${ROOTPWA}/amplitude/particleDataTable.txt"


# find mass bin directories
MASS_BINS=$(find ${ROOTPWA_DATA_DIR} -type d -regex '.*/[0-9]+.[0-9]+' -printf '%f\n' | sort -n)
if [[ -z "${MASS_BINS}" ]]
then
    echo "!!! error: cannot find any mass bins in ${ROOTPWA_DATA_DIR}"
    exit 1
fi

# convert real data files
for MASS_BIN in ${MASS_BINS}
do
    DIR=${ROOTPWA_DATA_DIR}/${MASS_BIN}
    EVT_FILE=${DIR}/${MASS_BIN}.evt
    ROOT_FILE=${DIR}/${MASS_BIN}.root
    echo ">>> info: converting '${EVT_FILE}' to '${ROOT_FILE}'"
    CMD="root -b -q rootlogon.C convertEvtToTree.C+\(\\\"${EVT_FILE}\\\",\\\"${ROOT_FILE}\\\"\)"
    echo ${CMD}
    time eval ${CMD}
    echo
done

# convert MC data files
for MASS_BIN in ${MASS_BINS}
do
    #DIR=${ROOTPWA_DATA_DIR}/${MASS_BIN}/MC
    DIR=${ROOTPWA_DATA_DIR}/${MASS_BIN}
    EVT_FILE=${DIR}/${MASS_BIN}.ps.evt
    ROOT_FILE=${DIR}/${MASS_BIN}.ps.root
    echo ">>> info: converting '${EVT_FILE}' to '${ROOT_FILE}'"
    CMD="root -b -q rootlogon.C convertEvtToTree.C+\(\\\"${EVT_FILE}\\\",\\\"${ROOT_FILE}\\\"\)"
    echo ${CMD}
    time eval ${CMD}
    echo
done

# PWA2000 uses the four-vectors as they are given in the .evt files.
# this may lead to inconsistencies in the mass values used in the .evt
# files and those used for the isobars. it also causes troubles when
# comparing PWA2000 amplitudes with those generated with the new
# framework. enabling the regeneration of the .evt files rewrites them
# with energies according to the mass in the particle data table and
# leaving the three-momenta essentially untouched
if [[ "${REGENERATE_EVT_FILES}" == "true" ]]
then

    # rewrite real data files
    for MASS_BIN in ${MASS_BINS}
    do
				DIR=${ROOTPWA_DATA_DIR}/${MASS_BIN}
				EVT_FILE=${DIR}/${MASS_BIN}.evt
				ROOT_FILE=${DIR}/${MASS_BIN}.root
				echo ">>> info: rewriting '${EVT_FILE}' based on '${ROOT_FILE}'"
				CMD="root -b -q rootlogon.C convertTreeToEvt.C+\(\\\"${ROOT_FILE}\\\",\\\"${EVT_FILE}\\\",\\\"${PDG_TABLE}\\\"\)"
				echo ${CMD}
				time eval ${CMD}
				echo
    done

    # rewrite MC data files
    for MASS_BIN in ${MASS_BINS}
    do
				DIR=${ROOTPWA_DATA_DIR}/${MASS_BIN}/MC
				EVT_FILE=${DIR}/${MASS_BIN}.ps.evt
				ROOT_FILE=${DIR}/${MASS_BIN}.ps.root
				echo ">>> info: rewriting '${EVT_FILE}' based on '${ROOT_FILE}'"
				CMD="root -b -q rootlogon.C convertTreeToEvt.C+\(\\\"${ROOT_FILE}\\\",\\\"${EVT_FILE}\\\",\\\"${PDG_TABLE}\\\"\)"
				echo ${CMD}
				time eval ${CMD}
				echo
    done

fi

echo ">>> info: ${0} finished on $(date)"
exit 0
