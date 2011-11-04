#!/bin/bash
#///////////////////////////////////////////////////////////////////////////
#//
#//    Copyright 2009 Sebastian Neubert
#//
#//    This file is part of rootpwa
#//
#//    rootpwa is free software: you can redistribute it and/or modify
#//    it under the terms of the GNU General Public License as published by
#//    the Free Software Foundation, either version 3 of the License, or
#//    (at your option) any later version.
#//
#//    rootpwa is distributed in the hope that it will be useful,
#//    but WITHOUT ANY WARRANTY; without even the implied warranty of
#//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#//    GNU General Public License for more details.
#//
#//    You should have received a copy of the GNU General Public License
#//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
#//
#///////////////////////////////////////////////////////////////////////////

NMB_EVENTS=20000
GEN_CONFIG_FILE=${ROOTPWA}/generators/3pinDiffractive2008.conf

# construct array of mass bins in ascending order
MASS_BINS=( $(find ${ROOTPWA_DATA_DIR} -type d -regex '.*/[0-9]+.[0-9]+' -printf '%f\n' | sort -n) )
MASS_BIN_IDX=$(( MASS_BIN_IDX - 1 ))  # SGE task IDs are 1 based

for (( i=0; i<"${#MASS_BINS[@]}"; i++ ))
do
    # get mass range for bin
    MASS_BIN=${MASS_BINS[i]}
    MASS_BIN_M_MIN=${MASS_BIN%.*}
    MASS_BIN_M_MAX=${MASS_BIN#*.}
    MASS_BIN_WIDTH=$(( MASS_BIN_M_MAX-MASS_BIN_M_MIN ))
    MASS_BIN_MC_DIR="${ROOTPWA_DATA_DIR}/${MASS_BIN}/MC"
    if [[ ! -d "${MASS_BIN_MC_DIR}" ]]
    then
				mkdir --parents --verbose ${MASS_BIN_MC_DIR}
    fi
    echo ">>> info: generating phase space data for mass bin ${MASS_BIN_M_MIN}.${MASS_BIN_M_MAX}"
    CMD="genpw -n ${NMB_EVENTS} -o ${MASS_BIN_MC_DIR}/${MASS_BIN_M_MIN}.${MASS_BIN_M_MAX}.genbod.root -r ${GEN_CONFIG_FILE} -M ${MASS_BIN_M_MIN} -B ${MASS_BIN_WIDTH}"
    echo ${CMD}
    #time eval ${CMD}

    # convert .evt file into ROOTPWA .root format
    EVT_FILE=${MASS_BIN_MC_DIR}/${MASS_BIN_M_MIN}.${MASS_BIN_M_MAX}.genbod.evt
    ROOT_FILE=${ROOTPWA_DATA_DIR}/${MASS_BIN}/${MASS_BIN_M_MIN}.${MASS_BIN_M_MAX}.ps.root
    echo ">>> info: converting '${EVT_FILE}' to '${ROOT_FILE}'"
    CMD="root -b -q ${ROOTPWA}/amplitude/rootlogon.C ${ROOTPWA}/amplitude/convertEvtToTree.C+\(\\\"${EVT_FILE}\\\",\\\"${ROOT_FILE}\\\"\)"
    echo ${CMD}
    time eval ${CMD}
    echo

done
