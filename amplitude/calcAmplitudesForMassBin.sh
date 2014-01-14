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
#      calculates amplitudes for a single mass bin specified by
#      (1-based) index MASS_BIN_INDEX
#
#      runs locally as well as on a batch system
#
#      job arguments may be passed via command line or environment
#      variable
#
#      a typical command line to run this script locally looks like this:
#      for i in $(seq 1 50);\
#        do calcAmplitudesForMassBin.sh ${i} &> ${ROOTPWA_LOGS_DIR}/calcAmplitudes.${i}.log;\
#      done
#
#      uses ROOTPWA environment variables
#      ROOTPWA_ENV_SET
#      ROOTPWA_BIN
#      ROOTPWA_KEYS_DIR
#      ROOTPWA_DATA_DIR
#
#
# Author List:
#      Boris Grube          TUM            (original author)
#
#
#-------------------------------------------------------------------------


# file layout parameters
# subdirectory names for amplitudes
DT_AMP_DIR_NAME=AMPS
PS_AMP_DIR_NAME=PSPAMPS
PSACC_AMP_DIR_NAME=ACCAMPS
# extensions of event data files; naming scheme is <mass bin><extension>
DT_FILE_EXT=.evt
PS_FILE_EXT=.genbod.evt
#PS_FILE_EXT=.ps.evt
PSACC_FILE_EXT=.acc.evt
#DT_FILE_EXT=.root
#PS_FILE_EXT=.ps.root
#PSACC_FILE_EXT=.acc.root
# amplitude file format
AMP_FILE_EXT=.amp
#AMP_FILE_EXT=.root
# integral file name
INT_FILE_NAME=norm.int
#INT_FILE_NAME=norm.root


# SGE setup
# if SGE is used as batch system, job arrays can be conveniently used
# to run this script; MASS_BIN_INDEX is set to SGE_TASK_ID
# define queue
##$ -l medium=TRUE,h_vmem=200M
# path for stdout files
##$ -o /afs/e18.ph.tum.de/user/bgrube/compass/pwa/4PionMuoProdPwa/logs
# merge sterr into stdout
#$ -j yes
# send mail when job is aborted, rescheduled, or suspended
##$ -M bgrube
##$ -m as
# enable path aliasing
#$ -cwd
# export all environment variables to job
#$ -V


echo ">>> info: called ${0} ${*}"


function usage {
    cat << EOF
usage: $(basename ${0}) [OPTIONS]... [MASS BIN INDEX]

creates amplitude and integral files for real data, phase-space MC,
and accepted phase-space MC in mass bin directory given by 1-based
index [MASS BIN INDEX] by running 'calcAmplitudes' and 'calcIntegrals'

options:
  -${KEY_PATTERN_OPT} '<PATTERN>'    glob pattern that defines set of key files to be processed
  -${MASS_BINS_DIR_OPT} <DIR>          path to directory with mass bins
  -${PDG_TABLE_OPT} <FILE>         path to PDG table file
  -${HELP_OPT}                display this help and exit
EOF
echo
printPar
if [[ "${#}" -eq 1 ]]
then
		exit ${1}
else
		exit 1
fi
}


function printPar {
    cat << EOF
parameter values:
    glob pattern that defines set of key files ... -${KEY_PATTERN_OPT} '${KEY_PATTERN}'
    path to directory with mass bins ............. -${MASS_BINS_DIR_OPT} '${MASS_BINS_DIR}'
    path to PDG table file ....................... -${PDG_TABLE_OPT} '${PDG_TABLE}'

    mass bin index ............................... ${MASS_BIN_INDEX}
EOF
}


function runCalcAmplitudes {
    local _DT_FILE="${1}"
    local _AMP_DIR="${2}"
    local _NMB_OF_KEYS=$(eval ls -1 ${KEY_PATTERN} | wc -l)
    # process all key files
    declare -i _COUNT_KEY=0
    for _KEY_FILE in ${KEY_PATTERN}
    do
				local _KEY_NAME=$(basename ${_KEY_FILE})
				local _AMP_FILE=${_AMP_DIR}/${_KEY_NAME/.key/${AMP_FILE_EXT}}
				(( ++_COUNT_KEY ))
				echo
				echo ">>> info: generating amplitudes for '${_KEY_NAME}' [${_COUNT_KEY}/${_NMB_OF_KEYS}]"
        # don't overwrite existing files
				if [[ -s "${_AMP_FILE}" ]]
				then
						echo "??? warning: file '${_AMP_FILE}' already exists. skipping."
				else
            local _CMD="${ROOTPWA_BIN}/calcAmplitudes -k ${_KEY_FILE} -p ${PDG_TABLE} -o ${_AMP_FILE} ${_DT_FILE}"
						echo "${_CMD}"
						time eval ${_CMD}
				fi
    done
}


function runInt {
    local _AMP_DIR="${1}"
    echo
    echo ">>> info: generating integrals for '${_AMP_DIR}'"
    local _CURRENT_DIR=$(pwd)
    # avoid problems with full path names in pwa2000
    cd "${_AMP_DIR}"
    local _CMD="${ROOTPWA_BIN}/int *.amp > ${INT_FILE_NAME}"  # works only with .amp and .int files
    echo "${_CMD}"
    time eval ${_CMD}
    cd ${_CURRENT_DIR}
}


function runCalcIntegrals {
    local _AMP_DIR="${1}"
    echo
    echo ">>> info: generating integrals for '${_AMP_DIR}'"
    local _CMD="${ROOTPWA_BIN}/calcIntegrals -o ${_AMP_DIR}/${INT_FILE_NAME} ${_AMP_DIR}/*${AMP_FILE_EXT}"
    echo "${_CMD}"
    time eval ${_CMD}
}


function runPhaseSpaceGen {
    local _NMB_EVENTS=50000
    local _SEED=1234567890
    local _PS_FILE="${1}"
    local _MASS_BIN_NAME=$(basename "${_PS_FILE}" "${PS_FILE_EXT}")
    local _BIN_M_MIN=${_MASS_BIN_NAME%.*}
    local _BIN_M_MAX=${_MASS_BIN_NAME#*.}
    echo
    echo ">>> info: generating ${_NMB_EVENTS} phase space events for mass bin [${_BIN_M_MIN}, ${_BIN_M_MAX}] MeV/c^2"
    local _CMD="root -b -q -l \"runPhaseSpaceGen.C(${_NMB_EVENTS}, \\\"${_PS_FILE}\\\", ${_BIN_M_MIN}, ${_BIN_M_MAX}, ${_SEED})\""
    echo "${_CMD}"
    time eval ${_CMD}
}


# take arguments either from command line, environment, or use default
# (in this order of priority)
# option variables
KEY_PATTERN_OPT="k"
MASS_BINS_DIR_OPT="m"
PDG_TABLE_OPT="p"
HELP_OPT="h"
# use default values, if variables are not defined in environment
if [[ -z "${KEY_PATTERN}" && ! -z "${ROOTPWA_KEYS_DIR}" ]]
then
    KEY_PATTERN="${ROOTPWA_KEYS_DIR}/*.key"
fi
if [[ -z "${MASS_BIN_INDEX}" ]]
then
    MASS_BIN_INDEX="1"
fi
# override MASS_BIN_INDEX with SGE_TASK_ID
if [[ -n "${SGE_TASK_ID}" ]]
then
    MASS_BIN_INDEX="${SGE_TASK_ID}"
fi
if [[ -z "${MASS_BINS_DIR}" ]]
then
    MASS_BINS_DIR="${ROOTPWA_DATA_DIR}"
fi
if [[ -z "${PDG_TABLE}" ]]
then
    PDG_TABLE="${ROOTPWA}/amplitude/particleDataTable.txt"
fi
# parse command line options
while getopts "${KEY_PATTERN_OPT}:${MASS_BINS_DIR_OPT}:${PDG_TABLE_OPT}:${HELP_OPT}" OPTION
do
    case ${OPTION} in
        ${KEY_PATTERN_OPT})   KEY_PATTERN=${OPTARG};;
        ${MASS_BINS_DIR_OPT}) MASS_BINS_DIR=${OPTARG};;
        ${PDG_TABLE_OPT})     PDG_TABLE=${OPTARG};;
        ${HELP_OPT})          usage 0;;
    esac
done
shift $((${OPTIND} - 1))  # just leave remaining arguments in $*
if [[ -n "${1}" ]]
then
    MASS_BIN_INDEX="${1}"
fi
if [[ -z "${KEY_PATTERN}" || -z "${MASS_BINS_DIR}" || -z "${PDG_TABLE}" || -z "${MASS_BIN_INDEX}" ]]
then
    usage 1
fi
if [[ "${ROOTPWA_ENV_SET}" != "true" ]]
then
    echo "!!! error: ROOTPWA environment is not setup. please source the ROOTPWA setup script first."
    exit 1
fi


# convert all input paths to absolute paths
KEY_PATTERN=$(readlink --canonicalize-missing "${KEY_PATTERN}")
MASS_BINS_DIR=$(readlink --canonicalize-missing "${MASS_BINS_DIR}")
PDG_TABLE=$(readlink --canonicalize-missing "${PDG_TABLE}")


echo ">>> info: ${0} started on $(date)"
printPar && echo


# construct array of mass bins in ascending order
MASS_BINS=( $(find ${MASS_BINS_DIR} -type d -regex '.*/[0-9]+.[0-9]+' -printf '%f\n' | sort -n) )
if [[ -z "${MASS_BINS}" ]]
then
    echo "!!! error: cannot find any mass bins in '${MASS_BINS_DIR}'"
    exit 1
fi
# find right mass bin
MASS_BIN_DIR=${MASS_BINS[$(expr ${MASS_BIN_INDEX} - 1)]}
if [[ -z "${MASS_BIN_DIR}" ]]
then
    echo "!!! error: cannot find mass bin with index ${MASS_BIN_INDEX}"
    exit 1
fi
MASS_BIN_DIR=${MASS_BINS_DIR}/${MASS_BIN_DIR}
if [[ ! -d "${MASS_BIN_DIR}" ]]
then
    echo "!!! error: mass bin directory '${MASS_BIN_DIR}' does not exist"
    exit 1
fi


# define mass bin file structure
# name of mass bin
MASS_BIN_NAME=$(basename ${MASS_BIN_DIR})
# directory for amplitude files from real data
DT_AMP_DIR=${MASS_BIN_DIR}/${DT_AMP_DIR_NAME}
# path of real data file
DT_FILE=${MASS_BIN_DIR}/${MASS_BIN_NAME}${DT_FILE_EXT}
if [[ ! -e "${DT_FILE}" ]]
then
    echo "!!! error: data file '${DT_FILE}' does not exist. exiting."
    exit 1
fi
# directory for amplitude files from Monte-Carlo
PS_AMP_DIR=${MASS_BIN_DIR}/${PS_AMP_DIR_NAME}
# path of Monte-Carlo data file
PS_FILE=${MASS_BIN_DIR}/${MASS_BIN_NAME}${PS_FILE_EXT}
# directory for amplitude files from Monte-Carlo
PSACC_AMP_DIR=${MASS_BIN_DIR}/${PSACC_AMP_DIR_NAME}
# path of Monte-Carlo data file
PSACC_FILE=${MASS_BIN_DIR}/${MASS_BIN_NAME}${PSACC_FILE_EXT}


# check whether necessary file and directories are there
if [[ "$(eval ls -1 ${KEY_PATTERN} | wc -l)" == "0" ]]
then
    echo "!!! error: no key files match glob pattern '${KEY_PATTERN}'. exiting."
    exit 1
fi
if [[ ! -s "${PDG_TABLE}" ]]
then
    echo "!!! error: PDG table file '${PDG_TABLE}' does not exist. exiting."
    exit 1
fi


# calculate amplitudes
echo ">>> info: starting amplitude calculation for mass bin ${MASS_BIN_INDEX} in '${MASS_BIN_DIR}' using waveset(s) '${KEY_PATTERN}'"
echo ">>> info: $(eval ls -1 ${KEY_PATTERN} | wc -l) key files:"
eval ls -l ${KEY_PATTERN}
echo
echo ">>> info: ${0} started on $(date)"


echo "------------------------------------------------------------"
echo ">>> info: processing real data for mass bin '${MASS_BIN_DIR}'"
# create directory if necessary
if [[ ! -d "${DT_AMP_DIR}" ]]
then
    mkdir --parents --verbose "${DT_AMP_DIR}"
fi
# generate amplitude files for real data
runCalcAmplitudes "${DT_FILE}" "${DT_AMP_DIR}"
echo


echo "------------------------------------------------------------"
echo ">>> info: processing phase-space Monte Carlo data for mass bin '${MASS_BIN_DIR}'"
# generate amplitude files for phase-space Monte Carlo data
if [[ ! -s "${PS_FILE}" ]]
then
    echo "??? warning: phase-space MC data file '${PS_FILE}' does not exist"
    # generate phase space for normalization integral
    #runPhaseSpaceGen "${PS_FILE}"
fi
# create directory if necessary
if [[ ! -d "${PS_AMP_DIR}" ]]
then
    mkdir --parents --verbose "${PS_AMP_DIR}"
fi
# generate amplitude files for phase-space MC data
runCalcAmplitudes "${PS_FILE}" "${PS_AMP_DIR}"
# perform integration for phase-space MC data
#runInt "${PS_AMP_DIR}"
runCalcIntegrals "${PS_AMP_DIR}"
echo


echo "------------------------------------------------------------"
echo ">>> info: processing accepted phase-space Monte-Carlo data for mass bin '${MASS_BIN_DIR}'"
# generate amplitude files for accepted phase-space Monte-Carlo data
if [[ ! -s "${PSACC_FILE}" ]]
then
    echo "??? warning: accepted phase-space MC data file '${PSACC_FILE}' does not exist"
else
    # create directory if necessary
		if [[ ! -d "${PSACC_AMP_DIR}" ]]
		then
				mkdir --parents --verbose "${PSACC_AMP_DIR}"
		fi
    # generate amplitude files for accepted phase-space MC data
		runCalcAmplitudes "${PSACC_FILE}" "${PSACC_AMP_DIR}"
    # perform integration for accepted phase-space MC data
    #runInt "${PSACC_AMP_DIR}"
		runCalcIntegrals "${PSACC_AMP_DIR}"
fi
echo


echo ">>> info: ${0} successfully finished on $(date)"

exit 0
