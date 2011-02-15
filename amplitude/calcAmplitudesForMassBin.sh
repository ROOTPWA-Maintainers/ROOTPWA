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
# File and Version Information:
# $Rev::                             $: revision of last commit
# $Author::                          $: author of last commit
# $Date::                            $: date of last commit
#
# Description:
#      calculates amplitudes for a single mass bin specified by
#      (1-based) index MASS_BIN_INDEX; runs locally as well as a
#      batch system; job arguments may be passed via command line or
#      environment variable
#
#      a typical command line to run this script locally looks like this:
#      for i in $(seq 1 50);\
#        do calcAmplitudesForMassBin.sh ${i} &> ${PWA_LOGS_DIR}/calcAmplitudes.${i}.log;\
#      done
#
#      !NOTE! this script does not yet handle accepted MC data
#
#      uses ROOTPWA environment variables
#      PWA_ENV_SET
#      ROOTPWA_BIN
#      PWA_KEYS_DIR
#      PWA_DATA_DIR
#
#
# Author List:
#      Boris Grube          TUM            (original author)
#
#
#-------------------------------------------------------------------------


# SGE setup
# if SGE is used as batch system, job arrays can be conveniently used
# to run this script; MASS_BIN_INDEX is set to SGE_TASK_ID
# define queue
#$ -l medium=TRUE,h_vmem=200M
# path for stdout files
#$ -o /afs/e18.ph.tum.de/user/bgrube/compass/pwa/4PionMuoProdPwa/logs
# merge sterr into stdout
#$ -j yes
# send mail when job is aborted, rescheduled, or suspended
#$ -M bgrube
#$ -m as
# enable path aliasing
#$ -cwd
# export all environment variables to job
#$ -V


echo ">>> info: called ${0} ${*}"


function usage {
    cat << EOF
usage: $(basename ${0}) [OPTIONS]... [MASS BIN INDEX]

creates amplitude and normalization files for real and MC data in MASS BIN INDEX
by running 'calcAmplitudes' and 'int'

options:
  -${KEY_PATTERN_OPT} '<PATTERN>'    glob pattern that defines set of key files to be processed
  -${MASS_BINS_DIR_OPT} <DIR>          path to directory with mass bins
  -${PDG_TABLE_OPT} <FILE>         path to PDG table file
  -${SYM_LIST_OPT} <FILE>         path to list with amplitudes to symmetrize
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
    path to list with amplitudes to symmetrize ... -${SYM_LIST_OPT} '${SYM_LIST}'

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
				local _AMP_FILE=${_AMP_DIR}/${_KEY_NAME/.key/.amp}
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


function runSymmetrization {
    local _AMP_DIR="${1}"
    local _SYM_LIST="${2}"
    local _SYM_DIR="${_AMP_DIR}/SYM"
    if [[ ! -d "${_SYM_DIR}" ]]
    then
				mkdir --parents --verbose "${_SYM_DIR}"
    fi
    echo
    echo ">>> info: linking all amplitude files in '${_AMP_DIR}' to '${_SYM_DIR}'"
    local _CURRENT_DIR=$(pwd)
    cd ${_SYM_DIR}
    for _AMP in ${_AMP_DIR}/*.amp
    do
				ln --symbolic --verbose "../$(basename ${_AMP})"
    done
    cd ${_CURRENT_DIR}
    if [[ -n "${_SYM_LIST}" && -s "${_SYM_LIST}" ]]
    then
				echo
				echo ">>> info: symmetrizing amplitudes in '${_AMP_DIR}' using '${_SYM_LIST}'"
				echo "    writing results to '${_SYM_DIR}'"
	# read list into arrays
				IFS_OLD=${IFS}
				IFS=$'\n'
				local _LINES=( $(cat ${_SYM_LIST}) )
				IFS=${IFS_OLD}
				local _AMPS1=()       # array of input amplitude 1 names
				local _AMPS2=()       # array of input amplitude 2 names
				local _SYM_AMPS=()    # array of output amplitude names
				local _PHASES=()      # array of relative phases
				local _AMP_RATIOS=()  # array of amplitude ratios
				for (( IDX=0; IDX<"${#_LINES[@]}"; IDX+=5 ))
				do
						_AMPS1=( ${_AMPS1[@]} ${_LINES[IDX]} )
						_AMPS2=( ${_AMPS2[@]} ${_LINES[IDX+1]} )
						_SYM_AMPS=( ${_SYM_AMPS[@]} ${_LINES[IDX+2]} )
						_PHASES=( ${_PHASES[@]} ${_LINES[IDX+3]} )
						_AMP_RATIOS=( ${_AMP_RATIOS[@]} ${_LINES[IDX+4]} )
				done
				_CURRENT_DIR=$(pwd)
				cd "${_SYM_DIR}"
				for (( IDX=0; IDX<"${#_SYM_AMPS[@]}"; IDX++ ))
				do
						local _CMD="${ROOTPWA_BIN}/addamp ${_AMPS1[IDX]} ${_AMPS2[IDX]} ${_SYM_AMPS[IDX]} ${_PHASES[IDX]} ${_AMP_RATIOS[IDX]}"
						echo "${_CMD}"
            eval ${_CMD}
						rm ${_AMPS1[IDX]} ${_AMPS2[IDX]}
				done
				cd ${_CURRENT_DIR}
    fi
    # deference third argument
    eval "${3}=${_SYM_DIR}"
}


function runInt {
    local _AMP_DIR="${1}"
    echo
    echo ">>> info: generating integrals for '${_AMP_DIR}'"
    local _CURRENT_DIR=$(pwd)
    # avoid problems with full path names in pwa2000
    cd "${_AMP_DIR}"
    local _CMD="${ROOTPWA_BIN}/int *.amp > norm.int"
    echo "${_CMD}"
    time eval ${_CMD}
    cd ${_CURRENT_DIR}
}


function runCalcIntegrals {
    local _AMP_DIR="${1}"
    echo
    echo ">>> info: generating integrals for '${_AMP_DIR}'"
    local _CMD="${ROOTPWA_BIN}/calcIntegrals -o ${_AMP_DIR}/norm.root ${_AMP_DIR}/*.amp"
    echo "${_CMD}"
    time eval ${_CMD}
}


function runPhaseSpaceGen {
    local _NMB_EVENTS=50000
    local _SEED=1234567890
    local _MC_FILE="${1}"
    local _MASS_BIN_NAME=$(basename "${_MC_FILE}" ".ps.evt")
    local _BIN_M_MIN=${_MASS_BIN_NAME%.*}
    local _BIN_M_MAX=${_MASS_BIN_NAME#*.}
    echo
    echo ">>> info: generating ${_NMB_EVENTS} phase space events for mass bin [${_BIN_M_MIN}, ${_BIN_M_MAX}] MeV/c^2"
    local _CMD="root -b -q -l \"runPhaseSpaceGen.C(${_NMB_EVENTS}, \\\"${_MC_FILE}\\\", ${_BIN_M_MIN}, ${_BIN_M_MAX}, ${_SEED})\""
    echo "${_CMD}"
    time eval ${_CMD}
}


# take arguments either from command line, environment, or use default
# (in this order of priority)
# option variables
KEY_PATTERN_OPT="k"
MASS_BINS_DIR_OPT="m"
PDG_TABLE_OPT="p"
SYM_LIST_OPT="s"
HELP_OPT="h"
# use default values, if variables are not defined in environment
if [[ -z "${KEY_PATTERN}" && ! -z "${PWA_KEYS_DIR}" ]]
then
    KEY_PATTERN="${PWA_KEYS_DIR}/*.key"
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
    MASS_BINS_DIR="${PWA_DATA_DIR}"
fi
if [[ -z "${PDG_TABLE}" ]]
then
    PDG_TABLE="${ROOTPWA}/amplitude/particleDataTable.txt"
fi
if [[ -z "${SYM_LIST}" ]]
then
    SYM_LIST=""
fi
# parse command line options
while getopts "${KEY_PATTERN_OPT}:${MASS_BINS_DIR_OPT}:${PDG_TABLE_OPT}:${HELP_OPT}" OPTION
do
    case ${OPTION} in
        ${KEY_PATTERN_OPT})   KEY_PATTERN=${OPTARG};;
        ${MASS_BINS_DIR_OPT}) MASS_BINS_DIR=${OPTARG};;
        ${PDG_TABLE_OPT})     PDG_TABLE=${OPTARG};;
        ${SYM_LIST_OPT})      SYM_LIST=${OPTARG};;
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
if [[ -z "${PWA_ENV_SET}" ]]
then
    echo "!!! error: PWA environment is not setup. source setupThisProgDir.sh."
    exit 1
fi


# convert all input paths to absolute paths
KEY_PATTERN=$(readlink --canonicalize-missing "${KEY_PATTERN}")
MASS_BINS_DIR=$(readlink --canonicalize-missing "${MASS_BINS_DIR}")
PDG_TABLE=$(readlink --canonicalize-missing "${PDG_TABLE}")
SYM_LIST=$(readlink --canonicalize-missing "${SYM_LIST}")


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
# directory for .amp files from real data
DT_AMP_DIR=${MASS_BIN_DIR}/AMPS
# path of real data file
DT_FILE=${MASS_BIN_DIR}/${MASS_BIN_NAME}.root
if [[ ! -s "${DT_FILE}" ]]
then
    echo "!!! error: data file '${DT_FILE}' does not exist. exiting."
    exit 1
fi
# directory for .amp files from Monte-Carlo
MC_AMP_DIR=${MASS_BIN_DIR}/PSPAMPS
# path of Monte-Carlo data file
MC_FILE=${MASS_BIN_DIR}/${MASS_BIN_NAME}.ps.root


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
CURRENT_DIR=$(pwd)


echo "------------------------------------------------------------"
echo ">>> info: processing real data for mass bin '${MASS_BIN_DIR}'"
# create directory if necessary
if [[ ! -d "${DT_AMP_DIR}" ]]
then
    mkdir --parents --verbose "${DT_AMP_DIR}"
fi
# generate .amp files for real data
runCalcAmplitudes "${DT_FILE}" "${DT_AMP_DIR}"
# perform symmetrization for real data
runSymmetrization "${DT_AMP_DIR}" "${SYM_LIST}" DUMMY
echo


echo "------------------------------------------------------------"
echo ">>> info: processing Monte-Carlo data for mass bin '${MASS_BIN_DIR}'"
# generate .amp files for Monte-Carlo data
if [[ ! -s "${MC_FILE}" ]]
then
    # generate phase space for normalization
    echo "??? warning: MC data file '${MC_FILE}' does not exist"
    #runPhaseSpaceGen "${MC_FILE}"
fi
# create directory if necessary
if [[ ! -d "${MC_AMP_DIR}" ]]
then
    mkdir --parents --verbose "${MC_AMP_DIR}"
fi
# generate .amp files for MC data
runCalcAmplitudes "${MC_FILE}" "${MC_AMP_DIR}"
# perform symmetrization for MC data
runSymmetrization "${MC_AMP_DIR}" "${SYM_LIST}" MC_SYM_AMP_DIR
# perform integration for MC data
#runInt "${MC_AMP_DIR}"
runCalcIntegrals "${MC_AMP_DIR}"
cd "${MC_SYM_AMP_DIR}"
if [[ -f ../norm.int ]]
then
		ln -s ../norm.int norm.int
fi
if [[ -f ../norm.root ]]
then
		ln -s ../norm.root norm.root
fi
cd -
echo


echo ">>> info: ${0} successfully finished on $(date)"
cd ${CURRENT_DIR}

exit 0
