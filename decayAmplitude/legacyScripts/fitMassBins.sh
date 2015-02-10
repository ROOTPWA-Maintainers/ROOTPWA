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
#      performs PWA fit of a number of consecutive mass bins specified
#      by (1-based) start index MASS_BIN_INDEX and the number of bins
#      NMB_BINS_TO_COMBINE
#
#      runs locally as well as on a batch system
#
#      job arguments may be passed via command line or environment
#      variable
#
#      a typical command line to run this script locally looks like this:
#      for i in $(seq 1 50);\
#        do fitMassBins.sh -r 2 -N ${i} &> ${ROOTPWA_LOGS_DIR}/fitMassBins.${i}.log;\
#      done
#
#      uses PWA environment variable(s)
#      ROOTPWA_ENV_SET
#      ROOTPWA_BIN
#      ROOTPWA_WAVE_LIST
#      ROOTPWA_DATA_DIR
#      ROOTPWA_FITS_DIR
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
# amplitude file format
AMP_FILE_EXT=.amp
#AMP_FILE_EXT=.root
# integral file name
INT_FILE_NAME=norm.int
#INT_FILE_NAME=norm.root
# fit result file extension
FIT_FILE_EXT=.result.root


# SGE setup
# if SGE is used as batch system, job arrays can be conveniently used
# to run this script; MASS_BIN_INDEX is set to SGE_TASK_ID
# define queue
##$ -l medium=TRUE,h_vmem=200M
# path for stdout files
##$ -o /afs/e18.ph.tum.de/user/bgrube/compass/pwa/4PionMuoProdPwa/logs
# merge sterr into stdout
#$ -j yes
# send mail when job aborted, rescheduled, or suspended
##$ -M bgrube
##$ -m as
# enable path aliasing
#$ -cwd
# export all environment variables to job
#$ -V


echo ">>> info: called ${0} ${*}"


function usage {
    cat << EOF
usage: $(basename ${0}) [OPTIONS]... [BATCH INDEX]

performs PWA fit for a number of consecutive mass bins, by running pwafit of the ROOTPWA package

options:
  -${WAVE_LIST_OPT} <FILE>         path to wave list file
  -${RANK_OPT} <RANK>         rank of spin density matrix
  -${NMB_BINS_TO_COMBINE_OPT} <#>            number of mass bins processed by one fit job
  -${MASS_BINS_DIR_OPT} <DIR>          path to directory with mass bins
  -${FITS_DIR_OPT} <DIR>          path to output directory for fit results
  -${NORM_OPT}                if set, amplitudes are normalized
  -${VERBOSE_OPT}                if set, debug output is printed
  -${TOLERANCE_OPT} <#>            minimizer tolerance
  -${START_VAL_OPT}                if set, start values are fixed to ${DEFAULT_START_VAL}
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
    if [[ -z "${NORM}" || "${NORM}" != "-N" ]]
    then
				NORM_STATE="false"
    else
				NORM_STATE="true"
    fi
    if [[ "${VERBOSE}" = "-q" ]]
    then
				VERBOSE_STATE="false"
    else
				VERBOSE_STATE="true"
    fi
    if [[ -z "${FIX_START_VAL}" || "${FIX_START_VAL}" != "-x${DEFAULT_START_VAL}" ]]
    then
				START_VAL_STATE="false"
    else
				START_VAL_STATE="true"
    fi
    cat << EOF
parameter values:
    path to wave list file ..................... -${WAVE_LIST_OPT} ${WAVE_LIST}
    rank of spin density matrix ................ -${RANK_OPT} ${RANK}
    number of mass bins to fit ................. -${NMB_BINS_TO_COMBINE_OPT} ${NMB_BINS_TO_COMBINE}
    path to directory with mass bins ........... -${MASS_BINS_DIR_OPT} ${MASS_BINS_DIR}
    path to output directory for fit results ... -${FITS_DIR_OPT} ${FITS_DIR}
    normalize amplitudes ....................... ${NORM_STATE}
    print debug output ......................... ${VERBOSE_STATE}
    minimizer tolerance ........................ ${TOLERANCE}
    fix start values ........................... ${START_VAL_STATE}

    mass bin batch index ....................... ${SGE_TASK_ID}
EOF
}


# take arguments either from command line, environment, or use default
# command line has highest priority
# option variables
WAVE_LIST_OPT="w"
RANK_OPT="r"
NMB_BINS_TO_COMBINE_OPT="n"
MASS_BINS_DIR_OPT="m"
FITS_DIR_OPT="o"
NORM_OPT="N"
VERBOSE_OPT="v"
TOLERANCE_OPT="t"
START_VAL_OPT="x"
HELP_OPT="h"
# use default values, if variables are not defined in environment
if [[ -z "${WAVE_LIST}" ]]
then
    WAVE_LIST="${ROOTPWA_WAVE_LIST}"
fi
if [[ -z "${RANK}" ]]
then
    RANK=1
fi
if [[ -z "${NMB_BINS_TO_COMBINE}" ]]
then
    NMB_BINS_TO_COMBINE=1
fi
if [[ -z "${MASS_BINS_DIR}" ]]
then
    MASS_BINS_DIR="${ROOTPWA_DATA_DIR}"
fi
if [[ -z "${FITS_DIR}" ]]
then
    FITS_DIR="${ROOTPWA_FITS_DIR}"
fi
if [[ -z "${NORM}" || "${NORM}" != "-N" ]]
then
    NORM=
fi
if [[ -z "${TOLERANCE}" ]]
then
    TOLERANCE="1e-10"
fi
if [[ -z "${DEFAULT_START_VAL}" ]]
then
    DEFAULT_START_VAL="0.01"
fi
VERBOSE="-q"
if [[ -z "${SGE_TASK_ID}" ]]
then
    SGE_TASK_ID=1
fi
# parse command line options
while getopts "${WAVE_LIST_OPT}:${RANK_OPT}:${NMB_BINS_TO_COMBINE_OPT}:${MASS_BINS_DIR_OPT}:${FITS_DIR_OPT}:${NORM_OPT}${VERBOSE_OPT}${TOLERANCE_OPT}:${START_VAL_OPT}${HELP_OPT}" OPTION
do
    case ${OPTION} in
        ${WAVE_LIST_OPT})           WAVE_LIST=${OPTARG};;
        ${RANK_OPT})                RANK=${OPTARG};;
        ${NMB_BINS_TO_COMBINE_OPT}) NMB_BINS_TO_COMBINE=${OPTARG};;
        ${MASS_BINS_DIR_OPT})       MASS_BINS_DIR=${OPTARG};;
        ${FITS_DIR_OPT})            FITS_DIR=${OPTARG};;
        ${NORM_OPT})                NORM="-N";;
        ${VERBOSE_OPT})             VERBOSE="";;
        ${TOLERANCE_OPT})           TOLERANCE=${OPTARG};;
        ${START_VAL_OPT})           FIX_START_VAL="-x${DEFAULT_START_VAL}";;
        ${HELP_OPT})                usage 0;;
    esac
done
shift $((${OPTIND} - 1))  # remove used options and leave remaining arguments in $*
if [[ -n "${1}" ]]
then
    SGE_TASK_ID="${1}"
fi
if [[ -z "${WAVE_LIST}" || -z "${RANK}" || -z "${NMB_BINS_TO_COMBINE}" || -z "${MASS_BINS_DIR}" || -z "${FITS_DIR}" ]]
then
    usage 1
fi
if [[ "${ROOTPWA_ENV_SET}" != "true" ]]
then
    echo "!!! error: ROOTPWA environment is not setup. please source the ROOTPWA setup script first."
    exit 1
fi


echo ">>> info: ${0} started on $(date)"
printPar && echo


# construct array of mass bins in ascending order
MASS_BINS=( $(find ${MASS_BINS_DIR} -type d -regex '.*/[0-9]+.[0-9]+' -printf '%f\n' | sort -n) )
# get index and mass range for bin sequence
declare -i BIN_SEQ_IDX_MIN=0
declare -i BIN_SEQ_IDX_MAX=0
(( BIN_SEQ_IDX_MIN = (SGE_TASK_ID - 1) * NMB_BINS_TO_COMBINE ))
(( BIN_SEQ_IDX_MAX = SGE_TASK_ID * NMB_BINS_TO_COMBINE - 1 ))
if (( BIN_SEQ_IDX_MIN > ${#MASS_BINS[@]} - 1 ))
then
    echo "!!! error: lower bin sequence index $(expr ${BIN_SEQ_IDX_MAX} + 1) out of range. exiting."
    exit 1
fi
if (( BIN_SEQ_IDX_MAX > ${#MASS_BINS[@]} - 1 ))
then
    echo "??? warning: upper bin sequence index $(expr ${BIN_SEQ_IDX_MAX} + 1) out of range. set to maximum possible value."
    (( BIN_SEQ_IDX_MAX = ${#MASS_BINS[@]} - 1 ))
fi
if (( BIN_SEQ_IDX_MIN > BIN_SEQ_IDX_MAX ))
then
    echo "!!! error: lower bin sequence index $(expr ${BIN_SEQ_IDX_MIN} + 1) larger than upper $(expr ${BIN_SEQ_IDX_MAX} + 1). exiting."
    exit 1
fi
BIN_SEQ_M_MIN=${MASS_BINS[${BIN_SEQ_IDX_MIN}]%.*}
BIN_SEQ_M_MAX=${MASS_BINS[${BIN_SEQ_IDX_MAX}]#*.}
# path of output ROOT file
OUT_FILE=${FITS_DIR}/${BIN_SEQ_M_MIN}.${BIN_SEQ_M_MAX}${FIT_FILE_EXT}
# check whether wave list file exists
if [[ ! -s "${WAVE_LIST}" ]]
then
    echo "!!! error: wave list file ${WAVE_LIST} does not exist! exiting."
    exit 1
fi


# loop over bins in batch
echo "------------------------------------------------------------"
echo ">>> info: performing PWA fit for ${NMB_BINS_TO_COMBINE} mass bins: index range = [$(expr ${BIN_SEQ_IDX_MIN} + 1), $(expr ${BIN_SEQ_IDX_MAX} + 1)], mass range = [${BIN_SEQ_M_MIN}, ${BIN_SEQ_M_MAX}]"
echo
for (( IDX=BIN_SEQ_IDX_MIN; IDX<=BIN_SEQ_IDX_MAX; IDX++ ))
do
    MASS_BIN_NAME=${MASS_BINS[${IDX}]}
    # path to mass bin data
    AMP_DIR="${MASS_BINS_DIR}/${MASS_BIN_NAME}/${DT_AMP_DIR_NAME}"
    BIN_M_MIN=${MASS_BIN_NAME%.*}
    BIN_M_MAX=${MASS_BIN_NAME#*.}
    # check whether amplitude directory exists
    if [[ ! -d "${AMP_DIR}" ]]
    then
				echo "!!! error: amplitude directory '${AMP_DIR}' does not exist! exiting."
				exit 1
    fi
    PS_AMP_DIR="${MASS_BINS_DIR}/${MASS_BIN_NAME}/${PS_AMP_DIR_NAME}"
    # check whether directory with phase-space amplitudes exists
    if [[ ! -d "${PS_AMP_DIR}" ]]
    then
				echo "!!! error: directory '${PS_AMP_DIR}' with phase-space amplitudes does not exist! exiting."
				exit 1
    fi
    PSACC_AMP_DIR="${MASS_BINS_DIR}/${MASS_BIN_NAME}/${PSACC_AMP_DIR_NAME}"
    # check whether directory with accepted phase space amplitudes exists
    if [[ ! -d "${PSACC_AMP_DIR}" ]]
    then
				echo "??? warning: directory '${PSACC_AMP_DIR}' with accepted phase-space amplitudes does not exist. using phase space integrals in '${PS_AMP_DIR}'."
				PSACC_AMP_DIR=${PS_AMP_DIR}
    fi
    # set path to phase-space integrals
    PS_INT_FILE="${PS_AMP_DIR}/${INT_FILE_NAME}"
    PSACC_INT_FILE="${PSACC_AMP_DIR}/${INT_FILE_NAME}"
    # run pwafit
    NMB_NORM_EVENTS=$(sed '2q;d' ${PS_INT_FILE})
    echo "............................................................"
    echo ">>> info: starting PWA fit for bin $(expr ${IDX} + 1) with mass = [${BIN_M_MIN}, ${BIN_M_MAX}]"
    echo ">>> info: using amplitude files in ${AMP_DIR}"
    #ls -l --dereference ${AMP_DIR}/*.${AMP_FILE_EXT}
    echo ">>> info: using wavelist ${WAVE_LIST}"
    CMD="${ROOTPWA_BIN}/pwafit -l ${BIN_M_MIN} -u ${BIN_M_MAX} -w ${WAVE_LIST} -d ${AMP_DIR} -o ${OUT_FILE} ${NORM} -n ${PS_INT_FILE} -A ${NMB_NORM_EVENTS} -a ${PSACC_INT_FILE} -r ${RANK} ${VERBOSE} -t ${TOLERANCE} ${FIX_START_VAL}"
    echo "${CMD}"
    time eval ${CMD}
    echo
done


echo ">>> info: ${0} successfully finished on $(date)"
exit 0
