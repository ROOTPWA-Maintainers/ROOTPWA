#!/bin/bash

# first argument specifies number of jobs to run in parallel
# default is 1
NMB_JOBS=${1:-1}
# if second argument is given and not "TRUE", all output is written to stdout instead of log file
WRITE_LOG_FILE=${2:-"TRUE"}


# check for environment variable that specifies where to find Boost
if [[ ! "${BOOST_ROOT}" ]]
then
	echo '!!! error: ${BOOST_ROOT} environment variable is not set or empty. cannot compile BOOST libraries.'
	exit 1
fi


# check whether Boost directory exists and define log file location
BOOST_ROOT_DIR=$(readlink --canonicalize ${BOOST_ROOT})
LOG_FILE=$(basename ${0})
LOG_FILE=${BOOST_ROOT_DIR}/${LOG_FILE%.*}.log
if [[ -d "${BOOST_ROOT_DIR}" ]]
then
	echo ">>> compiling Boost libraries in '${BOOST_ROOT_DIR}'"
	cd "${BOOST_ROOT_DIR}"
else
	echo "!!! error: '${BOOST_ROOT_DIR}' does not exist"
	exit 1
fi


# make sure Boost MPI module is included in build
JAM_CONFIG_FILE_SOURCE="${BOOST_ROOT_DIR}/tools/build/example/user-config.jam"
JAM_CONFIG_FILE_TARGET="${BOOST_ROOT_DIR}/tools/build/src/user-config.jam"
if [[ ! -e "${JAM_CONFIG_FILE_SOURCE}" ]]
then
	echo "!!! error: '${JAM_CONFIG_FILE_SOURCE}' does not exist"
	exit 1
fi
if [[ ! -e "${JAM_CONFIG_FILE_TARGET}" ]]
then
	cp -v "${JAM_CONFIG_FILE_SOURCE}" "${JAM_CONFIG_FILE_TARGET}"
fi
LINE='using mpi ;'
if ! grep --quiet --line-regexp "${LINE}" "${JAM_CONFIG_FILE_TARGET}"
then
	echo "${LINE}" >> "${JAM_CONFIG_FILE_TARGET}"
fi


# helper function that checks return code of preceding command
function checkFail() {
	RETURN_CODE=${?}
	if [[ ${RETURN_CODE} != 0 ]]
	then
		echo "!!! error: command failure. ${1}"
		exit ${RETURN_CODE}
	fi
}


# compile libraries
# use file descriptor 3 for stdout and stderr
ERR_MSG_SUFFIX=""
if [[ ${WRITE_LOG_FILE} == "TRUE" ]]
then
	# connect file descriptor 3 to log file, append
	echo "    view result in ${LOG_FILE}"
	ERR_MSG_SUFFIX="; see ${LOG_FILE} for details"
	exec 3>> "${LOG_FILE}"
else
	# connect file descriptor 3 to stdout
	exec 3>&1
fi
echo "        ... bootstrapping"
./bootstrap.sh --prefix="$(pwd -P)" >&3 2>&1
checkFail "bootstrap failed"${ERR_MSG_SUFFIX}
echo "        ... linking headers"  # needed for Modular Boost repositories
./b2 -j ${NMB_JOBS} headers >&3 2>&1
checkFail "linking headers failed"${ERR_MSG_SUFFIX}
echo "        ... building libraries; running ${NMB_JOBS} jobs in parallel"
./b2 -j ${NMB_JOBS} -a >&3 2>&1
checkFail "building libraries failed"${ERR_MSG_SUFFIX}
# close file descriptor 3
exec 3>&-


echo "*** success: Boost libraries were successfully built"
exit 0
