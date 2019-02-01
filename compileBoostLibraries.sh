#!/bin/bash

# first argument specifies number of jobs to run in parallel
# default is 1
NMB_JOBS=${1:-1}


if [[ ! "${BOOST_ROOT}" ]]
then
	echo '!!! error: ${BOOST_ROOT} environment variable is not set or empty. cannot compile BOOST libraries.'
	exit 1
fi


BOOST_ROOT_DIR=$(readlink --canonicalize ${BOOST_ROOT})
LOG_FILE=$(basename ${0})
LOG_FILE=${BOOST_ROOT_DIR}/${LOG_FILE%.*}.log
if [[ -d ${BOOST_ROOT_DIR} ]]
then
	echo ">>> compiling Boost libraries in '${BOOST_ROOT_DIR}'"
	cd ${BOOST_ROOT_DIR}
else
	echo "!!! error: '${BOOST_ROOT_DIR}' does not exist"
	exit 1
fi


JAM_CONFIG_FILE_SOURCE="${BOOST_ROOT_DIR}/tools/build/example/user-config.jam"
JAM_CONFIG_FILE_TARGET="${BOOST_ROOT_DIR}/tools/build/src/user-config.jam"
if [[ ! -e ${JAM_CONFIG_FILE_SOURCE} ]]
then
	echo "!!! error: '${JAM_CONFIG_FILE_SOURCE}' does not exist"
	exit 1
fi
if [[ ! -e ${JAM_CONFIG_FILE_TARGET} ]]
then
	cp ${JAM_CONFIG_FILE_SOURCE} ${JAM_CONFIG_FILE_TARGET}
fi
LINE='using mpi ;'
if ! grep --quiet --line-regexp "${LINE}" ${JAM_CONFIG_FILE_TARGET}
then
	echo "${LINE}" >> ${JAM_CONFIG_FILE_TARGET}
fi


function checkFail() {
	RETURN_CODE=${?}
	if [[ ${RETURN_CODE} != 0 ]]
	then
		echo "!!! error: command failure. ${1}"
		exit ${RETURN_CODE}
	fi
}

echo "    view result in ${LOG_FILE}"
echo "        ... bootstrapping"
./bootstrap.sh --prefix=$(pwd -P) >> ${LOG_FILE}  2>&1
checkFail "bootstrap failed; see ${LOG_FILE} for details"
echo "        ... linking headers"  # needed for Modular Boost repositories
./b2 -j ${NMB_JOBS} headers >> ${LOG_FILE} 2>&1
checkFail "linking headers failed; see ${LOG_FILE} for details"
echo "        ... building libraries; running ${NMB_JOBS} jobs in parallel"
./b2 -j ${NMB_JOBS} -a >> ${LOG_FILE} 2>&1
checkFail "building libraries failed; see ${LOG_FILE} for details"


if grep --quiet --line-regexp "The Boost C++ Libraries were successfully built!" ${LOG_FILE}
then
	echo "*** success: Boost libraries were successfully built"
	exit 0
else
	echo "!!! error: there were problems building the libraries; see ${LOG_FILE} for details"
	exit 1
fi
