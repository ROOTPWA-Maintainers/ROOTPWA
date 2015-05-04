#!/bin/bash


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

echo "    view result in ${LOG_FILE}"
./bootstrap.sh --with-libraries=mpi,python,timer --prefix=. &> ${LOG_FILE}
./b2 -a >> ${LOG_FILE} 2>&1


exit 0
