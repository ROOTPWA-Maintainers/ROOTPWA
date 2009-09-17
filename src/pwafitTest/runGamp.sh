#!/bin/bash

#
# generates amplitude files for all key files that match KEY_PATTERN
# reads DATA_FILE and writes results to AMP_DIR
# 
# existing files are _not_ regenerated
#
# assumes that gamp and int are in path
#

# parameters
KEY_PATTERN="./keyfiles/SET?/*.key"
DATA_FILE="./2100.2140.genbod.evt"
AMP_DIR="./amplitudes"
PDG_TABLE="./pdgTable.txt"

CURRENT_DIR=$(pwd)
NMB_OF_KEYS=$(ls -1 ${KEY_PATTERN} | wc -l)

# get absolute paths
DATA_FILE=$(readlink -f "${DATA_FILE}")
AMP_DIR=$(readlink -f "${AMP_DIR}")
PDG_TABLE=$(readlink -f "${PDG_TABLE}")
mkdir -p "${AMP_DIR}"

# process all key files
echo ">>> ${0} started on $(date)"
declare -i COUNT_KEY=0
for KEY_FILE in ${KEY_PATTERN}
do
    KEY_NAME=$(basename "${KEY_FILE}")
    AMP_FILE=${AMP_DIR}/${KEY_NAME/key/amp}
    (( ++COUNT_KEY ))
    # don't overwrite existing files
    if [[ -s ${AMP_FILE} ]]
    then
	#echo "??? warning: file ${AMP_FILE} already exists. skipping."
	: #noop
    else
	echo "............................................................"
	echo ">>> processing ${KEY_NAME} (${COUNT_KEY}/${NMB_OF_KEYS})"
        # avoid problems with absolute amplitude path names in gamp
	cd $(dirname ${KEY_FILE})
        COMMAND="gamp -P ${PDG_TABLE} ${KEY_NAME} < ${DATA_FILE} > ${AMP_FILE}"
	echo "${COMMAND}"
	time eval ${COMMAND}
	cd ${CURRENT_DIR}
	echo ">>> gamp finished $(date)"
	echo
    fi
done

# avoid problems with absolute path names in int
cd ${AMP_DIR}
# don't overwrite existing files
if [[ -s "norm.int" ]]
then
    echo "??? warning: file norm.int already exists. skipping."
    : #noop
else
    echo "------------------------------------------------------------"
    echo ">>> running integrator"
    COMMAND="int *.amp > norm.int"
    echo "${COMMAND}"
    time eval ${COMMAND}
    cd ${CURRENT_DIR}
fi
echo

echo ">>> ${0} finished on $(date)"
exit 0
