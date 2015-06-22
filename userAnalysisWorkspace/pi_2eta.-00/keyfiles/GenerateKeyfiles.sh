#!/bin/bash

[[ -z ${DESTINATION_DIR} ]] && DESTINATION_DIR="keyfiles"
[[ -z ${PARTICLE_DATA_TABLE} ]] && PARTICLE_DATA_TABLE="../../../particleData/particleDataTable.txt"
[[ -z ${WAVESET_FILES} ]] && WAVESET_FILES=""

TEMPLATE_KEY_FILES="etaeta.template.key pi-eta.template.key"

# if WAVESET_FILES is not empty, only keep those keyfiles actually used in one
# of the wavesets.

# check and create directory to contain the keyfiles
if [[ ! -d ${DESTINATION_DIR} ]]
then
	mkdir -p ${DESTINATION_DIR}
fi

# generate the keyfiles
for TEMPLATE_KEY_FILE in ${TEMPLATE_KEY_FILES}
do
	generateWaveSet -p ${PARTICLE_DATA_TABLE} -o ${DESTINATION_DIR} -k ${TEMPLATE_KEY_FILE}
done

if [[ ! -z "${WAVESET_FILES}" ]]
then
	# copy wavesets to destination dir
	ALL_WAVESET_FILES=
	for WAVESET_FILE in ${WAVESET_FILES}
	do
		if [[ ! -e ${WAVESET_FILE} ]]
		then
			echo "Waveset file '${WAVESET_FILE}' does not exist."
		else
			if [[ ! -e ${DESTINATION_DIR}/${WAVESET_FILE} ]]
			then
				cp ${WAVESET_FILE} ${DESTINATION_DIR}/${WAVESET_FILE}
				ALL_WAVESET_FILES="${ALL_WAVESET_FILES} ${DESTINATION_DIR}/${WAVESET_FILE}"
			else
				echo "Waveset file '${WAVESET_FILE}' already exists in '${DESTINATION_DIR}'. Check manually that this file is correct."
				ALL_WAVESET_FILES="${ALL_WAVESET_FILES} ${DESTINATION_DIR}/${WAVESET_FILE}"
			fi
		fi
	done
	# create list of all waves in wavesets (removing the thresholds)
	awk '{print $1".key"}' ${ALL_WAVESET_FILES} | sort -u > temp.waves.keep
	# create list of all keyfiles just created
	for i in `ls -1 ${DESTINATION_DIR}/*.key` ; do basename $i ; done | sort -u > temp.waves.all
	if [[ `diff temp.waves.all temp.waves.keep | grep "^>" | wc -l` -gt 0 ]]
	then
		echo "The wavesets contain waves for which no keyfiles were created."
		diff temp.waves.all temp.waves.keep | grep "^>" | awk '{print $2}'
		echo "No keyfiles will be removed!"
	else
		rm -rf `diff temp.waves.all temp.waves.keep | grep "^<" | awk "{print \"${DESTINATION_DIR}/\"\\\$2}"`
	fi
	rm -rf temp.waves.all temp.waves.keep
fi
