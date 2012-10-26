#!/bin/bash
#-------------------------------------------------------------------------
# File and Version Information:
# $Rev::                             $: revision of last commit
# $Author::                          $: author of last commit
# $Date::                            $: date of last commit
#
# Description:
#      generates complete wave set from template key files
#
#
# Author List:
#      Boris Grube          TUM            (original author)
#
#
#-------------------------------------------------------------------------


PARTICLE_DATA_TABLE="../../amplitude/particleDataTable.txt"
PARTICLE_DECAY_LIST="./5piParticleDecays.txt"
TEMPLATE_KEY_FILES="5pi3pi2pi.template.key 5pi4pi2pi2pi.template.key 5pi4pi3pi.template.key"


for TEMPLATE_KEY_FILE in ${TEMPLATE_KEY_FILES}
do

		DEST_DIR=${TEMPLATE_KEY_FILE%.template.key}
		if [[ ! -d ${DEST_DIR} ]]
				then 
				mkdir --verbose ${DEST_DIR}
		fi
		LOG_FILE=${DEST_DIR}/${TEMPLATE_KEY_FILE}.log
		CMD="generateWaveSet -p \"${PARTICLE_DATA_TABLE}\" -d \"${PARTICLE_DECAY_LIST}\" -f -o \"${DEST_DIR}\" -n -k ${TEMPLATE_KEY_FILE} &> ${LOG_FILE}"
		echo ${CMD}		
		eval ${CMD}
		# postprocess key files
		for KEY_FILE in ${DEST_DIR}/*.key
		do
				./5piPostProcessKeyFiles.sed "${KEY_FILE}"
		done
		# create backup copy of template key file
		cp --verbose --force "${TEMPLATE_KEY_FILE}" "${DEST_DIR}"
		
done



exit 0
