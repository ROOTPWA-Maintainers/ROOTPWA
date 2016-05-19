#!/bin/bash

[[ -z ${DESTINATION_DIR} ]] && DESTINATION_DIR="keyfiles"
[[ -z ${PARTICLE_DATA_TABLE} ]] && PARTICLE_DATA_TABLE="../../../particleData/particleDataTable.txt"
[[ -z ${TEMPLATE_KEY_FILES} ]] && TEMPLATE_KEY_FILES="template.key"
[[ -z ${WAVESET_FILES} ]] && WAVESET_FILES=""

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

# post-process the keyfiles
# * change mass dependence for sigma isobars
for KEY_FILE in ${DESTINATION_DIR}/*sigma0=*.key
do
	if [[ ! -f ${KEY_FILE} ]]
	then
		continue
	fi
	awk '
BEGIN {
	idx = 0
}
{
	if (pos = index($0, "name = \"sigma0\";")) {
		depths[++idx] = pos
	}
	if (idx > 0) {
		if (index($0, "} );") == (depths[idx]-2)) {
			pre = ""
			while (length(pre) < (depths[idx]-2)) pre = pre"  ";
			print pre"massDep : "
			print pre"{"
			print pre"  name = \"piPiSWaveAuMorganPenningtonKachaev\";"
			print pre"};"
			--idx
		}
	}
	print $0
}
' ${KEY_FILE} > ${KEY_FILE}.new
	mv -f ${KEY_FILE}.new ${KEY_FILE}
done

# post-process the keyfiles
# * f0(980) isobars with three different mass dependences:
#   usual Breit-Wigner, special Breit-Wigner and Flatte
for KEY_FILE in ${DESTINATION_DIR}/*f0_980_0=*.key
do
	if [[ ! -f ${KEY_FILE} ]]
	then
		continue
	fi
	awk '
BEGIN {
	idx = 0
}
{
	if (pos = index($0, "name = \"f0(980)0\";")) {
		depths[++idx] = pos
	}
	if (idx > 0) {
		if (index($0, "} );") == (depths[idx]-2)) {
			pre = ""
			while (length(pre) < (depths[idx]-2)) pre = pre"  ";
			print pre"massDep : "
			print pre"{"
			print pre"  name = \"f_0(980)\";"
			print pre"};"
			--idx
		}
	}
	print $0
}
' ${KEY_FILE} > `echo ${KEY_FILE} | sed -e 's/f0_980_0/f0_980_0bw/g'`
done
for KEY_FILE in ${DESTINATION_DIR}/*f0_980_0=*.key
do
	if [[ ! -f ${KEY_FILE} ]]
	then
		continue
	fi
	awk '
BEGIN {
	idx = 0
}
{
	if (pos = index($0, "name = \"f0(980)0\";")) {
		depths[++idx] = pos
	}
	if (idx > 0) {
		if (index($0, "} );") == (depths[idx]-2)) {
			pre = ""
			while (length(pre) < (depths[idx]-2)) pre = pre"  ";
			print pre"massDep : "
			print pre"{"
			print pre"  name = \"f_0(980)Flatte\";"
			print pre"};"
			--idx
		}
	}
	print $0
}
' ${KEY_FILE} > `echo ${KEY_FILE} | sed -e 's/f0_980_0/f0_980_0fl/g'`
done

if [[ ! -z "${WAVESET_FILES}" ]]
then
	# copy wavesets to destination dir, and create copies for the various
	# f0(980) mass dependences
	ALL_WAVESET_FILES=
	for WAVESET_FILE in ${WAVESET_FILES}
	do
		if [[ ! -e ${WAVESET_FILE} ]]
		then
			echo "Waveset file '${WAVESET_FILE}' does not exist."
		else
			WAVESET_FILE_BASENAME=`basename ${WAVESET_FILE}`
			if [[ ! -e ${DESTINATION_DIR}/${WAVESET_FILE_BASENAME} ]]
			then
				cp ${WAVESET_FILE} ${DESTINATION_DIR}/${WAVESET_FILE_BASENAME}
				ALL_WAVESET_FILES="${ALL_WAVESET_FILES} ${DESTINATION_DIR}/${WAVESET_FILE_BASENAME}"
			else
				echo "Waveset file '${WAVESET_FILE_BASENAME}' already exists in '${DESTINATION_DIR}'. Check manually that this file is correct."
				ALL_WAVESET_FILES="${ALL_WAVESET_FILES} ${DESTINATION_DIR}/${WAVESET_FILE_BASENAME}"
			fi
		fi
	done
	for WAVESET_FILE in ${WAVESET_FILES}
	do
		if [[ -e ${WAVESET_FILE} ]]
		then
			NEW_WAVESET_FILE="`basename ${WAVESET_FILE}`.f0980bw"
			if [[ ! -e ${DESTINATION_DIR}/${NEW_WAVESET_FILE} ]]
			then
				sed -e 's/f0_980_0=/f0_980_0bw=/g' ${WAVESET_FILE} > ${DESTINATION_DIR}/${NEW_WAVESET_FILE}
				if cmp -s ${WAVESET_FILE} ${DESTINATION_DIR}/${NEW_WAVESET_FILE}
				then
					rm -rf ${DESTINATION_DIR}/${NEW_WAVESET_FILE}
				else
					ALL_WAVESET_FILES="${ALL_WAVESET_FILES} ${DESTINATION_DIR}/${NEW_WAVESET_FILE}"
				fi
			else
				echo "Waveset file '${NEW_WAVESET_FILE}' already exists in '${DESTINATION_DIR}'. Check manually that this file is correct."
				ALL_WAVESET_FILES="${ALL_WAVESET_FILES} ${DESTINATION_DIR}/${NEW_WAVESET_FILE}"
			fi
			NEW_WAVESET_FILE="`basename ${WAVESET_FILE}`.f0980fl"
			if [[ ! -e ${DESTINATION_DIR}/${NEW_WAVESET_FILE} ]]
			then
				sed -e 's/f0_980_0=/f0_980_0fl=/g' ${WAVESET_FILE} > ${DESTINATION_DIR}/${NEW_WAVESET_FILE}
				if cmp -s ${WAVESET_FILE} ${DESTINATION_DIR}/${NEW_WAVESET_FILE}
				then
					rm -rf ${DESTINATION_DIR}/${NEW_WAVESET_FILE}
				else
					ALL_WAVESET_FILES="${ALL_WAVESET_FILES} ${DESTINATION_DIR}/${NEW_WAVESET_FILE}"
				fi
			else
				echo "Waveset file '${NEW_WAVESET_FILE}' already exists in '${DESTINATION_DIR}'. Check manually that this file is correct."
				ALL_WAVESET_FILES="${ALL_WAVESET_FILES} ${DESTINATION_DIR}/${NEW_WAVESET_FILE}"
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
