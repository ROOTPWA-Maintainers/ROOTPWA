# this script tests the status of the Kpipi PWA workflow chain
# the variables in set_workspace_var.sh must be set correctly
# author: Promme@web.de ; jasinski@kph.uni-mainz.de
# worked with SVN version 290

#!/bin/bash

# test if ${ROOTPWA} is set (set it in .bashrc for example)
if [ -d ${ROOTPWA} ]; 
then
	echo -e "\n ROOTPWA variable is set to ${ROOTPWA} \n"
else
	echo -e "\E[37;31m \n ROOTPWA variable is not set correctly! "; tput sgr0
	echo -e "\E[37;31m Please set it first! "; tput sgr0
	echo -e "\n Aborting this script"; tput sgr0
	return 0;
fi

# load the variables (even if already done)
source ${ROOTPWA}/scripts/pwa_Kpipi_example/set_workspace_var.sh

# array with variables 1/0 for each step to hold
export STEP_NAME=(\
 "             getting files              "\
 "          setting up workspace          "\
 "     Filling flat phasespace events     "\
 "    MC acceptance of flat phasespace    "\
 "  Filtering RAW and MC data into bins   "\
 "  Generation of partial wave key files  "\
 " Calculation of Amplitudes and Integrals"\
 "      Specifying amplitudes for fit     "\
 "        Fitting the partial waves       "\
 "      Visualization of fit results      "\
 )
export STEP_DONE=( 0 0 0 0 0 0 0 0 0 0 );

# perform some test giving hints for the current PWA analysis status

# step 1: is the raw file in place?
if [ -f ${KPIPI_RAW_FILE_DIR}/${KPIPI_RAW_FILE} ]
then
	if [ -f ${KPIPI_MC_FILE_DIR}/${KPIPI_MC_FILE} ]
	then
		STEP_DONE[0]=1
	fi
fi

# step 6: search for at least one .key file
# we need the number of keyfiles for further tests first
let KEYFILE_COUNTER=0
if [ -d ${KPIPI_KEYFILE_DIR} ] 
then
	for KEYFILE in ${KPIPI_KEYFILE_DIR}/* 
	do
		EXT=${KEYFILE#*.}
		if [ ${EXT} == "key" ]
		then
			((KEYFILE_COUNTER++))
			STEP_DONE[5]=1
		fi
	done	
fi

# step 2: Is the workspace directory set correctly? Are there some bins?
# step 3: Searching for .genbod.evt files
STEP_DONE[2]=1
FILTERED_MC_EXACT_EVENTS_EXIST=1
# step 4: MC data checks
STEP_DONE[3]=1
# step 5: Searching for .evt files and .acc.evt files
STEP_DONE[4]=1
FILTERED_RAW_EVENTS_EXIST=1
FILTERED_MC_ACC_EVENTS_EXIST=1
# step 7: count the .amp files in the AMPS folder
STEP_DONE[6]=1
if [ -d ${KPIPI_WORK_DIR} ]
then
	# count valid directories
	let COUNTBINS=0;
	for BIN in ${KPIPI_WORK_DIR}/*
	do
		_BIN=$(basename ${BIN}) # get the directory name
		BINHIGH=${_BIN#*.}		# the number behind the "." is the hight bound 
		BINLOW=${_BIN%.*}		# the number in front of the "." is the low bound
		# not everything in ./* is always a valid folder. Check the name to have numbers	
		if echo ${BINLOW} | grep "^[0-9]*$">/tmp/aux
		then
			((COUNTBINS++));
			if [ ! -f ${BIN}/${BINLOW}.${BINHIGH}.genbod.evt ]
			then
				STEP_DONE[2]=0 # only one genbod file must be missing to turn the result false
				FILTERED_MC_EXACT_EVENTS_EXIST=0
			fi	
			if [ ! -f ${KPIPI_MC_FILE_DIR}/${KPIPI_MC_FILE} ]
			then
				STEP_DONE[3]=0 # this check is sensless here since the file is obtained by the web. To be changed.
			fi
			if [ ! -f ${BIN}/${BINLOW}.${BINHIGH}.evt ]
			then
				STEP_DONE[4]=0 # only one event file must be missing to turn the result false
				FILTERED_RAW_EVENTS_EXIST=0
			fi
			if [ ! -f ${BIN}/${BINLOW}.${BINHIGH}.acc.evt ]
			then
				STEP_DONE[4]=0 # only one acceptance file must be missing to turn the result false
				FILTERED_MC_ACC_EVENTS_EXIST=0
			fi
			# check for .amp and .int files in AMPS folder
			let AMPFILE_COUNTER=0
			for AMPFILE in cat $(ls ${BIN}/AMPS/)
			do
				EXT=${AMPFILE#*.}
				if [ ${EXT} == "amp" ]
				then
					((AMPFILE_COUNTER++))
				fi
			done
			if [ ! ${AMPFILE_COUNTER} == ${KEYFILE_COUNTER} ]
			then
				STEP_DONE[6]=0	
			fi
			if [ ! -f ${BIN}/AMPS/norm.int ]
			then
				STEP_DONE[6]=0
			fi
			# check for .amp and .int files in PSPAMPS folder
			let AMPFILE_COUNTER=0
			for AMPFILE in cat $(ls ${BIN}/PSPAMPS/)
			do
				EXT=${AMPFILE#*.}
				if [ ${EXT} == "amp" ]
				then
					((AMPFILE_COUNTER++))
				fi
			done
			if [ ! ${AMPFILE_COUNTER} == ${KEYFILE_COUNTER} ]
			then
				STEP_DONE[6]=0	
			fi
			if [ ! -f ${BIN}/PSPAMPS/norm.int ]
			then
				STEP_DONE[6]=0
			fi
			# check for .amp and .int files in ACCAMPS folder
			let AMPFILE_COUNTER=0
			for AMPFILE in cat $(ls ${BIN}/ACCAMPS/)
			do
				EXT=${AMPFILE#*.}
				if [ ${EXT} == "amp" ]
				then
					((AMPFILE_COUNTER++))
				fi
			done
			if [ ! ${AMPFILE_COUNTER} == ${KEYFILE_COUNTER} ]
			then
				STEP_DONE[6]=0	
			fi
			if [ ! -f ${BIN}/ACCAMPS/norm.int ]
			then
				STEP_DONE[6]=0
			fi
		fi
	done
	if [ ${COUNTBINS} == ${KPIPI_NBINS} ]
	then
		STEP_DONE[1]=1
	else
		# negate also some other steps in this case
		STEP_DONE[2]=0
		STEP_DONE[3]=0
		STEP_DONE[4]=0
		STEP_DONE[6]=0
	fi	
else
	# negate also some other steps in this case
	STEP_DONE[2]=0
	STEP_DONE[3]=0
	STEP_DONE[4]=0
	STEP_DONE[6]=0
fi

# step 8: search for wavelist
if [ -f ${KPIPI_WAVE_LIST} ] 
then
	STEP_DONE[7]=1
fi

# step 9: number of fit files should be equal to number of bins
let RESULTFILE_COUNTER=0
for RESULTFILE in ${KPIPI_FIT_DIR}/*.result.root
do
	((RESULTFILE_COUNTER++))
done
if [ ${RESULTFILE_COUNTER} == ${KPIPI_NBINS} ]
then
	STEP_DONE[8]=1
fi

# step 10: Is there already a file available showing the Intensities?
if [ -f ${KPIPI_FIT_DIR}/waveintensities.ps ]
then
	STEP_DONE[9]=1
fi

# show the results of the performed tests
echo -e "\n The current status of the PWA analysis chain is: \n"
for (( I=0; I<10; I++ )) # in
do
	if [ ${STEP_DONE[$I]} == "1" ]
	then
		echo -e " ${STEP_NAME[$I]}\E[37;32m completed"; tput sgr0
	else
		echo -e " ${STEP_NAME[$I]}\E[37;31m not complete "; tput sgr0
	fi
done
echo -e "\n Warning: the consistency of data is not beeing checked! \n"
