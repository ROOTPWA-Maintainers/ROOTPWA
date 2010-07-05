# author: P.Jasinski Promme@web.de , jasinski@kph.uni-mainz.de
# example for the usage of rootpwa [working with SVN version 278]

#!/bin/bash

_PWD=$PWD # store the path the script was started from

echo
echo -e '\E[37;31m Example for K pi pi PWA chain in rootpwa by Promme@web.de'; tput sgr0
date
echo

# ************************ variables to set ! **************************
	export KPIPI_RAW_FILE_DIR="/mnt/disk2/temp" # the directory where the data file will be copied to (must exist)
	export KPIPI_RAW_FILE="hist.root"			# the name of the file
	export KPIPI_WEB_FILE="http://www.staff.uni-mainz.de/jasinsk/demo-data/hist.root" # where to get the file from
	export KPIPI_FILE_TREE="events/tree_events" # path to the tree in the rootfile
	export KPIPI_NPARTICLES=4 					# number of particles in the tree per event

	export KPIPI_WORK_DIR="/mnt/disk2/analysis/Kp_Kppipi/PWA2/SET_test/" # work directory containing binned data and calculated amplitudes (path to directory must exist)
	export KPIPI_BIN_MIN=800  # lowest energy bin of invariant Kpipi mass spectrum to cut
	export KPIPI_BIN_MAX=3000 # highest energy bin of invariant Kpipi mass spectrum to cut
	export KPIPI_NBINS=110	  # number of bins to create -> number of folders that will be created

	export KPIPI_NPHASESPACE_MC_EVENTS=10000 # number of flat phasespace events to be MC generated per bin
	export KPIPI_MC_CONFIG_FILE_FLAT=${ROOTPWA}/generators/KpipiDiffractive2008.conf # configuration file for flat diffractive phase space

	export KPIPI_TREE_READER=${ROOTPWA}/keyfiles/keyKpipi/scripts/Tree_to_BNL_event.C # path to the tree converter that fills the bins

	export KPIPI_KEYFILE_DIR=${ROOTPWA}/keyfiles/keyKpipi/SETWA3_test 	# path to keyfile directory
	export KPIPI_KEYFILE_GENERATOR=${ROOTPWA}/keygen/genKey_Kpipi_WA3.C # path to keyfile generator
	export KPIPI_GENERATOR_PATH=${ROOTPWA}/keygen						# currently the path where the generator must be compiled in

	export KPIPI_FIT_DIR="/mnt/disk2/analysis/Kp_Kppipi/PWA2/FIT_test/" # path to the fit directory

# ************************ end variables *******************************

echo -e "\n **************** part 1 ******************"
echo      " *             getting files              *"
echo -e   " ******************************************\n"

if [ -f ${KPIPI_RAW_FILE_DIR}/${KPIPI_RAW_FILE} ]
then
	echo -e "\n root file ${KPIPI_RAW_FILE_DIR}/${KPIPI_RAW_FILE} "
	echo    " seems to exist "
	echo -e " skipping copy procedure! \n"
else
	if which wget >/dev/null; then
    	echo -e "\n wget found"
	else
		echo -e '\E[37;31m \n wget not found: Please download'
		echo -e "\E[37;34m ${KPIPI_WEB_FILE}" 
		echo -e "\E[37;31m with a webbrowser to "
		echo -e "\E[37;34m ${KPIPI_RAW_FILE_DIR}/${KPIPI_RAW_FILE} "
		echo -e "\n Aborting this script"; tput sgr0
		return 0
	fi

	echo " copying file from ${KPIPI_WEB_FILE} to ${KPIPI_RAW_FILE_DIR}/${KPIPI_RAW_FILE}"

	wget ${KPIPI_WEB_FILE} -O ${KPIPI_RAW_FILE_DIR}/${KPIPI_RAW_FILE}
	#return 0
fi

echo -e "\n **************** part 2 ******************"
echo      " *          setting up workspace          *"
echo -e   " ******************************************\n"

if [ -e ${ROOTPWA} ]; 
then
	echo -e "\n ROOTPWA variable is set to ${ROOTPWA} . \n"
else
	echo -e "\E[37;31m \n ROOTPWA variable is not set correctly! "; tput sgr0
	echo -e "\E[37;31m Please set it first! "; tput sgr0
	echo -e "\n Aborting this script"; tput sgr0
	return 0;
fi

if [ -e ${KPIPI_WORK_DIR} ]; 
then
	echo -e "\E[37;31m \n The directory "; tput sgr0
	echo -e "\E[37;34m ${KPIPI_WORK_DIR}"; tput sgr0
	echo -e "\E[37;31m exists already. Remove the contents (Y/N)? "; tput sgr0
	read CONTINUE
	if [ ${CONTINUE} == "Y" ];
	then
		echo -e "\n cleaning up ${KPIPI_WORK_DIR}\n"
		rm -rf ${KPIPI_WORK_DIR}
		mkdir ${KPIPI_WORK_DIR}
	fi
fi

let BINWIDTH=(${KPIPI_BIN_MAX}-${KPIPI_BIN_MIN})/${KPIPI_NBINS}
echo -e "\n creating ${KPIPI_NBINS} folder in ${BINWIDTH} (MeV) steps "
for (( I=${KPIPI_BIN_MIN}; $I<${KPIPI_BIN_MAX}; I+=${BINWIDTH} )) # in
do
	BINLOW=$I
	let BINHIGH=$I+$BINWIDTH
	FOLDERNAME=${BINLOW}.${BINHIGH}
	if [ -e ${KPIPI_WORK_DIR}/${FOLDERNAME} ]
	then
		echo -e "\E[37;31m \n ommiting existing directory "; tput sgr0
		echo -e "\E[37;34m ${FOLDERNAME} \n"; tput sgr0
	else
		echo -e "\n creating folder ${FOLDERNAME}"
		mkdir ${KPIPI_WORK_DIR}/${FOLDERNAME}
		cd ${KPIPI_WORK_DIR}/${FOLDERNAME}
		mkdir "ACCAMPS"
		mkdir "AMPS"
		mkdir "MC"
		mkdir "PSPAMPS" 
		cd -
	fi
done



echo -e "\n **************** part 3 ******************"
echo      " *     Filling flat phasespace events     *"
echo -e   " ******************************************\n"

if which genpw >/dev/null; then
	echo -e "\n genpw found"
else
	echo -e "\E[37;31m \n phasespace generator genpw not found:"; tput sgr0
	echo -e "\E[37;31m Please compile rootpwa according to "
	echo -e "\E[37;34m ${ROOTPWA}/INSTALL "
	echo -e "\E[37;31m instructions "; tput sgr0
	echo -e "\n Aborting this script"; tput sgr0
	return 0
fi

FILTERED_RAWEVENTS=""
for BIN in ${KPIPI_WORK_DIR}/*
do
	_BIN=$(basename ${BIN}) # get the directory name
	BINHIGH=${_BIN#*.}		# the number after the "." is the hight bound 
	BINLOW=${_BIN%.*}		# the number before the "." is the low bound
	# not everything in ./* is always a valid folder. Check the name to have numbers	
	if echo ${BINLOW} | grep "^[0-9]*$">/tmp/aux
	then
		echo -e "\n attempting to generate flat phasespace MC events in "
		echo -e "${BIN}"
		cd ${BIN}
		#echo "${BINHIGH} ${BINLOW} ${BINWIDTH}"
		if [ -e ${BINLOW}.${BINHIGH}.genbod.evt ]
		then
			echo -e "\E[37;31m ommiting: There is allready a file existing! "; tput sgr0	
		else
			#echo " doing it "
			genpw -n ${KPIPI_NPHASESPACE_MC_EVENTS} -M ${BINLOW} -B ${BINWIDTH} -c -r ${KPIPI_MC_CONFIG_FILE_FLAT}
		fi
		# by the way: do we have filtered RAW events in here?
		if [ -e ${BINLOW}.${BINHIGH}.evt ]
		then
			echo " FILE EXISTS "
			FILTERED_RAWEVENTS=1
		fi	
		cd -  		
	else
  		echo -e "\n skipping ${BIN}"
	fi
	rm /tmp/aux
done

echo -e "\n **************** part 4 ******************"
echo      " *      Filtering RAW data into bins      *"
echo -e   " ******************************************\n"

if [ "${FILTERED_RAWEVENTS}" == '1' ] # it was previously checked if events were already filtered
then
	echo -e "\E[37;31m ommiting: There are allready filtered files existing! "; tput sgr0 
else
	echo -e "\n compiling tree reader "
	g++ `root-config --glibs` `root-config --cflags` -O0 -o /tmp/readtree ${KPIPI_TREE_READER}
	if [ -f /tmp/readtree ]
	then
		echo -e " compilation succeeded, filling events. " 
		/tmp/readtree ${KPIPI_RAW_FILE_DIR}/${KPIPI_RAW_FILE} ${KPIPI_FILE_TREE} \
		${KPIPI_WORK_DIR} ${KPIPI_NBINS} ${KPIPI_BIN_MIN} ${KPIPI_BIN_MAX} ${KPIPI_NPARTICLES}
		rm /tmp/readtree
	else
		echo -e "\E[37;31m \n Error compiling to /tmp/readtree "; tput sgr0
		echo -e "\E[37;31m Please check! "; tput sgr0
		echo -e "\n Aborting this script"; tput sgr0
		return 0
	fi
fi

echo -e "\n **************** part 5 ******************"
echo      " *    MC acceptance of flat phasespace    *"
echo -e   " ******************************************\n"

echo -e "\E[37;31m \n This part is not implemented here and needs special treatment."
echo -e "\E[37;31m Please ask COMGEANT/CORAL experts to retrieve instructions."
echo -e "\E[37;31m We skip this part and perform PWA without MC acceptance correction."; tput sgr0

echo -e "\n **************** part 6 ******************"
echo      " *  Generation of partial wave key files  *"
echo -e   " ******************************************\n"

if [ -e ${KPIPI_KEYFILE_DIR} ]; 
then
	echo -e "\E[37;31m \n The directory "; tput sgr0
	echo -e "\E[37;34m ${KPIPI_KEYFILE_DIR}"; tput sgr0
	echo -e "\E[37;31m exists already. Remove the contents (Y/N)? "; tput sgr0
	read CONTINUE
	if [ ${CONTINUE} == "Y" ];
	then
		echo -e "\n cleaning up ${KPIPI_KEYFILE_DIR}\n"
		rm -rf ${KPIPI_KEYFILE_DIR}
	else
		echo -e "\E[37;31m \n Please specify a new key file folder! "; tput sgr0
		echo -e "\n Aborting this script"; tput sgr0
		return 0
	fi	
fi
mkdir ${KPIPI_KEYFILE_DIR}

# currently the key file generator must be compiled in the path containing the generators
cd ${KPIPI_GENERATOR_PATH}
mkdir backup # backup existing keyfiles in this folder (that should anyhow not last there)
mv ./*.key ./backup/
root -l -q -b ${KPIPI_KEYFILE_GENERATOR}+\(true,\"../keyfiles/keyKpipi/testEventsKpipi.evt\",\"./pdgTable.txt\",\"${KPIPI_KEYFILE_DIR}\"\)
echo " moving key files to destination folder if not done already "
mv ./*.key ${KPIPI_KEYFILE_DIR}/
rm -f ./*.key.C	
cd -

echo -e "\n **************** part 7 ******************"
echo      " * Calculation of Amplitudes and Integrals*"
echo -e   " ******************************************\n"

cd ${ROOTPWA}/scripts # there are other scripts used here
export KEYDIR=${KPIPI_KEYFILE_DIR} # needed for doamps.sh script 

for BIN in ${KPIPI_WORK_DIR}/*
do
	_BIN=$(basename ${BIN}) # get the directory name
	BINHIGH=${_BIN#*.}		# the number after the "." is the hight bound 
	BINLOW=${_BIN%.*}		# the number before the "." is the low bound
	# not everything in ./* is always a valid folder. Check the name to have numbers	
	if echo ${BINLOW} | grep "^[0-9]*$">/tmp/aux
	then
		echo -e "\n attempting to calculate amplitudes in "
		echo " PROCESSING MASS BIN: ${BIN}";
  		echo " Starting time: ";
  		date;
  		export AMPDIR=${BIN}/AMPS;
  		export FILE="${BIN}/${BINLOW}.${BINHIGH}.evt";
  		echo " input: ${FILE}";
  		./doamps.sh ${FILE};

  		# Do monte carlo
  		export FILE=${FILE/evt/genbod.evt};
  		echo "---- mc-input: ${FILE}";
  		export AMPDIR=${BIN}/PSPAMPS;
  		./doamps.sh ${FILE};

  		# Do integration
  		cd ${AMPDIR}
  		int *.amp > norm.int;
  		cd -;		
	else
  		echo -e "\n skipping ${BIN}"
	fi
	rm /tmp/aux
done

cd ${_PWD}

echo -e "\n **************** part 8 ******************"
echo      " *        Fitting the partial waves       *"
echo -e   " ******************************************\n"

if [ -e ${KPIPI_FIT_DIR} ]; 
then
	echo -e "\E[37;31m \n The directory "; tput sgr0
	echo -e "\E[37;34m ${KPIPI_FIT_DIR}"; tput sgr0
	echo -e "\E[37;31m exists already. Remove the contents (Y/N)? "; tput sgr0
	read CONTINUE
	if [ ${CONTINUE} == "Y" ];
	then
		echo -e "\n cleaning up ${KPIPI_FIT_DIR}\n"
		rm -rf ${KPIPI_FIT_DIR}
	else
		echo -e "\E[37;31m \n Please specify a new fit file folder! "; tput sgr0
		echo -e "\n Aborting this script"; tput sgr0
		return 0
	fi
fi
mkdir ${KPIPI_FIT_DIR}

cd ${_PWD}

echo -e "\n **************** part 9 ******************"
echo      " *     Visualization of fit results       *"
echo -e   " ******************************************\n"

echo -e "\n finished! \n"

return 0
