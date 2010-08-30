# this script generates flat phasespace events into the bins
# should be called by run_X_PWA_analysis.sh script
# author: P.Jasinski Promme@web.de , jasinski@kph.uni-mainz.de

source ${ROOTPWA}/scripts/pwa_Kpipi_example/set_workspace_var.sh

echo -e "\n **************** part 3 ******************"
echo      " *${STEP_NAME[2]}*"
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
	BINHIGH=${_BIN#*.}		# the number behind the "." is the hight bound 
	BINLOW=${_BIN%.*}		# the number in front of the "." is the low bound
    let BINWIDTH=(${KPIPI_BIN_MAX}-${KPIPI_BIN_MIN})/${KPIPI_NBINS}
	# not everything in ./* is always a valid folder. Check the name to have numbers	
	if echo ${BINLOW} | grep "^[0-9]*$">/tmp/aux
	then
		echo -e "\n attempting to generate flat phasespace MC events in "
		echo -e "${BIN}"
		cd ${BIN}
		#echo "${BINHIGH} ${BINLOW} ${BINWIDTH}"
		if [ -e ${BINLOW}.${BINHIGH}.genbod.evt ]
		then
			echo -e "\E[37;31m omitting: There is allready a file existing! "; tput sgr0	
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
