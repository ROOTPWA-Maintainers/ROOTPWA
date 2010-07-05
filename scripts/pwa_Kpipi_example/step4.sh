# this script filters data from a root file to binned BNL events
# should be called by run_X_PWA_analysis.sh script
# author: P.Jasinski Promme@web.de , jasinski@kph.uni-mainz.de

source ${ROOTPWA}/scripts/pwa_Kpipi_example/set_workspace_var.sh

echo -e "\n **************** part 4 ******************"
echo      " *${STEP_NAME[3]}*"
echo -e   " ******************************************\n"

if [ "${FILTERED_RAWEVENTS}" == '1' ] # it was previously checked if events were already filtered
then
	echo -e "\E[37;31m omitting: There are allready filtered files existing! "; tput sgr0 
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
