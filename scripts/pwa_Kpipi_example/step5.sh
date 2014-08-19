# this script filters data from a root file to binned BNL events
# should be called by run_X_PWA_analysis.sh script
# test_workspace.sh must be run before to set some additional variables
# author: P.Jasinski Promme@web.de , jasinski@kph.uni-mainz.de

source ${ROOTPWA}/scripts/pwa_Kpipi_example/set_workspace_var.sh

echo -e "\n **************** part 5 ******************"
echo      " *${STEP_NAME[4]}*"
echo -e   " ******************************************\n"

echo -e "\n compiling tree reader "
g++ `root-config --glibs` `root-config --cflags` -O0 -o /tmp/readtree ${KPIPI_TREE_READER}
if [ -f /tmp/readtree ]
then

	echo -e " compilation succeeded"
	if [ "${FILTERED_RAW_EVENTS_EXIST}" == '1' ]
	then
		echo -e "\E[37;31m RAW events seem to exist already. skipping! "; tput sgr0
	else
		echo -e " filling RAW events... "
		/tmp/readtree ${KPIPI_RAW_FILE_DIR}/${KPIPI_RAW_FILE} ${KPIPI_FILE_TREE} \
		${KPIPI_WORK_DIR} ${KPIPI_NBINS} ${KPIPI_BIN_MIN} ${KPIPI_BIN_MAX} ${KPIPI_NPARTICLES}
	fi

	if [ "${FILTERED_MC_ACC_EVENTS_EXIST}" == '1' ]
	then
		echo -e "\E[37;31m MC acceptance events seem to exist already. skipping! "; tput sgr0
	else
		echo -e " filling MC acceptance events... "
		/tmp/readtree ${KPIPI_MC_FILE_DIR}/${KPIPI_MC_FILE} ${KPIPI_FILE_TREE} \
		${KPIPI_WORK_DIR} ${KPIPI_NBINS} ${KPIPI_BIN_MIN} ${KPIPI_BIN_MAX} ${KPIPI_NPARTICLES}
	fi
	rm /tmp/readtree
else
	echo -e "\E[37;31m \n Error in compilation to /tmp/readtree "; tput sgr0
	echo -e "\E[37;31m Please check! "; tput sgr0
	echo -e "\n Aborting this script"; tput sgr0
	return 0
fi
