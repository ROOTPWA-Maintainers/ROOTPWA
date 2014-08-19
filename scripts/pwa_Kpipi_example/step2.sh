# this script sets the folder structure needed for pwa analysis
# should be called by run_X_PWA_analysis.sh script
# author: P.Jasinski Promme@web.de , jasinski@kph.uni-mainz.de

source ${ROOTPWA}/scripts/pwa_Kpipi_example/set_workspace_var.sh

echo -e "\n **************** part 2 ******************"
echo      " *${STEP_NAME[1]}*"
echo -e   " ******************************************\n"

if [ -d ${KPIPI_WORK_DIR} ];
then
	echo -e "\E[37;31m \n The directory "; tput sgr0
	echo -e "\E[37;34m ${KPIPI_WORK_DIR}"; tput sgr0
	echo -e "\E[37;31m exists already. Remove the contents (Y/N)? "; tput sgr0
	read CONTINUE
	if [ ${CONTINUE} == "Y" ];
	then
		echo -e "\n cleaning up ${KPIPI_WORK_DIR}\n"
		rm -rf ${KPIPI_WORK_DIR}
	fi
fi
mkdir ${KPIPI_WORK_DIR}

let BINWIDTH=(${KPIPI_BIN_MAX}-${KPIPI_BIN_MIN})/${KPIPI_NBINS}
echo -e "\n creating ${KPIPI_NBINS} folder in ${BINWIDTH} (MeV) steps "
for (( I=${KPIPI_BIN_MIN}; $I<${KPIPI_BIN_MAX}; I+=${BINWIDTH} )) # in
do
	BINLOW=$I
	let BINHIGH=$I+$BINWIDTH
	FOLDERNAME=${BINLOW}.${BINHIGH}
	if [ -e ${KPIPI_WORK_DIR}/${FOLDERNAME} ]
	then
		echo -e "\E[37;31m \n omitting existing directory "; tput sgr0
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

return 0
