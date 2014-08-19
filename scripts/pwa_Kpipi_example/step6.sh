# this script generates the key files
# should be called by run_X_PWA_analysis.sh script
# author: P.Jasinski Promme@web.de , jasinski@kph.uni-mainz.de

source ${ROOTPWA}/scripts/pwa_Kpipi_example/set_workspace_var.sh

echo -e "\n **************** part 6 ******************"
echo      " *${STEP_NAME[5]}*"
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
