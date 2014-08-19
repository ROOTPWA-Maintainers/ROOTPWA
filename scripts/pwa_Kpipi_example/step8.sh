# this script filters creates a list of amplitudes available
# by using the available keys
# should be called by run_X_PWA_analysis.sh script
# author: P.Jasinski Promme@web.de , jasinski@kph.uni-mainz.de

source ${ROOTPWA}/scripts/pwa_Kpipi_example/set_workspace_var.sh

echo -e "\n **************** part 8 ******************"
echo      " *${STEP_NAME[7]}*"
echo -e   " ******************************************\n"

# put the list with amplitudes to fit in the fit directory
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

# here we take all waves for all bins, you may change it binwise
# by applying different wavelists

ls ${KPIPI_KEYFILE_DIR} > "${KPIPI_WAVE_LIST}_temp"
# cat "${KPIPI_WAVE_LIST}_temp"
for AMP in $(cat ${KPIPI_WAVE_LIST}_temp); # keep only .key files
do
	EXT=${AMP#*.}
	FNAME=${AMP%.*}
	if [ ${EXT} == "key" ];
	then
		# save with the extension to .amp
		echo "${FNAME}.amp" >> ${KPIPI_WAVE_LIST}
	fi
done
rm "${KPIPI_WAVE_LIST}_temp"
echo -e "\E[37;31m \n the following amplitudes will be fitted: \n"; tput sgr0
cat ${KPIPI_WAVE_LIST}
