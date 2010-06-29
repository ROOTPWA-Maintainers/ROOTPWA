echo -e "\n **************** part 1 ******************"
echo      " *${STEP_NAME[0]}*"
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
