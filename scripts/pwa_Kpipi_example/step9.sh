# this script calls the fitting procedure
# should be called by run_X_PWA_analysis.sh script
# author: P.Jasinski Promme@web.de , jasinski@kph.uni-mainz.de

source ${ROOTPWA}/scripts/pwa_Kpipi_example/set_workspace_var.sh

echo -e "\n **************** part 9 ******************"
echo      " *${STEP_NAME[8]}*"
echo -e   " ******************************************\n"

if which pwafit >/dev/null; then
	echo -e "\n pwafit found"
else
	echo -e "\E[37;31m \n pwa fitter pwafit not found:"; tput sgr0
	echo -e "\E[37;31m Please compile rootpwa according to "
	echo -e "\E[37;34m ${ROOTPWA}/INSTALL "
	echo -e "\E[37;31m instructions "; tput sgr0
	echo -e "\n Aborting this script"; tput sgr0
	return 0
fi

RANK=1 # 1,2 or 3

OVERRIDE_FIT="N"
for BIN in ${KPIPI_WORK_DIR}/*
do
	_BIN=$(basename ${BIN}) # get the directory name
	BINHIGH=${_BIN#*.}		# the number after the "." is the high bound
	BINLOW=${_BIN%.*}		# the number in front of the "." is the low bound
	# not everything in ./* is always a valid folder. Check the name to have numbers
	if echo ${BINLOW} | grep "^[0-9]*$">/tmp/aux
	then
		OUTFILE=${KPIPI_FIT_DIR}/${_BIN}.result.root
		if [ -e ${OUTFILE} ];
		then
			if [ ${OVERRIDE_FIT}=="N" ]; # ask once to override if file exists already
			then
				echo -e "\E[37;31m \n fit result"
				echo -e "\E[37;34m ${OUTFILE}"
				echo -e "\E[37;31m already exists! Override all? (Y/N) "; tput sgr0
				read OVERRIDE_FIT
			fi
			rm ${OUTFILE}
		fi

		if [ -e ${OUTFILE} ]; # perform only if the file does not exist or was removed previously
		then
			echo -e "\E[37;31m \n skipping existing bin fit! "
			echo -e "\E[37;34m ${OUTFILE}"; tput sgr0
		else
			echo -e " \n fitting ${BIN}"
			echo -e " into ${OUTFILE}"
			echo -e " Starting time: ";
			date;
			cd ${BIN}/AMPS;
			NORM="../PSPAMPS/norm.int"
			if [ -e norm.int ];
			then
				echo -e " Default normalization integral already existent. "
			else
				echo -e " using ${NORM} for normalization ->\E[37;31m no acceptance correction! "; tput sgr0
				cp ${NORM} ./
			fi
			pwafit -q -w ${KPIPI_WAVE_LIST} -o ${OUTFILE} -r ${RANK} -l ${BINLOW} -N -u ${BINHIGH} -n ${NORM}
		fi
	else
		echo -e "\n skipping ${BIN}"
	fi
	rm /tmp/aux
done

cd ${_PWD}
