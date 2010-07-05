# this script calculates the amplitudes and integrals
# should be called by run_X_PWA_analysis.sh script
# author: P.Jasinski Promme@web.de , jasinski@kph.uni-mainz.de

source ${ROOTPWA}/scripts/pwa_Kpipi_example/set_workspace_var.sh

echo -e "\n **************** part 7 ******************"
echo      " *${STEP_NAME[6]}*"
echo -e   " ******************************************\n"

cd ${ROOTPWA}/scripts # there are other scripts used here
export KEYDIR=${KPIPI_KEYFILE_DIR} # needed in doamps.sh script 

for BIN in ${KPIPI_WORK_DIR}/*
do
	_BIN=$(basename ${BIN}) # get the directory name
	BINHIGH=${_BIN#*.}		# the number behind the "." is the high bound 
	BINLOW=${_BIN%.*}		# the number in front of the "." is the low bound
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
