# this script creates and runs/submits scripts
# for the calculation of the amplitudes and integrals
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
  		# echo " Starting time: ";
  		# date;

		if [ "${KPIPI_CALC_AMP_FARM}" == 'mainz' ] 
		then
			# write and submit a scrip for data amplitude calculation
			export AMPDIR=${BIN}/AMPS;
	  		export FILE="${BIN}/${BINLOW}.${BINHIGH}.evt";
	  		echo " input: ${FILE}";
			export BATCHFILE="${BIN}/calcamp_batch_${KPIPI_CALC_AMP_FARM}.evt.sh"
			cat > ${BATCHFILE} << EOF
#!/bin/bash
# 
#PBS -N calcampJob
#PBS -j oe
#PBS -o ${BATCHFILE}.out
#PBS -V
#PBS -l nodes=1:x86_64

export PATH=$PBS_O_PATH
cd $PBS_O_WORKDIR
./doamps.sh ${FILE}
			EOF
			echo " submitting batchscript ${BATCHFILE}"
			qsub ${BATCHFILE}

	  		# write and submit a script for monte carlo amplitude calculation
	  		export FILE=${FILE/evt/genbod.evt};
	  		export AMPDIR=${BIN}/PSPAMPS;
			echo " mc-input: ${FILE}";
	  		export BATCHFILE="${BIN}/calcamp_batch_${KPIPI_CALC_AMP_FARM}.genbod.evt.sh"
			cat > ${BATCHFILE} << EOF
#!/bin/bash
# 
#PBS -N calcmcampJob
#PBS -j oe
#PBS -o ${BATCHFILE}.out
#PBS -V
#PBS -l nodes=1:x86_64

export PATH=$PBS_O_PATH
cd $PBS_O_WORKDIR
./doamps.sh ${FILE}
cd ${AMPDIR}
int *.amp > norm.int;
			EOF
			echo " submitting batchscript ${BATCHFILE}"
			qsub ${BATCHFILE}

			# write and submit a script for monte carlo accepted calculation
	  		export FILE=${FILE/genbod.evt/acc.evt};
	  		export AMPDIR=${BIN}/ACCAMPS;
			echo " acc-input: ${FILE}";
	  		export BATCHFILE="${BIN}/calcamp_batch_${KPIPI_CALC_AMP_FARM}.acc.evt.sh"
			cat > ${BATCHFILE} << EOF
#!/bin/bash
# 
#PBS -N calcmcaccampJob
#PBS -j oe
#PBS -o ${BATCHFILE}.out
#PBS -V
#PBS -l nodes=1:x86_64

export PATH=$PBS_O_PATH
cd $PBS_O_WORKDIR
./doamps.sh ${FILE}
cd ${AMPDIR}
int *.amp > norm.int;
			EOF
			echo " submitting batchscript ${BATCHFILE}"
			qsub ${BATCHFILE}

		else # if mainz
		if [ "${KPIPI_CALC_AMP_FARM}" == 'local' ]
		then
			# data amplitude calculation
			export AMPDIR=${BIN}/AMPS;
	  		export FILE="${BIN}/${BINLOW}.${BINHIGH}.evt";
	  		echo " input: ${FILE}";
	  		./doamps.sh ${FILE};

	  		# monte carlo amplitude calculation
	  		export FILE=${FILE/evt/genbod.evt};
	  		echo "---- mc-input: ${FILE}";
	  		export AMPDIR=${BIN}/PSPAMPS;
	  		./doamps.sh ${FILE};

	  		# Do integration
	  		cd ${AMPDIR}
	  		int *.amp > norm.int;
			cd -;

			# monte carlo accepted calculation
	  		export FILE=${FILE/genbod.evt/acc.evt};
	  		echo "---- acc-input: ${FILE}";
	  		export AMPDIR=${BIN}/ACCAMPS;
	  		./doamps.sh ${FILE};

	  		# Do integration
	  		cd ${AMPDIR}
	  		int *.amp > norm.int;		
	  		cd -;
		else # if local
			echo " amplitude analysis for ${KPIPI_CALC_AMP_FARM} not implemented yet! "
		fi # if local
		fi # if mainz		
	else
  		echo -e "\n skipping ${BIN}"
	fi
	rm /tmp/aux
done

cd ${_PWD}
