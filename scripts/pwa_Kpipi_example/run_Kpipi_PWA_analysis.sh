# author: P.Jasinski Promme@web.de , jasinski@kph.uni-mainz.de
# example for the usage of rootpwa [working with SVN version 278]
# ${ROOTPWA} must be set to rootpwa base directory
# please check also paths in set_workspace_var.sh first!

#!/bin/bash

export _PWD=$PWD # store the path the script was started from
cd /tmp/ # if something goes wrong usually plenty of data is still created

echo
echo -e '\E[37;31m Example for K pi pi PWA chain in rootpwa by Promme@web.de'; tput sgr0
date
echo

# load the variables and paths needed to run this file
source ${ROOTPWA}/scripts/pwa_Kpipi_example/set_workspace_var.sh
source ${ROOTPWA}/scripts/pwa_Kpipi_example/test_workspace.sh # exports also ${STEP_NAME[]} variable

RUNSCRIPT=( "N" "N" "N" "N" "N" "N" "N" "N" "N" "N" ) # setup an array to determine which scripts to run

echo -e "\n Perform all unfinished steps? (Y/N)"
read PERFORM_ALL_UNFINISHED
for (( STEP=0; STEP<10; STEP++ )) # in
do
	if [ ${STEP_DONE[${STEP}]} == "0" ]
	then
		if [ ! ${PERFORM_ALL_UNFINISHED} == "Y" ]
		then
			echo -e " Do you want to run step ${STEP}: ${STEP_NAME[${STEP}]}? (Y/N)"
			read RUNSCRIPT[${STEP}]
		else
			RUNSCRIPT[${STEP}]="Y"
		fi
	else
		RUNSCRIPT[${STEP}]="N"
	fi
done	

# scripts to run
RUNSCRIPT_FILE=(\
 "${ROOTPWA}/scripts/pwa_Kpipi_example/step1.sh"\
 "${ROOTPWA}/scripts/pwa_Kpipi_example/step2.sh"\
 "${ROOTPWA}/scripts/pwa_Kpipi_example/step3.sh"\
 "${ROOTPWA}/scripts/pwa_Kpipi_example/step4.sh"\
 "${ROOTPWA}/scripts/pwa_Kpipi_example/step5.sh"\
 "${ROOTPWA}/scripts/pwa_Kpipi_example/step6.sh"\
 "${ROOTPWA}/scripts/pwa_Kpipi_example/step7.sh"\
 "${ROOTPWA}/scripts/pwa_Kpipi_example/step8.sh"\
 "${ROOTPWA}/scripts/pwa_Kpipi_example/step9.sh"\
 "${ROOTPWA}/scripts/pwa_Kpipi_example/step10.sh"\
)

for (( STEP=0; STEP<10; STEP++ )) # in
do 
	# echo "${STEP} is ${RUNSCRIPT[${STEP}]}"
	if [ ${RUNSCRIPT[${STEP}]} == "Y" ]
	then
		source ${RUNSCRIPT_FILE[${STEP}]}
	fi
done

cd ${_PWD}

echo -e "\n finished! \n"

return 0
