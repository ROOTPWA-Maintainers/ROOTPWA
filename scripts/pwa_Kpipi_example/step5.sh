# this script does nothing but printing
# should be called by run_X_PWA_analysis.sh script
# author: P.Jasinski Promme@web.de , jasinski@kph.uni-mainz.de

source ${ROOTPWA}/scripts/pwa_Kpipi_example/set_workspace_var.sh

echo -e "\n **************** part 5 ******************"
echo      " *${STEP_NAME[4]}*"
echo -e   " ******************************************\n"

echo -e "\E[37;31m \n This part is not implemented here and needs special treatment."
echo -e "\E[37;31m Please ask COMGEANT/CORAL experts to retrieve instructions."
echo -e "\E[37;31m We skip this part and perform PWA without MC acceptance correction."; tput sgr0
