# this script takes the MC acceptance events and puts them
# into the bins. The MC analysis part is not included here.
# should be called by run_X_PWA_analysis.sh script
# author: P.Jasinski Promme@web.de , jasinski@kph.uni-mainz.de

source ${ROOTPWA}/scripts/pwa_Kpipi_example/set_workspace_var.sh

echo -e "\n **************** part 4 ******************"
echo      " *${STEP_NAME[3]}*"
echo -e   " ******************************************\n"

echo -e "\E[37;31m \n This part is partialy implemented here and needs special treatment."
echo -e "\E[37;31m Please ask COMGEANT/CORAL experts to retrieve instructions."
echo -e "\E[37;31m We skip the COMGEANT/CORAL/PHAST chain and use directly the results."; tput sgr0
