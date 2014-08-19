# this script demonstrates how to visualize the fitting results
# currently Intensities only are shown
# should be called by run_X_PWA_analysis.sh script
# author: P.Jasinski Promme@web.de , jasinski@kph.uni-mainz.de

source ${ROOTPWA}/scripts/pwa_Kpipi_example/set_workspace_var.sh

echo -e "\n **************** part 10 *****************"
echo      " *${STEP_NAME[9]}*"
echo -e   " ******************************************\n"

# visualization must be performed in the rootscript folder containing all needed (logon) scripts
cd ${ROOTPWA}/src/rootscripts
# write a macro to load into root
rm /tmp/rootmacro.C
echo "void rootmacro(){" >> /tmp/rootmacro.C
echo "TTree* tree = loadFitResult(\"${KPIPI_FIT_DIR}/*.result.root\");" >> /tmp/rootmacro.C
echo "plotAllIntensities(tree, true);" >> /tmp/rootmacro.C
echo "};" >> /tmp/rootmacro.C
# run the lines of code
root -l -q /tmp/rootmacro.C
rm /tmp/rootmacro.C

# hopefully a ps file was created containing all intensities
mv ${ROOTPWA}/src/rootscripts/waveIntensities.ps ${KPIPI_FIT_DIR}
# ps2pdf ${KPIPI_FIT_DIR}/waveintensities.ps

cd ${_PWD}
