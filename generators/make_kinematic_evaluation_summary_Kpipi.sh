#this script has to be run after calling predict in the rootpwa_gui in the fit folder it self
# currently only working for Kpipi case correctly

#!/bin/bash
echo "removing old summy file"
rm kinevalplots.root
echo "commulating summary plots"
hadd kinevalplots.root kineval_plots.*.root
echo "creating overview summary plots"
root -l -b -q "${ROOTPWA}/generators/rootpalettes.C" "${ROOTPWA}/generators/rootlogon.C" "${ROOTPWA}/generators/plotGlobalWeightedEvts_Kpipi.C+(\"kinevalplots.root\", \"kinevalplots_summary.root\")"