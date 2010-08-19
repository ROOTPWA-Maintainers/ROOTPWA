# this script sets all variables needed to run the Kpipi analysis
# author: Promme@web.de ; jasinski@kph.uni-mainz.de
# worked with SVN version 258

#!/bin/bash

# ************************ variables to set ! **************************
	export KPIPI_WORKSPACE="/mnt/disk2/analysis/Kp_Kppipi/PWA_20MeV" # if you specify here an exisiting directory with ~5GB space you may leave the rest unchanged

	export KPIPI_RAW_FILE_DIR=${KPIPI_WORKSPACE} # the directory where the data file will be copied to (must exist)
	export KPIPI_RAW_FILE="hist.root"			 # the name of the file
	export KPIPI_MC_FILE_DIR=${KPIPI_WORKSPACE}  # the directory where the mc-acceptance data file will be copied to (must exist)
    export KPIPI_MC_FILE="hist_mc_acc.root";     # the name of the file
	export KPIPI_WEB_FILE="http://www.staff.uni-mainz.de/jasinsk/demo-data/hist.root" # where to get the file from
	export KPIPI_WEB_MC_FILE="http://www.staff.uni-mainz.de/jasinsk/demo-data/hist_mc_acc.root" # where to get the MC acceptance file from		
	export KPIPI_FILE_TREE="events/tree_events_v2" # path to the tree in the rootfile
	export KPIPI_NPARTICLES=0 					 # number of particles in the tree per event (since v2 of trees not important anymore)

	export KPIPI_WORK_DIR="${KPIPI_WORKSPACE}/SET_20MeV/" # work directory containing binned data and calculated amplitudes
	export KPIPI_BIN_MIN=800  # lowest energy bin of invariant Kpipi mass spectrum to cut
	export KPIPI_BIN_MAX=3000 # highest energy bin of invariant Kpipi mass spectrum to cut
	export KPIPI_NBINS=110	  # number of bins to create -> number of folders that will be created

	export KPIPI_NPHASESPACE_MC_EVENTS=100000 # number of flat phasespace events to be MC generated per bin
	export KPIPI_MC_CONFIG_FILE_FLAT=${ROOTPWA}/generators/KpipiDiffractive2008.conf # configuration file for flat diffractive phase space

	export KPIPI_TREE_READER=${ROOTPWA}/keyfiles/keyKpipi/scripts/Tree_to_BNL_event.C # path to the tree converter that fills the bins

	export KPIPI_KEYFILE_DIR=${ROOTPWA}/keyfiles/keyKpipi/SETWA3_full 	# path to keyfile directory
	export KPIPI_KEYFILE_GENERATOR=${ROOTPWA}/keygen/genKey_Kpipi_WA3.C # path to keyfile generator
	export KPIPI_GENERATOR_PATH=${ROOTPWA}/keygen						# currently the path where the generator must be compiled in

	export KPIPI_FIT_DIR="${KPIPI_WORKSPACE}/FIT_20MeV/" # path to the fit directory

	export KPIPI_WAVE_LIST="${KPIPI_FIT_DIR}/wavelist" # path to wave list containing the amplitude names to fit
# ************************ end variables *******************************
