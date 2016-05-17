#!/bin/bash

### BEGIN CONFIGURATION ###

	### BEGIN PHASE SPACE ###

	#MASS BIN (MeV)
	MASS=1800
	BINWIDTH=20

	NMB_PS_EVENTS=500	

	### END PHASE SPACE ###

	PARTICLE_DATA_TABLE="${ROOTPWA}/particleData/particleDataTable.txt"


	### BEGIN PSEUDO DATA ###

	NMB_PSEUDO_EVENTS=500

	### END PSEUDO DATA ###


### END CONFIGURATION ###

if [ -d "${ROOTPWA}" ]; then
	echo "RootPwa installation found in: ${ROOTPWA}"
	echo "Starting MC test ..."
else
	echo "Error no installation of RootPwa found! Aborting."
	exit 1
fi

TESTDIR=${PWD}/ROOTPWA_MC_TEST


### check if directories exist!

echo "Creating directories ..."

mkdir "${TESTDIR}"
mkdir "${TESTDIR}/data"
mkdir "${TESTDIR}/weighted_mc_data"
mkdir "${TESTDIR}/amps"
mkdir "${TESTDIR}/ints"
mkdir "${TESTDIR}/keyfiles"
mkdir "${TESTDIR}/fits"

cd "${TESTDIR}"

echo "Copying necessary files to test directory ..."

cp "${ROOTPWA}/userAnalysisWorkspace/3pi.--+/generator.conf" "./"
cp "${ROOTPWA}/userAnalysisWorkspace/3pi.--+/keyfiles/GenerateKeyfiles.sh" "./"
cp "${ROOTPWA}/userAnalysisWorkspace/3pi.--+/keyfiles/template.key" "./"
cp "${ROOTPWA}/userAnalysisWorkspace/3pi.--+/keyfiles/wavelist.compass.2008.88waves" "./"

cp "${ROOTPWA}/pyInterface/rootpwa.config" "./"


echo "Generating phase space ..."

### create variables 

${ROOTPWA}/pyInterface/genpw.py -n $NMB_PS_EVENTS -p "${PARTICLE_DATA_TABLE}" -M $MASS -B $BINWIDTH "./generator.conf" -o "./data/phase_space_MASS_${MASS}-$((MASS+BINWIDTH))_N_${NMB_PS_EVENTS}.root"

DESTINATION_DIR="${TESTDIR}/keyfiles"
export PARTICLE_DATA_TABLE
export DESTINATION_DIR
./GenerateKeyfiles.sh 

rm ${TESTDIR}/keyfiles/*[f0_980_0bw=[*.key
rm ${TESTDIR}/keyfiles/*[f0_980_0=[*.key

${ROOTPWA}/pyInterface/userInterface/createFileManager.py


${ROOTPWA}/pyInterface/userInterface/calcAmplitudes.py -e generated

${ROOTPWA}/pyInterface/userInterface/calcIntegrals.py -e generated

${ROOTPWA}/pyInterface/genPseudoData.py "./generator.conf" "${ROOTPWA}/pyInterface/mcTest/reference_fit/bin65_c2pap_bestfits_converged_MASS_1800_1820_N45340.root" "./ints/integral_binID-0_2.root" "./weighted_mc_data/weighted_pseudoData_MASS_${MASS}-$((MASS+BINWIDTH))_N_${NMB_PSEUDO_EVENTS}.root" -s 654321 -n ${NMB_PSEUDO_EVENTS} -M ${MASS} -B ${BINWIDTH} 


${ROOTPWA}/pyInterface/deWeight.py "./weighted_mc_data/weighted_pseudoData_MASS_${MASS}-$((MASS+BINWIDTH))_N_${NMB_PSEUDO_EVENTS}.root" "./data/pseudoData_MASS_${MASS}-$((MASS+BINWIDTH))_N_${NMB_PSEUDO_EVENTS}.root"

mv ./fileManager.pkl ./fileManager.pkl.old

${ROOTPWA}/pyInterface/userInterface/createFileManager.py

${ROOTPWA}/pyInterface/userInterface/calcAmplitudes.py -e real

${ROOTPWA}/pyInterface/userInterface/calcIntegrals.py -e real
