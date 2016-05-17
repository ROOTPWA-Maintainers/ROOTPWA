#!/bin/bash

### BEGIN CONFIGURATION ###

	### BEGIN PHASE SPACE ###

		# mass bin (MeV)
		MASS=1800
		BINWIDTH=20

		# number of phase space events to generate 
		NMB_PS_EVENTS=500	

	### END PHASE SPACE ###

	### BEGIN PARTICLE DATA ###

		# particle data table
		PARTICLE_DATA_TABLE="${ROOTPWA}/particleData/particleDataTable.txt"

	### END PARTICLE DATA ###

	### BEGIN PSEUDO DATA ###

		# number of weighted pseudo data events to generate
		NMB_PSEUDO_EVENTS=500

	### END PSEUDO DATA ###


#-- END CONFIGURATION --#


### BEGIN PREPARATION ###

if [ -d "${ROOTPWA}" ]; then
	echo "RootPwa installation found in: ${ROOTPWA}"
	echo "Starting MC test ..."
else
	echo "Error no installation of RootPwa found! Aborting."
	exit 1
fi

TESTDIR=${PWD}/ROOTPWA_MC_TEST

if [ -d "${TESTDIR}" ]; then
	echo "Test directory ${TESTDIR} already exists! Aborting."
	exit 1
fi

echo "Creating directories ..."

mkdir "${TESTDIR}"
mkdir "${TESTDIR}/data"
mkdir "${TESTDIR}/weighted_mc_data"
mkdir "${TESTDIR}/amps"
mkdir "${TESTDIR}/ints"
mkdir "${TESTDIR}/keyfiles"
mkdir "${TESTDIR}/reference_fit"
mkdir "${TESTDIR}/fits"
mkdir "${TESTDIR}/log"

echo "Test directory and subfolders created."
echo "Changing to test directory."

cd "${TESTDIR}"

exec > >(tee -i ./log/logfile.txt)
exec 2>&1

echo "Copying necessary files to test directory ..."

cp "${ROOTPWA}/userAnalysisWorkspace/3pi.--+/generator_noBeamSimulation.conf" "./"
cp "${ROOTPWA}/userAnalysisWorkspace/3pi.--+/keyfiles/GenerateKeyfiles.sh" "./"
cp "${ROOTPWA}/userAnalysisWorkspace/3pi.--+/keyfiles/template.key" "./"
cp "${ROOTPWA}/userAnalysisWorkspace/3pi.--+/keyfiles/wavelist.compass.2008.88waves" "./"
cp "${ROOTPWA}/pyInterface/rootpwa.config" "./"
cp "${ROOTPWA}/pyInterface/mcTest/reference_fit/bin65_c2pap_bestfits_converged_MASS_1800_1820_N45340.root.example" "./reference_fit/bin65_c2pap_bestfits_converged_MASS_1800_1820_N45340.root"

#-- END PREPARATION --#

### BEGIN MONTE CARLO GENERATION ###

# generate phase space data
echo "Generating phase space ..."
${ROOTPWA}/pyInterface/genpw.py -n $NMB_PS_EVENTS -p "${PARTICLE_DATA_TABLE}" -M $MASS -B $BINWIDTH "./generator_noBeamSimulation.conf" -o "./data/phase_space_MASS_${MASS}-$((MASS+BINWIDTH))_N_${NMB_PS_EVENTS}.root"

# generate the keyfiles
echo "Generating keyfiles ..."
DESTINATION_DIR="${TESTDIR}/keyfiles"
export PARTICLE_DATA_TABLE
export DESTINATION_DIR
./GenerateKeyfiles.sh 

# remove some parameterizations
rm ${TESTDIR}/keyfiles/*[f0_980_0bw=[*.key
rm ${TESTDIR}/keyfiles/*[f0_980_0=[*.key

# create first file manager
echo "Create first file manager ..."
${ROOTPWA}/pyInterface/userInterface/createFileManager.py

# calculate amplitudes for generated data
echo "Calculate amplitudes for generated data ..."
${ROOTPWA}/pyInterface/userInterface/calcAmplitudes.py -e generated

# calculate integrals for generated data
echo "Calculate integrals for generated data ..."
${ROOTPWA}/pyInterface/userInterface/calcIntegrals.py -e generated

# generate weighted MC pseudo data
echo "Generate weighted MC pseudo data ..."
${ROOTPWA}/pyInterface/genPseudoData.py "./generator_noBeamSimulation.conf" "${TESTDIR}/reference_fit/bin65_c2pap_bestfits_converged_MASS_1800_1820_N45340.root" "./ints/integral_binID-0_2.root" "./weighted_mc_data/weighted_pseudoData_MASS_${MASS}-$((MASS+BINWIDTH))_N_${NMB_PSEUDO_EVENTS}.root" -s 654321 -n ${NMB_PSEUDO_EVENTS} -M ${MASS} -B ${BINWIDTH} 

# deweight the pseudo data
echo "Deweight pseudo data ..."
${ROOTPWA}/pyInterface/deWeight.py "./weighted_mc_data/weighted_pseudoData_MASS_${MASS}-$((MASS+BINWIDTH))_N_${NMB_PSEUDO_EVENTS}.root" "./data/pseudoData_MASS_${MASS}-$((MASS+BINWIDTH))_N_${NMB_PSEUDO_EVENTS}.root"

# backup old file manager ...
echo "Move first file manager (backup) ..."
mv ./fileManager.pkl ./fileManager.pkl.old

# ... create a new file manager
echo "Create second file manager ..."
${ROOTPWA}/pyInterface/userInterface/createFileManager.py

# calculate amplitudes for 'real' (= MC) data
echo "Calculate amplitudes for 'real' (= MC) data ..."
${ROOTPWA}/pyInterface/userInterface/calcAmplitudes.py -e real

# 'fake' accepeted integrals -> acceptance = 1
echo "'Fake' accepted integrals ..."
ln -s "${TESTDIR}/ints/integral_binID-0_2.root" "${TESTDIR}/ints/integral_binID-0_3.root"

#-- END MONTE CARLO GENERATION --#

### BEGIN FIT TEST ###

echo "Test pwaFit.py without prior ..."
${ROOTPWA}/pyInterface/userInterface/pwaFit.py "./fits/pwaTest_NONLOPT_NOPRIOR.root" -w wavelist.compass.2008.88waves

echo "Test pwaFit.py with prior ..."
${ROOTPWA}/pyInterface/userInterface/pwaFit.py "./fits/pwaTest_NONLOPT_CAUCHY_PRIOR_WIDTH_0.5.root" -w wavelist.compass.2008.88waves -C -P 0.5

echo "Test pwaNloptFit.py without prior ..."
${ROOTPWA}/pyInterface/userInterface/pwaNloptFit.py "./fits/pwaTest_NLOPT_NOPRIOR.root" -w wavelist.compass.2008.88waves

echo "Test pwaNloptFit.py with prior ..."
${ROOTPWA}/pyInterface/userInterface/pwaNloptFit.py "./fits/pwaTest_NLOPT_CAUCHY_PRIOR_WIDTH_0.5.root" -w wavelist.compass.2008.88waves -C -P 0.5

#-- END FIT TEST --#

exit 0
