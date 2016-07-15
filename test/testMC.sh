#!/bin/bash

### BEGIN CONFIGURATION ###

		TESTDIR=$(pwd)/ROOTPWA_MC_TEST

	### BEGIN PHASE SPACE ###

		# mass bin [MeV/c^2]
		MASS=1800
		BINWIDTH=20

		# number of phase space events to generate
		NMB_PS_EVENTS=500

		SEED_PS=123456

	#-- END PHASE SPACE --#

	### BEGIN PARTICLE DATA ###

		# particle data table
		PARTICLE_DATA_TABLE="${ROOTPWA}/particleData/particleDataTable.txt"

	#-- END PARTICLE DATA --#

	### BEGIN PSEUDO DATA ###

		# number of weighted pseudo data events to generate
		NMB_PSEUDO_EVENTS=500

		SEED_PSEUDO=654321

	#-- END PSEUDO DATA --#

	### BEGIN DEWEIGHT ###

		SEED_DEWEIGHT=789012

	#-- END DEWEIGHT --#

	### BEGIN FIT ###

		SEED_FIT=210987

	#-- END FIT --#

#-- END CONFIGURATION --#


### SOME CONVENIENCE FUNCTIONS

THIS_SCRIPT=$(basename ${0})

function printInfo {
		echo ">>> ${THIS_SCRIPT}: info: ${1}"
}

function printErr {
		echo "!!! ${THIS_SCRIPT}: ${BASH_LINENO[0]}: error: ${1}"
		exit 1
}


### BEGIN PREPARATION ###

if [ -d "${ROOTPWA}" ]; then
	printInfo "Installation of ROOTPWA found in: ${ROOTPWA}"
	echo "    Starting MC test ..."
else
	printErr "No installation of ROOTPWA found. Aborting."
fi

if [ -d "${TESTDIR}" ]; then
	printErr "Test directory ${TESTDIR} already exists. Aborting."
fi

printInfo "Creating directories ..."

mkdir --verbose "${TESTDIR}"
mkdir --verbose "${TESTDIR}/data"
mkdir --verbose "${TESTDIR}/weighted_mc_data"
mkdir --verbose "${TESTDIR}/keyfiles"
mkdir --verbose "${TESTDIR}/reference_fit"
mkdir --verbose "${TESTDIR}/fits"
mkdir --verbose "${TESTDIR}/log"

echo "    Test directory and subfolders created."
printInfo "Changing to test directory ${TESTDIR}."

cd "${TESTDIR}"

exec > >(tee -i ./log/logfile.txt)
exec 2>&1

printInfo "Copying necessary files to test directory ..."

cp --verbose "${ROOTPWA}/userAnalysisWorkspace/3pi.--+/generator_noBeamSimulation.conf" "./"
cp --verbose "${ROOTPWA}/userAnalysisWorkspace/3pi.--+/keyfiles/wavelist.compass.2008.88waves" "./"
cp --verbose "${ROOTPWA}/rootpwa.config" "./"
cp --verbose "${ROOTPWA}/test/mcTest/reference_fit/bin65_c2pap_bestfits_converged_MASS_1800_1820_N45340.root.example" "./reference_fit/bin65_c2pap_bestfits_converged_MASS_1800_1820_N45340.root"

# make sure the integral binning in the config file is set correctly
sed -i 's/^integralBinning.*$/integralBinning                        = [ { "mass": (1.8, 1.82) } ]/' rootpwa.config

#-- END PREPARATION --#

### BEGIN MONTE CARLO GENERATION ###

# generate phase space data
printInfo "Generating phase space ..."
if ! ${ROOTPWA}/build/bin/genpw -s ${SEED_PS} -n ${NMB_PS_EVENTS} -p "${PARTICLE_DATA_TABLE}" -M ${MASS} -B ${BINWIDTH} "./generator_noBeamSimulation.conf" -o "./data/phase_space_MASS_${MASS}-$((MASS+BINWIDTH))_N_${NMB_PS_EVENTS}.root"; then
	printErr "Generation of phase space was not successful. Aborting."
fi

# generate the keyfiles
printInfo "Generating keyfiles ..."
DESTINATION_DIR="${TESTDIR}/keyfiles"
export PARTICLE_DATA_TABLE
export DESTINATION_DIR
export WAVESET_FILES="wavelist.compass.2008.88waves"
export TEMPLATE_KEY_FILES="${ROOTPWA}/userAnalysisWorkspace/3pi.--+/keyfiles/template.key"
if ! ${ROOTPWA}/userAnalysisWorkspace/3pi.--+/keyfiles/GenerateKeyfiles.sh; then
	printErr "Generation of keyfiles was not successful. Aborting."
fi

# remove some parameterizations
rm ${TESTDIR}/keyfiles/*[f0_980_0bw=[*.key
rm ${TESTDIR}/keyfiles/*[f0_980_0=[*.key

# create first file manager
printInfo "Create first file manager ..."
if ! ${ROOTPWA}/build/bin/createFileManager; then
	printErr "Creation of first file manager was not successful. Aborting."
fi

# calculate amplitudes for generated data
printInfo "Calculate amplitudes for generated data ..."
if ! ${ROOTPWA}/build/bin/calcAmplitudes -e generated; then
	printErr "Calculation of amplitudes for generated data was not successful. Aborting."
fi

# calculate integrals for generated data
printInfo "Calculate integrals for generated data ..."
if ! ${ROOTPWA}/build/bin/calcIntegrals -e generated; then
	printErr "Calculation of integrals for generated was not successful. Aborting."
fi

# generate weighted MC pseudo data
printInfo "Generate weighted MC pseudo data ..."
if ! ${ROOTPWA}/build/bin/genPseudoData "./generator_noBeamSimulation.conf" "${TESTDIR}/reference_fit/bin65_c2pap_bestfits_converged_MASS_1800_1820_N45340.root" "./ints/integral_binID-0_2.root" "./weighted_mc_data/weighted_pseudoData_MASS_${MASS}-$((MASS+BINWIDTH))_N_${NMB_PSEUDO_EVENTS}.root" -s ${SEED_PSEUDO} -n ${NMB_PSEUDO_EVENTS} -M ${MASS} -B ${BINWIDTH}; then
	printErr "Generation of weighted MC pseudo data was not successful. Aborting."
fi

# deweight the pseudo data
printInfo "Deweight pseudo data ..."
if ! ${ROOTPWA}/build/bin/deWeight "./weighted_mc_data/weighted_pseudoData_MASS_${MASS}-$((MASS+BINWIDTH))_N_${NMB_PSEUDO_EVENTS}.root" "./data/pseudoData_MASS_${MASS}-$((MASS+BINWIDTH))_N_${NMB_PSEUDO_EVENTS}.root" -s ${SEED_DEWEIGHT}; then
	printErr "Deweighting of weighted MC pseudo data was not successful. Aborting."
fi

# backup old file manager ...
printInfo "Move first file manager (backup) ..."
mv ./fileManager.pkl ./fileManager.pkl.old

# ... create a new file manager
printInfo "Create second file manager ..."
if ! ${ROOTPWA}/build/bin/createFileManager; then
	printErr "Creation of second file manager was not successful. Aborting."
fi

# calculate amplitudes for 'real' (= MC) data
printInfo "Calculate amplitudes for 'real' (= MC) data ..."
if ! ${ROOTPWA}/build/bin/calcAmplitudes -e real; then
	printErr "Calculation of amplitudes for 'real' (= MC) data  was not successful. Aborting."
fi

#-- END MONTE CARLO GENERATION --#

### BEGIN FIT TEST ###

printInfo "Test pwaFit without prior ..."
if ! ${ROOTPWA}/build/bin/pwaFit "./fits/pwaTest_NONLOPT_NOPRIOR.root" --noAcceptance -w wavelist.compass.2008.88waves -s ${SEED_FIT}; then
	printErr "pwaFit was not successful. Aborting."
fi

printInfo "Test pwaFit with prior ..."
if ! ${ROOTPWA}/build/bin/pwaFit "./fits/pwaTest_NONLOPT_CAUCHY_PRIOR_WIDTH_0.5.root" --noAcceptance -w wavelist.compass.2008.88waves -C -P 0.5 -s ${SEED_FIT}; then
	printErr "pwaFit with prior was not successful. Aborting."
fi

printInfo "Test pwaNloptFit without prior ..."
if ! ${ROOTPWA}/build/bin/pwaNloptFit "./fits/pwaTest_NLOPT_NOPRIOR.root" --noAcceptance -w wavelist.compass.2008.88waves -s ${SEED_FIT}; then
	printErr "pwaNloptFit was not successful. Aborting."
fi

printInfo "Test pwaNloptFit with prior ..."
if ! ${ROOTPWA}/build/bin/pwaNloptFit "./fits/pwaTest_NLOPT_CAUCHY_PRIOR_WIDTH_0.5.root" --noAcceptance -w wavelist.compass.2008.88waves -C -P 0.5 -s ${SEED_FIT}; then
	printErr "pwaNloptFit with prior was not successful. Aborting."
fi

#-- END FIT TEST --#

exit 0
