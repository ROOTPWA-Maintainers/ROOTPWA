#!/bin/bash

### BEGIN CONFIGURATION ###

	### BEGIN PHASE SPACE ###

		# mass bin (MeV)
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


### BEGIN PREPARATION ###

if [ -d "${ROOTPWA}" ]; then
	echo "Installation of RootPwa found in: ${ROOTPWA}"
	echo "Starting MC test ..."
else
	echo "No installation of RootPwa found. Aborting..."
	exit 1
fi

TESTDIR=${PWD}/ROOTPWA_MC_TEST

if [ -d "${TESTDIR}" ]; then
	echo "Test directory ${TESTDIR} already exists. Aborting..."
	exit 1
fi

echo "Creating directories ..."

mkdir "${TESTDIR}"
mkdir "${TESTDIR}/data"
mkdir "${TESTDIR}/weighted_mc_data"
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
cp "${ROOTPWA}/userAnalysisWorkspace/3pi.--+/keyfiles/wavelist.compass.2008.88waves" "./"
cp "${ROOTPWA}/pyInterface/rootpwa.config" "./"
cp "${ROOTPWA}/pyInterface/mcTest/reference_fit/bin65_c2pap_bestfits_converged_MASS_1800_1820_N45340.root.example" "./reference_fit/bin65_c2pap_bestfits_converged_MASS_1800_1820_N45340.root"

# make sure the integral binning in the config file is set correctly
sed -i.bak 's/^integralBinning.*$/integralBinning                        = [ { "mass": (1.8, 1.82) } ]/' rootpwa.config
rm -f rootpwa.config.bak

#-- END PREPARATION --#

### BEGIN MONTE CARLO GENERATION ###

# generate phase space data
echo "Generating phase space ..."
if ! ${ROOTPWA}/build/bin/genpw.py -s $SEED_PS -n $NMB_PS_EVENTS -p "${PARTICLE_DATA_TABLE}" -M $MASS -B $BINWIDTH "./generator_noBeamSimulation.conf" -o "./data/phase_space_MASS_${MASS}-$((MASS+BINWIDTH))_N_${NMB_PS_EVENTS}.root"; then
	echo "Generation of phase space was not successful. Aborting..."
	exit 1
fi

# generate the keyfiles
echo "Generating keyfiles ..."
DESTINATION_DIR="${TESTDIR}/keyfiles"
export PARTICLE_DATA_TABLE
export DESTINATION_DIR
export WAVESET_FILES="wavelist.compass.2008.88waves"
export TEMPLATE_KEY_FILES="${ROOTPWA}/userAnalysisWorkspace/3pi.--+/keyfiles/template.key"
if ! ${ROOTPWA}/userAnalysisWorkspace/3pi.--+/keyfiles/GenerateKeyfiles.sh; then
	echo "Generation of keyfiles was not successful. Aborting..."
	exit 1
fi

# remove some parameterizations
rm ${TESTDIR}/keyfiles/*[f0_980_0bw=[*.key
rm ${TESTDIR}/keyfiles/*[f0_980_0=[*.key

# create first file manager
echo "Create first file manager ..."
if ! ${ROOTPWA}/build/bin/createFileManager.py; then
	echo "Creation of first file manager was not successful. Aborting..."
	exit 1
fi

# calculate amplitudes for generated data
echo "Calculate amplitudes for generated data ..."
if ! ${ROOTPWA}/build/bin/calcAmplitudes.py -e generated; then
	echo "Calculation of amplitudes for generated data was not successful. Aborting..."
	exit 1
fi

# calculate integrals for generated data
echo "Calculate integrals for generated data ..."
if ! ${ROOTPWA}/build/bin/calcIntegrals.py -e generated; then
	echo "Calculation of integrals for generated was not successful. Aborting..."
	exit 1
fi

# generate weighted MC pseudo data
echo "Generate weighted MC pseudo data ..."
if ! ${ROOTPWA}/build/bin/genPseudoData.py "./generator_noBeamSimulation.conf" "${TESTDIR}/reference_fit/bin65_c2pap_bestfits_converged_MASS_1800_1820_N45340.root" "./ints/integral_binID-0_2.root" "./weighted_mc_data/weighted_pseudoData_MASS_${MASS}-$((MASS+BINWIDTH))_N_${NMB_PSEUDO_EVENTS}.root" -s $SEED_PSEUDO -n ${NMB_PSEUDO_EVENTS} -M ${MASS} -B ${BINWIDTH}; then 
	echo "Generation of weighted MC pseudo data was not successful. Aborting..."
	exit 1
fi

# deweight the pseudo data
echo "Deweight pseudo data ..."
if ! ${ROOTPWA}/build/bin/deWeight.py "./weighted_mc_data/weighted_pseudoData_MASS_${MASS}-$((MASS+BINWIDTH))_N_${NMB_PSEUDO_EVENTS}.root" "./data/pseudoData_MASS_${MASS}-$((MASS+BINWIDTH))_N_${NMB_PSEUDO_EVENTS}.root" -s $SEED_DEWEIGHT; then
	echo "Deweighting of weighted MC pseudo data was not successful. Aborting..."
	exit 1
fi

# backup old file manager ...
echo "Move first file manager (backup) ..."
mv ./fileManager.pkl ./fileManager.pkl.old

# ... create a new file manager
echo "Create second file manager ..."
if ! ${ROOTPWA}/build/bin/createFileManager.py; then
	echo "Creation of second file manager was not successful. Aborting..."
	exit 1
fi

# calculate amplitudes for 'real' (= MC) data
echo "Calculate amplitudes for 'real' (= MC) data ..."
if ! ${ROOTPWA}/build/bin/calcAmplitudes.py -e real; then
	echo "Calculation of amplitudes for 'real' (= MC) data  was not successful. Aborting..."
	exit 1
fi

#-- END MONTE CARLO GENERATION --#

### BEGIN FIT TEST ###

echo "Test pwaFit.py without prior ..."
if ! ${ROOTPWA}/build/bin/pwaFit.py "./fits/pwaTest_NONLOPT_NOPRIOR.root" --noAcceptance -w wavelist.compass.2008.88waves -s $SEED_FIT; then
	echo "pwaFit was not successful. Aborting..."
	exit 1
fi

echo "Test pwaFit.py with prior ..."
if ! ${ROOTPWA}/build/bin/pwaFit.py "./fits/pwaTest_NONLOPT_CAUCHY_PRIOR_WIDTH_0.5.root" --noAcceptance -w wavelist.compass.2008.88waves -C -P 0.5 -s $SEED_FIT; then
	echo "pwaFit with prior was not successful. Aborting..."
	exit 1
fi

echo "Test pwaNloptFit.py without prior ..."
if ! ${ROOTPWA}/build/bin/pwaNloptFit.py "./fits/pwaTest_NLOPT_NOPRIOR.root" --noAcceptance -w wavelist.compass.2008.88waves -s $SEED_FIT; then
	echo "pwaNloptFit was not successful. Aborting..."
	exit 1
fi

echo "Test pwaNloptFit.py with prior ..."
if ! ${ROOTPWA}/build/bin/pwaNloptFit.py "./fits/pwaTest_NLOPT_CAUCHY_PRIOR_WIDTH_0.5.root" --noAcceptance -w wavelist.compass.2008.88waves -C -P 0.5 -s $SEED_FIT; then
	echo "pwaNloptFit with prior was not successful. Aborting..."
	exit 1
fi

#-- END FIT TEST --#

exit 0
