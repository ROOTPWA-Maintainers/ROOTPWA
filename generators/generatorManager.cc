
#include <map>

#include <boost/assign/std/vector.hpp>
#include <libconfig.h++>

#include <TVector3.h>

#include "diffractivePhaseSpace.h"
#include "generator.hpp"
#include "generatorParameters.hpp"
#include "libConfigUtils.hpp"
#include "particleDataTable.h"
#include "reportingUtils.hpp"
#include "generatorPickerFunctions.hpp"
#include "generatorManager.h"


using namespace std;
using namespace rpwa;

bool generatorManager::_debug = false;

generatorManager::generatorManager()
	: _pickerFunction(NULL),
	  _reactionFileRead(false),
	  _generator(NULL) { };


generatorManager::~generatorManager() { };


bool generatorManager::readReactionFile(const string& fileName) {
	using namespace boost::assign;
	using namespace libconfig;

	if(_reactionFileRead) {
		printWarn << "reading reaction file twice." << endl;
	}

	printInfo << "reading reaction file '" << fileName << "'." << endl;

	Config configFile;
	if(not parseLibConfigFile(fileName, configFile, _debug)) {
		printErr << "could not read reaction file '" << fileName << "'." << endl;
		return false;
	}

	const Setting& configRoot = configFile.getRoot();

	// Read the beam settings.
	const Setting* configBeam = findLibConfigGroup(configRoot, "beam");
	if(not configBeam) {
		printErr << "'beam' section not found in reaction file." << endl;
		return false;
	} else {
		map<string, Setting::Type> mandatoryArguments;
		insert (mandatoryArguments)
			("name", Setting::TypeString)
			("momentum", Setting::TypeFloat)
			("sigma_momentum", Setting::TypeFloat)
			("DxDz", Setting::TypeFloat)
			("sigma_DxDz", Setting::TypeFloat)
			("DyDz", Setting::TypeFloat)
			("sigma_DyDz", Setting::TypeFloat);
		if(not checkIfAllVariablesAreThere(configBeam, mandatoryArguments)) {
			printErr << "'beam' section in reaction file contains errors." << endl;
			return false;
		}
		string beamParticleName;
		configBeam->lookupValue("name", beamParticleName);
		const particleProperties* beamParticle = particleDataTable::entry(beamParticleName);
		if(not beamParticle) {
			printErr << "invalid beam particle." << endl;
			return false;
		} else {
			_beam.particle = *beamParticle;
			_beam.momentum = (*configBeam)["momentum"];
			_beam.momentumSigma = (*configBeam)["sigma_momentum"];
			_beam.DxDz = (*configBeam)["DxDz"];
			_beam.DxDzSigma = (*configBeam)["sigma_DxDz"];
			_beam.DyDz = (*configBeam)["DyDz"];
			_beam.DyDzSigma = (*configBeam)["sigma_DyDz"];
			printSucc << "initialized beam parameters." << endl;
			_beam.print(printInfo);
		}
	} // Finished with the beam settings.

	// Read the target settings.
	const Setting* configTarget = findLibConfigGroup(configRoot, "target");
	if(not configTarget) {
		printErr << "no 'target' section found in reaction file." << endl;
		return false;
	} else {
		map<string, Setting::Type> mandatoryArguments;
		insert (mandatoryArguments)
			("targetParticleName", Setting::TypeString)
			("recoilParticleName", Setting::TypeString)
			("pos", Setting::TypeArray)
			("radius", Setting::TypeFloat)
			("length", Setting::TypeFloat);
		if(not checkIfAllVariablesAreThere(configTarget, mandatoryArguments)) {
			printErr << "'target' section in reaction file contains errors." << endl;
			return false;
		}
		const Setting& configTargetPos = (*configTarget)["pos"];
		if(configTargetPos.getLength() != 3) {
			printErr << "'pos' in 'target' section has to have a length of 3." << endl;
			return false;
		}
		if(not configTargetPos[0].isNumber()) {
			printErr << "'pos' in 'target' section has to be made up of numbers." << endl;
			return false;
		}
		string targetParticleName;
		string recoilParticleName;
		configTarget->lookupValue("targetParticleName", targetParticleName);
		configTarget->lookupValue("recoilParticleName", recoilParticleName);
		const particleProperties* targetParticle = particleDataTable::entry(targetParticleName);
		const particleProperties* recoilParticle = particleDataTable::entry(recoilParticleName);
		if(not (targetParticle && recoilParticle)) {
			printErr << "invalid target or recoil particle" << endl;
			return false;
		}
		_target.targetParticle = *targetParticle;
		_target.recoilParticle = *recoilParticle;
		_target.position.SetXYZ(configTargetPos[0], configTargetPos[1], configTargetPos[2]);
		_target.radius = (*configTarget)["radius"];
		_target.length = (*configTarget)["length"];
		printSucc << "initialized target parameters." << endl;
		_target.print(printInfo);
	}// Finished with the target settings.

	// Read the final state parameters.
	const Setting* configFinalState = findLibConfigGroup(configRoot, "finalstate");
	if(not configFinalState) {
		printErr << "'finalstate' section not found in reaction file." << endl;
		return false;
	} else {
		map<string, Setting::Type> mandatoryArguments;
		insert (mandatoryArguments)
/*			("mass_min", Setting::TypeFloat)
		    ("mass_max", Setting::TypeFloat)
		    ("t_slope", Setting::TypeFloat)
		    ("t_min", Setting::TypeFloat)
*/
		    ("particles", Setting::TypeArray);
		if(not checkIfAllVariablesAreThere(configFinalState, mandatoryArguments)) {
			printErr << "'finalstate' section in reaction file contains errors." << endl;
			return false;
		}
		const Setting& configFinalStateParticles = (*configFinalState)["particles"];
		if(configFinalStateParticles.getLength() < 1) {
			printErr << "'particles' in 'finalstate' section has to have at least one entry." << endl;
			return false;
		}
		if(configFinalStateParticles[0].getType() != Setting::TypeString) {
			printErr << "'particles' in 'finalstate' section has to be made up of strings." << endl;
			return false;
		}
/*
		_finalState->minimumMass = (*configFinalState)["mass_min"];
		_finalState->maximumMass = (*configFinalState)["mass_max"];
		_finalState->tPrimeSlope = (*configFinalState)["t_slope"];
		_finalState->minimumTPrime = (*configFinalState)["t_min"];
*/
		for(int i = 0; i < configFinalStateParticles.getLength(); ++i) {
			string particleName = configFinalStateParticles[i];
			const particleProperties* particle = particleDataTable::entry(particleName);
			if(not particle) {
				printErr << "invalid final state particle" << endl;
				return false;
			}
			_finalState.particles.push_back(*particle);
		}
		printSucc << "initialized final state parameters." << endl;
		_finalState.print(printInfo);
	} // Finished final state parameters.

	// Get t'- and m-dependence.
	const Setting* configTAndMDependence = findLibConfigGroup(configRoot, "t_and_m_dependence");
	if(not configTAndMDependence) {
		printErr << "'t_and_m_dependence' section not found in reaction file." << endl;
		return false;
	} else {
		map<string, Setting::Type> mandatoryArguments;
		insert (mandatoryArguments)
			("function", Setting::TypeString)
			("settings", Setting::TypeGroup);
		if(not checkIfAllVariablesAreThere(configTAndMDependence, mandatoryArguments)) {
			printErr << "'configTAndMDependence' section in reaction file contains errors." << endl;
			return false;
		}
		const Setting& settings = (*configTAndMDependence)["settings"];
		string functionName;
		configTAndMDependence->lookupValue("function", functionName);
		if(functionName == "uniformMassExponentialT") {
			_pickerFunction = new uniformMassExponentialTPicker();
		} else if(functionName == "blabla") {
			//nothing to see here, move along
		} else {
			printErr << "'function' name '" << functionName << "' unknown." << endl;
			return false;
		}
		if(not _pickerFunction->init(settings)) {
			printErr << "Could not initialize 'function' " << functionName << "." << endl;
			delete _pickerFunction;
			_pickerFunction = NULL;
			return false;
		}
		printSucc << "initialized t' and mass dependence '" << functionName << "'." << endl;
		_pickerFunction->print(printInfo);
	} // Finished with t'- and m-dependence.

	printSucc << "read reaction file '" << fileName << "'." << endl;
	_reactionFileRead = true;
	return true;

}


bool generatorManager::initializeGenerator() {

	printInfo << "initializing event generator." << endl;

	if(not _reactionFileRead) {
		printErr << "trying to initialize generator before reading reaction file." << endl;
		return false;
	}
	if(_generator != NULL) {
		printWarn << "generator already initialized. Overwriting old generator." << endl;
		delete _generator;
		_generator = NULL;
	}
	_generator = new diffractivePhaseSpace();
	_generator->setBeam(_beam);
	_generator->setTarget(_target);

	return false;

}
