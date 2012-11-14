
#include <map>

#include <boost/assign/std/vector.hpp>
#include <libconfig.h++>

#include <TVector3.h>

#include "libConfigUtils.hpp"
#include "reportingUtils.hpp"
#include "generatorManager.h"

using namespace std;
using namespace rpwa;

bool generatorManager::_debug = false;

generatorManager::generatorManager()
	: _beam(NULL),
	  _target(NULL),
	  _finalState(NULL) { };

generatorManager::~generatorManager() { };

bool generatorManager::readReactionFile(const string& fileName) {
	using namespace boost::assign;
	using namespace libconfig;

	printInfo << "reading reaction file '" << fileName << "'." << endl;

	Config configFile;
	if(not parseLibConfigFile(fileName, configFile, _debug)) {
		printErr << "could not read reaction file '" << fileName << "'." << endl;
		return false;
	}

	const Setting& configRoot = configFile.getRoot();

	// Read the beam settings.
	const Setting* configBeam = findLibConfigGroup(configRoot, "beam");
	if(configBeam) {
		bool active = true;
		int activeFlag;
		if(not configBeam->lookupValue("active", activeFlag)) {
			printWarn << "'active' not found in 'beam' section, assuming beam simulation is enabled." << endl;
		} else {
			if(activeFlag == 1) {
				printInfo << "beam simulation enabled." << endl;
			} else if (activeFlag == 0) {
				printInfo << "beam simulation disabled." << endl;
				active = false;
			} else {
				printWarn << "'active' flag in 'beam' section should be 0 or 1, found '"
				          << activeFlag << "' instead. Beam simulation disabled."
				          << endl;
				active = false;
			}
		}
		if(active) {
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
				printWarn << " Beam simulation disabled." << endl;
			} else {
				_beam = new Beam();
				configBeam->lookupValue("name", _beam->particleName);
				_beam->momentum = (*configBeam)["momentum"];
				_beam->momentumSigma = (*configBeam)["sigma_momentum"];
				_beam->DxDz = (*configBeam)["DxDz"];
				_beam->DxDzSigma = (*configBeam)["sigma_DxDz"];
				_beam->DyDz = (*configBeam)["DyDz"];
				_beam->DyDzSigma = (*configBeam)["sigma_DyDz"];
				printSucc << "initialized beam parameters." << endl;
				_beam->print(printInfo);
			}
		} // end if(active)
	} else {
		printWarn << "no 'beam' section found in configuration file. "
		          << "Beam simulation disabled." << endl;
	} // Finished with the beam settings.

	// Read the target settings.
	const Setting* configTarget = findLibConfigGroup(configRoot, "target");
	if(not configTarget) {
		printErr << "no 'target' section found in configuration file." << endl;
		return false;
	} else {
		map<string, Setting::Type> mandatoryArguments;
		insert (mandatoryArguments)
			("name", Setting::TypeString)
			("pos", Setting::TypeArray)
			("radius", Setting::TypeFloat)
			("length", Setting::TypeFloat);
		if(not checkIfAllVariablesAreThere(configTarget, mandatoryArguments)) {
			printErr << "'target' section in configuration file contains errors." << endl;
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
		_target = new Target();
		configTarget->lookupValue("name", _target->particleName);
		_target->position.SetXYZ(configTargetPos[0], configTargetPos[1], configTargetPos[2]);
		_target->radius = (*configTarget)["radius"];
		_target->length = (*configTarget)["length"];
		printSucc << "initialized target parameters." << endl;
		_target->print(printInfo);
	}// Finished with the target settings.

	// Read the final state parameters.
	const Setting* configFinalState = findLibConfigGroup(configRoot, "finalstate");
	if(not configFinalState) {
		printErr << "'finalstate' section not found in configuration file." << endl;
		return false;
	}
	map<string, Setting::Type> mandatoryArguments;
	insert (mandatoryArguments)
		("mass_min", Setting::TypeFloat)
		("mass_max", Setting::TypeFloat)
		("t_slope", Setting::TypeFloat)
		("t_min", Setting::TypeFloat)
		("particles", Setting::TypeArray);
	if(not checkIfAllVariablesAreThere(configFinalState, mandatoryArguments)) {
		printErr << "'finalstate' section in configuration file contains errors." << endl;
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
	_finalState = new FinalState();
	_finalState->minimumMass = (*configFinalState)["mass_min"];
	_finalState->maximumMass = (*configFinalState)["mass_max"];
	_finalState->tPrimeSlope = (*configFinalState)["t_slope"];
	_finalState->minimumTPrime = (*configFinalState)["t_min"];
	for(int i = 0; i < configFinalStateParticles.getLength(); ++i) {
		_finalState->particleNames.push_back(configFinalStateParticles[i]);
	}
	printSucc << "initialized final state parameters." << endl;
	_finalState->print(printInfo);
	// Finished final state parameters.

	printSucc << "read reaction file '" << fileName << "'." << endl;
	return true;

}
