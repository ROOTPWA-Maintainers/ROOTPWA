
#include <map>

#include <boost/assign/std/vector.hpp>
#include <libconfig.h++>

#include <TVector3.h>

#include "libConfigUtils.hpp"
#include "particleDataTable.h"
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
		map<string, Setting::Type> mandatoryArguments;
		insert (mandatoryArguments)
			("active", Setting::TypeBoolean)
			("name", Setting::TypeString)
			("momentum", Setting::TypeFloat)
			("sigma_momentum", Setting::TypeFloat)
			("DxDz", Setting::TypeFloat)
			("sigma_DxDz", Setting::TypeFloat)
			("DyDz", Setting::TypeFloat)
			("sigma_DyDz", Setting::TypeFloat);
		if(not checkIfAllVariablesAreThere(configBeam, mandatoryArguments)) {
			printWarn << " Beam simulation disabled." << endl;
		} else if((*configBeam)["active"]) {
			_beam = new Beam();
			string beamParticleName;
			configBeam->lookupValue("name", beamParticleName);
			const particleProperties* beamParticle = particleDataTable::entry(beamParticleName);
			if(not beamParticle) {
				printWarn << "invalid beam particle." << endl;
				delete _beam;
				_beam = NULL;
			} else {
				_beam->particle = *beamParticle;
				_beam->momentum = (*configBeam)["momentum"];
				_beam->momentumSigma = (*configBeam)["sigma_momentum"];
				_beam->DxDz = (*configBeam)["DxDz"];
				_beam->DxDzSigma = (*configBeam)["sigma_DxDz"];
				_beam->DyDz = (*configBeam)["DyDz"];
				_beam->DyDzSigma = (*configBeam)["sigma_DyDz"];
				printSucc << "initialized beam parameters." << endl;
				_beam->print(printInfo);
			}
		} else {
			printInfo << "beam simulation disabled by 'active' flag." << endl;
		}
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
		string targetParticleName;
		configTarget->lookupValue("name", targetParticleName);
		const particleProperties* targetParticle = particleDataTable::entry(targetParticleName);
		if(not targetParticle) {
			printErr << "invalid target particle" << endl;
			delete _target;
			_target = NULL;
			return false;
		}
		_target->particle = *targetParticle;
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
	} else {
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
			string particleName = configFinalStateParticles[i];
			const particleProperties* particle = particleDataTable::entry(particleName);
			if(not particle) {
				printErr << "invalid final state particle" << endl;
				delete _finalState;
				_finalState = NULL;
				return false;
			}
			_finalState->particles.push_back(*particle);
		}
		printSucc << "initialized final state parameters." << endl;
		_finalState->print(printInfo);
	} // Finished final state parameters.

	printSucc << "read reaction file '" << fileName << "'." << endl;
	return true;

}
