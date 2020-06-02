
#include <map>

#include <boost/assign/std/vector.hpp>
#include <libconfig.h++>

#include <TVector3.h>

#include "diffractivePhaseSpace.h"
#include "generator.h"
#include "generatorParameters.hpp"
#include "libConfigUtils.hpp"
#include "particleDataTable.h"
#include "beamAndVertexGenerator.h"
#include "reportingUtils.hpp"
#include "generatorPickerFunctions.h"
#include "generatorManager.h"


using namespace boost;
using namespace std;
using namespace rpwa;

bool generatorManager::_debug = false;

generatorManager::generatorManager()
	: _beamAndVertexGenerator(beamAndVertexGeneratorPtr(new beamAndVertexGenerator())),
	  _pickerFunction(massAndTPrimePickerPtr()),
	  _beamFileName(""),
	  _reactionFileRead(false),
	  _generator(NULL) { };


generatorManager::~generatorManager() {
	delete _generator;
}


#ifdef USE_BAT
rpwa::importanceSamplerPtr generatorManager::getImportanceSampler(rpwa::modelIntensityPtr model) {

	if(not _reactionFileRead) {
		printErr << "trying to initialize importance sampler before reading reaction file." << endl;
		return rpwa::importanceSamplerPtr();
	}

	rpwa::importanceSamplerPtr sampler(new importanceSampler(model,
	                                                         _beamAndVertexGenerator,
	                                                         _pickerFunction,
	                                                         _beam,
	                                                         _target,
	                                                         _finalState));

	return sampler;

}
#endif


unsigned int generatorManager::event() {
	if(not _reactionFileRead) {
		printErr << "cannot generate event before reading the reaction file." << endl;
		throw;
	}
	return _generator->event();
}


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
			("momentumSigma", Setting::TypeFloat)
			("DxDz", Setting::TypeFloat)
			("DxDzSigma", Setting::TypeFloat)
			("DyDz", Setting::TypeFloat)
			("DyDzSigma", Setting::TypeFloat);
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
			_beam.momentumSigma = (*configBeam)["momentumSigma"];
			_beam.DxDz = (*configBeam)["DxDz"];
			_beam.DxDzSigma = (*configBeam)["DxDzSigma"];
			_beam.DyDz = (*configBeam)["DyDz"];
			_beam.DyDzSigma = (*configBeam)["DyDzSigma"];
			printSucc << "initialized beam parameters." << endl;
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
			("length", Setting::TypeFloat)
			("interactionLength", Setting::TypeFloat);
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
		if(not (targetParticle and recoilParticle)) {
			printErr << "invalid target or recoil particle" << endl;
			return false;
		}
		double interactionLength;
		configTarget->lookupValue("interactionLength", interactionLength);
		if(interactionLength < 0) {
			printErr << "interaction length has to be >0 (set it to 0 to disable this feature)." << endl;
			return false;
		}
		_target.targetParticle = *targetParticle;
		_target.recoilParticle = *recoilParticle;
		_target.position.SetXYZ(configTargetPos[0], configTargetPos[1], configTargetPos[2]);
		_target.radius = (*configTarget)["radius"];
		_target.length = (*configTarget)["length"];
		_target.interactionLength = interactionLength;
		printSucc << "initialized target parameters." << endl;
	}// Finished with the target settings.

	// Read the settings for the beam simulation.
	const Setting* configBeamSimulation = findLibConfigGroup(configRoot, "beamSimulation");
	if(not configBeamSimulation) {
		printWarn << "no 'beamSimulation' section found in reaction file. Beam package disabled." << endl;
	} else {
		map<string, Setting::Type> mandatoryArguments;
		insert (mandatoryArguments)
			("active", Setting::TypeBoolean);
		if(not checkIfAllVariablesAreThere(configBeamSimulation, mandatoryArguments)) {
			printErr << "'beamSimulation' section in reaction file contains errors." << endl;
			return false;
		} else if((*configBeamSimulation)["active"]) {
			if(_beamFileName == "") {
				if(not configBeamSimulation->lookupValue("beamFile", _beamFileName)) {
					printErr << "'beamSimulation' section in reaction file is missing the 'beamFile' entry." << endl;
					return false;
				}
			} else {
				printInfo << "Beamfile command line override '" << _beamFileName << "' found." << endl;
			}
			double sigmaScalingFactor = 1.;
			if(not configBeamSimulation->exists("sigmaScalingFactor")) {
				printWarn << "'sigmaScalingFactor' is missing from 'beamSimulation' section." << endl;
			} else {
				if(not configBeamSimulation->lookupValue("sigmaScalingFactor", sigmaScalingFactor)) {
					printErr << "Could not read 'sigmaScalingFactor' in 'beamSimulation' section." << endl;
					return false;
				} else {
					if(sigmaScalingFactor < 0.) {
						printErr << "'sigmaScalingFactor' in 'beamSimulation' section must be positive." << endl;
						return false;
					}
				}
			}

			bool takeZpositionFromData = false;
			if(configBeamSimulation->exists("takeZpositionFromData")) {
				if(not configBeamSimulation->lookupValue("takeZpositionFromData", takeZpositionFromData)) {
					printErr << "Could not read 'takeZpositionFromData' in 'beamSimulation' section." << endl;
					return false;
				}
			}

			_beamAndVertexGenerator->setSigmaScalingFactor(sigmaScalingFactor);
			if(not _beamAndVertexGenerator->loadBeamFile(_beamFileName)) {
				printErr << "could not initialize beam and vertex generator." << endl;
				return false;
			}
			_beamAndVertexGenerator->setTakeZpositionFromData(takeZpositionFromData);
			printSucc << "initialized beam package." << endl;
		} else {
			printInfo << "beam package disabled." << endl;
		}
	} // Finished with the beam simulation settings.

	// Read the final state parameters.
	const Setting* configFinalState = findLibConfigGroup(configRoot, "finalstate");
	if(not configFinalState) {
		printErr << "'finalstate' section not found in reaction file." << endl;
		return false;
	} else {
		map<string, Setting::Type> mandatoryArguments;
		insert (mandatoryArguments)
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
			_pickerFunction = massAndTPrimePickerPtr(new uniformMassExponentialTPicker());
		} else if(functionName == "polynomialMassAndTPrime") {
			_pickerFunction = massAndTPrimePickerPtr(new polynomialMassAndTPrimeSlopePicker());
		} else if(functionName == "uniformMassAndT") {
			_pickerFunction = massAndTPrimePickerPtr(new uniformMassAndTPicker());
		} else {
			printErr << "'function' name '" << functionName << "' unknown." << endl;
			return false;
		}
		if(not _pickerFunction->init(settings)) {
			printErr << "Could not initialize 'function' " << functionName << "." << endl;
			_pickerFunction = massAndTPrimePickerPtr();
			return false;
		}
		printSucc << "initialized t' and mass dependence '" << functionName << "'." << endl;
	} // Finished with t'- and m-dependence.

	// Features which have been removed during clean-up.
	if(configRoot.exists("importance")) {
		printWarn << "importance sampling removed for the time being. "
		          << "If you want it back, check SVN revision 1072." << endl;
	}
	if(configRoot.exists("resonances")) {
		printWarn << "resonances weighing removed for the time being. "
		          << "If you want it back, check SVN revision 1072." << endl;
	}

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
	_generator->setTPrimeAndMassPicker(_pickerFunction);
	_generator->setPrimaryVertexGenerator(_beamAndVertexGenerator);
	_generator->setDecayProducts(_finalState.particles);

	printSucc << "event generator initialized" << endl;
	return true;

}


void generatorManager::overrideMassRange(double lowerLimit, double upperLimit) {

	if(not _reactionFileRead) {
		printErr << "reaction file has to have been read to override the mass range." << endl;
		throw;
	}
	_pickerFunction->overrideMassRange(lowerLimit, upperLimit);

}


void generatorManager::readBeamfileSequentially(bool readBeamfileSequentially) {

	if(not _reactionFileRead) {
		printErr << "reaction file has to have been read to set this option (readBeamfileSequentially)." << endl;
		throw;
	}
	if(not _beamAndVertexGenerator) {
		printErr << "beam and vertex package seems to be disabled, unable to read beamfile sequentially." << endl;
		throw;
	}
	_beamAndVertexGenerator->setBeamfileSequentialReading(readBeamfileSequentially);

}


void generatorManager::randomizeBeamfileStartingPosition() {

	if(not _reactionFileRead) {
		printErr << "reaction file has to have been read to set this option (randomizeBeamfileStartingPosition)." << endl;
		throw;
	}
	if(not _beamAndVertexGenerator) {
		printErr << "beam and vertex package seems to be disabled, unable to randomize beamfile starting position." << endl;
		throw;
	}
	_beamAndVertexGenerator->randomizeBeamfileStartingPosition();

}


ostream& generatorManager::print(ostream& out) const {

	out << "generatorManager parameter collection:" << endl;
	_beam.print(out);
	_target.print(out);
	if(_beamAndVertexGenerator) {
		_beamAndVertexGenerator->print(out);
	}
	_finalState.print(out);
	if(_pickerFunction) {
		_pickerFunction->print(out);
	}
	return out;

}
