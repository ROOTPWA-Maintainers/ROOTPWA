///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010-2012 Sebastian Neubert (TUM)
//    Copyright 2014-2016 Sebastian Uhl (TUM)
//
//    This file is part of ROOTPWA
//
//    ROOTPWA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ROOTPWA is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ROOTPWA.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      implementation of the master class of the resonance fit
//
//-------------------------------------------------------------------------


#include "resonanceFit.h"
#include "resonanceFitInternal.h"

#include <iostream>

#include <boost/assign/std/vector.hpp>

#include <yaml-cpp/yaml.h>

#include <TFormula.h>

#include <reportingUtils.hpp>
#include <yamlCppUtils.hpp>

#include "components.h"
#include "data.h"
#include "fsmd.h"
#include "information.h"
#include "model.h"
#include "parameter.h"
#include "resonanceFitHelper.h"


std::map<std::string, double>
rpwa::resonanceFit::readFitQuality(const YAML::Node& configRoot)
{
	if(rpwa::resonanceFit::debug()) {
		printDebug << "reading 'fitquality'." << std::endl;
	}

	const YAML::Node& configFitQuality = configRoot["fitquality"];

	std::map<std::string, double> fitQuality;
	if(not configFitQuality) {
		// it is perfectly okay for a config file to not
		// contain fit quality information
		return fitQuality;
	}

	if(not configFitQuality.IsMap()) {
		printErr << "'fitquality' is not a YAML map." << std::endl;
		throw;
	}

	for(YAML::const_iterator it = configFitQuality.begin(); it != configFitQuality.end(); ++it) {
		if(not checkVariableType(it->first, rpwa::YamlCppUtils::TypeString)) {
			printErr << "entries in 'fitquality' must be pairs of 'string' and 'double'." << std::endl;
			throw;
		}

		const std::string key = it->first.as<std::string>();

		double value;
		if(checkVariableType(it->second, rpwa::YamlCppUtils::TypeFloat)) {
			value = it->second.as<double>();
		} else if(checkVariableType(it->second, rpwa::YamlCppUtils::TypeString) and it->second.as<std::string>() == "nan") {
			// some systems cannot convert the string 'nan' to a
			// floating point number, so this is to be done
			// manually.
			value = std::numeric_limits<double>::has_quiet_NaN ? std::numeric_limits<double>::quiet_NaN() : 0.0;
		} else if(checkVariableType(it->second, rpwa::YamlCppUtils::TypeString) and it->second.as<std::string>() == "inf") {
			// some systems cannot convert the string 'inf' to a
			// floating point number, so this is to be done
			// manually.
			value = std::numeric_limits<double>::has_infinity ? std::numeric_limits<double>::infinity() : 0.0;
		} else {
			printErr << "entries in 'fitquality' must be pairs of 'string' and 'double'." << std::endl;
			throw;
		}

		if(fitQuality.count(key) != 0) {
			printErr << "variable '" << key << "' of 'fitquality' given multiple times." << std::endl;
			throw;
		}
		fitQuality[key] = value;

		if(rpwa::resonanceFit::debug()) {
			printDebug << "read key '" << key << "' with value '" << value << "'." << std::endl;
		}
	}

	return fitQuality;
}


std::vector<std::string>
rpwa::resonanceFit::readFreeParameters(const YAML::Node& configRoot)
{
	if(rpwa::resonanceFit::debug()) {
		printDebug << "reading 'freeparameters'." << std::endl;
	}

	const YAML::Node& configFreeParameters = configRoot["freeparameters"];

	std::vector<std::string> freeParameters;
	if(not configFreeParameters) {
		printWarn << "release order of parameters not specified in configuration file." << std::endl;
		return freeParameters;
	}

	if(not configFreeParameters.IsSequence()) {
		printErr << "'freeparameters' is not a YAML sequence." << std::endl;
		throw;
	}

	const size_t nrItems = configFreeParameters.size();
	if(nrItems == 0) {
		printErr << "'freeparameters' is an empty sequence, when defined it must at least contain one entry." << std::endl;
		throw;
	}

	for(size_t idxItem = 0; idxItem < nrItems; ++idxItem) {
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading of entry " << idxItem << " in 'freeparameters'." << std::endl;
		}

		if(not checkVariableType(configFreeParameters[idxItem], rpwa::YamlCppUtils::TypeString)) {
			printErr << "'freeparameters' entry at index " << idxItem << " is not a string." << std::endl;
			throw;
		}

		freeParameters.push_back(configFreeParameters[idxItem].as<std::string>());

		if(rpwa::resonanceFit::debug()) {
			printDebug << "read parameters to release: '" << freeParameters.back() << "'." << std::endl;
		}
	}

	return freeParameters;
}


namespace {


	std::vector<rpwa::resonanceFit::information::bin>
	readInformationFitResults(const YAML::Node& configInput)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading 'fitresults'." << std::endl;
		}

		const YAML::Node& configInputFitResults = configInput["fitresults"];

		if(not configInputFitResults) {
			printErr << "'fitresults' does not exist in 'input'." << std::endl;
			throw;
		}
		if(not configInputFitResults.IsSequence()) {
			printErr << "'fitresults' is not a YAML sequence." << std::endl;
			throw;
		}

		std::vector<rpwa::resonanceFit::information::bin> bins;

		const size_t nrFitResults = configInputFitResults.size();
		for(size_t idxFitResult = 0; idxFitResult < nrFitResults; ++idxFitResult) {
			if(rpwa::resonanceFit::debug()) {
				printDebug << "reading of entry " << idxFitResult << " in 'fitresults'." << std::endl;
			}

			const YAML::Node& configInputFitResult = configInputFitResults[idxFitResult];

			std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
			boost::assign::insert(mandatoryArguments)
			                     ("name", rpwa::YamlCppUtils::TypeString)
			                     ("tPrimeMean", rpwa::YamlCppUtils::TypeFloat);
			if(not checkIfAllVariablesAreThere(configInputFitResult, mandatoryArguments)) {
				printErr << "'fitresults' entry at index " << idxFitResult << " does not contain all required variables." << std::endl;
				throw;
			}

			const std::string fileName = configInputFitResult["name"].as<std::string>();
			const double tPrimeMean = configInputFitResult["tPrimeMean"].as<double>();

			double rescaleErrors = 1.0;
			if(configInputFitResult["rescaleErrors"]) {
				if(checkVariableType(configInputFitResult["rescaleErrors"], rpwa::YamlCppUtils::TypeFloat)) {
					rescaleErrors = configInputFitResult["rescaleErrors"].as<double>();
				} else {
					printErr << "variable 'rescaleErrors' of 'fitresults' entry at index " << idxFitResult << " is not a floating point number." << std::endl;
					throw;
				}
			}

			if(rpwa::resonanceFit::debug()) {
				printDebug << "read file name of fit results of mass-independent fit: '" << fileName << "'." << std::endl;
				printDebug << "read mean t' value: '" << tPrimeMean << "'." << std::endl;
				printDebug << "rescale errors by factor: '" << rescaleErrors << "'." << std::endl;
			}

			std::vector<std::string> sysFileNames;
			// get information for plotting of systematic error
			const YAML::Node& configInputFitResultSystematics = configInputFitResult["systematics"];
			if(configInputFitResultSystematics) {
				if(not configInputFitResultSystematics.IsSequence()) {
					printErr << "'systematics' is not a YAML sequence." << std::endl;
					throw;
				}

				const size_t nrSystematics = configInputFitResultSystematics.size();
				if(rpwa::resonanceFit::debug()) {
					printDebug << "going to read information for " << nrSystematics << " files containing information for systematic errors." << std::endl;
				}

				for(size_t idxSystematics = 0; idxSystematics < nrSystematics; ++idxSystematics) {
					if(not checkVariableType(configInputFitResultSystematics[idxSystematics], rpwa::YamlCppUtils::TypeString)) {
						printErr << "'systematics' entry at index " << idxSystematics << " is not a string." << std::endl;
						throw;
					}

					sysFileNames.push_back(configInputFitResultSystematics[idxSystematics].as<std::string>());
				}
			}

			bins.push_back(rpwa::resonanceFit::information::bin(fileName,
			                                                    tPrimeMean,
			                                                    rescaleErrors,
			                                                    sysFileNames));

			if(rpwa::resonanceFit::debug()) {
				printDebug << bins.back() << std::endl;
			}
		}

		return bins;
	}


	std::vector<rpwa::resonanceFit::information::wave>
	readInformationWaves(const YAML::Node& configInput)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading 'waves'." << std::endl;
		}

		const YAML::Node& configInputWaves = configInput["waves"];

		if(not configInputWaves) {
			printErr << "'waves' does not exist in 'input'." << std::endl;
			throw;
		}
		if(not configInputWaves.IsSequence()) {
			printErr << "'waves' is not a YAML sequence." << std::endl;
			throw;
		}

		std::vector<rpwa::resonanceFit::information::wave> waves;

		const size_t nrWaves = configInputWaves.size();
		for(size_t idxWave = 0; idxWave < nrWaves; ++idxWave) {
			if(rpwa::resonanceFit::debug()) {
				printDebug << "reading of entry " << idxWave << " in 'waves'." << std::endl;
			}

			const YAML::Node& configInputWave = configInputWaves[idxWave];

			std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
			boost::assign::insert(mandatoryArguments)
			                     ("name", rpwa::YamlCppUtils::TypeString);
			if(not checkIfAllVariablesAreThere(configInputWave, mandatoryArguments)) {
				printErr << "'waves' entry at index " << idxWave << " does not contain all required variables." << std::endl;
				throw;
			}

			const std::string waveName = configInputWave["name"].as<std::string>();

			// check that wave does not yet exist
			for(size_t idxCheckWave = 0; idxCheckWave < waves.size(); ++idxCheckWave) {
				if(waves[idxCheckWave].waveName() == waveName) {
					printErr << "wave '" << waveName << "' defined twice." << std::endl;
					throw;
				}
				if(std::find(waves[idxCheckWave].waveNameAlternatives().begin(), waves[idxCheckWave].waveNameAlternatives().end(), waveName) != waves[idxCheckWave].waveNameAlternatives().end()) {
					printErr << "wave '" << waveName << "' already defined as alternative name of wave '" << waves[idxCheckWave].waveName() << "'." << std::endl;
					throw;
				}
			}

			double massLower = -1.;
			if(configInputWave["massLower"]) {
				if(checkVariableType(configInputWave["massLower"], rpwa::YamlCppUtils::TypeFloat)) {
					massLower = configInputWave["massLower"].as<double>();
				} else {
					printErr << "variable 'massLower' of 'waves' entry at index " << idxWave << " is not a floating point number." << std::endl;
					throw;
				}
			}
			double massUpper = -1.;
			if(configInputWave["massUpper"]) {
				if(checkVariableType(configInputWave["massUpper"], rpwa::YamlCppUtils::TypeFloat)) {
					massUpper = configInputWave["massUpper"].as<double>();
				} else {
					printErr << "variable 'massUpper' of 'waves' entry at index " << idxWave << " is not a floating point number." << std::endl;
					throw;
				}
			}

			if(rpwa::resonanceFit::debug()) {
				printDebug << "read wave name: '" << waveName << "'." << std::endl;
				printDebug << "read mass range: '" << massLower << "' to '" << massUpper << "'." << std::endl;
			}

			std::vector<std::string> waveNameAlternatives;
			if(configInputWave["alternativeNames"]) {
				if(checkVariableType(configInputWave["alternativeNames"], rpwa::YamlCppUtils::TypeSequence)) {
					for(size_t idxAlt = 0; idxAlt < configInputWave["alternativeNames"].size(); ++idxAlt) {
						if(not checkVariableType(configInputWave["alternativeNames"][idxAlt], rpwa::YamlCppUtils::TypeString)) {
							printErr << "element " << idxAlt << " of variable 'alternativeNames' of 'waves' entry at index " << idxWave << " is not a string." << std::endl;
							throw;
						}
						const std::string waveNameAlternative = configInputWave["alternativeNames"][idxAlt].as<std::string>();

						// check that the alternative name does not yet exist
						if(waveNameAlternative == waveName) {
							printErr << "alternative name '" << waveNameAlternative << "' is equal to name of wave '" << waveName << "'." << std::endl;
							throw;
						}
						if(std::find(waveNameAlternatives.begin(), waveNameAlternatives.end(), waveNameAlternative) != waveNameAlternatives.end()) {
							printErr << "alternative name '" << waveNameAlternative << "' of wave '" << waveName << "' defined twice." << std::endl;
							throw;
						}
						for(size_t idxCheckWave = 0; idxCheckWave < waves.size(); ++idxCheckWave) {
							if(waves[idxCheckWave].waveName() == waveNameAlternative) {
								printErr << "alternative name '" << waveNameAlternative << "' of wave '" << waveName << "' already defined as separate wave." << std::endl;
								throw;
							}
							if(std::find(waves[idxCheckWave].waveNameAlternatives().begin(), waves[idxCheckWave].waveNameAlternatives().end(), waveNameAlternative) != waves[idxCheckWave].waveNameAlternatives().end()) {
								printErr << "alternative name '" << waveNameAlternative << "' of wave '" << waveName << "' already defined as alternative name of wave '" << waves[idxCheckWave].waveName() << "'." << std::endl;
								throw;
							}
						}

						waveNameAlternatives.push_back(waveNameAlternative);
					}
				} else {
					printErr << "variable 'alternativeNames' of 'waves' entry at index " << idxWave << " is not a sequence." << std::endl;
					throw;
				}
			}

			waves.push_back(rpwa::resonanceFit::information::wave(waveName,
			                                                      std::make_pair(massLower, massUpper),
			                                                      waveNameAlternatives));

			if(rpwa::resonanceFit::debug()) {
				printDebug << waves.back() << std::endl;
			}
		}

		return waves;
	}


}


rpwa::resonanceFit::informationConstPtr
rpwa::resonanceFit::readInformation(const YAML::Node& configRoot)
{
	if(rpwa::resonanceFit::debug()) {
		printDebug << "reading 'input'." << std::endl;
	}

	const YAML::Node& configInput = configRoot["input"];

	if(not configInput) {
		printErr << "'input' does not exist in configuration file." << std::endl;
		throw;
	}

	const std::vector<rpwa::resonanceFit::information::bin> bins = readInformationFitResults(configInput);
	const std::vector<rpwa::resonanceFit::information::wave> waves = readInformationWaves(configInput);

	return std::make_shared<rpwa::resonanceFit::information>(bins, waves);
}


namespace {


	void
	readModelAnchors(const YAML::Node& configModel,
	                 std::string& anchorWaveName,
	                 std::string& anchorComponentName)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading 'anchorwave'." << std::endl;
		}

		const YAML::Node& configAnchors = configModel["anchorwave"];
		if(not configAnchors) {
			printErr << "'anchorwave' is not a valid YAML node." << std::endl;
			throw;
		}
		if(not configAnchors.IsMap()) {
			printErr << "'anchorwave' is not a YAML map." << std::endl;
			throw;
		}

		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("name", rpwa::YamlCppUtils::TypeString)
		                     ("resonance", rpwa::YamlCppUtils::TypeString);
		if(not checkIfAllVariablesAreThere(configAnchors, mandatoryArguments)) {
			printErr << "'anchorwave' does not contain all required variables." << std::endl;
			throw;
		}

		anchorWaveName = configAnchors["name"].as<std::string>();
		anchorComponentName = configAnchors["resonance"].as<std::string>();
	}


	void
	readModelParameter(const YAML::Node& configComponent,
	                   const std::string& parameterName,
	                   rpwa::resonanceFit::parameter& parameter)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading parameter '" << parameterName << "'." << std::endl;
		}

		const YAML::Node& configParameter = configComponent[parameterName];
		if(not configParameter) {
			printErr << "final-state mass-dependence does not define parameter '" << parameterName << "'." << std::endl;
			throw;
		}

		parameter.setName(parameterName);

		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("val", rpwa::YamlCppUtils::TypeFloat)
		                     ("fix", rpwa::YamlCppUtils::TypeBoolean);
		if(not checkIfAllVariablesAreThere(configParameter, mandatoryArguments)) {
			printErr << "'" << parameterName << "' of final-state mass-dependence does not contain all required variables." << std::endl;
			throw;
		}

		const double startValue = configParameter["val"].as<double>();
		parameter.setStartValue(startValue);

		if(configParameter["error"]) {
			if(checkVariableType(configParameter["error"], rpwa::YamlCppUtils::TypeFloat)) {
				parameter.setStartError(configParameter["error"].as<double>());
			} else if(checkVariableType(configParameter["error"], rpwa::YamlCppUtils::TypeString) and configParameter["error"].as<std::string>() == "nan") {
				// some systems cannot convert the string 'nan'
				// to a floating point number, so this is to be
				// done manually. in any case this does not
				// really matter as the error is usually not
				// used.
				parameter.setStartError(std::numeric_limits<double>::has_quiet_NaN ? std::numeric_limits<double>::quiet_NaN() : 0.0);
			} else {
				printErr << "variable 'error' for parameter '" << parameterName << "' of final-state mass-dependence defined, but not a floating point number." << std::endl;
				throw;
			}
		}

		parameter.setFixed(configParameter["fix"].as<bool>());

		if(configParameter["lower"]) {
			if(checkVariableType(configParameter["lower"], rpwa::YamlCppUtils::TypeFloat)) {
				parameter.setLimitedLower(true);
				parameter.setLimitLower(configParameter["lower"].as<double>());
			} else {
				printErr << "variable 'lower' for parameter '" << parameterName << "' of final-state mass-dependence defined, but not a floating point number." << std::endl;
				throw;
			}
		} else {
			parameter.setLimitedLower(false);
		}
		if(configParameter["upper"]) {
			if(checkVariableType(configParameter["upper"], rpwa::YamlCppUtils::TypeFloat)) {
				parameter.setLimitedUpper(true);
				parameter.setLimitUpper(configParameter["upper"].as<double>());
			} else {
				printErr << "variable 'upper' for parameter '" << parameterName << "' of final-state mass-dependence defined, but not a floating point number." << std::endl;
				throw;
			}
		} else {
			parameter.setLimitedUpper(false);
		}

		if(configParameter["step"]) {
			if(checkVariableType(configParameter["step"], rpwa::YamlCppUtils::TypeFloat)) {
				parameter.setStep(configParameter["step"].as<double>());
			} else {
				printErr << "variable 'step' for parameter '" << parameterName << "' of final-state mass-dependence defined, but not a floating point number." << std::endl;
				throw;
			}
		}
	}


	double
	readModelComponentDecayChannelBranchingRatio(const YAML::Node& configDecayChannel)
	{
		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("branchingRatio", rpwa::YamlCppUtils::TypeFloat);
		if(not checkIfAllVariablesAreThere(configDecayChannel, mandatoryArguments)) {
			printErr << "decay channel does not contain all required variables." << std::endl;
			throw;
		}

		return configDecayChannel["branchingRatio"].as<double>();
	}


	double
	readModelComponentDecayChannelExponent(const YAML::Node& configDecayChannel)
	{
		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("exponent", rpwa::YamlCppUtils::TypeFloat);
		if(not checkIfAllVariablesAreThere(configDecayChannel, mandatoryArguments)) {
			printErr << "decay channel does not contain all required variables." << std::endl;
			throw;
		}

		return configDecayChannel["exponent"].as<double>();
	}


	double
	readModelComponentDecayChannelExponent(const YAML::Node& configDecayChannel,
	                                       const double defaultExponent)
	{
		if(configDecayChannel["exponent"]) {
			if(checkVariableType(configDecayChannel["exponent"], rpwa::YamlCppUtils::TypeFloat)) {
				return configDecayChannel["exponent"].as<double>();
			} else {
				printErr << "variable 'exponent' defined, but not a floating-point number." << std::endl;
				throw;
			}
		}

		printInfo << "variable 'exponent' not defined, using default value " << defaultExponent << "." << std::endl;
		return defaultExponent;
	}


	void
	readModelComponentDecayChannelIntegral(const YAML::Node& configDecayChannel,
	                                       std::vector<double>& masses,
	                                       std::vector<double>& values)
	{
		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("integral", rpwa::YamlCppUtils::TypeSequence);
		if(not checkIfAllVariablesAreThere(configDecayChannel, mandatoryArguments)) {
			printErr << "decay channel does not contain all required variables." << std::endl;
			throw;
		}

		const YAML::Node& integrals = configDecayChannel["integral"];

		const size_t nrValues = integrals.size();
		if(nrValues < 2) {
			printErr << "phase-space integral has to contain at least two points." << std::endl;
			throw;
		}

		masses.clear();
		values.clear();

		for(size_t idx = 0; idx < nrValues; ++idx) {
			const YAML::Node& integral = integrals[idx];
			if(not integral.IsSequence()) {
				printErr << "phase-space integral has to consist of arrays of two floating point numbers." << std::endl;
				throw;
			}
			if(integral.size() != 2) {
				printErr << "phase-space integral has to consist of arrays of two floating point numbers." << std::endl;
				throw;
			}
			if(not checkVariableType(integral[0], rpwa::YamlCppUtils::TypeFloat)) {
				printErr << "phase-space integral has to consist of arrays of two floating point numbers." << std::endl;
				throw;
			}
			if(not checkVariableType(integral[1], rpwa::YamlCppUtils::TypeFloat)) {
				printErr << "phase-space integral has to consist of arrays of two floating point numbers." << std::endl;
				throw;
			}

			const double mass = integral[0].as<double>();
			const double val = integral[1].as<double>();

			if(masses.size() > 0 and masses.back() > mass) {
				printErr << "masses of phase-space integral have to be strictly ordered." << std::endl;
				throw;
			}

			masses.push_back(mass);
			values.push_back(val);
		}
	}


	double
	readModelComponentDecayChannelMIsobar1(const YAML::Node& configDecayChannel)
	{
		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("mIsobar1", rpwa::YamlCppUtils::TypeFloat);
		if(not checkIfAllVariablesAreThere(configDecayChannel, mandatoryArguments)) {
			printErr << "decay channel does not contain all required variables." << std::endl;
			throw;
		}

		return configDecayChannel["mIsobar1"].as<double>();
	}


	double
	readModelComponentDecayChannelMIsobar2(const YAML::Node& configDecayChannel)
	{
		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("mIsobar2", rpwa::YamlCppUtils::TypeFloat);
		if(not checkIfAllVariablesAreThere(configDecayChannel, mandatoryArguments)) {
			printErr << "decay channel does not contain all required variables." << std::endl;
			throw;
		}

		return configDecayChannel["mIsobar2"].as<double>();
	}


	int
	readModelComponentDecayChannelRelAngularMom(const YAML::Node& configDecayChannel)
	{
		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("relAngularMom", rpwa::YamlCppUtils::TypeInt);
		if(not checkIfAllVariablesAreThere(configDecayChannel, mandatoryArguments)) {
			printErr << "decay channel does not contain all required variables." << std::endl;
			throw;
		}

		return configDecayChannel["relAngularMom"].as<int>();
	}


	int
	readModelComponentDecayChannelRelAngularMom(const YAML::Node& configDecayChannel,
	                                            const int defaultRelAngularMom)
	{
		if(configDecayChannel["relAngularMom"]) {
			if(checkVariableType(configDecayChannel["relAngularMom"], rpwa::YamlCppUtils::TypeInt)) {
				return configDecayChannel["relAngularMom"].as<int>();
			} else {
				printErr << "variable 'relAngularMom' defined, but not an integer." << std::endl;
				throw;
			}
		}

		printInfo << "variable 'relAngularMom' not defined, using default value " << defaultRelAngularMom << "." << std::endl;
		return defaultRelAngularMom;
	}


	std::vector<rpwa::resonanceFit::component::channel>
	readModelComponentDecayChannels(const YAML::Node& configComponent,
	                                const std::string& componentName,
	                                const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                                const rpwa::resonanceFit::baseDataConstPtr& fitData)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading 'decaychannels'." << std::endl;
		}

		std::map<std::string, size_t> waveIndices;
		for(size_t idxWave = 0; idxWave < fitInformation->nrWaves(); ++idxWave) {
			const rpwa::resonanceFit::information::wave& wave = fitInformation->getWave(idxWave);

			waveIndices[wave.waveName()] = idxWave;
			for(size_t idxAlt = 0; idxAlt < wave.waveNameAlternatives().size(); ++idxAlt) {
				waveIndices[wave.waveNameAlternatives()[idxAlt]] = idxWave;
			}
		}

		std::map<std::string, std::vector<size_t> > waveBins;
		for(size_t idxBin = 0; idxBin < fitInformation->nrBins(); ++idxBin) {
			for(size_t idxWave = 0; idxWave < fitInformation->nrWaves(); ++idxWave) {
				waveBins[fitData->waveNames()[idxBin][idxWave]].push_back(idxBin);
			}
		}

		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("decaychannels", rpwa::YamlCppUtils::TypeSequence);
		if(not checkIfAllVariablesAreThere(configComponent, mandatoryArguments)) {
			printErr << "'components' entry of component '" << componentName << "' does not contain all required variables." << std::endl;
			throw;
		}

		const YAML::Node& configDecayChannels = configComponent["decaychannels"];

		std::vector<rpwa::resonanceFit::component::channel> decayChannels;
		const size_t nrDecayChannels = configDecayChannels.size();
		for(size_t idxDecayChannel = 0; idxDecayChannel < nrDecayChannels; ++idxDecayChannel) {
			const YAML::Node& configDecayChannel = configDecayChannels[idxDecayChannel];

			std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
			boost::assign::insert(mandatoryArguments)
			                     ("amp", rpwa::YamlCppUtils::TypeString);
			if(not checkIfAllVariablesAreThere(configDecayChannel, mandatoryArguments)) {
				printErr << "one of the decay channels of component '" << componentName << "' does not contain all required variables." << std::endl;
				throw;
			}

			const std::string waveName = configDecayChannel["amp"].as<std::string>();

			// check that a wave with this wave is not yet in the decay channels
			for(size_t idxChannel = 0; idxChannel < decayChannels.size(); ++idxChannel) {
				if(decayChannels[idxChannel].getWaveName() == waveName) {
					printErr << "wave '" << waveName << "' defined twice in the decay channels of component '" << componentName << "'." << std::endl;
					throw;
				}
			}

			// get index of wave in array for wave names, phase-space integrals, ...
			const std::map<std::string, size_t>::const_iterator it = waveIndices.find(waveName);
			if(it == waveIndices.end()) {
				printErr << "wave '" << waveName << "' not in fit, but used as decay channel in component '" << componentName << "'." << std::endl;
				throw;
			}
			const size_t waveIdx = it->second;

			// get list of bins this wave is defined in
			const std::map<std::string, std::vector<size_t> >::const_iterator it2 = waveBins.find(waveName);
			if(it2 == waveBins.end()) {
				printErr << "wave '" << waveName << "' not in fit, but used as decay channel in component '" << componentName << "'." << std::endl;
				throw;
			}
			const std::vector<size_t>& binsForWave = it2->second;

			// get a view for the current wave for all bins and all mass bins
			boost::multi_array<double, 3>::const_array_view<2>::type view = fitData->phaseSpaceIntegrals()[boost::indices[boost::multi_array<double, 3>::index_range()][boost::multi_array<double, 3>::index_range()][waveIdx]];
			decayChannels.push_back(rpwa::resonanceFit::component::channel(waveIdx,
			                                                               waveName,
			                                                               binsForWave,
			                                                               fitData->nrMassBins(),
			                                                               fitData->massBinCenters(),
			                                                               view));
		}

		return decayChannels;
	}


	template<typename T>
	rpwa::resonanceFit::componentPtr
	readModelComponent(const YAML::Node& /*configComponent*/,
	                   const rpwa::resonanceFit::informationConstPtr& /*fitInformation*/,
	                   const rpwa::resonanceFit::baseDataConstPtr& fitData,
	                   const size_t id,
	                   const std::string& name,
	                   const std::vector<rpwa::resonanceFit::parameter>& parameters,
	                   const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
	                   const bool useBranchings)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading component of non-specialized type." << std::endl;
		}

		return std::make_shared<T>(id,
		                           name,
		                           parameters,
		                           decayChannels,
		                           fitData->nrMassBins(),
		                           fitData->massBinCenters(),
		                           useBranchings);
	}


	template<>
	rpwa::resonanceFit::componentPtr
	readModelComponent<rpwa::resonanceFit::dynamicWidthBreitWigner>(const YAML::Node& configComponent,
	                                                                const rpwa::resonanceFit::informationConstPtr& /*fitInformation*/,
	                                                                const rpwa::resonanceFit::baseDataConstPtr& fitData,
	                                                                const size_t id,
	                                                                const std::string& name,
	                                                                const std::vector<rpwa::resonanceFit::parameter>& parameters,
	                                                                const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
	                                                                const bool useBranchings)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading component of type 'dynamicWidthBreitWigner'." << std::endl;
		}

		std::vector<double> branchingRatio;
		std::vector<int> relAngularMom;
		std::vector<double> mIsobar1;
		std::vector<double> mIsobar2;

		const YAML::Node& configDecayChannels = configComponent["decaychannels"];
		if(not configDecayChannels) {
			printErr << "a component of type 'dynamicWidthBreitWigner' has no decay channels." << std::endl;
			throw;
		}

		const size_t nrDecayChannels = configDecayChannels.size();
		for(size_t idxDecayChannel = 0; idxDecayChannel < nrDecayChannels; ++idxDecayChannel) {
			const YAML::Node& configDecayChannel = configDecayChannels[idxDecayChannel];

			branchingRatio.push_back(readModelComponentDecayChannelBranchingRatio(configDecayChannel));
			relAngularMom.push_back(readModelComponentDecayChannelRelAngularMom(configDecayChannel));
			mIsobar1.push_back(readModelComponentDecayChannelMIsobar1(configDecayChannel));
			mIsobar2.push_back(readModelComponentDecayChannelMIsobar2(configDecayChannel));
		}

		const YAML::Node& configExtraDecayChannels = configComponent["extradecaychannels"];
		const size_t nrExtraDecayChannels = configExtraDecayChannels ? configExtraDecayChannels.size() : 0;
		for(size_t idxExtraDecayChannel = 0; idxExtraDecayChannel < nrExtraDecayChannels; ++idxExtraDecayChannel) {
			const YAML::Node& configExtraDecayChannel = configExtraDecayChannels[idxExtraDecayChannel];

			branchingRatio.push_back(readModelComponentDecayChannelBranchingRatio(configExtraDecayChannel));
			relAngularMom.push_back(readModelComponentDecayChannelRelAngularMom(configExtraDecayChannel));
			mIsobar1.push_back(readModelComponentDecayChannelMIsobar1(configExtraDecayChannel));
			mIsobar2.push_back(readModelComponentDecayChannelMIsobar2(configExtraDecayChannel));
		}

		return std::make_shared<rpwa::resonanceFit::dynamicWidthBreitWigner>(id,
		                                                                     name,
		                                                                     parameters,
		                                                                     decayChannels,
		                                                                     fitData->nrMassBins(),
		                                                                     fitData->massBinCenters(),
		                                                                     useBranchings,
		                                                                     branchingRatio,
		                                                                     relAngularMom,
		                                                                     mIsobar1,
		                                                                     mIsobar2);
	}


	template<>
	rpwa::resonanceFit::componentPtr
	readModelComponent<rpwa::resonanceFit::integralWidthBreitWigner>(const YAML::Node& configComponent,
	                                                                 const rpwa::resonanceFit::informationConstPtr& /*fitInformation*/,
	                                                                 const rpwa::resonanceFit::baseDataConstPtr& fitData,
	                                                                 const size_t id,
	                                                                 const std::string& name,
	                                                                 const std::vector<rpwa::resonanceFit::parameter>& parameters,
	                                                                 const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
	                                                                 const bool useBranchings)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading component of type 'integralWidthBreitWigner'." << std::endl;
		}

		std::vector<double> branchingRatio;
		std::vector<std::vector<double> > masses;
		std::vector<std::vector<double> > values;

		const YAML::Node& configDecayChannels = configComponent["decaychannels"];
		if(not configDecayChannels) {
			printErr << "a component of type 'integralWidthBreitWigner' has no decay channels." << std::endl;
			throw;
		}

		const size_t nrDecayChannels = configDecayChannels.size();
		for(size_t idxDecayChannel = 0; idxDecayChannel < nrDecayChannels; ++idxDecayChannel) {
			const YAML::Node& configDecayChannel = configDecayChannels[idxDecayChannel];

			branchingRatio.push_back(readModelComponentDecayChannelBranchingRatio(configDecayChannel));

			std::vector<double> tempMasses;
			std::vector<double> tempValues;
			readModelComponentDecayChannelIntegral(configDecayChannel, tempMasses, tempValues);
			masses.push_back(tempMasses);
			values.push_back(tempValues);
		}

		const YAML::Node& configExtraDecayChannels = configComponent["extradecaychannels"];
		const size_t nrExtraDecayChannels = configExtraDecayChannels ? configExtraDecayChannels.size() : 0;
		for(size_t idxExtraDecayChannel = 0; idxExtraDecayChannel < nrExtraDecayChannels; ++idxExtraDecayChannel) {
			const YAML::Node& configExtraDecayChannel = configExtraDecayChannels[idxExtraDecayChannel];

			branchingRatio.push_back(readModelComponentDecayChannelBranchingRatio(configExtraDecayChannel));

			std::vector<double> tempMasses;
			std::vector<double> tempValues;
			readModelComponentDecayChannelIntegral(configExtraDecayChannel, tempMasses, tempValues);
			masses.push_back(tempMasses);
			values.push_back(tempValues);
		}

		return std::make_shared<rpwa::resonanceFit::integralWidthBreitWigner>(id,
		                                                                      name,
		                                                                      parameters,
		                                                                      decayChannels,
		                                                                      fitData->nrMassBins(),
		                                                                      fitData->massBinCenters(),
		                                                                      useBranchings,
		                                                                      branchingRatio,
		                                                                      masses,
		                                                                      values);
	}


	template<>
	rpwa::resonanceFit::componentPtr
	readModelComponent<rpwa::resonanceFit::exponentialBackground>(const YAML::Node& configComponent,
	                                                              const rpwa::resonanceFit::informationConstPtr& /*fitInformation*/,
	                                                              const rpwa::resonanceFit::baseDataConstPtr& fitData,
	                                                              const size_t id,
	                                                              const std::string& name,
	                                                              const std::vector<rpwa::resonanceFit::parameter>& parameters,
	                                                              const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
	                                                              const bool useBranchings)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading component of type 'exponentialBackground'." << std::endl;
		}

		const int relAngularMom = readModelComponentDecayChannelRelAngularMom(configComponent, 0);
		const double mIsobar1 = readModelComponentDecayChannelMIsobar1(configComponent);
		const double mIsobar2 = readModelComponentDecayChannelMIsobar2(configComponent);
		const double exponent = readModelComponentDecayChannelExponent(configComponent, 2.0);

		return std::make_shared<rpwa::resonanceFit::exponentialBackground>(id,
		                                                                   name,
		                                                                   parameters,
		                                                                   decayChannels,
		                                                                   fitData->nrMassBins(),
		                                                                   fitData->massBinCenters(),
		                                                                   useBranchings,
		                                                                   relAngularMom,
		                                                                   mIsobar1,
		                                                                   mIsobar2,
		                                                                   exponent);
	}


	template<>
	rpwa::resonanceFit::componentPtr
	readModelComponent<rpwa::resonanceFit::tPrimeDependentBackground>(const YAML::Node& configComponent,
	                                                                  const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                                                                  const rpwa::resonanceFit::baseDataConstPtr& fitData,
	                                                                  const size_t id,
	                                                                  const std::string& name,
	                                                                  const std::vector<rpwa::resonanceFit::parameter>& parameters,
	                                                                  const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
	                                                                  const bool useBranchings)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading component of type 'tPrimeDependentBackground'." << std::endl;
		}

		std::vector<double> tPrimeMeans;
		for(size_t idxBin = 0; idxBin < fitInformation->nrBins(); ++idxBin) {
			tPrimeMeans.push_back(fitInformation->getBin(idxBin).tPrimeMean());
		}

		const int relAngularMom = readModelComponentDecayChannelRelAngularMom(configComponent, 0);
		const double mIsobar1 = readModelComponentDecayChannelMIsobar1(configComponent);
		const double mIsobar2 = readModelComponentDecayChannelMIsobar2(configComponent);
		const double exponent = readModelComponentDecayChannelExponent(configComponent, 2.0);

		return std::make_shared<rpwa::resonanceFit::tPrimeDependentBackground>(id,
		                                                                       name,
		                                                                       parameters,
		                                                                       decayChannels,
		                                                                       fitData->nrMassBins(),
		                                                                       fitData->massBinCenters(),
		                                                                       useBranchings,
		                                                                       tPrimeMeans,
		                                                                       relAngularMom,
		                                                                       mIsobar1,
		                                                                       mIsobar2,
		                                                                       exponent);
	}


	template<>
	rpwa::resonanceFit::componentPtr
	readModelComponent<rpwa::resonanceFit::exponentialBackgroundIntegral>(const YAML::Node& configComponent,
	                                                                      const rpwa::resonanceFit::informationConstPtr& /*fitInformation*/,
	                                                                      const rpwa::resonanceFit::baseDataConstPtr& fitData,
	                                                                      const size_t id,
	                                                                      const std::string& name,
	                                                                      const std::vector<rpwa::resonanceFit::parameter>& parameters,
	                                                                      const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
	                                                                      const bool useBranchings)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading component of type 'exponentialBackgroundIntegral'." << std::endl;
		}

		std::vector<double> masses;
		std::vector<double> values;
		readModelComponentDecayChannelIntegral(configComponent, masses, values);
		const double exponent = readModelComponentDecayChannelExponent(configComponent);

		return std::make_shared<rpwa::resonanceFit::exponentialBackgroundIntegral>(id,
		                                                                           name,
		                                                                           parameters,
		                                                                           decayChannels,
		                                                                           fitData->nrMassBins(),
		                                                                           fitData->massBinCenters(),
		                                                                           useBranchings,
		                                                                           masses,
		                                                                           values,
		                                                                           exponent);
	}


	template<>
	rpwa::resonanceFit::componentPtr
	readModelComponent<rpwa::resonanceFit::tPrimeDependentBackgroundIntegral>(const YAML::Node& configComponent,
	                                                                          const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                                                                          const rpwa::resonanceFit::baseDataConstPtr& fitData,
	                                                                          const size_t id,
	                                                                          const std::string& name,
	                                                                          const std::vector<rpwa::resonanceFit::parameter>& parameters,
	                                                                          const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
	                                                                          const bool useBranchings)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading component of type 'tPrimeDependentBackgroundIntegral'." << std::endl;
		}

		std::vector<double> tPrimeMeans;
		for(size_t idxBin = 0; idxBin < fitInformation->nrBins(); ++idxBin) {
			tPrimeMeans.push_back(fitInformation->getBin(idxBin).tPrimeMean());
		}

		std::vector<double> masses;
		std::vector<double> values;
		readModelComponentDecayChannelIntegral(configComponent, masses, values);
		const double exponent = readModelComponentDecayChannelExponent(configComponent);

		return std::make_shared<rpwa::resonanceFit::tPrimeDependentBackgroundIntegral>(id,
		                                                                               name,
		                                                                               parameters,
		                                                                               decayChannels,
		                                                                               fitData->nrMassBins(),
		                                                                               fitData->massBinCenters(),
		                                                                               useBranchings,
		                                                                               tPrimeMeans,
		                                                                               masses,
		                                                                               values,
		                                                                               exponent);
	}


	template<typename T>
	rpwa::resonanceFit::componentPtr
	readModelComponent(const YAML::Node& configComponent,
	                   const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                   const rpwa::resonanceFit::baseDataConstPtr& fitData,
	                   const size_t id,
	                   const bool useBranchings)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading component of type 'component'." << std::endl;
		}

		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("name", rpwa::YamlCppUtils::TypeString);
		if(not checkIfAllVariablesAreThere(configComponent, mandatoryArguments)) {
			printErr << "'components' entry does not contain all required variables." << std::endl;
			throw;
		}

		const std::string name = configComponent["name"].as<std::string>();

		std::vector<rpwa::resonanceFit::parameter> parameters = T::getDefaultParameters();
		for(size_t idxParameter = 0; idxParameter < parameters.size(); ++idxParameter) {
			const std::string parameterName = parameters[idxParameter].name();
			readModelParameter(configComponent,
			                   parameterName,
			                   parameters[idxParameter]);
		}

		const std::vector<rpwa::resonanceFit::component::channel> decayChannels = readModelComponentDecayChannels(configComponent,
		                                                                                                          name,
		                                                                                                          fitInformation,
		                                                                                                          fitData);

		return readModelComponent<T>(configComponent,
		                             fitInformation,
		                             fitData,
		                             id,
		                             name,
		                             parameters,
		                             decayChannels,
		                             useBranchings);
	}


	rpwa::resonanceFit::componentPtr
	readModelComponent(const YAML::Node& configComponent,
	                   const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                   const rpwa::resonanceFit::baseDataConstPtr& fitData,
	                   const size_t id,
	                   const bool useBranchings)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading component." << std::endl;
		}

		std::string type = "fixedWidthBreitWigner";
		if(configComponent["type"]) {
			if(not checkVariableType(configComponent["type"], rpwa::YamlCppUtils::TypeString)) {
				printErr << "a component has a type that is not a string." << std::endl;
				throw;
			}
			type = configComponent["type"].as<std::string>();
		}

		if(rpwa::resonanceFit::debug()) {
			printDebug << "found component of type '" << type << "'." << std::endl;
		}

		rpwa::resonanceFit::componentPtr component;
		if(type == "fixedWidthBreitWigner") {
			component = readModelComponent<rpwa::resonanceFit::fixedWidthBreitWigner>(configComponent,
			                                                                          fitInformation,
			                                                                          fitData,
			                                                                          id,
			                                                                          useBranchings);
		} else if(type == "dynamicWidthBreitWigner") {
			component = readModelComponent<rpwa::resonanceFit::dynamicWidthBreitWigner>(configComponent,
			                                                                            fitInformation,
			                                                                            fitData,
			                                                                            id,
			                                                                            useBranchings);
		} else if(type == "integralWidthBreitWigner") {
			component = readModelComponent<rpwa::resonanceFit::integralWidthBreitWigner>(configComponent,
			                                                                             fitInformation,
			                                                                             fitData,
			                                                                             id,
			                                                                             useBranchings);
		} else if(type == "constantBackground") {
			component = readModelComponent<rpwa::resonanceFit::constantBackground>(configComponent,
			                                                                       fitInformation,
			                                                                       fitData,
			                                                                       id,
			                                                                       useBranchings);
		} else if(type == "exponentialBackground") {
			component = readModelComponent<rpwa::resonanceFit::exponentialBackground>(configComponent,
			                                                                          fitInformation,
			                                                                          fitData,
			                                                                          id,
			                                                                          useBranchings);
		} else if(type == "tPrimeDependentBackground") {
			component = readModelComponent<rpwa::resonanceFit::tPrimeDependentBackground>(configComponent,
			                                                                              fitInformation,
			                                                                              fitData,
			                                                                              id,
			                                                                              useBranchings);
		} else if(type == "exponentialBackgroundIntegral") {
			component = readModelComponent<rpwa::resonanceFit::exponentialBackgroundIntegral>(configComponent,
			                                                                                  fitInformation,
			                                                                                  fitData,
			                                                                                  id,
			                                                                                  useBranchings);
		} else if(type == "tPrimeDependentBackgroundIntegral") {
			component = readModelComponent<rpwa::resonanceFit::tPrimeDependentBackgroundIntegral>(configComponent,
			                                                                                      fitInformation,
			                                                                                      fitData,
			                                                                                      id,
			                                                                                      useBranchings);
		} else {
			printErr << "unknown type '" << type << "'." << std::endl;
			throw;
		}

		if(rpwa::resonanceFit::debug()) {
			component->print(printDebug);
		}

		return component;
	}


	rpwa::resonanceFit::componentPtr
	readModelComponent(const YAML::Node& configComponent,
	                   const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                   const rpwa::resonanceFit::baseDataConstPtr& fitData,
	                   rpwa::resonanceFit::parameters& fitParameters,
	                   rpwa::resonanceFit::parameters& fitParametersError,
	                   const size_t id,
	                   const bool useBranchings)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading component and its parameters." << std::endl;
		}

		const rpwa::resonanceFit::componentPtr& component = readModelComponent(configComponent,
		                                                                       fitInformation,
		                                                                       fitData,
		                                                                       id,
		                                                                       useBranchings);

		fitParameters.resize(id+1, component->getNrChannels(), component->getNrParameters(), fitInformation->nrBins());
		fitParametersError.resize(id+1, component->getNrChannels(), component->getNrParameters(), fitInformation->nrBins());

		for(size_t idxParameter = 0; idxParameter < component->getNrParameters(); ++idxParameter) {
			fitParameters.setParameter(id, idxParameter, component->getParameter(idxParameter).startValue());
			fitParametersError.setParameter(id, idxParameter, component->getParameter(idxParameter).startError());
		}

		for(size_t idxDecayChannel = 0; idxDecayChannel < component->getNrChannels(); ++idxDecayChannel) {
			const rpwa::resonanceFit::component::channel& channel = component->getChannel(idxDecayChannel);
			const YAML::Node& configDecayChannel = configComponent["decaychannels"][idxDecayChannel];

			if(component->mapCouplingToMasterChannel(component->mapChannelToCoupling(idxDecayChannel)) == idxDecayChannel) {
				const YAML::Node& configCouplings = configDecayChannel["couplings"];
				if(not configCouplings) {
					printErr << "decay channel '" << channel.getWaveName() << "' of component '" << component->getName() << "' has no couplings." << std::endl;
					throw;
				}
				if(not configCouplings.IsSequence()) {
					printErr << "decay channel '" << channel.getWaveName() << "' of component '" << component->getName() << "' has no couplings." << std::endl;
					throw;
				}

				const size_t nrCouplings = configCouplings.size();
				if(nrCouplings != channel.getBins().size()) {
					printErr << "decay channel '" << channel.getWaveName() << "' of component '" << component->getName() << "' has " << nrCouplings << " couplings, not " << channel.getBins().size() << "." << std::endl;
					throw;
				}

				for(size_t idxCoupling = 0; idxCoupling < nrCouplings; ++idxCoupling) {
					const YAML::Node& configCoupling = configCouplings[idxCoupling];
					if(not checkVariableType(configCoupling, rpwa::YamlCppUtils::TypeSequence)) {
						printErr << "one of the couplings of the decay channel '" << channel.getWaveName() << "' of the component '" << component->getName() << "' is not a YAML sequence." << std::endl;
						throw;
					}
					if(configCoupling.size() != 2) {
						printErr << "one of the couplings of the decay channel '" << channel.getWaveName() << "' of the component '" << component->getName() << "' does not contain exactly two entries." << std::endl;
						throw;
					}

					if(not checkVariableType(configCoupling[0], rpwa::YamlCppUtils::TypeFloat)) {
						printErr << "real part of one of the couplings of the decay channel '" << channel.getWaveName() << "' of the component '" << component->getName() << "' is not a floating point number." << std::endl;
						throw;
					}
					if(not checkVariableType(configCoupling[1], rpwa::YamlCppUtils::TypeFloat)) {
						printErr << "imaginary part of one of the couplings of the decay channel '" << channel.getWaveName() << "' of the component '" << component->getName() << "' is not a floating point number." << std::endl;
						throw;
					}

					const double couplingReal = configCoupling[0].as<double>();
					const double couplingImag = configCoupling[1].as<double>();

					const std::complex<double> coupling(couplingReal, couplingImag);
					fitParameters.setCoupling(id, component->mapChannelToCoupling(idxDecayChannel), channel.getBins()[idxCoupling], coupling);
				}
			}

			if(component->getNrBranchings() > 1) {
				if(component->mapBranchingToMasterChannel(component->mapChannelToBranching(idxDecayChannel)) == idxDecayChannel) {
					const YAML::Node& configBranching = configDecayChannel["branching"];
					if(not configBranching) {
						printErr << "decay channel '" << channel.getWaveName() << "' of component '" << component->getName() << "' has no branching." << std::endl;
						throw;
					}
					if(not configBranching.IsSequence()) {
						printErr << "decay channel '" << channel.getWaveName() << "' of component '" << component->getName() << "' has no branching." << std::endl;
						throw;
					}

					if(configBranching.size() != 2) {
						printErr << "branching of the decay channel '" << channel.getWaveName() << "' of the component '" << component->getName() << "' does not contain exactly two entries." << std::endl;
						throw;
					}

					if(not checkVariableType(configBranching[0], rpwa::YamlCppUtils::TypeFloat)) {
						printErr << "real part of branching of the decay channel '" << channel.getWaveName() << "' of the component '" << component->getName() << "' is not a floating point number." << std::endl;
						throw;
					}
					if(not checkVariableType(configBranching[1], rpwa::YamlCppUtils::TypeFloat)) {
						printErr << "imaginary part of branching of the decay channel '" << channel.getWaveName() << "' of the component '" << component->getName() << "' is not a floating point number." << std::endl;
						throw;
					}

					double branchingReal = configBranching[0].as<double>();
					double branchingImag = configBranching[1].as<double>();

					if(component->isBranchingFixed(component->mapChannelToBranching(idxDecayChannel))) {
						// the first branching should always be 1.
						if(branchingReal != 1.0 or branchingImag != 0.0) {
							printWarn << "branching of the decay channel '" << channel.getWaveName() << "' of the component '" << component->getName() << "' forced to 1." << std::endl;
							branchingReal = 1.0;
							branchingImag = 0.0;
						}
					}

					const std::complex<double> branching(branchingReal, branchingImag);
					fitParameters.setBranching(id, component->mapChannelToBranching(idxDecayChannel), branching);
				}
			} else {
				const std::complex<double> branching(1.0, 0.0);
				fitParameters.setBranching(id, 0, branching);
			}
		}

		return component;
	}


	std::vector<rpwa::resonanceFit::componentPtr>
	readModelComponents(const YAML::Node& configModel,
	                    const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                    const rpwa::resonanceFit::baseDataConstPtr& fitData,
	                    rpwa::resonanceFit::parameters& fitParameters,
	                    rpwa::resonanceFit::parameters& fitParametersError,
	                    const bool useBranchings)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading 'components'." << std::endl;
		}

		const YAML::Node& configComponents = configModel["components"];
		if(not configComponents) {
			printErr << "'components' is not a valid YAML node." << std::endl;
			throw;
		}
		if(not configComponents.IsSequence()) {
			printErr << "'components' is not a YAML sequence." << std::endl;
			throw;
		}

		const size_t nrComponents = configComponents.size();
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading " << nrComponents << " components from configuration file." << std::endl;
		}

		std::vector<rpwa::resonanceFit::componentPtr> components;
		for(size_t idxComponent = 0; idxComponent < nrComponents; ++idxComponent) {
			const YAML::Node& configComponent = configComponents[idxComponent];

			const rpwa::resonanceFit::componentPtr& component = readModelComponent(configComponent,
			                                                                       fitInformation,
			                                                                       fitData,
			                                                                       fitParameters,
			                                                                       fitParametersError,
			                                                                       components.size(),
			                                                                       useBranchings);

			for(size_t idx = 0; idx < components.size(); ++idx) {
				if(components[idx]->getName() == component->getName()) {
					printErr << "component '" << component->getName() << "' defined twice." << std::endl;
					throw;
				}
			}

			components.push_back(component);
		}

		return components;
	}


	void
	readModelFsmdBin(const YAML::Node& configFsmd,
	                 std::shared_ptr<TFormula>& function,
	                 boost::multi_array<rpwa::resonanceFit::parameter, 1>& parameters)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading 'finalStateMassDependence' for an individual bin." << std::endl;
		}

		if(not configFsmd.IsMap()) {
			printErr << "'finalStateMassDependence' for an individual bin is not a YAML map." << std::endl;
			throw;
		}

		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("formula", rpwa::YamlCppUtils::TypeString);
		if(not checkIfAllVariablesAreThere(configFsmd, mandatoryArguments)) {
			printErr << "'finalStateMassDependence' for an individual bin does not contain all required variables." << std::endl;
			return throw;
		}

		const std::string formula = configFsmd["formula"].as<std::string>();
		function.reset(new TFormula("finalStateMassDependence", formula.c_str()));

		const size_t nrParameters = function->GetNpar();
		parameters.resize(boost::extents[nrParameters]);

		for(size_t idxParameter = 0; idxParameter < nrParameters; ++idxParameter) {
			const std::string parameterName = function->GetParName(idxParameter);

			// set default value for step size
			parameters[idxParameter].setStep(0.0001);

			readModelParameter(configFsmd,
			                   parameterName,
			                   parameters[idxParameter]);
		}
	}


	rpwa::resonanceFit::fsmdPtr
	readModelFsmd(const YAML::Node& configModel,
	              const rpwa::resonanceFit::baseDataConstPtr& fitData,
	              rpwa::resonanceFit::parameters& fitParameters,
	              rpwa::resonanceFit::parameters& fitParametersError,
	              const size_t id)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading 'finalStateMassDependence'." << std::endl;
		}

		const YAML::Node& configFsmd = configModel["finalStateMassDependence"];
		if(not configFsmd) {
			// final-state mass-dependence might not be specified
			return rpwa::resonanceFit::fsmdPtr();
		}

		if(not configFsmd.IsMap() and not configFsmd.IsSequence()) {
			printErr << "'finalStateMassDependence' is not a YAML map or sequence." << std::endl;
			throw;
		}

		rpwa::resonanceFit::fsmdPtr fsmd;
		if(configFsmd.IsMap()) {
			// a single final-state mass-dependence is given
			std::shared_ptr<TFormula> function;
			boost::multi_array<rpwa::resonanceFit::parameter, 1> parameters;

			readModelFsmdBin(configFsmd,
			                 function,
			                 parameters);

			const size_t nrParameters = function->GetNpar();
			fitParameters.resize(id+1, 0, nrParameters, 0);
			fitParametersError.resize(id+1, 0, nrParameters, 0);

			for(size_t idxParameter = 0; idxParameter < nrParameters; ++idxParameter) {
				fitParameters.setParameter(id, idxParameter, parameters[idxParameter].startValue());
				fitParametersError.setParameter(id, idxParameter, parameters[idxParameter].startError());
			}

			fsmd.reset(new rpwa::resonanceFit::fsmd(id,
			                                        fitData->nrMassBins(),
			                                        fitData->massBinCenters(),
			                                        function,
			                                        parameters));
		} else {
			// a final-state mass-dependence for each bin is given
			std::vector<std::shared_ptr<TFormula> > functions;
			boost::multi_array<rpwa::resonanceFit::parameter, 2> parameters;

			size_t nrParametersSum = 0;
			const size_t nrBins = configFsmd.size();
			for(size_t idxBin = 0; idxBin < nrBins; ++idxBin) {
				std::shared_ptr<TFormula> tempFunction;
				boost::multi_array<rpwa::resonanceFit::parameter, 1> tempParameters;

				readModelFsmdBin(configFsmd[idxBin],
				                 tempFunction,
				                 tempParameters);

				rpwa::resonanceFit::adjustSizeAndSet(functions, idxBin, tempFunction);
				rpwa::resonanceFit::adjustSizeAndSet(parameters, idxBin, tempParameters);

				const size_t idxParameterFirst = nrParametersSum;
				const size_t nrParameters = tempFunction->GetNpar();
				nrParametersSum += nrParameters;

				fitParameters.resize(id+1, 0, nrParametersSum, nrBins);
				fitParametersError.resize(id+1, 0, nrParametersSum, nrBins);

				for(size_t idxParameter = 0; idxParameter < nrParameters; ++idxParameter) {
					fitParameters.setParameter(id, idxParameterFirst+idxParameter, tempParameters[idxParameter].startValue());
					fitParametersError.setParameter(id, idxParameterFirst+idxParameter, tempParameters[idxParameter].startError());
				}
			}

			fsmd.reset(new rpwa::resonanceFit::fsmd(id,
			                                        fitData->nrMassBins(),
			                                        fitData->massBinCenters(),
			                                        functions,
			                                        parameters));
		}

		if(rpwa::resonanceFit::debug()) {
			printDebug << *fsmd << std::endl;
		}

		return fsmd;
	}


}


rpwa::resonanceFit::modelConstPtr
rpwa::resonanceFit::readModel(const YAML::Node& configRoot,
                              const rpwa::resonanceFit::informationConstPtr& fitInformation,
                              const rpwa::resonanceFit::baseDataConstPtr& fitData,
                              rpwa::resonanceFit::parameters& fitParameters,
                              rpwa::resonanceFit::parameters& fitParametersError,
                              const bool useBranchings)
{
	if(rpwa::resonanceFit::debug()) {
		printDebug << "reading 'model'." << std::endl;
	}

	const YAML::Node& configModel = configRoot["model"];

	if(not configModel) {
		printErr << "'model' does not exist in configuration file." << std::endl;
		throw;
	}

	std::string anchorWaveName;
	std::string anchorComponentName;
	readModelAnchors(configModel,
	                 anchorWaveName,
	                 anchorComponentName);

	const std::vector<rpwa::resonanceFit::componentPtr>& components = readModelComponents(configModel,
	                                                                                      fitInformation,
	                                                                                      fitData,
	                                                                                      fitParameters,
	                                                                                      fitParametersError,
	                                                                                      useBranchings);

	// get information for creating the final-state mass-dependence
	const rpwa::resonanceFit::fsmdPtr& fsmd = readModelFsmd(configModel,
	                                                        fitData,
	                                                        fitParameters,
	                                                        fitParametersError,
	                                                        components.size());

	return std::make_shared<rpwa::resonanceFit::model>(fitInformation,
	                                                   components,
	                                                   fsmd,
	                                                   anchorWaveName,
	                                                   anchorComponentName);
}
