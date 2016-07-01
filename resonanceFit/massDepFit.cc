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


#include "massDepFit.h"


#include <algorithm>
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>
#include <string>

#include <boost/assign/std/vector.hpp>
#include <boost/tokenizer.hpp>

#include <yaml-cpp/yaml.h>

#include <TTree.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TString.h>
#include <TComplex.h>
#include <TRandom3.h>
#include <TStopwatch.h>

#include "fileUtils.hpp"
#include "fitResult.h"
#include "massDepFitComponents.h"
#include "massDepFitFsmd.h"
#include "massDepFitFunction.h"
#include "massDepFitModel.h"
#include "reportingUtils.hpp"
#include "yamlCppUtils.hpp"


bool rpwa::massDepFit::massDepFit::_debug = false;


rpwa::massDepFit::massDepFit::massDepFit()
	: _sysPlotting(false),
	  _sameMassBinning(true),
	  _maxMassBins(0),
	  _nrSystematics(0),
	  _nrWaves(0)
{
}


bool
rpwa::massDepFit::massDepFit::readConfig(const YAML::Node& configRoot,
                                         rpwa::massDepFit::model& fitModel,
                                         rpwa::massDepFit::parameters& fitParameters,
                                         rpwa::massDepFit::parameters& fitParametersError,
                                         int& minStatus,
                                         double& chi2,
                                         unsigned int& ndf,
                                         const std::string& valTreeName,
                                         const std::string& valBranchName)
{
	// fit result information
	const YAML::Node& configFitquality = configRoot["fitquality"];
	if(configFitquality) {
		if(not readConfigFitquality(configFitquality, minStatus, chi2, ndf)) {
			printErr << "error while reading 'fitquality' in configuration file." << std::endl;
			return false;
		}
	} else {
		minStatus = 0;
		chi2 = 0.;
		ndf = 0;
	}

	// input section
	const YAML::Node& configInput = configRoot["input"];
	if(not configInput) {
		printErr << "'input' does not exist in configuration file." << std::endl;
		return false;
	}
	if(not readConfigInput(configInput)) {
		printErr << "error while reading 'input' in configuration file." << std::endl;
		return false;
	}

	// extract information from fit results
	if(not readInFiles(valTreeName, valBranchName)) {
		printErr << "error while reading fit result." << std::endl;
		return false;
	}

	// extract information for systematic errors
	if(not readSystematicsFiles(valTreeName, valBranchName)) {
		printErr << "error while reading fit results for systematic errors." << std::endl;
		return false;
	}

	// prepare mass limits
	if(not prepareMassLimits()) {
		printErr << "error while determining which bins to use in the fit." << std::endl;
		return false;
	}

	// set-up fit model (resonances, background, final-state mass-dependence)
	const YAML::Node& configModel = configRoot["model"];
	if(not configModel) {
		printErr << "'model' does not exist in configuration file." << std::endl;
		return false;
	}
	if(not readConfigModel(configModel, fitModel, fitParameters, fitParametersError)) {
		printErr << "error while reading 'model' in configuration file." << std::endl;
		return false;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigFitquality(const YAML::Node& configFitquality,
                                                   int& minStatus,
                                                   double& chi2,
                                                   unsigned int& ndf) const
{
	if(not configFitquality) {
		printErr << "'configFitquality' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configFitquality.IsMap()) {
		printErr << "'fitquality' is not a YAML map." << std::endl;
		return false;
	}

	if(configFitquality["minStatus"]) {
		if(checkVariableType(configFitquality["minStatus"], YamlCppUtils::TypeInt)) {
			minStatus = configFitquality["minStatus"].as<int>();
		} else {
			printErr << "variable 'minStatus' of 'fitquality' is not a floating point number." << std::endl;
			return false;
		}
	} else {
		minStatus = 0;
	}

	if(configFitquality["chi2"]) {
		if(checkVariableType(configFitquality["chi2"], YamlCppUtils::TypeFloat)) {
			chi2 = configFitquality["chi2"].as<double>();
		} else {
			printErr << "variable 'chi2' of 'fitquality' is not a floating point number." << std::endl;
			return false;
		}
	} else {
		chi2 = 0.;
	}

	if(configFitquality["ndf"]) {
		if(checkVariableType(configFitquality["ndf"], YamlCppUtils::TypeInt)) {
			ndf = configFitquality["ndf"].as<unsigned int>();
		} else {
			printErr << "variable 'ndf' of 'fitquality' is not a floating point number." << std::endl;
			return false;
		}
	} else {
		ndf = 0;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigInput(const YAML::Node& configInput)
{
	if(not configInput) {
		printErr << "'configInput' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configInput.IsMap()) {
		printErr << "'input' is not a YAML map." << std::endl;
		return false;
	}

	// get information about fit results from mass-independent
	const YAML::Node& configInputFitResults = configInput["fitresults"];
	if(not configInputFitResults) {
		printErr << "'fitresults' does not exist in 'input'." << std::endl;
		return false;
	}
	if(not readConfigInputFitResults(configInputFitResults)) {
		printErr << "error while reading 'fitresults' in 'input'." << std::endl;
		return false;
	}

	// get information about waves to be used in the fit
	const YAML::Node& configInputWaves = configInput["waves"];
	if(not configInputWaves) {
		printErr << "'waves' does not exist in 'input'." << std::endl;
		return false;
	}
	if(not readConfigInputWaves(configInputWaves)) {
		printErr << "error while reading 'waves' in 'input'." << std::endl;
		return false;
	}

	// get information for plotting of systematic error
	const YAML::Node& configInputSystematics = configInput["systematics"];
	if(configInputSystematics) {
		if(not readConfigInputSystematics(configInputSystematics)) {
			printErr << "error while reading 'systematics' in 'input'." << std::endl;
			return false;
		}
	} else {
		_sysPlotting = false;
	}

	// get information for which parameters to release in which order
	const YAML::Node& configInputFreeParameters = configInput["freeparameters"];
	if(configInputFreeParameters) {
		if(not readConfigInputFreeParameters(configInputFreeParameters)) {
			printErr << "error while reading 'freeparameters' in 'input'." << std::endl;
			return false;
		}
	} else {
		_freeParameters.clear();
		_freeParameters.push_back("coupling branching");
		_freeParameters.push_back("coupling branching mass m0");
		_freeParameters.push_back("*");
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigInputFitResults(const YAML::Node& configInputFitResults)
{
	if(not configInputFitResults) {
		printErr << "'configInputFitResults' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configInputFitResults.IsSequence()) {
		printErr << "'fitresults' is not a YAML sequence." << std::endl;
		return false;
	}

	_inFileName.clear();

	const size_t nrFitResults = configInputFitResults.size();
	for(size_t idxFitResult=0; idxFitResult<nrFitResults; ++idxFitResult) {
		if(_debug) {
			printDebug << "reading of entry " << idxFitResult << " in 'fitresults'." << std::endl;
		}

		const YAML::Node& configInputFitResult = configInputFitResults[idxFitResult];

		std::map<std::string, YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("name", YamlCppUtils::TypeString)
		                     ("tPrimeMean", YamlCppUtils::TypeFloat);
		if(not checkIfAllVariablesAreThere(configInputFitResult, mandatoryArguments)) {
			printErr << "'fitresults' entry at index " << idxFitResult << " does not contain all required variables." << std::endl;
			return false;
		}

		const std::string fileName = configInputFitResult["name"].as<std::string>();
		_inFileName.push_back(fileName);

		if(_debug) {
			printDebug << "read file name of fit results of mass-independent fit: '" << fileName << "'." << std::endl;
		}

		const double tPrimeMean = configInputFitResult["tPrimeMean"].as<double>();
		_tPrimeMeans.push_back(tPrimeMean);

		if(_debug) {
			printDebug << "read mean t' value: '" << tPrimeMean << "'." << std::endl;
		}

		if(configInputFitResult["rescaleErrors"]) {
			if(checkVariableType(configInputFitResult["rescaleErrors"], YamlCppUtils::TypeFloat)) {
				_rescaleErrors.push_back(configInputFitResult["rescaleErrors"].as<double>());
			} else {
				printErr << "variable 'rescaleErrors' of 'fitresults' entry at index " << idxFitResult << " is not a floating point number." << std::endl;
				return false;
			}
		} else {
			_rescaleErrors.push_back(1.);
		}

		if(_debug) {
			printDebug << "rescale errors by factor: '" << _rescaleErrors.back() << "'." << std::endl;
		}
	}

	_nrBins = _inFileName.size();

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigInputWaves(const YAML::Node& configInputWaves)
{
	if(not configInputWaves) {
		printErr << "'configInputWaves' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configInputWaves.IsSequence()) {
		printErr << "'waves' is not a YAML sequence." << std::endl;
		return false;
	}

	_nrWaves = configInputWaves.size();
	if(_debug) {
		printDebug << "going to read information of " << _nrWaves << " waves to be used in the fit." << std::endl;
	}

	for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
		const YAML::Node& configInputWave = configInputWaves[idxWave];

		std::map<std::string, YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("name", YamlCppUtils::TypeString);
		if(not checkIfAllVariablesAreThere(configInputWave, mandatoryArguments)) {
			printErr << "'waves' entry at index " << idxWave << " does not contain all required variables." << std::endl;
			return false;
		}

		const std::string name = configInputWave["name"].as<std::string>();

		double massLower = -1.;
		if(configInputWave["massLower"]) {
			if(checkVariableType(configInputWave["massLower"], YamlCppUtils::TypeFloat)) {
				massLower = configInputWave["massLower"].as<double>();
			} else {
				printErr << "variable 'massLower' of 'waves' entry at index " << idxWave << " is not a floating point number." << std::endl;
				return false;
			}
		}
		double massUpper = -1.;
		if(configInputWave["massUpper"]) {
			if(checkVariableType(configInputWave["massUpper"], YamlCppUtils::TypeFloat)) {
				massUpper = configInputWave["massUpper"].as<double>();
			} else {
				printErr << "variable 'massUpper' of 'waves' entry at index " << idxWave << " is not a floating point number." << std::endl;
				return false;
			}
		}

		std::vector<std::string> alternativeNames;
		if(configInputWave["alternativeNames"]) {
			if(checkVariableType(configInputWave["alternativeNames"], YamlCppUtils::TypeSequence)) {
				for(size_t idxAlt = 0; idxAlt < configInputWave["alternativeNames"].size(); ++idxAlt) {
					if(not checkVariableType(configInputWave["alternativeNames"][idxAlt], YamlCppUtils::TypeString)) {
						printErr << "element " << idxAlt << " of variable 'alternativeNames' of 'waves' entry at index " << idxWave << " is not a string." << std::endl;
						return false;
					}
					const std::string alternativeName = configInputWave["alternativeNames"][idxAlt].as<std::string>();

					// check that the alternative name does not yet exist
					if(alternativeName == name) {
						printErr << "alternative name '" << alternativeName << "' is equal to name of wave '" << name << "'." << std::endl;
						return false;
					}
					if(find(alternativeNames.begin(), alternativeNames.end(), alternativeName) != alternativeNames.end()) {
						printErr << "alternative name '" << alternativeName << "' of wave '" << name << "' defined twice." << std::endl;
						return false;
					}
					if(find(_waveNames.begin(), _waveNames.end(), alternativeName) != _waveNames.end()) {
						printErr << "alternative name '" << alternativeName << "' of wave '" << name << "' already defined as separate wave." << std::endl;
						return false;
					}
					for(size_t i = 0; i < _waveNameAlternatives.size(); ++i) {
						if(find(_waveNameAlternatives[i].begin(), _waveNameAlternatives[i].end(), alternativeName) != _waveNameAlternatives[i].end()) {
							printErr << "alternative name '" << alternativeName << "' of wave '" << name << "' already defined as alternative name of wave '" << _waveNames[i] << "'." << std::endl;
							return false;
						}
					}

					alternativeNames.push_back(alternativeName);
				}
			} else {
				printErr << "variable 'alternativeNames' of 'waves' entry at index " << idxWave << " is not a sequence." << std::endl;
				return false;
			}
		}

		// check that wave does not yet exist
		if(find(_waveNames.begin(), _waveNames.end(), name) != _waveNames.end()) {
			printErr << "wave '" << name << "' defined twice." << std::endl;
			return false;
		}
		for(size_t i = 0; i < _waveNameAlternatives.size(); ++i) {
			if(find(_waveNameAlternatives[i].begin(), _waveNameAlternatives[i].end(), name) != _waveNameAlternatives[i].end()) {
				printErr << "wave '" << name << "' already defined as alternative name of wave '" << _waveNames[i] << "'." << std::endl;
				return false;
			}
		}

		_waveNames.push_back(name);
		_waveNameAlternatives.push_back(alternativeNames);
		_waveIndices[name] = _waveNames.size() - 1;
		for(size_t idxAlt = 0; idxAlt < alternativeNames.size(); ++idxAlt) {
			_waveIndices[alternativeNames[idxAlt]] = _waveNames.size() - 1;
		}
		_waveMassLimits.push_back(std::make_pair(massLower, massUpper));

		if(_debug) {
			std::ostringstream alternatives;
			if(alternativeNames.size() > 0) {
				alternatives << ", alternative names: ";
				for(size_t idxAlt = 0; idxAlt < alternativeNames.size(); ++idxAlt) {
					if(idxAlt != 0) {
						alternatives << ", ";
					}
					alternatives << "'" << alternativeNames[idxAlt] << "'";
				}
			}

			printDebug << idxWave << ": " << name << " (mass range: " << massLower << "-" << massUpper << " GeV/c^2, index: " << _waveIndices[name] << alternatives.str() << ")" << std::endl;
		}
	}

	if(_debug) {
		printDebug << "read " << _nrWaves << " in total." << std::endl;
	}

	std::ostringstream output;
	for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
		output << "    " << _waveNames[idxWave] << std::endl;
	}
	printInfo << _nrWaves << " waves to be used in fit:" << std::endl
	          << output.str();

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigInputSystematics(const YAML::Node& configInputSystematics)
{
	if(not configInputSystematics) {
		printErr << "'configInputSystematics' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configInputSystematics.IsSequence()) {
		printErr << "'systematics' is not a YAML sequence." << std::endl;
		return false;
	}

	_nrSystematics = configInputSystematics.size();
	if(_debug) {
		printDebug << "going to read information for " << _nrSystematics << " files containing information for systematic errors." << std::endl;
	}

	if(_nrSystematics > 0) {
		_sysPlotting = true;
	}

	for(size_t idxSystematics=0; idxSystematics<_nrSystematics; ++idxSystematics) {
		if(not checkVariableType(configInputSystematics[idxSystematics], YamlCppUtils::TypeString)) {
			printErr << "'systematics' entry at index " << idxSystematics << " is not a string." << std::endl;
			return false;
		}

		const std::string fileName = configInputSystematics[idxSystematics].as<std::string>();
		if(_debug) {
			printDebug << "'" << fileName << "' will be read to get information for systematic errors." << std::endl;
		}
		_sysFileNames.push_back(fileName);
	}

	// this also includes the main fitresult file
	++_nrSystematics;
	if(_debug) {
		printDebug << "in total " << _nrSystematics << " files to be read to get information for systematic errors." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigInputFreeParameters(const YAML::Node& configInputFreeParameters)
{
	if(not configInputFreeParameters) {
		printErr << "'configInputFreeParameters' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configInputFreeParameters.IsSequence()) {
		printErr << "'freeparameters' is not a YAML sequence." << std::endl;
		return false;
	}

	_freeParameters.clear();

	const size_t nrItems = configInputFreeParameters.size();
	if(nrItems == 0) {
		printErr << "'freeparameters' is an empty sequence, when defined it must at least contain one entry." << std::endl;
		return false;
	}

	if(_debug) {
		printDebug << "going to extract " << nrItems << " items from 'freeparameters'." << std::endl;
	}

	for(size_t idxItem=0; idxItem<nrItems; ++idxItem) {
		if(not checkVariableType(configInputFreeParameters[idxItem], YamlCppUtils::TypeString)) {
			printErr << "'freeparameters' entry at index " << idxItem << " is not a string." << std::endl;
			return false;
		}

		const std::string name = configInputFreeParameters[idxItem].as<std::string>();
		if(_debug) {
			printDebug << idxItem << ": '" << name << "'." << std::endl;
		}
		_freeParameters.push_back(name);
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigModel(const YAML::Node& configModel,
                                              rpwa::massDepFit::model& fitModel,
                                              rpwa::massDepFit::parameters& fitParameters,
                                              rpwa::massDepFit::parameters& fitParametersError)
{
	if(not configModel) {
		printErr << "'configModel' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configModel.IsMap()) {
		printErr << "'model' is not a YAML map." << std::endl;
		return false;
	}

	if(_debug) {
		printDebug << "reading fit model from configuration file." << std::endl;
	}

	// get information about anchor wave
	const YAML::Node& configAnchorWave = configModel["anchorwave"];
	if(not configAnchorWave) {
		printErr << "'anchorwave' does not exist in 'model'." << std::endl;
		return false;
	}
	if(not readConfigModelAnchorWave(configAnchorWave)) {
		printErr << "error while reading 'anchorwave' in 'model'." << std::endl;
		return false;
	}

	// read information for the individual components
	const YAML::Node& configComponents = configModel["components"];
	if(not configComponents) {
		printErr << "'components' does not exist in 'model'." << std::endl;
		return false;
	}
	if(not readConfigModelComponents(configComponents, fitModel, fitParameters, fitParametersError)) {
		printErr << "error while reading 'components' in 'model'." << std::endl;
		return false;
	}

	// get information for creating the final-state mass-dependence
	const YAML::Node& configFsmd = configModel["finalStateMassDependence"];
	if(configFsmd) {
		if(not readConfigModelFsmd(configFsmd, fitModel, fitParameters, fitParametersError)) {
			printErr << "error while reading 'finalStateMassDependence' in 'model'." << std::endl;
			return false;
		}
	} else {
		printInfo << "not using final-state mass-dependence." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigModelAnchorWave(const YAML::Node& configAnchorWave)
{
	if(not configAnchorWave) {
		printErr << "'configInputAnchorWave' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configAnchorWave.IsMap()) {
		printErr << "'anchorwave' is not a YAML map." << std::endl;
		return false;
	}

	std::map<std::string, YamlCppUtils::Type> mandatoryArguments;
	boost::assign::insert(mandatoryArguments)
	                     ("name", YamlCppUtils::TypeString)
	                     ("resonance", YamlCppUtils::TypeString);
	if(not checkIfAllVariablesAreThere(configAnchorWave, mandatoryArguments)) {
		printErr << "'anchorwave' does not contain all required variables." << std::endl;
		return false;
	}

	_anchorWaveName = configAnchorWave["name"].as<std::string>();
	_anchorComponentName = configAnchorWave["resonance"].as<std::string>();

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigModelComponents(const YAML::Node& configComponents,
                                                        rpwa::massDepFit::model& fitModel,
                                                        rpwa::massDepFit::parameters& fitParameters,
                                                        rpwa::massDepFit::parameters& fitParametersError) const
{
	if(not configComponents) {
		printErr << "'configComponents' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configComponents.IsSequence()) {
		printErr << "'components' is not a YAML sequence." << std::endl;
		return false;
	}

	const size_t nrComponents = configComponents.size();

	if(_debug) {
		printDebug << "reading " << nrComponents << " components from configuration file." << std::endl;
	}

	for(size_t idxComponent=0; idxComponent<nrComponents; ++idxComponent) {
		const YAML::Node& configComponent = configComponents[idxComponent];

		std::map<std::string, YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("name", YamlCppUtils::TypeString);
		if(not checkIfAllVariablesAreThere(configComponent, mandatoryArguments)) {
			printErr << "'components' entry at index " << idxComponent << " does not contain all required variables." << std::endl;
			return false;
		}

		const std::string name = configComponent["name"].as<std::string>();

		for(size_t idx=0; idx<fitModel.getNrComponents(); ++idx) {
			if(fitModel.getComponent(idx)->getName() == name) {
				printErr << "component '" << name << "' defined twice." << std::endl;
				return false;
			}
		}

		std::string type;
		if(configComponent["type"]) {
			if(not checkVariableType(configComponent["type"], YamlCppUtils::TypeString)) {
				printErr << "component '" << name << "' has a type that is not a string." << std::endl;
				return false;
			}
			type = configComponent["type"].as<std::string>();
		} else {
			if(_debug) {
				printDebug << "component '" << name << "' has no type, use 'fixedWidthBreitWigner'." << std::endl;
			}
			type = "fixedWidthBreitWigner";
		}

		if(_debug) {
			printDebug << "found component '" << name << "' with type '" << type << "'." << std::endl;
		}

		rpwa::massDepFit::componentPtr component;
		if(type == "fixedWidthBreitWigner") {
			component = std::make_shared<rpwa::massDepFit::fixedWidthBreitWigner>(fitModel.getNrComponents(), name);
		} else if(type == "dynamicWidthBreitWigner") {
			component = std::make_shared<rpwa::massDepFit::dynamicWidthBreitWigner>(fitModel.getNrComponents(), name);
		} else if(type == "integralWidthBreitWigner") {
			component = std::make_shared<rpwa::massDepFit::integralWidthBreitWigner>(fitModel.getNrComponents(), name);
		} else if(type == "constantBackground") {
			component = std::make_shared<rpwa::massDepFit::constantBackground>(fitModel.getNrComponents(), name);
		} else if(type == "exponentialBackground") {
			component = std::make_shared<rpwa::massDepFit::exponentialBackground>(fitModel.getNrComponents(), name);
		} else if(type == "tPrimeDependentBackground") {
			component = std::make_shared<rpwa::massDepFit::tPrimeDependentBackground>(fitModel.getNrComponents(), name);
			std::dynamic_pointer_cast<tPrimeDependentBackground>(component)->setTPrimeMeans(_tPrimeMeans);
		} else if(type == "exponentialBackgroundIntegral") {
			component = std::make_shared<rpwa::massDepFit::exponentialBackgroundIntegral>(fitModel.getNrComponents(), name);
		} else if(type == "tPrimeDependentBackgroundIntegral") {
			component = std::make_shared<rpwa::massDepFit::tPrimeDependentBackgroundIntegral>(fitModel.getNrComponents(), name);
			std::dynamic_pointer_cast<tPrimeDependentBackgroundIntegral>(component)->setTPrimeMeans(_tPrimeMeans);
		} else {
			printErr << "unknown type '" << type << "' for component '" << name << "'." << std::endl;
			return false;
		}

		if(not component->init(configComponent, fitParameters, fitParametersError, _nrMassBins, _massBinCenters, _waveIndices, _waveBins, _inPhaseSpaceIntegrals, fitModel.useBranchings(), _debug)) {
			printErr << "error while initializing component '" << name << "' of type '" << type << "'." << std::endl;
			return false;
		}

		fitModel.add(component);
	}

	std::ostringstream output;
	for(size_t idxComponent=0; idxComponent<fitModel.getNrComponents(); ++idxComponent) {
		output << "    " << fitModel.getComponent(idxComponent)->getName() << std::endl;
	}
	printInfo << "fitting " << fitModel.getNrComponents() << " components to the data:" << std::endl
	          << output.str();

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigModelFsmd(const YAML::Node& configFsmd,
                                                  rpwa::massDepFit::model& fitModel,
                                                  rpwa::massDepFit::parameters& fitParameters,
                                                  rpwa::massDepFit::parameters& fitParametersError) const
{
	if(not configFsmd) {
		printErr << "'configFsmd' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configFsmd.IsMap() and not configFsmd.IsSequence()) {
		printErr << "'finalStateMassDependence' is not a YAML map or sequence." << std::endl;
		return false;
	}

	if(_debug) {
		printDebug << "reading final-state mass-dependence from configuration file." << std::endl;
	}

	rpwa::massDepFit::fsmdPtr fsmd(new rpwa::massDepFit::fsmd(fitModel.getNrComponents()));
	if(not fsmd->init(configFsmd, fitParameters, fitParametersError, _nrBins, _sameMassBinning, _debug)) {
		printErr << "error while initializing final-state mass-dependence." << std::endl;
		return false;
	}
	fitModel.setFsmd(fsmd);

	printInfo << "using final-state mass-dependence as defined in the configuration file." << std::endl;

	return true;
}


bool
rpwa::massDepFit::massDepFit::init(rpwa::massDepFit::model& fitModel,
                                   rpwa::massDepFit::function& fitFunction)
{
	if(not fitModel.init(_nrBins,
	                     _waveNames,
	                     _waveNameAlternatives,
	                     _anchorWaveName,
	                     _anchorComponentName)) {
		printErr << "error while initializing the fit model." << std::endl;
		return false;
	}

	if(not fitFunction.init(&fitModel,
	                        _nrMassBins,
	                        _massBinCenters,
	                        _inProductionAmplitudes,
	                        _inProductionAmplitudesCovariance,
	                        _inSpinDensityMatrices,
	                        _inSpinDensityCovarianceMatrices,
	                        _wavePairMassBinLimits)) {
		printErr << "error while initializing the function to minimize." << std::endl;
		return false;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::writeConfig(std::ostream& output,
                                          const rpwa::massDepFit::model& fitModel,
                                          const rpwa::massDepFit::parameters& fitParameters,
                                          const rpwa::massDepFit::parameters& fitParametersError,
                                          const int minStatus,
                                          const double chi2,
                                          const unsigned int ndf) const
{
	if(_debug) {
		printDebug << "writing configuration file." << std::endl;
	}

	YAML::Emitter yamlOutput(output);
	yamlOutput << YAML::BeginMap;

	yamlOutput << YAML::Key << "fitquality";
	yamlOutput << YAML::Value;
	if(not writeConfigFitquality(yamlOutput, minStatus, chi2, ndf)) {
		printErr << "error while writing 'fitquality' to result file." << std::endl;
		return false;
	}

	yamlOutput << YAML::Key << "input";
	yamlOutput << YAML::Value;
	if(not writeConfigInput(yamlOutput)) {
		printErr << "error while writing 'input' to result file." << std::endl;
		return false;
	}

	yamlOutput << YAML::Key << "model";
	yamlOutput << YAML::Value;
	if(not writeConfigModel(yamlOutput, fitModel, fitParameters, fitParametersError)) {
		printErr << "error while writing 'model' to result file." << std::endl;
		return false;
	}

	yamlOutput << YAML::EndMap;

	// newline at end-of-file
	output << std::endl;

	return true;
}


bool
rpwa::massDepFit::massDepFit::writeConfigFitquality(YAML::Emitter& yamlOutput,
                                                    const int minStatus,
                                                    const double chi2,
                                                    const unsigned int ndf) const
{
	if(_debug) {
		printDebug << "writing 'fitquality'." << std::endl;
	}

	yamlOutput << YAML::BeginMap;
	yamlOutput << YAML::Key << "minStatus";
	yamlOutput << YAML::Value << minStatus;
	yamlOutput << YAML::Key << "chi2";
	yamlOutput << YAML::Value << chi2;
	yamlOutput << YAML::Key << "ndf";
	yamlOutput << YAML::Value << ndf;
	yamlOutput << YAML::Key << "redchi2";
	yamlOutput << YAML::Value << ((ndf>0) ? (chi2/(double)ndf) : 0.);
	yamlOutput << YAML::EndMap;

	return true;
}


bool
rpwa::massDepFit::massDepFit::writeConfigInput(YAML::Emitter& yamlOutput) const
{
	if(_debug) {
		printDebug << "writing 'input'." << std::endl;
	}

	yamlOutput << YAML::BeginMap;

	yamlOutput << YAML::Key << "fitresults";
	yamlOutput << YAML::Value;
	if(not writeConfigInputFitResults(yamlOutput)) {
		printErr << "error while writing 'fitresults' to result file." << std::endl;
		return false;
	}

	yamlOutput << YAML::Key << "waves";
	yamlOutput << YAML::Value;
	if(not writeConfigInputWaves(yamlOutput)) {
		printErr << "error while writing 'waves' to result file." << std::endl;
		return false;
	}

	if(_sysPlotting) {
		yamlOutput << YAML::Key << "systematics";
		yamlOutput << YAML::Value;
		if(not writeConfigInputSystematics(yamlOutput)) {
			printErr << "error while writing 'systematics' to result file." << std::endl;
			return false;
		}
	}

	yamlOutput << YAML::Key << "freeparameters";
	yamlOutput << YAML::Value;
	if(not writeConfigInputFreeParameters(yamlOutput)) {
		printErr << "error while writing 'freeparameters' to result file." << std::endl;
		return false;
	}

	yamlOutput << YAML::EndMap;

	return true;
}


bool
rpwa::massDepFit::massDepFit::writeConfigInputFitResults(YAML::Emitter& yamlOutput) const
{
	if(_debug) {
		printDebug << "writing 'fitresults'." << std::endl;
	}

	yamlOutput << YAML::BeginSeq;

	for (size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
		yamlOutput << YAML::BeginMap;

		yamlOutput << YAML::Key << "name";
		yamlOutput << YAML::Value << _inFileName[idxBin];

		yamlOutput << YAML::Key << "tPrimeMean";
		yamlOutput << YAML::Value << _tPrimeMeans[idxBin];

		if(_rescaleErrors[idxBin] != 1.) {
			yamlOutput << YAML::Key << "rescaleErrors";
			yamlOutput << YAML::Value << _rescaleErrors[idxBin];
		}

		yamlOutput << YAML::EndMap;
	}

	yamlOutput << YAML::EndSeq;

	return true;
}


bool
rpwa::massDepFit::massDepFit::writeConfigInputWaves(YAML::Emitter& yamlOutput) const
{
	if(_debug) {
		printDebug << "writing 'waves'." << std::endl;
	}

	yamlOutput << YAML::BeginSeq;

	for (size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
		yamlOutput << YAML::BeginMap;

		yamlOutput << YAML::Key << "name";
		yamlOutput << YAML::Value << _waveNames[idxWave];

		if (_waveMassLimits[idxWave].first >= 0) {
			yamlOutput << YAML::Key << "massLower";
			yamlOutput << YAML::Value << _waveMassLimits[idxWave].first;
		}
		if (_waveMassLimits[idxWave].second >= 0) {
			yamlOutput << YAML::Key << "massUpper";
			yamlOutput << YAML::Value << _waveMassLimits[idxWave].second;
		}

		if(_waveNameAlternatives[idxWave].size() > 0) {
			yamlOutput << YAML::Key << "alternativeNames";
			yamlOutput << YAML::Value;

			yamlOutput << YAML::Flow;
			yamlOutput << YAML::BeginSeq;

			for(size_t idxAlt = 0; idxAlt < _waveNameAlternatives[idxWave].size(); ++idxAlt) {
				yamlOutput << _waveNameAlternatives[idxWave][idxAlt];
			}

			yamlOutput << YAML::EndSeq;
			yamlOutput << YAML::Block;
		}

		yamlOutput << YAML::EndMap;
	}

	yamlOutput << YAML::EndSeq;

	return true;
}


bool
rpwa::massDepFit::massDepFit::writeConfigInputSystematics(YAML::Emitter& yamlOutput) const
{
	if(_debug) {
		printDebug << "writing 'systematics'." << std::endl;
	}

	yamlOutput << _sysFileNames;

	return true;
}


bool
rpwa::massDepFit::massDepFit::writeConfigInputFreeParameters(YAML::Emitter& yamlOutput) const
{
	if(_debug) {
		printDebug << "writing 'freeparameters'." << std::endl;
	}

	yamlOutput << _freeParameters;

	return true;
}


bool
rpwa::massDepFit::massDepFit::writeConfigModel(YAML::Emitter& yamlOutput,
                                               const rpwa::massDepFit::model& fitModel,
                                               const rpwa::massDepFit::parameters& fitParameters,
                                               const rpwa::massDepFit::parameters& fitParametersError) const
{
	if(_debug) {
		printDebug << "writing 'model'." << std::endl;
	}

	yamlOutput << YAML::BeginMap;

	yamlOutput << YAML::Key << "anchorwave";
	yamlOutput << YAML::Value;
	if(not writeConfigModelAnchorWave(yamlOutput)) {
		printErr << "error while writing 'anchorwave' to result file." << std::endl;
		return false;
	}

	yamlOutput << YAML::Key << "components";
	yamlOutput << YAML::Value;
	if(not writeConfigModelComponents(yamlOutput, fitModel, fitParameters, fitParametersError)) {
		printErr << "error while writing 'components' to result file." << std::endl;
		return false;
	}

	yamlOutput << YAML::Key << "finalStateMassDependence";
	yamlOutput << YAML::Value;
	if(fitModel.getFsmd()) {
		if(not writeConfigModelFsmd(yamlOutput, fitModel, fitParameters, fitParametersError)) {
			printErr << "error while writing 'finalStateMassDependence' to result file." << std::endl;
			return false;
		}
	}

	yamlOutput << YAML::EndMap;

	return true;
}


bool
rpwa::massDepFit::massDepFit::writeConfigModelAnchorWave(YAML::Emitter& yamlOutput) const
{
	if(_debug) {
		printDebug << "writing 'anchorwave'." << std::endl;
	}

	yamlOutput << YAML::BeginMap;

	yamlOutput << YAML::Key << "name";
	yamlOutput << YAML::Value << _anchorWaveName;

	yamlOutput << YAML::Key << "resonance";
	yamlOutput << YAML::Value << _anchorComponentName;

	yamlOutput << YAML::EndMap;

	return true;
}


bool
rpwa::massDepFit::massDepFit::writeConfigModelComponents(YAML::Emitter& yamlOutput,
                                                         const rpwa::massDepFit::model& fitModel,
                                                         const rpwa::massDepFit::parameters& fitParameters,
                                                         const rpwa::massDepFit::parameters& fitParametersError) const
{
	if(_debug) {
		printDebug << "writing 'components'." << std::endl;
	}

	yamlOutput << YAML::BeginSeq;

	const size_t nrComponents = fitModel.getNrComponents();
	for(size_t idxComponent=0; idxComponent<nrComponents; ++idxComponent) {
		if(not fitModel.getComponent(idxComponent)->write(yamlOutput, fitParameters, fitParametersError, fitModel.useBranchings(), _debug)) {
			printErr << "error while writing component at index " << idxComponent << " to result file." << std::endl;
			return false;
		}
	}

	yamlOutput << YAML::EndSeq;

	return true;
}


bool
rpwa::massDepFit::massDepFit::writeConfigModelFsmd(YAML::Emitter& yamlOutput,
                                                   const rpwa::massDepFit::model& fitModel,
                                                   const rpwa::massDepFit::parameters& fitParameters,
                                                   const rpwa::massDepFit::parameters& fitParametersError) const
{
	if(_debug) {
		printDebug << "writing 'finalStateMassDependence'." << std::endl;
	}

	if(not fitModel.getFsmd()) {
		printErr << "writing final-state mass-dependence requested, but there is no final-state mass-dependence." << std::endl;
		return false;
	}

	if(not fitModel.getFsmd()->write(yamlOutput, fitParameters, fitParametersError, _debug)) {
		printErr << "error while writing final-state mass-dependence to result file." << std::endl;
		return false;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readInFiles(const std::string& valTreeName,
                                          const std::string& valBranchName)
{
	for(size_t idxBin = 0; idxBin < _nrBins; ++idxBin) {
		if(not readInFile(idxBin, valTreeName, valBranchName)) {
			printErr << "error while reading file entry " << idxBin << "." << std::endl;
			return false;
		}
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readInFile(const size_t idxBin,
                                         const std::string& valTreeName,
                                         const std::string& valBranchName)
{
	if(_debug) {
		printDebug << "reading fit result from file '" << _inFileName[idxBin] << "'." << std::endl;
	}

	TFile* inFile = TFile::Open(_inFileName[idxBin].c_str());
	if(not inFile) {
		printErr << "input file '" << _inFileName[idxBin] << "' not found."<< std::endl;
		return false;
	}
	if(inFile->IsZombie()) {
		printErr << "error while reading input file '" << _inFileName[idxBin] << "'."<< std::endl;
		delete inFile;
		return false;
	}

	if(_debug) {
		printDebug << "searching for tree '" << valTreeName << "' in file '" << _inFileName[idxBin] << "'." << std::endl;
	}

	TTree* inTree;
	inFile->GetObject(valTreeName.c_str(), inTree);
	if(not inTree) {
		printErr << "input tree '" << valTreeName << "' not found in input file '" << _inFileName[idxBin] << "'."<< std::endl;
		delete inFile;
		return false;
	}

	if(_debug) {
		printDebug << "searching for branch '" << valBranchName << "' in tree '" << valTreeName << "'." << std::endl;
	}

	fitResult* inFit = NULL;
	if(inTree->SetBranchAddress(valBranchName.c_str(), &inFit)) {
		printErr << "branch '" << valBranchName << "' not found in input tree '" << valTreeName << "'." << std::endl;
		delete inFile;
		return false;
	}

	double massMin;
	double massMax;
	double massStep;
	size_t nrMassBins;
	boost::multi_array<double, 1> massBinCenters;
	if(not readFitResultMassBins(inTree, inFit, massMin, massMax, massStep, nrMassBins, massBinCenters)) {
		printErr << "could not extract mass bins from fit result tree in '" << _inFileName[idxBin] << "'." << std::endl;
		delete inFile;
		return false;
	}

	if(_maxMassBins < nrMassBins) {
		// take care to correctly copy the covariance matrices
		boost::multi_array<TMatrixT<double>, 2> tempProductionAmplitudesCovariance(std::vector<size_t>(_inProductionAmplitudesCovariance.shape(), _inProductionAmplitudesCovariance.shape()+_inProductionAmplitudesCovariance.num_dimensions()));
		boost::multi_array<TMatrixT<double>, 2> tempSpinDensityCovarianceMatrices(std::vector<size_t>(_inSpinDensityCovarianceMatrices.shape(), _inSpinDensityCovarianceMatrices.shape()+_inSpinDensityCovarianceMatrices.num_dimensions()));
		if(idxBin > 0) {
			for(size_t idx = 0; idx < _nrBins; ++idx) {
				for(size_t idxMass = 0; idxMass < _nrMassBins[idx]; ++idxMass) {
					tempProductionAmplitudesCovariance[idx][idxMass].ResizeTo(_inProductionAmplitudesCovariance[idx][idxMass]);
					tempProductionAmplitudesCovariance[idx][idxMass] = _inProductionAmplitudesCovariance[idx][idxMass];
					tempSpinDensityCovarianceMatrices[idx][idxMass].ResizeTo(_inSpinDensityCovarianceMatrices[idx][idxMass]);
					tempSpinDensityCovarianceMatrices[idx][idxMass] = _inSpinDensityCovarianceMatrices[idx][idxMass];
				}
			}
			_inProductionAmplitudesCovariance.resize(boost::extents[0][0]);
			_inSpinDensityCovarianceMatrices.resize(boost::extents[0][0]);
		}

		_maxMassBins = nrMassBins;

		// resize all array to store the information
		_massMaxs.resize(_nrBins);
		_massMins.resize(_nrBins);
		_massSteps.resize(_nrBins);
		_nrMassBins.resize(_nrBins);
		_massBinCenters.resize(boost::extents[_nrBins][_maxMassBins]);

		_inProductionAmplitudes.resize(boost::extents[_nrBins][_maxMassBins][_nrWaves]);
		_inProductionAmplitudesCovariance.resize(boost::extents[_nrBins][_maxMassBins]);
		_inSpinDensityMatrices.resize(boost::extents[_nrBins][_maxMassBins][_nrWaves][_nrWaves]);
		_inSpinDensityCovarianceMatrices.resize(boost::extents[_nrBins][_maxMassBins]);
		_inPhaseSpaceIntegrals.resize(boost::extents[_nrBins][_maxMassBins][_nrWaves]);
		_inIntensities.resize(boost::extents[_nrBins][_maxMassBins][_nrWaves][2]);
		_inPhases.resize(boost::extents[_nrBins][_maxMassBins][_nrWaves][_nrWaves][2]);

		if(idxBin > 0) {
			for(size_t idx = 0; idx < _nrBins; ++idx) {
				for(size_t idxMass = 0; idxMass < _nrMassBins[idx]; ++idxMass) {
					_inProductionAmplitudesCovariance[idx][idxMass].ResizeTo(tempProductionAmplitudesCovariance[idx][idxMass]);
					_inProductionAmplitudesCovariance[idx][idxMass] = tempProductionAmplitudesCovariance[idx][idxMass];
					_inSpinDensityCovarianceMatrices[idx][idxMass].ResizeTo(tempSpinDensityCovarianceMatrices[idx][idxMass]);
					_inSpinDensityCovarianceMatrices[idx][idxMass] = tempSpinDensityCovarianceMatrices[idx][idxMass];
				}
			}
		}
	}

	_massMaxs[idxBin] = massMax;
	_massMins[idxBin] = massMin;
	_massSteps[idxBin] = massStep;
	_nrMassBins[idxBin] = nrMassBins;
	_massBinCenters[idxBin] = massBinCenters;

	bool readMapping(false);
	std::vector<Long64_t> inMapping;
	if(_sameMassBinning and idxBin > 0) {
		if(checkFitResultMassBins(inTree, inFit, idxBin-1, inMapping)) {
			if(_debug) {
				printDebug << "bin " << idxBin << " has the same mass binning as bin " << idxBin-1 << "." << std::endl;
			}

			readMapping = true;
		} else {
			printInfo << "bin " << idxBin << " does not have the same mass binning as bin " << idxBin-1 << ", use individual mass binning for each bin." << std::endl;

			_sameMassBinning = false;
		}
	}
	if(not readMapping) {
		if(not checkFitResultMassBins(inTree, inFit, idxBin, inMapping)) {
			printErr << "error while checking and mapping mass bins from fit result tree in '" << _inFileName[idxBin] << "'." << std::endl;
			delete inFile;
			return false;
		}
	}

	std::vector<std::string> waveNames;
	boost::multi_array<std::complex<double>, 2> tempProductionAmplitudes;
	boost::multi_array<TMatrixT<double>, 1> tempProductionAmplitudesCovariance;
	boost::multi_array<std::complex<double>, 3> tempSpinDensityMatrices;
	boost::multi_array<TMatrixT<double>, 1> tempSpinDensityCovarianceMatrices;
	boost::multi_array<double, 3> tempIntensities;
	boost::multi_array<double, 4> tempPhases;
	if(not readFitResultMatrices(inTree, inFit, inMapping, _rescaleErrors[idxBin], waveNames, tempProductionAmplitudes, tempProductionAmplitudesCovariance,
	                             tempSpinDensityMatrices, tempSpinDensityCovarianceMatrices, tempIntensities, tempPhases)) {
		printErr << "error while reading spin-density matrix from fit result tree in '" << _inFileName[idxBin] << "'." << std::endl;
		delete inFile;
		return false;
	}
	for(std::vector<std::string>::const_iterator it = waveNames.begin(); it != waveNames.end(); ++it) {
		_waveBins[*it].push_back(idxBin);
	}
	_inProductionAmplitudes[idxBin] = tempProductionAmplitudes;
	_inSpinDensityMatrices[idxBin] = tempSpinDensityMatrices;
	_inIntensities[idxBin] = tempIntensities;
	_inPhases[idxBin] = tempPhases;
	for(size_t i = 0; i < nrMassBins; ++i) {
		_inProductionAmplitudesCovariance[idxBin][i].ResizeTo(tempProductionAmplitudesCovariance[i]);
		_inProductionAmplitudesCovariance[idxBin][i] = tempProductionAmplitudesCovariance[i];
		_inSpinDensityCovarianceMatrices[idxBin][i].ResizeTo(tempSpinDensityCovarianceMatrices[i]);
		_inSpinDensityCovarianceMatrices[idxBin][i] = tempSpinDensityCovarianceMatrices[i];
	}

	boost::multi_array<double, 2> tempPhaseSpaceIntegrals;
	if(not readFitResultIntegrals(inTree, inFit, inMapping, waveNames, tempPhaseSpaceIntegrals)) {
		printErr << "error while reading phase-space integrals from fit result tree in '" << _inFileName[idxBin] << "'." << std::endl;
		delete inFile;
		return false;
	}
	_inPhaseSpaceIntegrals[idxBin] = tempPhaseSpaceIntegrals;

	delete inFile;
	return true;
}


bool
rpwa::massDepFit::massDepFit::readSystematicsFiles(const std::string& valTreeName,
                                                   const std::string& valBranchName)
{
	if(not _sysPlotting) {
		return true;
	}

	// FIXME: make systematic errors work with multiple bins
	if(_nrBins > 1) {
		printErr << "systematic errors not yet supported for multiple bins." << std::endl;
		return false;
	}

	if(_debug) {
		printDebug << "reading fit results for systematic errors from " << _nrSystematics << " files." << std::endl;
	}

	_sysSpinDensityMatrices.resize(boost::extents[_nrSystematics][_maxMassBins][_nrWaves][_nrWaves]);
	_sysSpinDensityCovarianceMatrices.resize(boost::extents[_nrSystematics][_maxMassBins]);
	_sysIntensities.resize(boost::extents[_nrSystematics][_maxMassBins][_nrWaves][2]);
	_sysPhases.resize(boost::extents[_nrSystematics][_maxMassBins][_nrWaves][_nrWaves][2]);

	_sysSpinDensityMatrices[0] = _inSpinDensityMatrices[0];
	_sysIntensities[0] = _inIntensities[0];
	_sysPhases[0] = _inPhases[0];
	for(size_t i = 0; i < _maxMassBins; ++i) {
		_sysSpinDensityCovarianceMatrices[0][i].ResizeTo(_inSpinDensityCovarianceMatrices[0][i]);
		_sysSpinDensityCovarianceMatrices[0][i] = _inSpinDensityCovarianceMatrices[0][i];
	}

	for(size_t idxSystematics=1; idxSystematics<_nrSystematics; ++idxSystematics) {
		if(not readSystematicsFile(idxSystematics, valTreeName, valBranchName)) {;
			printErr << "error while reading fit results for systematic errors." << std::endl;
			return false;
		}
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readSystematicsFile(const size_t idxSystematics,
                                                  const std::string& valTreeName,
                                                  const std::string& valBranchName)
{
	if(_debug) {
		printDebug << "reading fit result for systematics for index " << idxSystematics << " from file '" << _sysFileNames[idxSystematics-1] << "'." << std::endl;
	}

	TFile* sysFile = TFile::Open(_sysFileNames[idxSystematics-1].c_str());
	if(not sysFile) {
		printErr << "input file '" << _sysFileNames[idxSystematics-1] << "' not found."<< std::endl;
		return false;
	}
	if(sysFile->IsZombie()) {
		printErr << "error while reading input file '" << _sysFileNames[idxSystematics-1] << "'."<< std::endl;
		delete sysFile;
		return false;
	}

	if(_debug) {
		printDebug << "searching for tree '" << valTreeName << "' in file '" << _sysFileNames[idxSystematics-1] << "'." << std::endl;
	}

	TTree* sysTree;
	sysFile->GetObject(valTreeName.c_str(), sysTree);
	if(not sysTree) {
		printErr << "input tree '" << valTreeName << "' not found in input file '" << _sysFileNames[idxSystematics-1] << "'."<< std::endl;
		delete sysFile;
		return false;
	}

	if(_debug) {
		printDebug << "searching for branch '" << valBranchName << "' in tree '" << valTreeName << "'." << std::endl;
	}

	fitResult* sysFit = NULL;
	if(sysTree->SetBranchAddress(valBranchName.c_str(), &sysFit)) {
		printErr << "branch '" << valBranchName << "' not found in input tree '" << valTreeName << "'." << std::endl;
		delete sysFile;
		return false;
	}

	std::vector<Long64_t> sysMapping;
	if(not checkFitResultMassBins(sysTree, sysFit, 0, sysMapping)) {
		printErr << "error while checking and mapping mass bins from fit result tree in '" << _sysFileNames[idxSystematics-1] << "'." << std::endl;
		delete sysFile;
		return false;
	}

	std::vector<std::string> waveNames;
	boost::multi_array<std::complex<double>, 2> tempProductionAmplitudes;
	boost::multi_array<TMatrixT<double>, 1> tempProductionAmplitudesCovariance;
	boost::multi_array<std::complex<double>, 3> tempSpinDensityMatrices;
	boost::multi_array<TMatrixT<double>, 1> tempSpinDensityCovarianceMatrices;
	boost::multi_array<double, 3> tempIntensities;
	boost::multi_array<double, 4> tempPhases;
	if(not readFitResultMatrices(sysTree, sysFit, sysMapping, _rescaleErrors[0], waveNames, tempProductionAmplitudes, tempProductionAmplitudesCovariance,
	                             tempSpinDensityMatrices, tempSpinDensityCovarianceMatrices, tempIntensities, tempPhases)) {
		printErr << "error while reading spin-density matrix from fit result tree in '" << _sysFileNames[idxSystematics-1] << "'." << std::endl;
		delete sysFile;
		return false;
	}
	_sysSpinDensityMatrices[idxSystematics] = tempSpinDensityMatrices;
	_sysIntensities[idxSystematics] = tempIntensities;
	_sysPhases[idxSystematics] = tempPhases;
	for(size_t i = 0; i < _nrMassBins[0]; ++i) {
		_sysSpinDensityCovarianceMatrices[idxSystematics][i].ResizeTo(tempSpinDensityCovarianceMatrices[i]);
		_sysSpinDensityCovarianceMatrices[idxSystematics][i] = tempSpinDensityCovarianceMatrices[i];
	}

	delete sysFile;
	return true;
}


bool
rpwa::massDepFit::massDepFit::checkFitResultMassBins(TTree* tree,
                                                     rpwa::fitResult* fit,
                                                     const size_t idxBin,
                                                     std::vector<Long64_t>& mapping) const
{
	if(not tree or not fit) {
		printErr << "'tree' or 'fit' is not a pointer to a valid object." << std::endl;
		return false;
	}

	// reset mapping
	mapping.assign(_nrMassBins[idxBin], std::numeric_limits<Long64_t>::max());

	// extract data from tree
	const Long64_t nrEntries = tree->GetEntries();

	if(_debug) {
		printDebug << "check that the centers of mass bins of " << nrEntries << " entries in tree are at a known place, "
		           << "and map the " << _nrMassBins[idxBin] << " mass bins to those entries." << std::endl;
	}

	for(Long64_t idx=0; idx<nrEntries; ++idx) {
		if(tree->GetEntry(idx) == 0) {
			printErr << "error while reading entry " << idx << " from tree." << std::endl;
			return false;
		}
		//FIXME: this would also be the place to select the best fit in case one file contains more than one fit result per mass bin
		const double mass = fit->massBinCenter();

		if(_debug) {
			printDebug << "entry " << idx << ": center of mass bin at " << mass << " GeV/c^2" << std::endl;
		}

		bool found = false;
		size_t idxMass=0;
		while(idxMass<_nrMassBins[idxBin]) {
			if(std::abs(_massBinCenters[idxBin][idxMass]-mass) < 1000.*std::numeric_limits<double>::epsilon()) {
				found = true;
				break;
			}
			++idxMass;
		}

		if(not found) {
			printErr << "could not map mass bin centered at " << mass << " GeV/c^2 to a known mass bin." << std::endl;
			return false;
		}

		if(mapping[idxMass] != std::numeric_limits<Long64_t>::max()) {
			printErr << "cannot map tree entry " << idx << " to mass bin " << idxMass << " (" << _massBinCenters[idxBin][idxMass] << " GeV/c^2)  "
			         << "which is already mapped to tree entry " << mapping[idxMass] << "." << std::endl;
			return false;
		}

		if(_debug) {
			printDebug << "mapping mass bin " << idxMass << " (" << _massBinCenters[idxBin][idxMass] << " GeV/c^2) to tree entry " << idx << "." << std::endl;
		}
		mapping[idxMass] = idx;
	} // end loop over entries in tree

	// check that all mass bins are mapped
	for(size_t idx=0; idx<mapping.size(); ++idx) {
		if(mapping[idx] == std::numeric_limits<Long64_t>::max()) {
			printErr << "mass bin " << idx << " (" << _massBinCenters[idxBin][idx] << " GeV/c^2) not mapped." << std::endl;
			return false;
		}
	}

	if(_debug) {
		std::ostringstream output;
		for(size_t idx=0; idx<mapping.size(); ++idx) {
			output << " " << idx << "->" << mapping[idx];
		}
		printDebug << "established mapping:" << output.str() << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readFitResultMassBins(TTree* tree,
                                                    rpwa::fitResult* fit,
                                                    double& massMin,
                                                    double& massMax,
                                                    double& massStep,
                                                    size_t& nrMassBins,
                                                    boost::multi_array<double, 1>& massBinCenters) const
{
	if(not tree or not fit) {
		printErr << "'tree' or 'fit' is not a pointer to a valid object." << std::endl;
		return false;
	}

	// extract data from tree
	const Long64_t nrEntries = tree->GetEntries();

	if(_debug) {
		printDebug << "getting center of mass bins from " << nrEntries << " entries in tree." << std::endl;
	}

	nrMassBins = 0;
	massBinCenters.resize(boost::extents[_maxMassBins]);
	for(Long64_t idx=0; idx<nrEntries; ++idx) {
		if(tree->GetEntry(idx) == 0) {
			printErr << "error while reading entry " << idx << " from tree." << std::endl;
			return false;
		}
		const double newMass = fit->massBinCenter();

		if(_debug) {
			printDebug << "entry " << idx << ": center of mass bin at " << newMass << " GeV/c^2" << std::endl;
		}

		bool found = false;
		for(size_t idxMass = 0; idxMass < nrMassBins; ++idxMass) {
			if(std::abs(massBinCenters[idxMass]-newMass) < 1000.*std::numeric_limits<double>::epsilon()) {
				found = true;
				if(_debug) {
					printDebug << "this center of mass bin already was encountered before." << std::endl;
				}
				break;
			}
		}

		if(not found) {
			if(nrMassBins >= _maxMassBins) {
				massBinCenters.resize(boost::extents[nrMassBins+1]);
			}

			massBinCenters[nrMassBins++] = newMass;
		}
	} // end loop over entries in tree

	// sort mass bins
	std::sort(massBinCenters.data(), massBinCenters.data() + nrMassBins);

	printInfo << "found " << nrMassBins << " mass bins, center of first and last mass bins: "
	          << massBinCenters[0] << " and " << massBinCenters[nrMassBins - 1] << " GeV/c^2." << std::endl;

	massStep = (massBinCenters[nrMassBins - 1] - massBinCenters[0]) / (nrMassBins - 1);
	for(size_t idxMass = 1; idxMass < nrMassBins; ++idxMass) {
		if(std::abs(massBinCenters[idxMass]-massBinCenters[idxMass-1] - massStep) > 1000.*std::numeric_limits<double>::epsilon()) {
			printErr << "mass distance between bins " << idxMass-1 << " (" << massBinCenters[idxMass-1] << " GeV/c^2) and "
			         << idxMass << " (" << massBinCenters[idxMass] << " GeV/c^2) does not agree with nominal distance "
			         << massStep << " GeV/c^2" << std::endl;
			return false;
		}
	}
	if(_debug) {
		printDebug << "distance between two mass bins is " << massStep << " GeV/c^2." << std::endl;
	}

	massMin = massBinCenters[0] - massStep / 2;
	massMax = massBinCenters[nrMassBins - 1] + massStep / 2;
	if(_debug) {
		printDebug << "mass bins cover the mass range from " << massMin << " to " << massMax << " GeV/c^2." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readFitResultMatrices(TTree* tree,
                                                    rpwa::fitResult* fit,
                                                    const std::vector<Long64_t>& mapping,
                                                    const double rescaleErrors,
                                                    std::vector<std::string>& waveNames,
                                                    boost::multi_array<std::complex<double>, 2>& productionAmplitudes,
                                                    boost::multi_array<TMatrixT<double>, 1>& productionAmplitudesCovariance,
                                                    boost::multi_array<std::complex<double>, 3>& spinDensityMatrices,
                                                    boost::multi_array<TMatrixT<double>, 1>& spinDensityCovarianceMatrices,
                                                    boost::multi_array<double, 3>& intensities,
                                                    boost::multi_array<double, 4>& phases) const
{
	if(not tree or not fit) {
		printErr << "'tree' or 'fit' is not a pointer to a valid object." << std::endl;
		return false;
	}

	if(_debug) {
		printDebug << "reading spin-density matrices for " << _nrWaves << " waves from fit result." << std::endl;
	}

	// read wave names from first fit result in tree
	waveNames.resize(_nrWaves);
	if(tree->GetEntry(0) == 0) {
		printErr << "error while reading entry " << 0 << " from tree." << std::endl;
		return false;
	}
	for(size_t idxWave = 0; idxWave < _nrWaves; ++idxWave) {
		int idx = fit->waveIndex(_waveNames[idxWave]);
		// try alternative wave names
		for(size_t idxAlt = 0; idxAlt < _waveNameAlternatives[idxWave].size(); ++idxAlt) {
			const int altIdx = fit->waveIndex(_waveNameAlternatives[idxWave][idxAlt]);
			if(altIdx != -1) {
				if(idx != -1) {
					printErr << "more than one wave name or alternative wave name is matching wave in fit result for wave '" << _waveNames[idxWave] << "'." << std::endl;
					return false;
				}
				idx = altIdx;
			}
		}
		if(idx == -1) {
			printErr << "wave '" << _waveNames[idxWave] << "' not in fit result." << std::endl;
			return false;
		}
		waveNames[idxWave] = fit->waveName(idx);
	}

	productionAmplitudes.resize(boost::extents[_maxMassBins][_nrWaves]);
	productionAmplitudesCovariance.resize(boost::extents[_maxMassBins]);

	spinDensityMatrices.resize(boost::extents[_maxMassBins][_nrWaves][_nrWaves]);
	spinDensityCovarianceMatrices.resize(boost::extents[_maxMassBins]);

	intensities.resize(boost::extents[_maxMassBins][_nrWaves][2]);
	phases.resize(boost::extents[_maxMassBins][_nrWaves][_nrWaves][2]);

	for(size_t idxMass = 0; idxMass < mapping.size(); ++idxMass) {
		if(_debug) {
			printDebug << "reading entry " << mapping[idxMass] << " for mass bin " << idxMass << " from tree." << std::endl;
		}
		// FIXME: in case of reading the fit result for a systematic tree this might happen, so this should be allowed in certain cases
		if(tree->GetEntry(mapping[idxMass]) == 0) {
			printErr << "error while reading entry " << mapping[idxMass] << " from tree." << std::endl;
			return false;
		}

		spinDensityCovarianceMatrices[idxMass].ResizeTo(_nrWaves * (_nrWaves+1), _nrWaves * (_nrWaves+1));
		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			const int idx = fit->waveIndex(waveNames[idxWave]);
			if(idx == -1) {
				printErr << "wave '" << _waveNames[idxWave] << "' not in fit result." << std::endl;
				return false;
			}

			intensities[idxMass][idxWave][0] = fit->intensity(idx);
			intensities[idxMass][idxWave][1] = fit->intensityErr(idx) * sqrt(rescaleErrors);

			for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
				const int jdx = fit->waveIndex(waveNames[jdxWave]);
				if(jdx == -1) {
					printErr << "wave '" << _waveNames[jdxWave] << "' not in fit result." << std::endl;
					return false;
				}

				phases[idxMass][idxWave][jdxWave][0] = fit->phase(idx, jdx);
				phases[idxMass][idxWave][jdxWave][1] = fit->phaseErr(idx, jdx) * sqrt(rescaleErrors);

				spinDensityMatrices[idxMass][idxWave][jdxWave] = fit->spinDensityMatrixElem(idx, jdx);

				if(jdxWave >= idxWave) {
					const TMatrixT<double> spinDensityMatrixElemCov = fit->spinDensityMatrixElemCov(idx, jdx) * rescaleErrors;

					const size_t idxCov = _nrWaves*(_nrWaves+1) - (_nrWaves-idxWave)*(_nrWaves-idxWave+1) + 2*(jdxWave-idxWave);
					spinDensityCovarianceMatrices[idxMass].SetSub(idxCov, idxCov, spinDensityMatrixElemCov);
				}
			}
		}

		// for the production amplitudes loop over the production
		// amplitudes of the fit result
		std::vector<unsigned int> prodAmpIndicesForCov(_nrWaves);
		for(unsigned int idxProdAmp=0; idxProdAmp < fit->nmbProdAmps(); ++idxProdAmp) {
			const std::string waveName = fit->waveNameForProdAmp(idxProdAmp);

			const std::map<std::string, size_t>::const_iterator it = _waveIndices.find(waveName);
			// most of the waves are ignored
			if(it == _waveIndices.end()) {
				continue;
			}
			size_t idxWave = it->second;

			int rank = fit->rankOfProdAmp(idxProdAmp);
			// TODO: multiple ranks, in that case also check that rank is not -1
			if(rank != 0) {
				printErr << "can only handle rank-1 fit (production amplitude '" << fit->prodAmpName(idxProdAmp)
				         << "' of wave '" << waveName << "' has rank " << rank << ")." << std::endl;
				return false;
			}

			productionAmplitudes[idxMass][idxWave] = fit->prodAmp(idxProdAmp);

			prodAmpIndicesForCov[idxWave] = idxProdAmp;
		}
		const TMatrixT<double> prodAmpCov = fit->prodAmpCov(prodAmpIndicesForCov) * rescaleErrors;
		productionAmplitudesCovariance[idxMass].ResizeTo(prodAmpCov);
		productionAmplitudesCovariance[idxMass] = prodAmpCov;

		if(_debug) {
			std::ostringstream outputProdAmp;
			std::ostringstream outputProdAmpCovariance;
			std::ostringstream output;
			std::ostringstream outputCovariance;

			outputProdAmp << " (";
			for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
				outputProdAmp << " " << productionAmplitudes[idxMass][idxWave];

				outputProdAmpCovariance << " (";
				output << " (";
				for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
					output << " " << spinDensityMatrices[idxMass][idxWave][jdxWave];

					outputProdAmpCovariance << " (";
					for(size_t idx=0; idx<2; ++idx) {
						outputProdAmpCovariance << " (";
						for(size_t jdx=0; jdx<2; ++jdx) {
							outputProdAmpCovariance << " " << productionAmplitudesCovariance[idxMass](idxWave+idx, jdxWave+jdx);
						}
						outputProdAmpCovariance << " )";
					}
					outputProdAmpCovariance << " )";

					if(jdxWave >= idxWave) {
						const size_t idxCov = _nrWaves*(_nrWaves+1) - (_nrWaves-idxWave)*(_nrWaves-idxWave+1) + 2*(jdxWave-idxWave);
						outputCovariance << " (";
						for(size_t idx = 0; idx < 2; ++idx) {
							outputCovariance << " (";
							for(size_t jdx = 0; jdx < 2; ++jdx) {
								outputCovariance << " " << spinDensityCovarianceMatrices[idxMass](idxCov+idx, idxCov+jdx);
							}
							outputCovariance << " )";
						}
						outputCovariance << " )";
					}
				}
				outputProdAmpCovariance << " )";
				output << " )";
			}
			outputProdAmp << " )";

			printDebug << "production amplitudes: " << outputProdAmp.str() << std::endl;
			printDebug << "production amplitudes covariances: " << outputProdAmpCovariance.str() << std::endl;
			printDebug << "spin-density matrix: " << output.str() << std::endl;
			printDebug << "spin-density covariance matrix (diagonal elements): " << outputCovariance.str() << std::endl;
		}
	} // end loop over mass bins

	return true;
}


bool
rpwa::massDepFit::massDepFit::readFitResultIntegrals(TTree* tree,
                                                     rpwa::fitResult* fit,
                                                     const std::vector<Long64_t>& mapping,
                                                     const std::vector<std::string>& waveNames,
                                                     boost::multi_array<double, 2>& phaseSpaceIntegrals) const
{
	if(not tree or not fit) {
		printErr << "'tree' or 'fit' is not a pointer to a valid object." << std::endl;
		return false;
	}

	phaseSpaceIntegrals.resize(boost::extents[_maxMassBins][_nrWaves]);

	if(_debug) {
		printDebug << "reading phase-space integrals for " << _nrWaves << " waves from fit result." << std::endl;
	}

	for(size_t idxMass = 0; idxMass < mapping.size(); ++idxMass) {
		if(_debug) {
			printDebug << "reading entry " << mapping[idxMass] << " for mass bin " << idxMass << " from tree." << std::endl;
		}
		if(tree->GetEntry(mapping[idxMass]) == 0) {
			printErr << "error while reading entry " << mapping[idxMass] << " from tree." << std::endl;
			return false;
		}

		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			const double ps = fit->phaseSpaceIntegral(waveNames[idxWave]);
			phaseSpaceIntegrals[idxMass][idxWave] = ps;
		}
	}

	if(_debug) {
		for(size_t idxWave = 0; idxWave < _nrWaves; ++idxWave) {
			std::ostringstream output;
			for(size_t idxMass = 0; idxMass < mapping.size(); ++idxMass) {
				output << " " << phaseSpaceIntegrals[idxMass][idxWave];
			}
			printDebug << "phase-space integrals for wave '" << _waveNames[idxWave] << "' (" << idxWave << "):" << output.str() << std::endl;
		}
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::prepareMassLimits()
{
	_waveMassBinLimits.resize(boost::extents[_nrBins][_nrWaves]);
	_wavePairMassBinLimits.resize(boost::extents[_nrBins][_nrWaves][_nrWaves]);

	for(size_t idxBin = 0; idxBin < _nrBins; ++idxBin) {
		if(not prepareMassLimit(idxBin)) {
			printErr << "error while determine bins to use in fit for bin " << idxBin << "." << std::endl;
			return false;
		}
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::prepareMassLimit(const size_t idxBin)
{
	if(_debug) {
		printDebug << "determine which mass bins to use in the fit for " << _nrMassBins[idxBin] << " mass bins, center of first and last mass bins: "
		           << _massBinCenters[idxBin][0] << " and " << _massBinCenters[idxBin][_nrMassBins[idxBin] - 1] << " GeV/c^2." << std::endl;
	}

	for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
		size_t binFirst = 0;
		size_t binLast = _nrMassBins[idxBin]-1;
		for(size_t idxMass = 0; idxMass < _nrMassBins[idxBin]; ++idxMass) {
			if(_massBinCenters[idxBin][idxMass] < _waveMassLimits[idxWave].first) {
				binFirst = idxMass+1;
			}
			if(_massBinCenters[idxBin][idxMass] == _waveMassLimits[idxWave].first) {
				binFirst = idxMass;
			}
			if(_massBinCenters[idxBin][idxMass] <= _waveMassLimits[idxWave].second) {
				binLast = idxMass;
			}
		}
		if(_waveMassLimits[idxWave].first < 0) {
			binFirst = 0;
		}
		if(_waveMassLimits[idxWave].second < 0) {
			binLast = _nrMassBins[idxBin]-1;
		}
		if(_debug) {
			printDebug << idxWave << ": " << _waveNames[idxWave] << ": "
			           << "mass range: " << (_waveMassLimits[idxWave].first<0. ? _massMins[idxBin] : _waveMassLimits[idxWave].first)
			           << "-" << (_waveMassLimits[idxWave].second<0. ? _massMaxs[idxBin] : _waveMassLimits[idxWave].second) << " GeV/c^2, "
			           << "bin range " << binFirst << "-" << binLast << std::endl;
		}
		_waveMassBinLimits[idxBin][idxWave] = std::make_pair(binFirst, binLast);
	}

	for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
		for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
			_wavePairMassBinLimits[idxBin][idxWave][jdxWave] = std::make_pair(std::max(_waveMassBinLimits[idxBin][idxWave].first,  _waveMassBinLimits[idxBin][jdxWave].first),
			                                                                  std::min(_waveMassBinLimits[idxBin][idxWave].second, _waveMassBinLimits[idxBin][jdxWave].second));
		}
	}

	if(_debug) {
		printDebug << "waves and mass limits:" << std::endl;
		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			std::ostringstream output;
			for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
				output << _wavePairMassBinLimits[idxBin][idxWave][jdxWave].first << "-" << _wavePairMassBinLimits[idxBin][idxWave][jdxWave].second << " ";
			}
			printDebug << _waveNames[idxWave] << " " << _waveMassBinLimits[idxBin][idxWave].first << "-" << _waveMassBinLimits[idxBin][idxWave].second
			           << ": " << output.str() << std::endl;
		}
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::createPlots(const rpwa::massDepFit::model& fitModel,
                                          const rpwa::massDepFit::parameters& fitParameters,
                                          rpwa::massDepFit::cache& cache,
                                          TFile* outFile,
                                          const bool rangePlotting,
                                          const size_t extraBinning) const
{
	if(_debug) {
		printDebug << "start creating plots." << std::endl;
	}

	for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
		TDirectory* outDirectory = NULL;
		if(_nrBins == 1) {
			outDirectory = outFile;
		} else {
			std::ostringstream name;
			name << "bin" << idxBin;
			outDirectory = outFile->mkdir(name.str().c_str());
		}

		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			if(not createPlotsWave(fitModel, fitParameters, cache, outDirectory, rangePlotting, extraBinning, idxWave, idxBin)) {
				printErr << "error while creating intensity plots for wave '" << _waveNames[idxWave] << "' in bin " << idxBin << "." << std::endl;
				return false;
			}
		}

		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			for(size_t jdxWave=idxWave+1; jdxWave<_nrWaves; ++jdxWave) {
				if(not createPlotsWavePair(fitModel, fitParameters, cache, outDirectory, rangePlotting, extraBinning, idxWave, jdxWave, idxBin)) {
					printErr << "error while creating intensity plots for wave pair '" << _waveNames[idxWave] << "' and '" << _waveNames[jdxWave] << "' in bin " << idxBin << "." << std::endl;
					return false;
				}
			}
		}

		if(fitModel.getFsmd() and (fitModel.getFsmd()->getNrBins() != 1 or not _sameMassBinning)) {
			if(not createPlotsFsmd(fitModel, fitParameters, cache, outDirectory, rangePlotting, extraBinning, idxBin)) {
				printErr << "error while creating plots for final-state mass-dependence in bin " << idxBin << "." << std::endl;
				return false;
			}
		}
	}

	if(_nrBins != 1 and _sameMassBinning and fitModel.isMappingEqualInAllBins()) {
		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			if(not createPlotsWaveSum(fitModel, fitParameters, cache, outFile, rangePlotting, extraBinning, idxWave)) {
				printErr << "error while creating intensity plots for wave '" << _waveNames[idxWave] << "' for sum over all bins." << std::endl;
				return false;
			}
		}
	}

	if(fitModel.getFsmd() and (fitModel.getFsmd()->getNrBins() == 1 and _sameMassBinning)) {
		if(not createPlotsFsmd(fitModel, fitParameters, cache, outFile, rangePlotting, extraBinning, 0)) {
			printErr << "error while creating plots for final-state mass-dependence." << std::endl;
			return false;
		}
	}

	if(_debug) {
		printDebug << "finished creating plots." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::createPlotsWave(const rpwa::massDepFit::model& fitModel,
                                              const rpwa::massDepFit::parameters& fitParameters,
                                              rpwa::massDepFit::cache& cache,
                                              TDirectory* outDirectory,
                                              const bool rangePlotting,
                                              const size_t extraBinning,
                                              const size_t idxWave,
                                              const size_t idxBin) const
{
	if(_debug) {
		printDebug << "start creating plots for wave '" << _waveNames[idxWave] << "' in bin " << idxBin << "." << std::endl;
	}

	TMultiGraph graphs;
	graphs.SetName(_waveNames[idxWave].c_str());
	graphs.SetTitle(_waveNames[idxWave].c_str());

	TGraphErrors* systematics = NULL;
	if(_sysPlotting) {
		systematics = new TGraphErrors;
		systematics->SetName((_waveNames[idxWave] + "__sys").c_str());
		systematics->SetTitle((_waveNames[idxWave] + "__sys").c_str());
		systematics->SetLineColor(kAzure-9);
		systematics->SetFillColor(kAzure-9);
		graphs.Add(systematics, "2");
	}

	TGraphErrors* data = new TGraphErrors;
	data->SetName((_waveNames[idxWave] + "__data").c_str());
	data->SetTitle((_waveNames[idxWave] + "__data").c_str());
	graphs.Add(data, "P");

	TGraph* fit = new TGraph;
	fit->SetName((_waveNames[idxWave] + "__fit").c_str());
	fit->SetTitle((_waveNames[idxWave] + "__fit").c_str());
	fit->SetLineColor(kRed);
	fit->SetLineWidth(2);
	fit->SetMarkerColor(kRed);
	graphs.Add(fit, "L");

	TGraph* phaseSpace = new TGraph;
	phaseSpace->SetName((_waveNames[idxWave] + "__ps").c_str());
	phaseSpace->SetTitle((_waveNames[idxWave] + "__ps").c_str());
	graphs.Add(phaseSpace, "L");

	const std::vector<std::pair<size_t, size_t> >& compChannel = fitModel.getComponentChannel(idxBin, idxWave);
	std::vector<TGraph*> components;
	for(size_t idxComponents=0; idxComponents<compChannel.size(); ++idxComponents) {
		const size_t idxComponent = compChannel[idxComponents].first;
		TGraph* component = new TGraph;
		component->SetName((_waveNames[idxWave] + "__" + fitModel.getComponent(idxComponent)->getName()).c_str());
		component->SetTitle((_waveNames[idxWave] + "__" + fitModel.getComponent(idxComponent)->getName()).c_str());

		Color_t color = kBlue;
		if(fitModel.getComponent(idxComponent)->getName().find("bkg") != std::string::npos) {
			color = kMagenta;
		}
		component->SetLineColor(color);
		component->SetMarkerColor(color);

		graphs.Add(component, "L");
		components.push_back(component);
	}

	// plot data
	double maxIE = -std::numeric_limits<double>::max();
	for(size_t point = 0; point <= (_nrMassBins[idxBin]-1); ++point) {
		const size_t idxMass = point;
		const double mass = _massBinCenters[idxBin][idxMass];
		const double halfBin = _massSteps[idxBin]/2.;

		data->SetPoint(point, mass, _inIntensities[idxBin][idxMass][idxWave][0]);
		data->SetPointError(point, halfBin, _inIntensities[idxBin][idxMass][idxWave][1]);
		maxIE = std::max(maxIE, _inIntensities[idxBin][idxMass][idxWave][0]+_inIntensities[idxBin][idxMass][idxWave][1]);

		if(_sysPlotting) {
			double maxSI = -std::numeric_limits<double>::max();
			double minSI = std::numeric_limits<double>::max();
			for(size_t idxSystematics=0; idxSystematics<_nrSystematics; ++idxSystematics) {
				maxSI = std::max(maxSI, _sysIntensities[idxSystematics][idxMass][idxWave][0]);
				minSI = std::min(minSI, _sysIntensities[idxSystematics][idxMass][idxWave][0]);
			}
			systematics->SetPoint(point, mass, (maxSI+minSI)/2.);
			systematics->SetPointError(point, halfBin, (maxSI-minSI)/2.);
			maxIE = std::max(maxIE, maxSI);
		}
	}

	// plot fit, either over full or limited mass range
	const size_t firstPoint = rangePlotting ? (extraBinning*_waveMassBinLimits[idxBin][idxWave].first) : 0;
	const size_t lastPoint = rangePlotting ? (extraBinning*_waveMassBinLimits[idxBin][idxWave].second) : (extraBinning*(_nrMassBins[idxBin]-1));
	for(size_t point=firstPoint; point<=lastPoint; ++point) {
		const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
		const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? _massBinCenters[idxBin][idxMass] : (_massBinCenters[idxBin][point/extraBinning] + (point%extraBinning) * _massSteps[idxBin]/extraBinning);

		const double intensity = fitModel.intensity(fitParameters, cache, idxWave, idxBin, mass, idxMass);
		fit->SetPoint(point-firstPoint, mass, intensity);
		maxIE = std::max(maxIE, intensity);

		for(size_t idxComponents=0; idxComponents<compChannel.size(); ++idxComponents) {
			const size_t idxComponent = compChannel[idxComponents].first;
			const size_t idxChannel = compChannel[idxComponents].second;

			std::complex<double> prodAmp = fitModel.getComponent(idxComponent)->val(fitParameters, cache, idxBin, mass, idxMass);
			prodAmp *= fitModel.getComponent(idxComponent)->getCouplingPhaseSpace(fitParameters, cache, idxChannel, idxBin, mass, idxMass);
			if(fitModel.getFsmd()) {
				prodAmp *= fitModel.getFsmd()->val(fitParameters, cache, idxBin, mass, idxMass);
			}

			components[idxComponents]->SetPoint(point-firstPoint, mass, norm(prodAmp));
			maxIE = std::max(maxIE, norm(prodAmp));
		}
	}

	boost::multi_array<double, 2>::const_array_view<1>::type viewM = _massBinCenters[boost::indices[idxBin][boost::multi_array<double, 2>::index_range(0, _nrMassBins[idxBin])]];
	boost::multi_array<double, 3>::const_array_view<1>::type viewInt = _inPhaseSpaceIntegrals[boost::indices[idxBin][boost::multi_array<double, 3>::index_range(0, _nrMassBins[idxBin])][idxWave]];
	ROOT::Math::Interpolator phaseSpaceInterpolator(std::vector<double>(viewM.begin(), viewM.end()), std::vector<double>(viewInt.begin(), viewInt.end()), ROOT::Math::Interpolation::kLINEAR);

	// plot phase-space
	double maxP = -std::numeric_limits<double>::max();
	for(size_t point = 0; point <= (extraBinning*(_nrMassBins[idxBin]-1)); ++point) {
		const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
		const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? _massBinCenters[idxBin][idxMass] : (_massBinCenters[idxBin][point/extraBinning] + (point%extraBinning) * _massSteps[idxBin]/extraBinning);

		double ps = pow((idxMass != std::numeric_limits<size_t>::max()) ? _inPhaseSpaceIntegrals[idxBin][idxMass][idxWave] : phaseSpaceInterpolator.Eval(mass), 2);
		if(fitModel.getFsmd()) {
			ps *= std::norm(fitModel.getFsmd()->val(fitParameters, cache, idxBin, mass, idxMass));
		}
		phaseSpace->SetPoint(point, mass, ps);
		maxP = std::max(maxP, ps);
	}

	// scale phase-space graph to half-height of intensity graphs
	for(Int_t idx=0; idx<phaseSpace->GetN(); ++idx) {
		double x, y;
		phaseSpace->GetPoint(idx, x, y);
		phaseSpace->SetPoint(idx, x, y * 0.5 * maxIE/maxP);
	}

	outDirectory->cd();
	graphs.Write();

	return true;
}


bool
rpwa::massDepFit::massDepFit::createPlotsWaveSum(const rpwa::massDepFit::model& fitModel,
                                                 const rpwa::massDepFit::parameters& fitParameters,
                                                 rpwa::massDepFit::cache& cache,
                                                 TDirectory* outDirectory,
                                                 const bool rangePlotting,
                                                 const size_t extraBinning,
                                                 const size_t idxWave) const
{
	if(_debug) {
		printDebug << "start creating plots for wave '" << _waveNames[idxWave] << "' for sum over all bins." << std::endl;
	}

	// all mass binnings must be the same to be able to create the sum plots
	if(not _sameMassBinning or not fitModel.isMappingEqualInAllBins()) {
		printErr << "cannot create plots for wave '" << _waveNames[idxWave] << "' for sum over all bins if the bins used different mass binnings." << std::endl;
		return false;
	}
	const size_t idxBin = 0;

	TMultiGraph graphs;
	graphs.SetName(_waveNames[idxWave].c_str());
	graphs.SetTitle(_waveNames[idxWave].c_str());

	TGraphErrors* data = new TGraphErrors;
	data->SetName((_waveNames[idxWave] + "__data").c_str());
	data->SetTitle((_waveNames[idxWave] + "__data").c_str());
	graphs.Add(data, "P");

	TGraph* fit = new TGraph;
	fit->SetName((_waveNames[idxWave] + "__fit").c_str());
	fit->SetTitle((_waveNames[idxWave] + "__fit").c_str());
	fit->SetLineColor(kRed);
	fit->SetLineWidth(2);
	fit->SetMarkerColor(kRed);
	graphs.Add(fit, "L");

	const std::vector<std::pair<size_t, size_t> >& compChannel = fitModel.getComponentChannel(idxBin, idxWave);
	std::vector<TGraph*> components;
	for(size_t idxComponents=0; idxComponents<compChannel.size(); ++idxComponents) {
		const size_t idxComponent = compChannel[idxComponents].first;
		TGraph* component = new TGraph;
		component->SetName((_waveNames[idxWave] + "__" + fitModel.getComponent(idxComponent)->getName()).c_str());
		component->SetTitle((_waveNames[idxWave] + "__" + fitModel.getComponent(idxComponent)->getName()).c_str());

		Color_t color = kBlue;
		if(fitModel.getComponent(idxComponent)->getName().find("bkg") != std::string::npos) {
			color = kMagenta;
		}
		component->SetLineColor(color);
		component->SetMarkerColor(color);

		graphs.Add(component, "L");
		components.push_back(component);
	}

	// plot data
	for(size_t point = 0; point <= (_nrMassBins[idxBin]-1); ++point) {
		const size_t idxMass = point;
		const double mass = _massBinCenters[idxBin][idxMass];
		const double halfBin = _massSteps[idxBin]/2.;

		double sum = 0.;
		double error2 = 0.;
		for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
			sum += _inIntensities[idxBin][idxMass][idxWave][0];
			error2 += std::pow(_inIntensities[idxBin][idxMass][idxWave][1], 2);
		}
		data->SetPoint(point, mass, sum);
		data->SetPointError(point, halfBin, sqrt(error2));
	}

	// plot fit, either over full or limited mass range
	const size_t firstPoint = rangePlotting ? (extraBinning*_waveMassBinLimits[idxBin][idxWave].first) : 0;
	const size_t lastPoint = rangePlotting ? (extraBinning*_waveMassBinLimits[idxBin][idxWave].second) : (extraBinning*(_nrMassBins[idxBin]-1));
	for(size_t point=firstPoint; point<=lastPoint; ++point) {
		const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
		const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? _massBinCenters[idxBin][idxMass] : (_massBinCenters[idxBin][point/extraBinning] + (point%extraBinning) * _massSteps[idxBin]/extraBinning);

		double sum = 0.;
		for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
			sum += fitModel.intensity(fitParameters, cache, idxWave, idxBin, mass, idxMass);
		}
		fit->SetPoint(point-firstPoint, mass, sum);

		for(size_t idxComponents=0; idxComponents<compChannel.size(); ++idxComponents) {
			const size_t idxComponent = compChannel[idxComponents].first;
			const size_t idxChannel = compChannel[idxComponents].second;

			double sum = 0.;
			for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
				std::complex<double> prodAmp = fitModel.getComponent(idxComponent)->val(fitParameters, cache, idxBin, mass, idxMass);
				prodAmp *= fitModel.getComponent(idxComponent)->getCouplingPhaseSpace(fitParameters, cache, idxChannel, idxBin, mass, idxMass);
				if(fitModel.getFsmd()) {
					prodAmp *= fitModel.getFsmd()->val(fitParameters, cache, idxBin, mass, idxMass);
				}
				sum += norm(prodAmp);
			}
			components[idxComponents]->SetPoint(point-firstPoint, mass, sum);
		}
	}

	outDirectory->cd();
	graphs.Write();

	return true;
}


bool
rpwa::massDepFit::massDepFit::createPlotsWavePair(const rpwa::massDepFit::model& fitModel,
                                                  const rpwa::massDepFit::parameters& fitParameters,
                                                  rpwa::massDepFit::cache& cache,
                                                  TDirectory* outDirectory,
                                                  const bool rangePlotting,
                                                  const size_t extraBinning,
                                                  const size_t idxWave,
                                                  const size_t jdxWave,
                                                  const size_t idxBin) const
{
	if(_debug) {
		printDebug << "start creating plots for wave pair '" << _waveNames[idxWave] << "' and '" << _waveNames[jdxWave] << "' in bin " << idxBin << "." << std::endl;
	}

	const std::string phaseName = _waveNames[idxWave] + "__" + _waveNames[jdxWave] + "__phase";
	const std::string realName = _waveNames[idxWave] + "__" + _waveNames[jdxWave] + "__real";
	const std::string imagName = _waveNames[idxWave] + "__" + _waveNames[jdxWave] + "__imag";

	TMultiGraph phase;
	phase.SetName(phaseName.c_str());
	phase.SetTitle(phaseName.c_str());

	TMultiGraph real;
	real.SetName(realName.c_str());
	real.SetTitle(realName.c_str());

	TMultiGraph imag;
	imag.SetName(imagName.c_str());
	imag.SetTitle(imagName.c_str());

	TGraphErrors* phaseSystematics = NULL;
	TGraphErrors* realSystematics = NULL;
	TGraphErrors* imagSystematics = NULL;
	if(_sysPlotting) {
		phaseSystematics = new TGraphErrors;
		phaseSystematics->SetName((phaseName + "__sys").c_str());
		phaseSystematics->SetTitle((phaseName + "__sys").c_str());
		phaseSystematics->SetLineColor(kAzure-9);
		phaseSystematics->SetFillColor(kAzure-9);
		phase.Add(phaseSystematics, "2");

		realSystematics = new TGraphErrors;
		realSystematics->SetName((realName + "__sys").c_str());
		realSystematics->SetTitle((realName + "__sys").c_str());
		realSystematics->SetLineColor(kAzure-9);
		realSystematics->SetFillColor(kAzure-9);
		real.Add(realSystematics, "2");

		imagSystematics = new TGraphErrors;
		imagSystematics->SetName((imagName + "__sys").c_str());
		imagSystematics->SetTitle((imagName + "__sys").c_str());
		imagSystematics->SetLineColor(kAzure-9);
		imagSystematics->SetFillColor(kAzure-9);
		imag.Add(imagSystematics, "2");
	}

	TGraphErrors* phaseData = new TGraphErrors;
	phaseData->SetName((phaseName + "__data").c_str());
	phaseData->SetTitle((phaseName + "__data").c_str());
	phase.Add(phaseData, "P");

	TGraphErrors* realData = new TGraphErrors;
	realData->SetName((realName + "__data").c_str());
	realData->SetTitle((realName + "__data").c_str());
	real.Add(realData, "P");

	TGraphErrors* imagData = new TGraphErrors;
	imagData->SetName((imagName + "__data").c_str());
	imagData->SetTitle((imagName + "__data").c_str());
	imag.Add(imagData, "P");

	TGraph* phaseFit = new TGraph;
	phaseFit->SetName((phaseName + "__fit").c_str());
	phaseFit->SetTitle((phaseName + "__fit").c_str());
	phaseFit->SetLineColor(kRed);
	phaseFit->SetLineWidth(2);
	phaseFit->SetMarkerColor(kRed);
	phase.Add(phaseFit, "L");

	TGraph* realFit = new TGraph;
	realFit->SetName((realName + "__fit").c_str());
	realFit->SetTitle((realName + "__fit").c_str());
	realFit->SetLineColor(kRed);
	realFit->SetLineWidth(2);
	realFit->SetMarkerColor(kRed);
	real.Add(realFit, "L");

	TGraph* imagFit = new TGraph;
	imagFit->SetName((imagName + "__fit").c_str());
	imagFit->SetTitle((imagName + "__fit").c_str());
	imagFit->SetLineColor(kRed);
	imagFit->SetLineWidth(2);
	imagFit->SetMarkerColor(kRed);
	imag.Add(imagFit, "L");

	// plot data
	for(size_t point = 0; point <= (_nrMassBins[idxBin]-1); ++point) {
		const size_t idxMass = point;
		const double mass = _massBinCenters[idxBin][idxMass];
		const double halfBin = _massSteps[idxBin]/2.;

		phaseData->SetPoint(point, mass, _inPhases[idxBin][idxMass][idxWave][jdxWave][0]);
		phaseData->SetPointError(point, halfBin, _inPhases[idxBin][idxMass][idxWave][jdxWave][1]);

		const size_t idxCov = _nrWaves*(_nrWaves+1) - ((jdxWave >= idxWave) ? ((_nrWaves-idxWave)*(_nrWaves-idxWave+1) - 2*(jdxWave-idxWave)) : ((_nrWaves-jdxWave)*(_nrWaves-jdxWave+1) - 2*(idxWave-jdxWave)));

		realData->SetPoint(point, mass, _inSpinDensityMatrices[idxBin][idxMass][idxWave][jdxWave].real());
		realData->SetPointError(point, halfBin, sqrt(_inSpinDensityCovarianceMatrices[idxBin][idxMass](idxCov, idxCov)));

		imagData->SetPoint(point, mass, _inSpinDensityMatrices[idxBin][idxMass][idxWave][jdxWave].imag());
		imagData->SetPointError(point, halfBin, sqrt(_inSpinDensityCovarianceMatrices[idxBin][idxMass](idxCov+1, idxCov+1)));

		if(_sysPlotting) {
			const double dataP = _inPhases[idxBin][idxMass][idxWave][jdxWave][0];
			double maxSP = -std::numeric_limits<double>::max();
			double minSP = std::numeric_limits<double>::max();

			double maxSR = -std::numeric_limits<double>::max();
			double minSR = std::numeric_limits<double>::max();

			double maxSI = -std::numeric_limits<double>::max();
			double minSI = std::numeric_limits<double>::max();

			for(size_t idxSystematics=0; idxSystematics<_nrSystematics; ++idxSystematics) {
				double sysP = _sysPhases[idxSystematics][idxMass][idxWave][jdxWave][0];
				if(std::abs(sysP+360.-dataP) < std::abs(sysP-dataP)) {
					sysP = sysP+360;
				} else if(std::abs(sysP-360.-dataP) < std::abs(sysP-dataP)) {
					sysP = sysP-360;
				}
				maxSP = std::max(maxSP, sysP);
				minSP = std::min(minSP, sysP);

				maxSR = std::max(maxSR, _sysSpinDensityMatrices[idxSystematics][idxMass][idxWave][jdxWave].real());
				minSR = std::min(minSR, _sysSpinDensityMatrices[idxSystematics][idxMass][idxWave][jdxWave].real());

				maxSI = std::max(maxSI, _sysSpinDensityMatrices[idxSystematics][idxMass][idxWave][jdxWave].imag());
				minSI = std::min(minSI, _sysSpinDensityMatrices[idxSystematics][idxMass][idxWave][jdxWave].imag());
			}
			phaseSystematics->SetPoint(point, mass, (maxSP+minSP)/2.);
			phaseSystematics->SetPointError(point, halfBin, (maxSP-minSP)/2.);

			realSystematics->SetPoint(point, mass, (maxSR+minSR)/2.);
			realSystematics->SetPointError(point, halfBin, (maxSR-minSR)/2.);

			imagSystematics->SetPoint(point, mass, (maxSI+minSI)/2.);
			imagSystematics->SetPointError(point, halfBin, (maxSI-minSI)/2.);
		}
	}

	// plot fit, either over full or limited mass range
	const size_t firstPoint = rangePlotting ? (extraBinning*_wavePairMassBinLimits[idxBin][idxWave][jdxWave].first) : 0;
	const size_t lastPoint = rangePlotting ? (extraBinning*_wavePairMassBinLimits[idxBin][idxWave][jdxWave].second) : (extraBinning*(_nrMassBins[idxBin]-1));
	for(size_t point=firstPoint; point<=lastPoint; ++point) {
		const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
		const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? _massBinCenters[idxBin][idxMass] : (_massBinCenters[idxBin][point/extraBinning] + (point%extraBinning) * _massSteps[idxBin]/extraBinning);

		const std::complex<double> element = fitModel.spinDensityMatrix(fitParameters, cache, idxWave, jdxWave, idxBin, mass, idxMass);
		realFit->SetPoint(point-firstPoint, mass, element.real());
		imagFit->SetPoint(point-firstPoint, mass, element.imag());
	}

	// keep track of phase over full mass range
	TGraph phaseFitAll;
	for(size_t point = 0; point <= (extraBinning*(_nrMassBins[idxBin]-1)); ++point) {
		const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
		const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? _massBinCenters[idxBin][idxMass] : (_massBinCenters[idxBin][point/extraBinning] + (point%extraBinning) * _massSteps[idxBin]/extraBinning);

		const double phase = fitModel.phase(fitParameters, cache, idxWave, jdxWave, idxBin, mass, idxMass) * TMath::RadToDeg();

		if(point != 0) {
			int bestOffs = 0;
			double bestDiff = std::numeric_limits<double>::max();

			double x;
			double prev;
			phaseFitAll.GetPoint(point-1, x, prev);
			for(int offs=-5; offs<6; ++offs) {
				if(std::abs(phase + offs*360. - prev) < bestDiff) {
					bestDiff = std::abs(phase + offs*360. - prev);
					bestOffs = offs;
				}
			}

			phaseFitAll.SetPoint(point, mass, phase + bestOffs*360.);
		} else {
			phaseFitAll.SetPoint(point, mass, phase);
		}
	}

	// rectify phase graphs
	for(size_t point = 0; point <= (extraBinning*(_nrMassBins[idxBin]-1)); ++point) {
		const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
		const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? _massBinCenters[idxBin][idxMass] : (_massBinCenters[idxBin][point/extraBinning] + (point%extraBinning) * _massSteps[idxBin]/extraBinning);

		double x;
		double valueFit;
		phaseFitAll.GetPoint(point, x, valueFit);

		if (idxMass != std::numeric_limits<size_t>::max()) {
			int bestOffs = 0;
			double bestDiff = std::numeric_limits<double>::max();

			double data;
			phaseData->GetPoint(idxMass, x, data);
			for(int offs=-5; offs<6; ++offs) {
				if(std::abs(data + offs*360. - valueFit) < bestDiff) {
					bestDiff = std::abs(data + offs*360. - valueFit);
					bestOffs = offs;
				}
			}

			phaseData->SetPoint(idxMass, x, data + bestOffs*360.);
			if(_sysPlotting) {
				phaseSystematics->GetPoint(idxMass, x, data);
				phaseSystematics->SetPoint(idxMass, x, data + bestOffs*360.);
			}
		}

		// check that this mass bin should be taken into account for this
		// combination of waves
		if(point < firstPoint || point > lastPoint) {
			continue;
		}

		phaseFit->SetPoint(point-firstPoint, mass, valueFit);
	}

	outDirectory->cd();
	phase.Write();
	real.Write();
	imag.Write();

	return true;
}


bool
rpwa::massDepFit::massDepFit::createPlotsFsmd(const rpwa::massDepFit::model& fitModel,
                                              const rpwa::massDepFit::parameters& fitParameters,
                                              rpwa::massDepFit::cache& cache,
                                              TDirectory* outDirectory,
                                              const bool /*rangePlotting*/,
                                              const size_t extraBinning,
                                              const size_t idxBin) const
{
	if(_debug) {
		printDebug << "start creating plots for final-state mass-dependence." << std::endl;
	}

	TGraph graph;
	graph.SetName("finalStateMassDependence");
	graph.SetTitle("finalStateMassDependence");

	for(size_t point = 0; point <= (extraBinning*(_nrMassBins[idxBin]-1)); ++point) {
		const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
		const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? _massBinCenters[idxBin][idxMass] : (_massBinCenters[idxBin][point/extraBinning] + (point%extraBinning) * _massSteps[idxBin]/extraBinning);

		graph.SetPoint(point, mass, std::norm(fitModel.getFsmd()->val(fitParameters, cache, idxBin, mass, idxMass)));
	}

	outDirectory->cd();
	graph.Write();

	return true;
}
