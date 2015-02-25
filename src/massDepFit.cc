///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010-2012 Sebastian Neubert (TUM)
//    Copyright 2014,2015 Sebastian Uhl (TUM)
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

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <TTree.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TMatrixD.h>
#include <TMultiGraph.h>
#include <TString.h>
#include <TComplex.h>
#include <TRandom3.h>
#include <TStopwatch.h>

#include <libconfig.h++>

#include "ampIntegralMatrix.h"
#include "fileUtils.hpp"
#include "fitResult.h"
#include "libConfigUtils.hpp"
#include "massDepFitComponents.h"
#include "massDepFitFsmd.h"
#include "massDepFitLikeli.h"
#include "massDepFitModel.h"
#include "reportingUtils.hpp"


bool rpwa::massDepFit::massDepFit::_debug = false;


rpwa::massDepFit::massDepFit::massDepFit()
	: _sysPlotting(false),
	  _nrMassBins(0),
	  _nrSystematics(0),
	  _nrWaves(0)
{
}


bool
rpwa::massDepFit::massDepFit::readConfig(const libconfig::Setting* configRoot,
                                         rpwa::massDepFit::model& fitModel,
                                         rpwa::massDepFit::parameters& fitParameters,
                                         const std::string& valTreeName,
                                         const std::string& valBranchName)
{
	// input section
	const libconfig::Setting* configInput = findLibConfigGroup(*configRoot, "input");
	if(not configInput) {
		printErr << "'input' section in configuration file does not exist." << std::endl;
		return false;
	}
	if(not readConfigInput(configInput)) {
		printErr << "error while reading 'input' section from configuration file." << std::endl;
		return false;
	}

	// extract information from fit results
	if(not readInFiles(valTreeName, valBranchName)) {
		printErr << "error while trying to read fit result." << std::endl;
		return false;
	}

	// extract information for systematic errors
	if(not readSystematicsFiles(valTreeName, valBranchName)) {
		printErr << "error while trying to read fit results for systematic errors." << std::endl;
		return false;
	}

	// prepare mass limits
	if(not prepareMassLimits()) {
		printErr << "error determine which bins to use in the fit." << std::endl;
		return false;
	}

	// set-up fit model (resonances, background, final-state mass dependence
	if(not readConfigModel(configRoot, fitModel, fitParameters)) {
		printErr << "error while reading fit model from configuration file." << std::endl;
		return false;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigInput(const libconfig::Setting* configInput)
{
	if(not configInput) {
		printErr << "'configInput' is not a pointer to a valid object." << std::endl;
		return false;
	}

	// get information about fit results from mass-independent
	const libconfig::Setting* configInputFitResults = findLibConfigList(*configInput, "fitresults");
	if(not configInputFitResults) {
		printErr << "'fitresults' list does not exist in section '" << configInput->getName() << "' in configuration file." << std::endl;
		return false;
	}
	if(not readConfigInputFitResults(configInputFitResults)) {
		printErr << "error while reading 'fitresults' in section '" << configInput->getName() << "' in configuration file." << std::endl;
		return false;
	}

	// get information about waves to be used in the fit
	const libconfig::Setting* configInputWaves = findLibConfigList(*configInput, "waves");
	if(not configInputWaves) {
		printErr << "'waves' list does not exist in section '" << configInput->getName() << "' in configuration file." << std::endl;
		return false;
	}
	if(not readConfigInputWaves(configInputWaves)) {
		printErr << "error while reading 'waves' in section '" << configInput->getName() << "' in configuration file." << std::endl;
		return false;
	}

	// get information for plotting of systematic error
	const libconfig::Setting* configInputSystematics = findLibConfigArray(*configInput, "systematics", false);
	if(not readConfigInputSystematics(configInputSystematics)) {
		printErr << "error while reading 'systematics' in section '" << configInput->getName() << "' in configuration file." << std::endl;
		return false;
	}

	// get information for which parameters to release in which order
	const libconfig::Setting* configInputFreeParameters = findLibConfigArray(*configInput, "freeparameters", false);
	if(not readConfigInputFreeParameters(configInputFreeParameters)) {
		printErr << "error while reading 'freeparameters' in section '" << configInput->getName() << "' in configuration file." << std::endl;
		return false;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigInputFitResults(const libconfig::Setting* configInputFitResults)
{
	if(not configInputFitResults) {
		printErr << "'configInputFitResults' is not a pointer to a valid object." << std::endl;
		return false;
	}

	_inFileName.clear();
	_inOverwritePhaseSpace.clear();

	const int nrFitResults = configInputFitResults->getLength();
	for(int idxFitResult=0; idxFitResult<nrFitResults; ++idxFitResult) {
		if(_debug) {
			printDebug << "reading of entry " << idxFitResult << " in 'fitresults'." << std::endl;
		}

		const libconfig::Setting* configInputFitResult = &((*configInputFitResults)[idxFitResult]);

		std::map<std::string, libconfig::Setting::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("name", libconfig::Setting::TypeString)
		                     ("tPrimeMean", libconfig::Setting::TypeFloat);
		if(not checkIfAllVariablesAreThere(configInputFitResult, mandatoryArguments)) {
			printErr << "'fitresults' list in 'input' section in configuration file contains errors." << std::endl;
			return false;
		}

		std::string fileName;
		configInputFitResult->lookupValue("name", fileName);
		_inFileName.push_back(fileName);

		if(_debug) {
			printDebug << "read file name of fit results of mass-independent fit: '" << fileName << "'." << std::endl;
		}

		double tPrimeMean;
		configInputFitResult->lookupValue("tPrimeMean", tPrimeMean);
		_tPrimeMeans.push_back(tPrimeMean);

		if(_debug) {
			printDebug << "read mean t' value: '" << tPrimeMean << "'." << std::endl;
		}

		std::vector<std::string> overwritePhaseSpace;

		const libconfig::Setting* configOverwrite = findLibConfigArray(*configInputFitResult, "overwritePhaseSpace", false);
		if(configOverwrite) {
			const int nrParts = configOverwrite->getLength();

			if(nrParts > 0 && (*configOverwrite)[0].getType() != libconfig::Setting::TypeString) {
				printErr << "contents of 'overwritePhaseSpace' array in 'input' needs to be strings." << std::endl;
				return false;
			}

			for(int idxPart=0; idxPart<nrParts; ++idxPart) {
				const std::string fileName = (*configOverwrite)[idxPart];
				overwritePhaseSpace.push_back(fileName);
			}

			if(_debug) {
				std::ostringstream output;
				output << "[ '";
				for(size_t idxPart=0; idxPart<overwritePhaseSpace.size(); ++idxPart) {
					if(idxPart > 0) {
						output << "', '";
					}
					output << overwritePhaseSpace[idxPart];
				}
				output << "' ]";
				printDebug << "phase-space integrals will be overwritten with files from this pattern: " << output.str() << " (" << overwritePhaseSpace.size() << " parts)" << std::endl;
			}
		}

		_inOverwritePhaseSpace.push_back(overwritePhaseSpace);
	}

	_nrBins = _inFileName.size();

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigInputWaves(const libconfig::Setting* configInputWaves)
{
	if(not configInputWaves) {
		printErr << "'configInputWaves' is not a pointer to a valid object." << std::endl;
		return false;
	}

	const int nrWaves = configInputWaves->getLength();
	if(_debug) {
		printDebug << "going to read information of " << nrWaves << " waves to be used in the fit." << std::endl;
	}

	for(int idxWave=0; idxWave<nrWaves; ++idxWave) {
		const libconfig::Setting* configInputWave = &((*configInputWaves)[idxWave]);

		std::map<std::string, libconfig::Setting::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("name", libconfig::Setting::TypeString);
		if(not checkIfAllVariablesAreThere(configInputWave, mandatoryArguments)) {
			printErr << "'waves' list in 'input' section in configuration file contains errors." << std::endl;
			return false;
		}

		std::string name;
		configInputWave->lookupValue("name", name);

		double massLower;
		if(not configInputWave->lookupValue("massLower", massLower)) {
			massLower = -1.;
		}
		double massUpper;
		if(not configInputWave->lookupValue("massUpper", massUpper)) {
			massUpper = -1.;
		}

		// check that wave does not yet exist
		if(find(_waveNames.begin(), _waveNames.end(), name) != _waveNames.end()) {
			printErr << "wave '" << name << "' defined twice." << std::endl;
			return false;
		}

		_waveNames.push_back(name);
		_waveIndices[name] = _waveNames.size() - 1;
		_waveMassLimits.push_back(std::make_pair(massLower, massUpper));

		if(_debug) {
			printDebug << idxWave << ": " << name << " (mass range: " << massLower << "-" << massUpper << " GeV/c^2, index: " << _waveIndices[name] << ")" << std::endl;
		}
	}

	_nrWaves = _waveNames.size();
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
rpwa::massDepFit::massDepFit::readConfigInputSystematics(const libconfig::Setting* configInputSystematics)
{
	// configInputSystematics might actually be a NULL pointer, in this
	// systematics is not plotted
	if(not configInputSystematics) {
		_sysPlotting = false;
		return true;
	}

	const int nrSystematics = configInputSystematics->getLength();
	if(_debug) {
		printDebug << "going to read information for " << nrSystematics << " files containing information for systematic errors." << std::endl;
	}

	if(nrSystematics > 0) {
		_sysPlotting = true;
	}

	if(nrSystematics > 0 && (*configInputSystematics)[0].getType() != libconfig::Setting::TypeString) {
		printErr << "contents of 'systematics' array in 'input' needs to be strings." << std::endl;
		return false;
	}

	for(int idxSystematics=0; idxSystematics<nrSystematics; ++idxSystematics) {
		const std::string fileName = (*configInputSystematics)[idxSystematics];
		if(_debug) {
			printDebug << "'" << fileName << "' will be read to get information for systematic errors." << std::endl;
		}
		_sysFileNames.push_back(fileName);
	}

	_nrSystematics = _sysFileNames.size() + 1;
	if(_debug) {
		printDebug << "in total " << _nrSystematics << " files to be read to get information for systematic errors." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigInputFreeParameters(const libconfig::Setting* configInputFreeParameters)
{
	if(not configInputFreeParameters) {
		_freeParameters.clear();
		_freeParameters.push_back("branching");
		_freeParameters.push_back("branching mass m0");
		_freeParameters.push_back("*");
		return true;
	}

	_freeParameters.clear();

	const int nrItems = configInputFreeParameters->getLength();

	if(_debug) {
		printDebug << "going to extract " << nrItems << " items from '" << configInputFreeParameters->getName() << "'." << std::endl;
	}

	for(int idxItem=0; idxItem<nrItems; ++idxItem) {
		const libconfig::Setting& item = (*configInputFreeParameters)[idxItem];

		if(item.getType() != libconfig::Setting::TypeString) {
			printErr << "'" << configInputFreeParameters->getName() << "' must be array of strings." << std::endl;
			return false;
		}

		const std::string name = item;

		if(_debug) {
			printDebug << idxItem << ": '" << name << "'." << std::endl;
		}

		_freeParameters.push_back(name);
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigModel(const libconfig::Setting* configRoot,
                                              rpwa::massDepFit::model& fitModel,
                                              rpwa::massDepFit::parameters& fitParameters)
{
	if(_debug) {
		printDebug << "reading fit model from configuration file." << std::endl;
	}

	// get components section of configuration file
	const libconfig::Setting* configModel = findLibConfigGroup(*configRoot, "model");
	if(not configModel) {
		printErr << "error while reading 'model' section in configuration file." << std::endl;
		return false;
	}

	// get information about anchor wave
	const libconfig::Setting* configAnchorWave = findLibConfigGroup(*configModel, "anchorwave");
	if(not readConfigModelAnchorWave(configAnchorWave)) {
		printErr << "error while reading 'anchorwave' in section '" << configModel->getName() << "' in configuration file." << std::endl;
		return false;
	}

	// read information for the individual components
	const libconfig::Setting* configComponents = findLibConfigList(*configModel, "components");
	if(not readConfigModelComponents(configComponents, fitModel, fitParameters)) {
		printErr << "error while reading 'components' in section '" << configModel->getName() << "' in configuration file." << std::endl;
		return false;
	}

	// get information for creating the final-state mass-dependence
	const libconfig::Setting* configFsmd = findLibConfigGroup(*configModel, "finalStateMassDependence", false);
	if(not readConfigModelFsmd(configFsmd, fitModel, fitParameters)) {
		printErr << "error while reading 'finalStateMassDependence' in section '" << configModel->getName() << "' in configuration file." << std::endl;
		return false;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigModelAnchorWave(const libconfig::Setting* configAnchorWave)
{
	if(not configAnchorWave) {
		printErr << "'configInputAnchorWave' is not a pointer to a valid object." << std::endl;
		return false;
	}

	std::map<std::string, libconfig::Setting::Type> mandatoryArguments;
	boost::assign::insert(mandatoryArguments)
	                     ("name", libconfig::Setting::TypeString)
	                     ("resonance", libconfig::Setting::TypeString);
	if(not checkIfAllVariablesAreThere(configAnchorWave, mandatoryArguments)) {
		printErr << "'anchorwave' list in 'input' section in configuration file contains errors." << std::endl;
		return false;
	}

	configAnchorWave->lookupValue("name", _anchorWaveName);
	configAnchorWave->lookupValue("resonance", _anchorComponentName);

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigModelComponents(const libconfig::Setting* configComponents,
                                                        rpwa::massDepFit::model& fitModel,
                                                        rpwa::massDepFit::parameters& fitParameters) const
{
	if(not configComponents) {
		printErr << "'configComponents' is not a pointer to a valid object." << std::endl;
		return false;
	}

	const int nrComponents = configComponents->getLength();

	if(_debug) {
		printDebug << "reading " << nrComponents << " components from configuration file." << std::endl;
	}

	for(int idxComponent=0; idxComponent<nrComponents; ++idxComponent) {
		const libconfig::Setting* configComponent = &((*configComponents)[idxComponent]);

		std::map<std::string, libconfig::Setting::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("name", libconfig::Setting::TypeString);
		if(not checkIfAllVariablesAreThere(configComponent, mandatoryArguments)) {
			printErr << "'components' list in 'model' section in configuration file contains errors." << std::endl;
			return false;
		}

		std::string name;
		configComponent->lookupValue("name", name);

		for(size_t idx=0; idx<fitModel.getNrComponents(); ++idx) {
			if(fitModel.getComponent(idx)->getName() == name) {
				printErr << "component '" << name << "' defined twice." << std::endl;
				return false;
			}
		}

		std::string type;
		if(not configComponent->lookupValue("type", type)) {
			if(_debug) {
				printDebug << "component '" << name << "' has no type, use 'fixedWidthBreitWigner'." << std::endl;
			}
			type = "fixedWidthBreitWigner";
		}

		if(_debug) {
			printDebug << "found component '" << name << "' with type '" << type << "'." << std::endl;
		}

		rpwa::massDepFit::component* component = NULL;
		if(type == "fixedWidthBreitWigner") {
			component = new fixedWidthBreitWigner(fitModel.getNrComponents(), name);
		} else if(type == "dynamicWidthBreitWigner") {
			component = new dynamicWidthBreitWigner(fitModel.getNrComponents(), name);
		} else if(type == "integralWidthBreitWigner") {
			component = new integralWidthBreitWigner(fitModel.getNrComponents(), name);
		} else if(type == "exponentialBackground") {
			component = new exponentialBackground(fitModel.getNrComponents(), name);
		} else if(type == "tPrimeDependentBackground") {
			component = new tPrimeDependentBackground(fitModel.getNrComponents(), name);
			((tPrimeDependentBackground*)component)->setTPrimeMeans(_tPrimeMeans);
		} else {
			printErr << "unknown type '" << type << "' for component '" << name << "'." << std::endl;
			return false;
		}

		if(not component->init(configComponent, fitParameters, _nrBins, _massBinCenters, _waveIndices, _inPhaseSpaceIntegrals, fitModel.useBranchings(), _debug)) {
			delete component;
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
rpwa::massDepFit::massDepFit::readConfigModelFsmd(const libconfig::Setting* configFsmd,
                                                  rpwa::massDepFit::model& fitModel,
                                                  rpwa::massDepFit::parameters& fitParameters) const
{
	// configFsmd might actually be a NULL pointer, in this the final-state
	// mass-dependence is not read
	if(not configFsmd) {
		printInfo << "not using final-state mass dependence." << std::endl;
		return true;
	}

	if(_debug) {
		printDebug << "reading final-state mass-dependence from configuration file." << std::endl;
	}

	rpwa::massDepFit::fsmd* fsmd = new rpwa::massDepFit::fsmd(fitModel.getNrComponents());
	if(not fsmd->init(configFsmd, fitParameters, _massBinCenters, _debug)) {
		delete fsmd;
		printErr << "error while initializing final-state mass-dependence." << std::endl;
		return false;
	}
	fitModel.setFsmd(fsmd);

	printInfo << "using final-state mass dependence as defined in the configuration file." << std::endl;

	return true;
}


bool
rpwa::massDepFit::massDepFit::init(rpwa::massDepFit::model& fitModel,
                                   rpwa::massDepFit::parameters& fitParameters,
                                   rpwa::massDepFit::likelihood& L)
{
	if(not fitModel.init(_waveNames,
	                     _anchorWaveName,
	                     _anchorComponentName)) {
		printErr << "error while initializing the fit model." << std::endl;
		return false;
	}

	if(not L.init(&fitModel,
	              _massBinCenters,
	              _inProductionAmplitudes,
	              _inProductionAmplitudesCovariance,
	              _inSpinDensityMatrices,
	              _inSpinDensityCovarianceMatrices,
	              _wavePairMassBinLimits)) {
		printErr << "error while initializing the likelihood calculator." << std::endl;
		return false;
	}

	return true;
}

bool
rpwa::massDepFit::massDepFit::updateConfig(libconfig::Setting* configRoot,
                                           const rpwa::massDepFit::model& fitModel,
                                           const rpwa::massDepFit::parameters& fitParameters,
                                           const ROOT::Math::Minimizer* minimizer,
                                           const double chi2,
                                           const int ndf,
                                           const double chi2red) const
{
	if(_debug) {
		printDebug << "updating configuration file." << std::endl;
	}

	const libconfig::Setting* configModel = findLibConfigGroup(*configRoot, "model");
	if(not updateConfigModel(configModel, fitModel, fitParameters, minimizer)) {
		printErr << "error while updating 'model' section of configuration file." << std::endl;
		return false;
	}

	libconfig::Setting* configFitquality = NULL;
	if(not configRoot->exists("fitquality")) {
		configFitquality = &(configRoot->add("fitquality", libconfig::Setting::TypeGroup));
	} else {
		configFitquality = &((*configRoot)["fitquality"]);
	}

	if(not configFitquality->exists("chi2")) {
		configFitquality->add("chi2", libconfig::Setting::TypeFloat);
	}
	(*configFitquality)["chi2"] = chi2;

	if(not configFitquality->exists("ndf")) {
		configFitquality->add("ndf", libconfig::Setting::TypeInt);
	}
	(*configFitquality)["ndf"] = ndf;

	if(not configFitquality->exists("redchi2")) {
		configFitquality->add("redchi2", libconfig::Setting::TypeFloat);
	}
	(*configFitquality)["redchi2"] = chi2red;

	return true;
}


bool
rpwa::massDepFit::massDepFit::updateConfigModel(const libconfig::Setting* configModel,
                                                const rpwa::massDepFit::model& fitModel,
                                                const rpwa::massDepFit::parameters& fitParameters,
                                                const ROOT::Math::Minimizer* minimizer) const
{
	if(_debug) {
		printDebug << "updating fit model in configuration file." << std::endl;
	}

	if(not configModel) {
		printErr << "error while updating 'model' section in configuration file." << std::endl;
		return false;
	}

	// update information of the individual components
	const libconfig::Setting* configComponents = findLibConfigList(*configModel, "components");
	if(not updateConfigModelComponents(configComponents, fitModel, fitParameters, minimizer)) {
		printErr << "error while updating 'components' in section '" << configModel->getName() << "' in configuration file." << std::endl;
		return false;
	}

	// update information of the final-state mass-dependence
	const libconfig::Setting* configFsmd = findLibConfigGroup(*configModel, "finalStateMassDependence", false);
	if(not updateConfigModelFsmd(configFsmd, fitModel, fitParameters, minimizer)) {
		printErr << "error while updating 'finalStateMassDependence' in section '" << configModel->getName() << "' in configuration file." << std::endl;
		return false;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::updateConfigModelComponents(const libconfig::Setting* configComponents,
                                                          const rpwa::massDepFit::model& fitModel,
                                                          const rpwa::massDepFit::parameters& fitParameters,
                                                          const ROOT::Math::Minimizer* minimizer) const
{
	if(not configComponents) {
		printErr << "'configComponents' is not a pointer to a valid object." << std::endl;
		return false;
	}

	const int nrComponents = configComponents->getLength();
	if(nrComponents < 0 || static_cast<size_t>(nrComponents) != fitModel.getNrComponents()) {
		printErr << "number of components in configuration file and fit model does not match." << std::endl;
		return false;
	}

	printInfo << "updating " << nrComponents << " components in configuration file." << std::endl;

	for(int idxComponent=0; idxComponent<nrComponents; ++idxComponent) {
		const libconfig::Setting* configComponent = &((*configComponents)[idxComponent]);

		std::map<std::string, libconfig::Setting::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("name", libconfig::Setting::TypeString);
		if(not checkIfAllVariablesAreThere(configComponent, mandatoryArguments)) {
			printErr << "'components' list in 'model' section in configuration file contains errors." << std::endl;
			return false;
		}

		std::string name;
		configComponent->lookupValue("name", name);

		const rpwa::massDepFit::component* component = NULL;
		for(size_t idx=0; idx<fitModel.getNrComponents(); ++idx) {
			if(fitModel.getComponent(idx)->getName() == name) {
				component = fitModel.getComponent(idx);
				break;
			}
		}
		if(not component) {
			printErr << "could not find component '" << name << "' in fit model." << std::endl;
			return false;
		}

		if(not component->update(configComponent, fitParameters, minimizer, fitModel.useBranchings(), _debug)) {
			printErr << "error while updating component '" << name << "'." << std::endl;
			return false;
		}
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::updateConfigModelFsmd(const libconfig::Setting* configFsmd,
                                                    const rpwa::massDepFit::model& fitModel,
                                                    const rpwa::massDepFit::parameters& fitParameters,
                                                    const ROOT::Math::Minimizer* minimizer) const
{
	// configFsmd might actually be a NULL pointer, in this the final-state
	// mass-dependence is not read
	if(not configFsmd) {
		if(fitModel.getFsmd() != NULL) {
			printErr << "no section 'finalStateMassDependence' in configuration file, but final-state mass-dependence exists." << std::endl;
			return false;
		}
		return true;
	}

	if(_debug) {
		printDebug << "updating final-state mass-dependence in configuration file." << std::endl;
	}

	if(not fitModel.getFsmd()->update(configFsmd, fitParameters, minimizer, _debug)) {
		printErr << "error while updating final-state mass-dependence." << std::endl;
		return false;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readInFiles(const std::string& valTreeName,
                                          const std::string& valBranchName)
{
	if(not readInFileFirst(valTreeName, valBranchName)) {
		printErr << "error while reading first file." << std::endl;
		return false;
	}

	for(size_t idxBin=1; idxBin<_nrBins; ++idxBin) {
		if(not readInFile(idxBin, valTreeName, valBranchName)) {
			printErr << "error while reading file entry " << idxBin << "." << std::endl;
			return false;
		}
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readInFileFirst(const std::string& valTreeName,
                                              const std::string& valBranchName)
{
	if(_debug) {
		printDebug << "reading fit result from file '" << _inFileName[0] << "'." << std::endl;
	}

	TFile* inFile = TFile::Open(_inFileName[0].c_str());
	if(not inFile) {
		printErr << "input file '" << _inFileName[0] << "' not found."<< std::endl;
		return false;
	}
	if(inFile->IsZombie()) {
		printErr << "error while reading input file '" << _inFileName[0] << "'."<< std::endl;
		delete inFile;
		return false;
	}

	if(_debug) {
		printDebug << "searching for tree '" << valTreeName << "' in file '" << _inFileName[0] << "'." << std::endl;
	}

	TTree* inTree;
	inFile->GetObject(valTreeName.c_str(), inTree);
	if(not inTree) {
		printErr << "input tree '" << valTreeName << "' not found in input file '" << _inFileName[0] << "'."<< std::endl;
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

	if(not readFitResultMassBins(inTree, inFit)) {
		printErr << "could not extract mass bins from fit result tree in '" << _inFileName[0] << "'." << std::endl;
		delete inFile;
		return false;
	}

	std::vector<Long64_t> inMapping;
	if(not checkFitResultMassBins(inTree, inFit, inMapping)) {
		printErr << "error while checking and mapping mass bins from fit result tree in '" << _inFileName[0] << "'." << std::endl;
		delete inFile;
		return false;
	}

	// resize all array to store the information
	_inProductionAmplitudes.resize(boost::extents[_nrBins][_nrMassBins][_nrWaves]);
	_inProductionAmplitudesCovariance.resize(boost::extents[_nrBins][_nrMassBins][_nrWaves][_nrWaves][2][2]);
	_inSpinDensityMatrices.resize(boost::extents[_nrBins][_nrMassBins][_nrWaves][_nrWaves]);
	_inSpinDensityCovarianceMatrices.resize(boost::extents[_nrBins][_nrMassBins][_nrWaves][_nrWaves][2][2]);
	_inPhaseSpaceIntegrals.resize(boost::extents[_nrBins][_nrMassBins][_nrWaves]);
	_inIntensities.resize(boost::extents[_nrBins][_nrMassBins][_nrWaves][2]);
	_inPhases.resize(boost::extents[_nrBins][_nrMassBins][_nrWaves][_nrWaves][2]);

	boost::multi_array<std::complex<double>, 2> tempProductionAmplitudes;
	boost::multi_array<double, 5> tempProductionAmplitudesCovariance;
	boost::multi_array<std::complex<double>, 3> tempSpinDensityMatrices;
	boost::multi_array<double, 5> tempSpinDensityCovarianceMatrices;
	boost::multi_array<double, 3> tempIntensities;
	boost::multi_array<double, 4> tempPhases;
	if(not readFitResultMatrices(inTree, inFit, inMapping, tempProductionAmplitudes, tempProductionAmplitudesCovariance,
	                             tempSpinDensityMatrices, tempSpinDensityCovarianceMatrices, tempIntensities, tempPhases)) {
		printErr << "error while reading spin-density matrix from fit result tree in '" << _inFileName[0] << "'." << std::endl;
		delete inFile;
		return false;
	}
	_inProductionAmplitudes[0] = tempProductionAmplitudes;
	_inProductionAmplitudesCovariance[0] = tempProductionAmplitudesCovariance;
	_inSpinDensityMatrices[0] = tempSpinDensityMatrices;
	_inSpinDensityCovarianceMatrices[0] = tempSpinDensityCovarianceMatrices;
	_inIntensities[0] = tempIntensities;
	_inPhases[0] = tempPhases;

	boost::multi_array<double, 2> tempPhaseSpaceIntegrals;
	if(not readFitResultIntegrals(inTree, inFit, inMapping, tempPhaseSpaceIntegrals)) {
		printErr << "error while reading phase-space integrals from fit result tree in '" << _inFileName[0] << "'." << std::endl;
		delete inFile;
		return false;
	}

	if(_inOverwritePhaseSpace[0].size() > 0) {
		if(not readPhaseSpaceIntegralMatrices(_inOverwritePhaseSpace[0], tempPhaseSpaceIntegrals)) {
			printErr << "error while reading phase-space integrals from integral matrices." << std::endl;
			delete inFile;
			return false;
		}
	}
	_inPhaseSpaceIntegrals[0] = tempPhaseSpaceIntegrals;

	delete inFile;
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

	std::vector<Long64_t> inMapping;
	if(not checkFitResultMassBins(inTree, inFit, inMapping)) {
		printErr << "error while checking and mapping mass bins from fit result tree in '" << _inFileName[idxBin] << "'." << std::endl;
		delete inFile;
		return false;
	}

	boost::multi_array<std::complex<double>, 2> tempProductionAmplitudes;
	boost::multi_array<double, 5> tempProductionAmplitudesCovariance;
	boost::multi_array<std::complex<double>, 3> tempSpinDensityMatrices;
	boost::multi_array<double, 5> tempSpinDensityCovarianceMatrices;
	boost::multi_array<double, 3> tempIntensities;
	boost::multi_array<double, 4> tempPhases;
	if(not readFitResultMatrices(inTree, inFit, inMapping, tempProductionAmplitudes, tempProductionAmplitudesCovariance,
	                             tempSpinDensityMatrices, tempSpinDensityCovarianceMatrices, tempIntensities, tempPhases)) {
		printErr << "error while reading spin-density matrix from fit result tree in '" << _inFileName[idxBin] << "'." << std::endl;
		delete inFile;
		return false;
	}
	_inProductionAmplitudes[idxBin] = tempProductionAmplitudes;
	_inProductionAmplitudesCovariance[idxBin] = tempProductionAmplitudesCovariance;
	_inSpinDensityMatrices[idxBin] = tempSpinDensityMatrices;
	_inSpinDensityCovarianceMatrices[idxBin] = tempSpinDensityCovarianceMatrices;
	_inIntensities[idxBin] = tempIntensities;
	_inPhases[idxBin] = tempPhases;

	boost::multi_array<double, 2> tempPhaseSpaceIntegrals;
	if(not readFitResultIntegrals(inTree, inFit, inMapping, tempPhaseSpaceIntegrals)) {
		printErr << "error while reading phase-space integrals from fit result tree in '" << _inFileName[idxBin] << "'." << std::endl;
		delete inFile;
		return false;
	}

	if(_inOverwritePhaseSpace[idxBin].size() > 0) {
		if(not readPhaseSpaceIntegralMatrices(_inOverwritePhaseSpace[idxBin], tempPhaseSpaceIntegrals)) {
			printErr << "error while reading phase-space integrals from integral matrices." << std::endl;
			delete inFile;
			return false;
		}
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

	_sysSpinDensityMatrices.resize(boost::extents[_nrSystematics][_nrMassBins][_nrWaves][_nrWaves]);
	_sysSpinDensityCovarianceMatrices.resize(boost::extents[_nrSystematics][_nrMassBins][_nrWaves][_nrWaves][2][2]);
	_sysIntensities.resize(boost::extents[_nrSystematics][_nrMassBins][_nrWaves][2]);
	_sysPhases.resize(boost::extents[_nrSystematics][_nrMassBins][_nrWaves][_nrWaves][2]);

	_sysSpinDensityMatrices[0] = _inSpinDensityMatrices[0];
	_sysSpinDensityCovarianceMatrices[0] = _inSpinDensityCovarianceMatrices[0];
	_sysIntensities[0] = _inIntensities[0];
	_sysPhases[0] = _inPhases[0];

	for(size_t idxSystematics=1; idxSystematics<_nrSystematics; ++idxSystematics) {
		readSystematicsFile(idxSystematics, valTreeName, valBranchName);
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
	if(not checkFitResultMassBins(sysTree, sysFit, sysMapping)) {
		printErr << "error while checking and mapping mass bins from fit result tree in '" << _sysFileNames[idxSystematics-1] << "'." << std::endl;
		delete sysFile;
		return false;
	}

	boost::multi_array<std::complex<double>, 2> tempProductionAmplitudes;
	boost::multi_array<double, 5> tempProductionAmplitudesCovariance;
	boost::multi_array<std::complex<double>, 3> tempSpinDensityMatrices;
	boost::multi_array<double, 5> tempSpinDensityCovarianceMatrices;
	boost::multi_array<double, 3> tempIntensities;
	boost::multi_array<double, 4> tempPhases;
	if(not readFitResultMatrices(sysTree, sysFit, sysMapping, tempProductionAmplitudes, tempProductionAmplitudesCovariance,
	                             tempSpinDensityMatrices, tempSpinDensityCovarianceMatrices, tempIntensities, tempPhases)) {
		printErr << "error while reading spin-density matrix from fit result tree in '" << _sysFileNames[idxSystematics-1] << "'." << std::endl;
		delete sysFile;
		return false;
	}
	_sysSpinDensityMatrices[idxSystematics] = tempSpinDensityMatrices;
	_sysSpinDensityCovarianceMatrices[idxSystematics] = tempSpinDensityCovarianceMatrices;
	_sysIntensities[idxSystematics] = tempIntensities;
	_sysPhases[idxSystematics] = tempPhases;

	delete sysFile;
	return true;
}


bool
rpwa::massDepFit::massDepFit::checkFitResultMassBins(TTree* tree,
                                                     rpwa::fitResult* fit,
                                                     std::vector<Long64_t>& mapping) const
{
	if(not tree or not fit) {
		printErr << "'tree' or 'fit' is not a pointer to a valid object." << std::endl;
		return false;
	}

	// reset mapping
	mapping.assign(_nrMassBins, std::numeric_limits<Long64_t>::max());

	// extract data from tree
	const Long64_t nrEntries = tree->GetEntries();

	if(_debug) {
		printDebug << "check that the centers of mass bins of " << nrEntries << " entries in tree are at a known place, "
		           << "and map the " << _nrMassBins << " mass bins to those entries." << std::endl;
	}

	for(Long64_t idx=0; idx<nrEntries; ++idx) {
		if(tree->GetEntry(idx) == 0) {
			printErr << "error while reading entry " << idx << " from tree." << std::endl;
			return false;
		}
		//FIXME: this would also be the place to select the best fit in case one file contains more than one fit result per mass bin
		const double mass = fit->massBinCenter() / 1000.;

		if(_debug) {
			printDebug << "entry " << idx << ": center of mass bin at " << mass << " GeV/c^2" << std::endl;
		}

		bool found = false;
		size_t idxMass=0;
		while(idxMass<_nrMassBins) {
			if(abs(_massBinCenters[idxMass]-mass) < 1000.*std::numeric_limits<double>::epsilon()) {
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
			printErr << "cannot map tree entry " << idx << " to mass bin " << idxMass << " (" << _massBinCenters[idxMass] << " GeV/c^2)  "
			         << "which is already mapped to tree entry " << mapping[idxMass] << "." << std::endl;
			return false;
		}

		if(_debug) {
			printDebug << "mapping mass bin " << idxMass << " (" << _massBinCenters[idxMass] << " GeV/c^2) to tree entry " << idx << "." << std::endl;
		}
		mapping[idxMass] = idx;
	} // end loop over entries in tree

	// check that all mass bins are mapped
	for(size_t idx=0; idx<mapping.size(); ++idx) {
		if(mapping[idx] == std::numeric_limits<Long64_t>::max()) {
			printErr << "mass bin " << idx << " (" << _massBinCenters[idx] << " GeV/c^2) not mapped." << std::endl;
			return false;
		}
	}

	if(_debug) {
		std::ostringstream output;
		for(size_t idx=0; idx<mapping.size(); ++idx) {
			output << " " << idx << "->" << mapping[idx];
		}
		printDebug << "etablished mapping:" << output.str() << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readFitResultMassBins(TTree* tree,
                                                    rpwa::fitResult* fit)
{
	if(not tree or not fit) {
		printErr << "'tree' or 'fit' is not a pointer to a valid object." << std::endl;
		return false;
	}

	// extract data from tree
	const Long64_t nrEntries = tree->GetEntries();
	_massBinCenters.clear();

	if(_debug) {
		printDebug << "getting center of mass bins from " << nrEntries << " entries in tree." << std::endl;
	}

	for(Long64_t idx=0; idx<nrEntries; ++idx) {
		if(tree->GetEntry(idx) == 0) {
			printErr << "error while reading entry " << idx << " from tree." << std::endl;
			return false;
		}
		const double newMass = fit->massBinCenter() / 1000.;

		if(_debug) {
			printDebug << "entry " << idx << ": center of mass bin at " << newMass << " GeV/c^2" << std::endl;
		}

		bool found = false;
		for(size_t idxMass=0; idxMass<_massBinCenters.size(); ++idxMass) {
			if(abs(_massBinCenters[idxMass]-newMass) < 1000.*std::numeric_limits<double>::epsilon()) {
				found = true;
				if(_debug) {
					printDebug << "this center of mass bin already was encountered before." << std::endl;
				}
				break;
			}
		}

		if(not found) {
			_massBinCenters.push_back(newMass);
		}
	} // end loop over entries in tree

	// sort mass bins
	sort(_massBinCenters.begin(), _massBinCenters.end());

	_nrMassBins = _massBinCenters.size();

	printInfo << "found " << _nrMassBins << " mass bins, center of first and last mass bins: "
	          << _massBinCenters[0] << " and " << _massBinCenters[_nrMassBins - 1] << " GeV/c^2." << std::endl;

	_massStep = (_massBinCenters[_nrMassBins - 1] - _massBinCenters[0]) / (_nrMassBins - 1);
	for(size_t idxMass=1; idxMass<_nrMassBins; ++idxMass) {
		if(abs(_massBinCenters[idxMass]-_massBinCenters[idxMass-1] - _massStep) > 1000.*std::numeric_limits<double>::epsilon()) {
			printErr << "mass distance between bins " << idxMass-1 << " (" << _massBinCenters[idxMass-1] << " GeV/c^2) and "
			         << idxMass << " (" << _massBinCenters[idxMass] << " GeV/c^2) does not agree with nominal distance "
			         << _massStep << " GeV/c^2" << std::endl;
			return false;
		}
	}
	if(_debug) {
		printDebug << "distance between two mass bins is " << _massStep << " GeV/c^2." << std::endl;
	}

	_massMin=_massBinCenters[0] - _massStep / 2;
	_massMax=_massBinCenters[_nrMassBins - 1] + _massStep / 2;
	if(_debug) {
		printDebug << "mass bins cover the mass range from " << _massMin << " to " << _massMax << " GeV/c^2." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readFitResultMatrices(TTree* tree,
                                                    rpwa::fitResult* fit,
                                                    const std::vector<Long64_t>& mapping,
                                                    boost::multi_array<std::complex<double>, 2>& productionAmplitudes,
                                                    boost::multi_array<double, 5>& productionAmplitudesCovariance,
                                                    boost::multi_array<std::complex<double>, 3>& spinDensityMatrices,
                                                    boost::multi_array<double, 5>& spinDensityCovarianceMatrices,
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

	productionAmplitudes.resize(boost::extents[_nrMassBins][_nrWaves]);
	productionAmplitudesCovariance.resize(boost::extents[_nrMassBins][_nrWaves][_nrWaves][2][2]);

	spinDensityMatrices.resize(boost::extents[_nrMassBins][_nrWaves][_nrWaves]);
	spinDensityCovarianceMatrices.resize(boost::extents[_nrMassBins][_nrWaves][_nrWaves][2][2]);

	intensities.resize(boost::extents[_nrMassBins][_nrWaves][2]);
	phases.resize(boost::extents[_nrMassBins][_nrWaves][_nrWaves][2]);

	for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
		if(_debug) {
			printDebug << "reading entry " << mapping[idxMass] << " for mass bin " << idxMass << " (" << _massBinCenters[idxMass] << " GeV/c^2) from tree." << std::endl;
		}
		// FIXME: in case of reading the fit result for a systematic tree this might happen, so this should be allowed in certain cases
		if(tree->GetEntry(mapping[idxMass]) == 0) {
			printErr << "error while reading entry " << mapping[idxMass] << " from tree." << std::endl;
			return false;
		}

		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			const int idx = fit->waveIndex(_waveNames[idxWave]);
			if(idx == -1) {
				printErr << "wave '" << _waveNames[idxWave] << "' not in fit result." << std::endl;
				return false;
			}

			intensities[idxMass][idxWave][0] = fit->intensity(idx);
			intensities[idxMass][idxWave][1] = fit->intensityErr(idx);

			for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
				const int jdx = fit->waveIndex(_waveNames[jdxWave]);
				if(jdx == -1) {
					printErr << "wave '" << _waveNames[jdxWave] << "' not in fit result." << std::endl;
					return false;
				}

				phases[idxMass][idxWave][jdxWave][0] = fit->phase(idx, jdx);
				phases[idxMass][idxWave][jdxWave][1] = fit->phaseErr(idx, jdx);

				spinDensityMatrices[idxMass][idxWave][jdxWave] = fit->spinDensityMatrixElem(idx, jdx);

				const TMatrixD covariance = fit->spinDensityMatrixElemCov(idx, jdx);
				spinDensityCovarianceMatrices[idxMass][idxWave][jdxWave][0][0] = covariance(0, 0);
				spinDensityCovarianceMatrices[idxMass][idxWave][jdxWave][0][1] = covariance(0, 1);
				spinDensityCovarianceMatrices[idxMass][idxWave][jdxWave][1][0] = covariance(1, 0);
				spinDensityCovarianceMatrices[idxMass][idxWave][jdxWave][1][1] = covariance(1, 1);
			}
		}

		// for the production amplitudes loop over the production
		// amplitudes of the fit result
		std::vector<unsigned int> prodAmpIndicesForCov(_nrWaves);
		for(unsigned int idxProdAmp=0; idxProdAmp < fit->nmbProdAmps(); ++idxProdAmp) {
			const std::string waveName = fit->waveNameForProdAmp(idxProdAmp).Data();

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

		const TMatrixD covariance = fit->prodAmpCov(prodAmpIndicesForCov);
		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
				productionAmplitudesCovariance[idxMass][idxWave][jdxWave][0][0] = covariance(2*idxWave + 0, 2*jdxWave + 0);
				productionAmplitudesCovariance[idxMass][idxWave][jdxWave][0][1] = covariance(2*idxWave + 0, 2*jdxWave + 1);
				productionAmplitudesCovariance[idxMass][idxWave][jdxWave][1][0] = covariance(2*idxWave + 1, 2*jdxWave + 0);
				productionAmplitudesCovariance[idxMass][idxWave][jdxWave][1][1] = covariance(2*idxWave + 1, 2*jdxWave + 1);
			}
		}

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
				outputCovariance << " (";
				for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
					output << " " << spinDensityMatrices[idxMass][idxWave][jdxWave];

					outputProdAmpCovariance << " (";
					outputCovariance << " (";
					for(size_t idx=0; idx<2; ++idx) {
						outputProdAmpCovariance << " (";
						outputCovariance << " (";
						for(size_t jdx=0; jdx<2; ++jdx) {
							outputProdAmpCovariance << " " << productionAmplitudesCovariance[idxMass][idxWave][jdxWave][idx][jdx];
							outputCovariance << " " << spinDensityCovarianceMatrices[idxMass][idxWave][jdxWave][idx][jdx];
						}
						outputProdAmpCovariance << " )";
						outputCovariance << " )";
					}
					outputProdAmpCovariance << " )";
					outputCovariance << " )";
				}
				outputProdAmpCovariance << " )";
				output << " )";
				outputCovariance << " )";
			}
			outputProdAmp << " )";

			printDebug << "production amplitudes: " << outputProdAmp.str() << std::endl;
			printDebug << "production amplitudes covariances: " << outputProdAmpCovariance.str() << std::endl;
			printDebug << "spin-density matrix: " << output.str() << std::endl;
			printDebug << "spin-density covariance matrix: " << outputCovariance.str() << std::endl;
		}
	} // end loop over mass bins

	return true;
}


bool
rpwa::massDepFit::massDepFit::readFitResultIntegrals(TTree* tree,
                                                     rpwa::fitResult* fit,
                                                     const std::vector<Long64_t>& mapping,
                                                     boost::multi_array<double, 2>& phaseSpaceIntegrals) const
{
	if(not tree or not fit) {
		printErr << "'tree' or 'fit' is not a pointer to a valid object." << std::endl;
		return false;
	}

	phaseSpaceIntegrals.resize(boost::extents[_nrMassBins][_nrWaves]);

	if(_debug) {
		printDebug << "reading phase-space integrals for " << _nrWaves << " waves from fit result." << std::endl;
	}

	for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
		if(_debug) {
			printDebug << "reading entry " << mapping[idxMass] << " for mass bin " << idxMass << " (" << _massBinCenters[idxMass] << " GeV/c^2) from tree." << std::endl;
		}
		if(tree->GetEntry(mapping[idxMass]) == 0) {
			printErr << "error while reading entry " << mapping[idxMass] << " from tree." << std::endl;
			return false;
		}

		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			const double ps = fit->phaseSpaceIntegral(_waveNames[idxWave]);
			phaseSpaceIntegrals[idxMass][idxWave] = ps;
		}
	}

	if(_debug) {
		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			std::ostringstream output;
			for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
				output << " " << phaseSpaceIntegrals[idxMass][idxWave];
			}
			printDebug << "phase-space integrals for wave '" << _waveNames[idxWave] << "' (" << idxWave << "):" << output.str() << std::endl;
		}
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readPhaseSpaceIntegralMatrices(const std::vector<std::string>& overwritePhaseSpace,
                                                             boost::multi_array<double, 2>& phaseSpaceIntegrals) const
{
	phaseSpaceIntegrals.resize(boost::extents[_nrMassBins][_nrWaves]);

	if(_debug) {
		printDebug << "reading phase-space integrals for " << _nrWaves << " waves from integral matrices." << std::endl;
	}

	for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
		std::ostringstream sFileName;
		for(size_t idxPart=0; idxPart<overwritePhaseSpace.size(); ++idxPart) {
			sFileName << overwritePhaseSpace[idxPart];
			if(idxPart == 0) {
				sFileName << (_massMin + idxMass*_massStep) * 1000.;
			} else if(idxPart == 1) {
				sFileName << (_massMin + (idxMass+1)*_massStep) * 1000.;
			}
		}
		const std::string fileName = sFileName.str();

		if(_debug) {
			printDebug << "reading phase-space integrals for mass bin " << idxMass << " (" << _massBinCenters[idxMass] << " GeV/c^2) from file '" << fileName << "'." << std::endl;
		}

		ampIntegralMatrix intMatrix;
		intMatrix.readAscii(fileName);

		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			const double ps = sqrt(abs(intMatrix.element(_waveNames[idxWave], _waveNames[idxWave])));
			phaseSpaceIntegrals[idxMass][idxWave] = ps;
		}
	}

	if(_debug) {
		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			std::ostringstream output;
			for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
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
	if(_debug) {
		printDebug << "determine which mass bins to use in the fit for " << _nrMassBins << " mass bins, center of first and last mass bins: "
		           << _massBinCenters[0] << " and " << _massBinCenters[_nrMassBins - 1] << " GeV/c^2." << std::endl;
	}

	_waveMassBinLimits.clear();
	for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
		size_t binFirst = 0;
		size_t binLast = _nrMassBins-1;
		for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
			if(_massBinCenters[idxMass] < _waveMassLimits[idxWave].first) {
				binFirst = idxMass+1;
			}
			if(_massBinCenters[idxMass] == _waveMassLimits[idxWave].first) {
				binFirst = idxMass;
			}
			if(_massBinCenters[idxMass] <= _waveMassLimits[idxWave].second) {
				binLast = idxMass;
			}
		}
		if(_waveMassLimits[idxWave].first < 0) {
			binFirst = 0;
		}
		if(_waveMassLimits[idxWave].second < 0) {
			binLast = _nrMassBins-1;
		}
		if(_debug) {
			printDebug << idxWave << ": " << _waveNames[idxWave] << ": "
			           << "mass range: " << (_waveMassLimits[idxWave].first<0. ? _massMin : _waveMassLimits[idxWave].first)
			           << "-" << (_waveMassLimits[idxWave].second<0. ? _massMax : _waveMassLimits[idxWave].second) << " GeV/c^2, "
			           << "bin range " << binFirst << "-" << binLast << std::endl;
		}
		_waveMassBinLimits.push_back(std::make_pair(binFirst, binLast));
	}

	_wavePairMassBinLimits.resize(boost::extents[_nrWaves][_nrWaves]);
	for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
		for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
			_wavePairMassBinLimits[idxWave][jdxWave] = std::make_pair(std::max(_waveMassBinLimits[idxWave].first,  _waveMassBinLimits[jdxWave].first),
			                                                          std::min(_waveMassBinLimits[idxWave].second, _waveMassBinLimits[jdxWave].second));
		}
	}

	if(_debug) {
		printDebug << "waves and mass limits:" << std::endl;
		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			std::ostringstream output;
			for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
				output << _wavePairMassBinLimits[idxWave][jdxWave].first << "-" << _wavePairMassBinLimits[idxWave][jdxWave].second << " ";
			}
			printDebug << _waveNames[idxWave] << " " << _waveMassBinLimits[idxWave].first << "-" << _waveMassBinLimits[idxWave].second
			           << ": " << output.str() << std::endl;
		}
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::createPlots(const rpwa::massDepFit::model& fitModel,
                                          const rpwa::massDepFit::parameters& fitParameters,
                                          TFile* outFile,
                                          const bool rangePlotting) const
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
			if(not createPlotsWave(fitModel, fitParameters, outDirectory, rangePlotting, idxWave, idxBin)) {
				printErr << "error while creating intensity plots for wave '" << _waveNames[idxWave] << "' in bin " << idxBin << "." << std::endl;
				return false;
			}
		}

		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			for(size_t jdxWave=idxWave+1; jdxWave<_nrWaves; ++jdxWave) {
				if(not createPlotsWavePair(fitModel, fitParameters, outDirectory, rangePlotting, idxWave, jdxWave, idxBin)) {
					printErr << "error while creating intensity plots for wave pair '" << _waveNames[idxWave] << "' and '" << _waveNames[jdxWave] << "' in bin " << idxBin << "." << std::endl;
					return false;
				}
			}
		}
	}

	if(fitModel.getFsmd() != NULL) {
		TGraph graph;
		graph.SetName("finalStateMassDependence");
		graph.SetTitle("finalStateMassDependence");
		graph.SetDrawOption("CP");

		Int_t point = -1;
		for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
			++point;
			const double mass = _massBinCenters[idxMass];

			graph.SetPoint(point, mass, pow(fitModel.getFsmd()->val(fitParameters, _massBinCenters[idxMass], idxMass), 2.));
		}

		outFile->cd();
		graph.Write();
	}

	if(_debug) {
		printDebug << "finished creating plots." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::createPlotsWave(const rpwa::massDepFit::model& fitModel,
                                              const rpwa::massDepFit::parameters& fitParameters,
                                              TDirectory* outDirectory,
                                              const bool rangePlotting,
                                              const size_t idxWave,
                                              const size_t idxBin) const
{
	if(_debug) {
		printDebug << "start creating plots for wave '" << _waveNames[idxWave] << "' in bin " << idxBin << "." << std::endl;
	}

	TMultiGraph graphs;
	graphs.SetName(_waveNames[idxWave].c_str());
	graphs.SetTitle(_waveNames[idxWave].c_str());
	graphs.SetDrawOption("AP");

	TGraphErrors* systematics = NULL;
	if(_sysPlotting) {
		systematics = new TGraphErrors;
		systematics->SetName((_waveNames[idxWave] + "__sys").c_str());
		systematics->SetTitle((_waveNames[idxWave] + "__sys").c_str());
		systematics->SetLineColor(kAzure-9);
		systematics->SetFillColor(kAzure-9);
		systematics->SetDrawOption("2");
		graphs.Add(systematics, "2");
	}

	TGraphErrors* data = new TGraphErrors;
	data->SetName((_waveNames[idxWave] + "__data").c_str());
	data->SetTitle((_waveNames[idxWave] + "__data").c_str());
	data->SetDrawOption("AP");
	graphs.Add(data, "P");

	TGraph* fit = new TGraph;
	fit->SetName((_waveNames[idxWave] + "__fit").c_str());
	fit->SetTitle((_waveNames[idxWave] + "__fit").c_str());
	fit->SetLineColor(kRed);
	fit->SetLineWidth(2);
	fit->SetMarkerColor(kRed);
	fit->SetDrawOption("AP");
	graphs.Add(fit, "CP");

	TGraph* phaseSpace = new TGraph;
	phaseSpace->SetName((_waveNames[idxWave] + "__ps").c_str());
	phaseSpace->SetTitle((_waveNames[idxWave] + "__ps").c_str());
	graphs.Add(phaseSpace, "CP");

	const std::vector<std::pair<size_t, size_t> >& compChannel = fitModel.getComponentChannel(idxWave);
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

		graphs.Add(component, "CP");
		components.push_back(component);
	}

	double maxP = -std::numeric_limits<double>::max();
	double maxIE = -std::numeric_limits<double>::max();
	Int_t point = -1;
	Int_t pointLimit = -1;
	for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
		++point;
		const double mass = _massBinCenters[idxMass];
		const double halfBin = _massStep/2.;

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

		double ps = pow(_inPhaseSpaceIntegrals[idxBin][idxMass][idxWave], 2);
		if(fitModel.getFsmd() != NULL) {
			ps *= pow(fitModel.getFsmd()->val(fitParameters, _massBinCenters[idxMass], idxMass), 2);
		}
		phaseSpace->SetPoint(point, mass, ps);
		maxP = std::max(maxP, ps);

		// check that this mass bin should be taken into account for this
		// combination of waves
		if(rangePlotting && (idxMass < _wavePairMassBinLimits[idxWave][idxWave].first || idxMass > _wavePairMassBinLimits[idxWave][idxWave].second)) {
			continue;
		}
		++pointLimit;

		const double intensity = fitModel.intensity(fitParameters, idxWave, idxBin, _massBinCenters[idxMass], idxMass);
		fit->SetPoint(pointLimit, mass, intensity);
		maxIE = std::max(maxIE, intensity);

		for(size_t idxComponents=0; idxComponents<compChannel.size(); ++idxComponents) {
			const size_t idxComponent = compChannel[idxComponents].first;
			const size_t idxChannel = compChannel[idxComponents].second;

			std::complex<double> prodAmp = fitModel.getComponent(idxComponent)->val(fitParameters, idxBin, _massBinCenters[idxMass]);
			prodAmp *= fitModel.getComponent(idxComponent)->getCouplingPhaseSpace(fitParameters, idxChannel, idxBin, _massBinCenters[idxMass], idxMass);
			if(fitModel.getFsmd() != NULL) {
				prodAmp *= fitModel.getFsmd()->val(fitParameters, _massBinCenters[idxMass], idxMass);
			}

			components[idxComponents]->SetPoint(pointLimit, mass, norm(prodAmp));
			maxIE = std::max(maxIE, norm(prodAmp));
		}
	}

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
rpwa::massDepFit::massDepFit::createPlotsWavePair(const rpwa::massDepFit::model& fitModel,
                                                  const rpwa::massDepFit::parameters& fitParameters,
                                                  TDirectory* outDirectory,
                                                  const bool rangePlotting,
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
	phase.SetDrawOption("AP");

	TMultiGraph real;
	real.SetName(realName.c_str());
	real.SetTitle(realName.c_str());
	real.SetDrawOption("AP");

	TMultiGraph imag;
	imag.SetName(imagName.c_str());
	imag.SetTitle(imagName.c_str());
	imag.SetDrawOption("AP");

	TGraphErrors* phaseSystematics = NULL;
	TGraphErrors* realSystematics = NULL;
	TGraphErrors* imagSystematics = NULL;
	if(_sysPlotting) {
		phaseSystematics = new TGraphErrors;
		phaseSystematics->SetName((phaseName + "__sys").c_str());
		phaseSystematics->SetTitle((phaseName + "__sys").c_str());
		phaseSystematics->SetLineColor(kAzure-9);
		phaseSystematics->SetFillColor(kAzure-9);
		phaseSystematics->SetDrawOption("2");
		phase.Add(phaseSystematics, "2");

		realSystematics = new TGraphErrors;
		realSystematics->SetName((realName + "__sys").c_str());
		realSystematics->SetTitle((realName + "__sys").c_str());
		realSystematics->SetLineColor(kAzure-9);
		realSystematics->SetFillColor(kAzure-9);
		realSystematics->SetDrawOption("2");
		real.Add(realSystematics, "2");

		imagSystematics = new TGraphErrors;
		imagSystematics->SetName((imagName + "__sys").c_str());
		imagSystematics->SetTitle((imagName + "__sys").c_str());
		imagSystematics->SetLineColor(kAzure-9);
		imagSystematics->SetFillColor(kAzure-9);
		imagSystematics->SetDrawOption("2");
		imag.Add(imagSystematics, "2");
	}

	TGraphErrors* phaseData = new TGraphErrors;
	phaseData->SetName((phaseName + "__data").c_str());
	phaseData->SetTitle((phaseName + "__data").c_str());
	phaseData->SetDrawOption("AP");
	phase.Add(phaseData, "P");

	TGraphErrors* realData = new TGraphErrors;
	realData->SetName((realName + "__data").c_str());
	realData->SetTitle((realName + "__data").c_str());
	realData->SetDrawOption("AP");
	real.Add(realData, "P");

	TGraphErrors* imagData = new TGraphErrors;
	imagData->SetName((imagName + "__data").c_str());
	imagData->SetTitle((imagName + "__data").c_str());
	imagData->SetDrawOption("AP");
	imag.Add(imagData, "P");

	TGraph* phaseFit = new TGraph;
	phaseFit->SetName((phaseName + "__fit").c_str());
	phaseFit->SetTitle((phaseName + "__fit").c_str());
	phaseFit->SetLineColor(kRed);
	phaseFit->SetLineWidth(2);
	phaseFit->SetMarkerColor(kRed);
	phaseFit->SetDrawOption("AP");
	phase.Add(phaseFit, "CP");

	TGraph* realFit = new TGraph;
	realFit->SetName((realName + "__fit").c_str());
	realFit->SetTitle((realName + "__fit").c_str());
	realFit->SetLineColor(kRed);
	realFit->SetLineWidth(2);
	realFit->SetMarkerColor(kRed);
	realFit->SetDrawOption("AP");
	real.Add(realFit, "CP");

	TGraph* imagFit = new TGraph;
	imagFit->SetName((imagName + "__fit").c_str());
	imagFit->SetTitle((imagName + "__fit").c_str());
	imagFit->SetLineColor(kRed);
	imagFit->SetLineWidth(2);
	imagFit->SetMarkerColor(kRed);
	imagFit->SetDrawOption("AP");
	imag.Add(imagFit, "CP");

	// keep track of phase over full mass range
	TGraph phaseFitAll;

	Int_t point = -1;
	Int_t pointLimit = -1;
	for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
		++point;
		const double mass = _massBinCenters[idxMass];
		const double halfBin = _massStep/2.;

		phaseData->SetPoint(point, mass, _inPhases[idxBin][idxMass][idxWave][jdxWave][0]);
		phaseData->SetPointError(point, halfBin, _inPhases[idxBin][idxMass][idxWave][jdxWave][1]);

		realData->SetPoint(point, mass, _inSpinDensityMatrices[idxBin][idxMass][idxWave][jdxWave].real());
		realData->SetPointError(point, halfBin, sqrt(_inSpinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][0]));

		imagData->SetPoint(point, mass, _inSpinDensityMatrices[idxBin][idxMass][idxWave][jdxWave].imag());
		imagData->SetPointError(point, halfBin, sqrt(_inSpinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][1]));

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
				if(abs(sysP+360.-dataP) < abs(sysP-dataP)) {
					sysP = sysP+360;
				} else if(abs(sysP-360.-dataP) < abs(sysP-dataP)) {
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

		phaseFitAll.SetPoint(point, mass, fitModel.phase(fitParameters, idxWave, jdxWave, idxBin, _massBinCenters[idxMass], idxMass) * TMath::RadToDeg());

		// check that this mass bin should be taken into account for this
		// combination of waves
		if(rangePlotting && (idxMass < _wavePairMassBinLimits[idxWave][jdxWave].first || idxMass > _wavePairMassBinLimits[idxWave][jdxWave].second)) {
			continue;
		}
		++pointLimit;

		const std::complex<double> element = fitModel.spinDensityMatrix(fitParameters, idxWave, jdxWave, idxBin, _massBinCenters[idxMass], idxMass);
		realFit->SetPoint(pointLimit, mass, element.real());
		imagFit->SetPoint(pointLimit, mass, element.imag());
	}

	// rectify phase graphs
	point = -1;
	pointLimit = -1;
	for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
		++point;
		const double mass = _massBinCenters[idxMass];

		double valueFit=0;
		if(point != 0) {
			int bestOffs = 0;
			double bestDiff = std::numeric_limits<double>::max();

			double x, prev, curr;
			phaseFitAll.GetPoint(point-1, x, prev);
			phaseFitAll.GetPoint(point, x, curr);
			for(int offs=-5; offs<6; ++offs) {
				if(abs(curr + offs*360. - prev) < bestDiff) {
					bestDiff = abs(curr + offs*360. - prev);
					bestOffs = offs;
				}
			}

			valueFit = curr + bestOffs*360.;
			phaseFitAll.SetPoint(point, x, valueFit);
		} else {
			double x;
			phaseFitAll.GetPoint(point, x, valueFit);
		}

		int bestOffs = 0;
		double bestDiff = std::numeric_limits<double>::max();

		double x, data;
		phaseData->GetPoint(point, x, data);
		for(int offs=-5; offs<6; ++offs) {
			if(abs(data + offs*360. - valueFit) < bestDiff) {
				bestDiff = abs(data + offs*360. - valueFit);
				bestOffs = offs;
			}
		}

		phaseData->SetPoint(point, x, data + bestOffs*360.);
		if(_sysPlotting) {
			phaseSystematics->GetPoint(point, x, data);
			phaseSystematics->SetPoint(point, x, data + bestOffs*360.);
		}

		// check that this mass bin should be taken into account for this
		// combination of waves
		if(rangePlotting && (idxMass < _wavePairMassBinLimits[idxWave][jdxWave].first || idxMass > _wavePairMassBinLimits[idxWave][jdxWave].second)) {
			continue;
		}
		++pointLimit;

		phaseFit->SetPoint(pointLimit, mass, valueFit);
	}

	outDirectory->cd();
	phase.Write();
	real.Write();
	imag.Write();

	return true;
}
