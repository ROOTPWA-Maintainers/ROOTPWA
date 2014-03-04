///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------
//
// Description:
//      fitting program for massdependent fit rootpwa
//      minimizes rpwa::massDepFit::likelihood function
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


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
#include <TF1.h>
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
#include "massDepFitLikeli.h"
#include "massDepFitModel.h"
#include "reportingUtils.hpp"
#include "reportingUtilsEnvironment.h"

#define MASSSCALE 0.001

using namespace std;
using namespace libconfig;
using namespace ROOT::Math;
using namespace rpwa;
using namespace boost;
using namespace boost::assign;


bool rpwa::massDepFit::massDepFit::_debug = false;


rpwa::massDepFit::massDepFit::massDepFit()
	: _sysPlotting(false),
	  _nrMassBins(0),
	  _nrSystematics(0),
	  _nrWaves(0)
{
}


bool
rpwa::massDepFit::massDepFit::readConfigInput(const Setting* configInput)
{
	if(not configInput) {
		printErr << "'configInput' is not a pointer to a valid object." << endl;
		return false;
	}

	// get information about fit results from mass-independent
	const Setting* configInputFitResults = findLibConfigList(*configInput, "fitresults");
	if(not configInputFitResults) {
		printErr << "'fitresults' list does not exist in section '" << configInput->getName() << "' in configuration file." << endl;
		return false;
	}
	if(not readConfigInputFitResults(configInputFitResults)) {
		printErr << "error while reading 'fitresults' in section '" << configInput->getName() << "' in configuration file." << endl;
		return false;
	}

	// get information about waves to be used in the fit
	const Setting* configInputWaves = findLibConfigList(*configInput, "waves");
	if(not configInputWaves) {
		printErr << "'waves' list does not exist in section '" << configInput->getName() << "' in configuration file." << endl;
		return false;
	}
	if(not readConfigInputWaves(configInputWaves)) {
		printErr << "error while reading 'waves' in section '" << configInput->getName() << "' in configuration file." << endl;
		return false;
	}

	// get information for plotting of systematic error
	const Setting* configInputSystematics = findLibConfigArray(*configInput, "systematics", false);
	if(not readConfigInputSystematics(configInputSystematics)) {
		printErr << "error while reading 'waves' in section '" << configInput->getName() << "' in configuration file." << endl;
		return false;
	}

	// get information for which parameters to release in which order
	const Setting* configInputFreeParameters = findLibConfigArray(*configInput, "freeparameters", false);
	if(not readConfigInputFreeParameters(configInputFreeParameters)) {
		printErr << "error while reading 'freeparameters' in section '" << configInput->getName() << "' in configuration file." << endl;
		return false;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigInputFitResults(const Setting* configInputFitResults)
{
	if(not configInputFitResults) {
		printErr << "'configInputFitResults' is not a pointer to a valid object." << endl;
		return false;
	}

	_inFileName.clear();
	_inOverwritePhaseSpace.clear();

	const int nrFitResults = configInputFitResults->getLength();
	for(int idxFitResult=0; idxFitResult<nrFitResults; ++idxFitResult) {
		if(_debug) {
			printDebug << "reading of entry " << idxFitResult << " in 'fitresults'." << endl;
		}

		const Setting* configInputFitResult = &((*configInputFitResults)[idxFitResult]);

		map<string, Setting::Type> mandatoryArguments;
		insert(mandatoryArguments)
		      ("name", Setting::TypeString);
		if(not checkIfAllVariablesAreThere(configInputFitResult, mandatoryArguments)) {
			printErr << "'fitresults' list in 'input' section in configuration file contains errors." << endl;
			return false;
		}

		string fileName;
		configInputFitResult->lookupValue("name", fileName);
		_inFileName.push_back(fileName);

		if(_debug) {
			printDebug << "read file name of fit results of mass-independent fit: '" << fileName << "'." << endl;
		}

		vector<string> overwritePhaseSpace;

		const Setting* configOverwrite = findLibConfigArray(*configInputFitResult, "overwritePhaseSpace", false);
		if(configOverwrite) {
			const int nrParts = configOverwrite->getLength();

			if(nrParts > 0 && (*configOverwrite)[0].getType() != Setting::TypeString) {
				printErr << "contents of 'overwritePhaseSpace' array in 'input' needs to be strings." << endl;
				return false;
			}

			for(int idxPart=0; idxPart<nrParts; ++idxPart) {
				const string fileName = (*configOverwrite)[idxPart];
				overwritePhaseSpace.push_back(fileName);
			}

			if(_debug) {
				ostringstream output;
				output << "[ '";
				for(size_t idxPart=0; idxPart<overwritePhaseSpace.size(); ++idxPart) {
					if(idxPart > 0) {
						output << "', '";
					}
					output << overwritePhaseSpace[idxPart];
				}
				output << "' ]";
				printDebug << "phase-space integrals will be overwritten with files from this pattern: " << output.str() << " (" << overwritePhaseSpace.size() << " parts)" << endl;
			}
		}

		_inOverwritePhaseSpace.push_back(overwritePhaseSpace);
	}

	_nrBins = _inFileName.size();

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigInputWaves(const Setting* configInputWaves)
{
	if(not configInputWaves) {
		printErr << "'configInputWaves' is not a pointer to a valid object." << endl;
		return false;
	}

	const int nrWaves = configInputWaves->getLength();
	if(_debug) {
		printDebug << "going to read information of " << nrWaves << " waves to be used in the fit." << endl;
	}

	for(int idxWave=0; idxWave<nrWaves; ++idxWave) {
		const Setting* configInputWave = &((*configInputWaves)[idxWave]);

		map<string, Setting::Type> mandatoryArguments;
		insert(mandatoryArguments)
		      ("name", Setting::TypeString);
		if(not checkIfAllVariablesAreThere(configInputWave, mandatoryArguments)) {
			printErr << "'waves' list in 'input' section in configuration file contains errors." << endl;
			return false;
		}

		string name;
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
			printErr << "wave '" << name << "' defined twice." << endl;
			return false;
		}

		_waveNames.push_back(name);
		_waveIndices[name] = _waveNames.size() - 1;
		_waveMassLimits.push_back(make_pair(massLower, massUpper));

		if(_debug) {
			printDebug << idxWave << ": " << name << " (mass range: " << massLower << "-" << massUpper << " MeV/c^2, index: " << _waveIndices[name] << ")" << endl;
		}
	}

	_nrWaves = _waveNames.size();
	if(_debug) {
		printDebug << "read " << _nrWaves << " in total." << endl;
	}

	ostringstream output;
	for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
		output << "    " << _waveNames[idxWave] << endl;
	}
	printInfo << _nrWaves << " waves to be used in fit:" << endl
	          << output.str();

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigInputSystematics(const Setting* configInputSystematics)
{
	// configInputSystematics might actually be a NULL pointer, in this
	// systematics is not plotted
	if(not configInputSystematics) {
		_sysPlotting = false;
		return true;
	}

	const int nrSystematics = configInputSystematics->getLength();
	if(_debug) {
		printDebug << "going to read information for " << nrSystematics << " files containing information for systematic errors." << endl;
	}

	if(nrSystematics > 0) {
		_sysPlotting = true;
	}

	if(nrSystematics > 0 && (*configInputSystematics)[0].getType() != Setting::TypeString) {
		printErr << "contents of 'systematics' array in 'input' needs to be strings." << endl;
		return false;
	}

	for(int idxSystematics=0; idxSystematics<nrSystematics; ++idxSystematics) {
		const string fileName = (*configInputSystematics)[idxSystematics];
		if(_debug) {
			printDebug << "'" << fileName << "' will be read to get information for systematic errors." << endl;
		}
		_sysFileNames.push_back(fileName);
	}

	_nrSystematics = _sysFileNames.size() + 1;
	if(_debug) {
		printDebug << "in total " << _nrSystematics << " files to be read to get information for systematic errors." << endl;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigInputFreeParameters(const Setting* configInputFreeParameters)
{
	if(not configInputFreeParameters) {
		_freeParameters.clear();
		_freeParameters.push_back("");
		_freeParameters.push_back("mass m0");
		_freeParameters.push_back("*");
		return true;
	}

	_freeParameters.clear();

	const int nrItems = configInputFreeParameters->getLength();

	if(_debug) {
		printDebug << "going to extract " << nrItems << " items from '" << configInputFreeParameters->getName() << "'." << endl;
	}

	for(int idxItem=0; idxItem<nrItems; ++idxItem) {
		const Setting& item = (*configInputFreeParameters)[idxItem];

		if(item.getType() != Setting::TypeString) {
			printErr << "'" << configInputFreeParameters->getName() << "' must be array of strings." << endl;
			return false;
		}

		const std::string name = item;

		if(_debug) {
			printDebug << idxItem << ": '" << name << "'." << endl;
		}

		_freeParameters.push_back(name);
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigModel(const Setting* configRoot,
                                              rpwa::massDepFit::model& fitModel)
{
	if(_debug) {
		printDebug << "reading fit model from configuration file." << endl;
	}

	// get components section of configuration file
	const Setting* configModel = findLibConfigGroup(*configRoot, "model");
	if(not configModel) {
		printErr << "error while reading 'components' section in configuration file." << endl;
		return false;
	}

	// get information about anchor wave
	const Setting* configAnchorWave = findLibConfigGroup(*configModel, "anchorwave");
	if(not readConfigModelAnchorWave(configAnchorWave)) {
		printErr << "error while reading 'anchorwave' in section '" << configModel->getName() << "' in configuration file." << endl;
		return false;
	}

	// read information for the individual components
	const Setting* configComponents = findLibConfigList(*configModel, "components");
	if(not readConfigModelComponents(configComponents, fitModel)) {
		printErr << "error while reading 'components' in section '" << configModel->getName() << "' in configuration file." << endl;
		return false;
	}

	// get information for creating the final-state mass-dependence
	const Setting* configFsmd = findLibConfigGroup(*configModel, "finalStateMassDependence", false);
	if(not readConfigModelFsmd(configFsmd, fitModel)) {
		printErr << "error while reading 'finalStateMassDependence' in section '" << configModel->getName() << "' in configuration file." << endl;
		return false;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigModelAnchorWave(const Setting* configAnchorWave)
{
	if(not configAnchorWave) {
		printErr << "'configInputAnchorWave' is not a pointer to a valid object." << endl;
		return false;
	}

	map<string, Setting::Type> mandatoryArguments;
	insert(mandatoryArguments)
	      ("name", Setting::TypeString)
	      ("resonance", Setting::TypeString);
	if(not checkIfAllVariablesAreThere(configAnchorWave, mandatoryArguments)) {
		printErr << "'anchorwave' list in 'input' section in configuration file contains errors." << endl;
		return false;
	}

	configAnchorWave->lookupValue("name", _anchorWaveName);
	configAnchorWave->lookupValue("resonance", _anchorComponentName);

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigModelComponents(const Setting* configComponents,
                                                        rpwa::massDepFit::model& fitModel) const
{
	if(not configComponents) {
		printErr << "'configComponents' is not a pointer to a valid object." << endl;
		return false;
	}

	const int nrComponents = configComponents->getLength();

	if(_debug) {
		printDebug << "reading " << nrComponents << " components from configuration file." << endl;
	}

	for(int idxComponent=0; idxComponent<nrComponents; ++idxComponent) {
		const Setting* configComponent = &((*configComponents)[idxComponent]);

		map<string, Setting::Type> mandatoryArguments;
		insert(mandatoryArguments)
		      ("name", Setting::TypeString);
		if(not checkIfAllVariablesAreThere(configComponent, mandatoryArguments)) {
			printErr << "'components' list in 'components' section in configuration file contains errors." << endl;
			return false;
		}

		string name;
		configComponent->lookupValue("name", name);

		for(size_t idx=0; idx<fitModel.getNrComponents(); ++idx) {
			if(fitModel.getComponent(idx)->getName() == name) {
				printErr << "component '" << name << "' defined twice." << endl;
				return false;
			}
		}

		string type;
		if(not configComponent->lookupValue("type", type)) {
			if(_debug) {
				printDebug << "component '" << name << "' has no type, use 'relativisticBreitWigner'." << endl;
			}
			type = "relativisticBreitWigner";
		}

		if(_debug) {
			printDebug << "found component '" << name << "' with type '" << type << "'." << endl;
		}

		rpwa::massDepFit::component* component = NULL;
		if(type == "relativisticBreitWigner") {
			component = new pwacomponent(name);
		} else if(type == "exponentialBackground") {
			component = new pwabkg(name);
		} else {
			printErr << "unknown type '" << type << "' for component '" << name << "'." << endl;
			return false;
		}

		if(not component->init(configComponent, _nrBins, _massBinCenters, _waveIndices, _inPhaseSpaceIntegrals, _debug)) {
			delete component;
			printErr << "error while initializing component '" << name << "' of type '" << type << "'." << endl;
			return false;
		}

		fitModel.add(component);
	}

	ostringstream output;
	for(size_t idx=0; idx<fitModel.getNrComponents(); ++idx) {
		output << "    " << fitModel.getComponent(idx)->getName() << endl;
	}
	printInfo << "fitting " << fitModel.getNrComponents() << " components to the data:" << endl
	          << output.str();

	return true;
}


bool
rpwa::massDepFit::massDepFit::readConfigModelFsmd(const Setting* configFsmd,
                                                  rpwa::massDepFit::model& fitModel) const
{
	// configFsmd might actually be a NULL pointer, in this the final-state
	// mass-dependence is not read
	if(not configFsmd) {
		printInfo << "not using final-state mass dependence." << endl;
		return true;
	}

	if(_debug) {
		printDebug << "reading final-state mass-dependence from configuration file." << endl;
	}

	map<string, Setting::Type> mandatoryArguments;
	insert(mandatoryArguments)
	      ("formula", Setting::TypeString)
	      ("val", Setting::TypeArray)
	      ("lower", Setting::TypeArray)
	      ("upper", Setting::TypeArray)
	      ("error", Setting::TypeArray)
	      ("fix", Setting::TypeArray);
	if(not checkIfAllVariablesAreThere(configFsmd, mandatoryArguments)) {
		printErr << "'finalStateMassDependence' section in configuration file contains errors." << endl;
		return false;
	}

	string formula;
	configFsmd->lookupValue("formula", formula);

	TF1* fsmdFunction = new TF1("finalStateMassDependence", formula.c_str(), _massMin, _massMax);
	const Int_t nrPar = fsmdFunction->GetNpar();

	const Setting& configFsmdValue = (*configFsmd)["val"];
	if(configFsmdValue.getLength() != nrPar) {
		printErr << "'val' in 'finalStateMassDependence' has to have a length of " << nrPar << "." << endl;
		return false;
	}
	if(configFsmdValue.getLength() > 0 and not configFsmdValue[0].isNumber()) {
		printErr << "'val' in 'finalStateMassDependence' has to be made up of numbers." << endl;
		return false;
	}

	const Setting& configFsmdLower = (*configFsmd)["lower"];
	if(configFsmdLower.getLength() != nrPar) {
		printErr << "'lower' in 'finalStateMassDependence' has to have a length of " << nrPar << "." << endl;
		return false;
	}
	if(configFsmdLower.getLength() > 0 and not configFsmdLower[0].isNumber()) {
		printErr << "'lower' in 'finalStateMassDependence' has to be made up of numbers." << endl;
		return false;
	}

	const Setting& configFsmdUpper = (*configFsmd)["upper"];
	if(configFsmdUpper.getLength() != nrPar) {
		printErr << "'upper' in 'finalStateMassDependence' has to have a length of " << nrPar << "." << endl;
		return false;
	}
	if(configFsmdUpper.getLength() > 0 and not configFsmdUpper[0].isNumber()) {
		printErr << "'upper' in 'finalStateMassDependence' has to be made up of numbers." << endl;
		return false;
	}

	const Setting& configFsmdError = (*configFsmd)["error"];
	if(configFsmdError.getLength() != nrPar) {
		printErr << "'error' in 'finalStateMassDependence' has to have a length of " << nrPar << "." << endl;
		return false;
	}
	if(configFsmdError.getLength() > 0 and not configFsmdError[0].isNumber()) {
		printErr << "'error' in 'finalStateMassDependence' has to be made up of numbers." << endl;
		return false;
	}

	const Setting& configFsmdFix = (*configFsmd)["fix"];
	if(configFsmdFix.getLength() != nrPar) {
		printErr << "'fix' in 'finalStateMassDependence' has to have a length of " << nrPar << "." << endl;
		return false;
	}
	if(configFsmdFix.getLength() > 0 and not configFsmdFix[0].isNumber()) {
		printErr << "'fix' in 'finalStateMassDependence' has to be made up of numbers." << endl;
		return false;
	}

	for(Int_t i=0; i<nrPar; ++i) {
		fsmdFunction->SetParameter(i, configFsmdValue[i]);
		fsmdFunction->SetParError(i, configFsmdError[i]);
		fsmdFunction->SetParLimits(i, configFsmdLower[i], configFsmdUpper[i]);

		if(((int)configFsmdFix[i]) != 0) {
			fsmdFunction->FixParameter(i, configFsmdValue[i]);
		}
	}

	fitModel.setFsmdFunction(fsmdFunction);

	printInfo << "using final-state mass dependence as defined in the configuration file." << endl;

	return true;
}


bool
rpwa::massDepFit::massDepFit::updateConfig(Setting* configRoot,
                                           const rpwa::massDepFit::model& fitModel,
                                           const Minimizer* minimizer,
                                           const double chi2,
                                           const int ndf,
                                           const double chi2red) const
{
	if(_debug) {
		printDebug << "updating configuration file." << endl;
	}

	const Setting* configModel = findLibConfigGroup(*configRoot, "model");
	if(not updateConfigModel(configModel, fitModel, minimizer)) {
		printErr << "error while updating 'model' section of configuration file." << endl;
		return false;
	}

	Setting* configFitquality = NULL;
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
rpwa::massDepFit::massDepFit::updateConfigModel(const Setting* configModel,
                                                const rpwa::massDepFit::model& fitModel,
                                                const Minimizer* minimizer) const
{
	if(_debug) {
		printDebug << "updating fit model in configuration file." << endl;
	}

	if(not configModel) {
		printErr << "error while updating 'components' section in configuration file." << endl;
		return false;
	}

	// update information of the individual components
	const Setting* configComponents = findLibConfigList(*configModel, "components");
	if(not updateConfigModelComponents(configComponents, fitModel, minimizer)) {
		printErr << "error while updating 'components' in section '" << configModel->getName() << "' in configuration file." << endl;
		return false;
	}

	// update information of the final-state mass-dependence
	const Setting* configFsmd = findLibConfigGroup(*configModel, "finalStateMassDependence", false);
	if(not updateConfigModelFsmd(configFsmd, fitModel, minimizer)) {
		printErr << "error while updating 'finalStateMassDependence' in section '" << configModel->getName() << "' in configuration file." << endl;
		return false;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::updateConfigModelComponents(const Setting* configComponents,
                                                          const rpwa::massDepFit::model& fitModel,
                                                          const Minimizer* minimizer) const
{
	if(not configComponents) {
		printErr << "'configComponents' is not a pointer to a valid object." << endl;
		return false;
	}

	const int nrComponents = configComponents->getLength();
	if(nrComponents < 0 || static_cast<size_t>(nrComponents) != fitModel.getNrComponents()) {
		printErr << "number of components in configuration file and fit model does not match." << endl;
		return false;
	}

	printInfo << "updating " << nrComponents << " components in configuration file." << endl;

	for(int idxComponent=0; idxComponent<nrComponents; ++idxComponent) {
		const Setting* configComponent = &((*configComponents)[idxComponent]);

		map<string, Setting::Type> mandatoryArguments;
		insert(mandatoryArguments)
		      ("name", Setting::TypeString);
		if(not checkIfAllVariablesAreThere(configComponent, mandatoryArguments)) {
			printErr << "'components' list in 'components' section in configuration file contains errors." << endl;
			return false;
		}

		string name;
		configComponent->lookupValue("name", name);

		const rpwa::massDepFit::component* component = NULL;
		for(size_t idx=0; fitModel.getNrComponents(); ++idx) {
			if(fitModel.getComponent(idx)->getName() == name) {
				component = fitModel.getComponent(idx);
				break;
			}
		}
		if(not component) {
			printErr << "could not find component '" << name << "' in fit model." << endl;
			return false;
		}

		if(not component->update(configComponent, minimizer, _debug)) {
			printErr << "error while updating component '" << name << "'." << endl;
			return false;
		}
	}


	return true;
}


bool
rpwa::massDepFit::massDepFit::updateConfigModelFsmd(const Setting* configFsmd,
                                                    const rpwa::massDepFit::model& fitModel,
                                                    const Minimizer* minimizer) const
{
	// configFsmd might actually be a NULL pointer, in this the final-state
	// mass-dependence is not read
	if(not configFsmd) {
		if(fitModel.getFsmdFunction()) {
			printErr << "no section 'finalStateMassDependence' in configuration file, but final-state mass-dependence exists." << endl;
			return false;
		}
		return true;
	}

	if(_debug) {
		printDebug << "updating final-state mass-dependence in configuration file." << endl;
	}

	const Setting& configFsmdValue = (*configFsmd)["val"];
	const Setting& configFsmdError = (*configFsmd)["error"];
	const Setting& configFsmdFix = (*configFsmd)["fix"];

	const unsigned int nrPar =fitModel.getFsmdFunction()->GetNpar();
	unsigned int iPar = 0;
	for(unsigned int i=0; i<nrPar; ++i) {
		if(((int)configFsmdFix[i]) == 0) {
			configFsmdValue[i] = fitModel.getFsmdParameter(iPar);

			ostringstream sName;
			sName << "PSP__" << iPar;
			const int varIndex = minimizer->VariableIndex(sName.str());
			if(varIndex == -1) {
				printErr << "variable '" << sName.str() << "' used to extract the error for one parameter "
				         << "of the final-state mass-dependence  not known to the minimizer." << endl;
				return false;
			}
			configFsmdError[i] = minimizer->Errors()[varIndex];

			++iPar;
		}
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readInFiles(const string& valTreeName,
                                          const string& valBranchName)
{
	if(not readInFileFirst(valTreeName, valBranchName)) {
		printErr << "error while reading first file." << endl;
		return false;
	}

	for(size_t idxBin=1; idxBin<_nrBins; ++idxBin) {
		if(not readInFile(idxBin, valTreeName, valBranchName)) {
			printErr << "error while reading file entry " << idxBin << "." << endl;
			return false;
		}
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readInFileFirst(const string& valTreeName,
                                              const string& valBranchName)
{
	if(_debug) {
		printDebug << "reading fit result from file '" << _inFileName[0] << "'." << endl;
	}

	TFile* inFile = TFile::Open(_inFileName[0].c_str());
	if(not inFile) {
		printErr << "input file '" << _inFileName[0] << "' not found."<< endl;
		return false;
	}
	if(inFile->IsZombie()) {
		printErr << "error while reading input file '" << _inFileName[0] << "'."<< endl;
		delete inFile;
		return false;
	}

	if(_debug) {
		printDebug << "searching for tree '" << valTreeName << "' in file '" << _inFileName[0] << "'." << endl;
	}

	TTree* inTree;
	inFile->GetObject(valTreeName.c_str(), inTree);
	if(not inTree) {
		printErr << "input tree '" << valTreeName << "' not found in input file '" << _inFileName[0] << "'."<< endl;
		delete inFile;
		return false;
	}

	if(_debug) {
		printDebug << "searching for branch '" << valBranchName << "' in tree '" << valTreeName << "'." << endl;
	}

	fitResult* inFit = NULL;
	if(inTree->SetBranchAddress(valBranchName.c_str(), &inFit)) {
		printErr << "branch '" << valBranchName << "' not found in input tree '" << valTreeName << "'." << endl;
		delete inFile;
		return false;
	}

	if(not readFitResultMassBins(inTree, inFit)) {
		printErr << "could not extract mass bins from fit result tree in '" << _inFileName[0] << "'." << endl;
		delete inFile;
		return false;
	}

	vector<Long64_t> inMapping;
	if(not checkFitResultMassBins(inTree, inFit, inMapping)) {
		printErr << "error while checking and mapping mass bins from fit result tree in '" << _inFileName[0] << "'." << endl;
		delete inFile;
		return false;
	}

	// resize all array to store the information
	_inSpinDensityMatrices.resize(extents[_nrBins][_nrMassBins][_nrWaves][_nrWaves]);
	_inSpinDensityCovarianceMatrices.resize(extents[_nrBins][_nrMassBins][_nrWaves][_nrWaves][2][2]);
	_inPhaseSpaceIntegrals.resize(extents[_nrBins][_nrMassBins][_nrWaves]);
	_inIntensities.resize(extents[_nrBins][_nrMassBins][_nrWaves][2]);
	_inPhases.resize(extents[_nrBins][_nrMassBins][_nrWaves][_nrWaves][2]);

	multi_array<complex<double>, 3> tempSpinDensityMatrices;
	multi_array<double, 5> tempSpinDensityCovarianceMatrices;
	multi_array<double, 3> tempIntensities;
	multi_array<double, 4> tempPhases;
	if(not readFitResultMatrices(inTree, inFit, inMapping, tempSpinDensityMatrices, tempSpinDensityCovarianceMatrices,
	                             tempIntensities, tempPhases)) {
		printErr << "error while reading spin-density matrix from fit result tree in '" << _inFileName[0] << "'." << endl;
		delete inFile;
		return false;
	}
	_inSpinDensityMatrices[0] = tempSpinDensityMatrices;
	_inSpinDensityCovarianceMatrices[0] = tempSpinDensityCovarianceMatrices;
	_inIntensities[0] = tempIntensities;
	_inPhases[0] = tempPhases;

	multi_array<double, 2> tempPhaseSpaceIntegrals;
	if(not readFitResultIntegrals(inTree, inFit, inMapping, tempPhaseSpaceIntegrals)) {
		printErr << "error while reading phase-space integrals from fit result tree in '" << _inFileName[0] << "'." << endl;
		delete inFile;
		return false;
	}

	if(_inOverwritePhaseSpace[0].size() > 0) {
		if(not readPhaseSpaceIntegralMatrices(_inOverwritePhaseSpace[0], tempPhaseSpaceIntegrals)) {
			printErr << "error while reading phase-space integrals from integral matrices." << endl;
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
                                         const string& valTreeName,
                                         const string& valBranchName)
{
	if(_debug) {
		printDebug << "reading fit result from file '" << _inFileName[idxBin] << "'." << endl;
	}

	TFile* inFile = TFile::Open(_inFileName[idxBin].c_str());
	if(not inFile) {
		printErr << "input file '" << _inFileName[idxBin] << "' not found."<< endl;
		return false;
	}
	if(inFile->IsZombie()) {
		printErr << "error while reading input file '" << _inFileName[idxBin] << "'."<< endl;
		delete inFile;
		return false;
	}

	if(_debug) {
		printDebug << "searching for tree '" << valTreeName << "' in file '" << _inFileName[idxBin] << "'." << endl;
	}

	TTree* inTree;
	inFile->GetObject(valTreeName.c_str(), inTree);
	if(not inTree) {
		printErr << "input tree '" << valTreeName << "' not found in input file '" << _inFileName[idxBin] << "'."<< endl;
		delete inFile;
		return false;
	}

	if(_debug) {
		printDebug << "searching for branch '" << valBranchName << "' in tree '" << valTreeName << "'." << endl;
	}

	fitResult* inFit = NULL;
	if(inTree->SetBranchAddress(valBranchName.c_str(), &inFit)) {
		printErr << "branch '" << valBranchName << "' not found in input tree '" << valTreeName << "'." << endl;
		delete inFile;
		return false;
	}

	vector<Long64_t> inMapping;
	if(not checkFitResultMassBins(inTree, inFit, inMapping)) {
		printErr << "error while checking and mapping mass bins from fit result tree in '" << _inFileName[idxBin] << "'." << endl;
		delete inFile;
		return false;
	}

	multi_array<complex<double>, 3> tempSpinDensityMatrices;
	multi_array<double, 5> tempSpinDensityCovarianceMatrices;
	multi_array<double, 3> tempIntensities;
	multi_array<double, 4> tempPhases;
	if(not readFitResultMatrices(inTree, inFit, inMapping, tempSpinDensityMatrices, tempSpinDensityCovarianceMatrices,
	                             tempIntensities, tempPhases)) {
		printErr << "error while reading spin-density matrix from fit result tree in '" << _inFileName[idxBin] << "'." << endl;
		delete inFile;
		return false;
	}
	_inSpinDensityMatrices[idxBin] = tempSpinDensityMatrices;
	_inSpinDensityCovarianceMatrices[idxBin] = tempSpinDensityCovarianceMatrices;
	_inIntensities[idxBin] = tempIntensities;
	_inPhases[idxBin] = tempPhases;

	multi_array<double, 2> tempPhaseSpaceIntegrals;
	if(not readFitResultIntegrals(inTree, inFit, inMapping, tempPhaseSpaceIntegrals)) {
		printErr << "error while reading phase-space integrals from fit result tree in '" << _inFileName[idxBin] << "'." << endl;
		delete inFile;
		return false;
	}

	if(_inOverwritePhaseSpace[idxBin].size() > 0) {
		if(not readPhaseSpaceIntegralMatrices(_inOverwritePhaseSpace[idxBin], tempPhaseSpaceIntegrals)) {
			printErr << "error while reading phase-space integrals from integral matrices." << endl;
			delete inFile;
			return false;
		}
	}
	_inPhaseSpaceIntegrals[idxBin] = tempPhaseSpaceIntegrals;

	delete inFile;
	return true;
}


bool
rpwa::massDepFit::massDepFit::readSystematicsFiles(const string& valTreeName,
                                                   const string& valBranchName)
{
	if(not _sysPlotting) {
		return true;
	}

	// FIXME: make systematic errors work with multiple bins
	if(_nrBins > 1) {
		printErr << "systematic errors not yet supported for multiple bins." << endl;
		return false;
	}

	if(_debug) {
		printDebug << "reading fit results for systematic errors from " << _nrSystematics << " files." << endl;
	}

	_sysSpinDensityMatrices.resize(extents[_nrSystematics][_nrMassBins][_nrWaves][_nrWaves]);
	_sysSpinDensityCovarianceMatrices.resize(extents[_nrSystematics][_nrMassBins][_nrWaves][_nrWaves][2][2]);
	_sysIntensities.resize(extents[_nrSystematics][_nrMassBins][_nrWaves][2]);
	_sysPhases.resize(extents[_nrSystematics][_nrMassBins][_nrWaves][_nrWaves][2]);

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
                                                  const string& valTreeName,
                                                  const string& valBranchName)
{
	if(_debug) {
		printDebug << "reading fit result for systematics for index " << idxSystematics << " from file '" << _sysFileNames[idxSystematics-1] << "'." << endl;
	}

	TFile* sysFile = TFile::Open(_sysFileNames[idxSystematics-1].c_str());
	if(not sysFile) {
		printErr << "input file '" << _sysFileNames[idxSystematics-1] << "' not found."<< endl;
		return false;
	}
	if(sysFile->IsZombie()) {
		printErr << "error while reading input file '" << _sysFileNames[idxSystematics-1] << "'."<< endl;
		delete sysFile;
		return false;
	}

	if(_debug) {
		printDebug << "searching for tree '" << valTreeName << "' in file '" << _sysFileNames[idxSystematics-1] << "'." << endl;
	}

	TTree* sysTree;
	sysFile->GetObject(valTreeName.c_str(), sysTree);
	if(not sysTree) {
		printErr << "input tree '" << valTreeName << "' not found in input file '" << _sysFileNames[idxSystematics-1] << "'."<< endl;
		delete sysFile;
		return false;
	}

	if(_debug) {
		printDebug << "searching for branch '" << valBranchName << "' in tree '" << valTreeName << "'." << endl;
	}

	fitResult* sysFit = NULL;
	if(sysTree->SetBranchAddress(valBranchName.c_str(), &sysFit)) {
		printErr << "branch '" << valBranchName << "' not found in input tree '" << valTreeName << "'." << endl;
		delete sysFile;
		return false;
	}

	vector<Long64_t> sysMapping;
	if(not checkFitResultMassBins(sysTree, sysFit, sysMapping)) {
		printErr << "error while checking and mapping mass bins from fit result tree in '" << _sysFileNames[idxSystematics-1] << "'." << endl;
		delete sysFile;
		return false;
	}

	multi_array<complex<double>, 3> tempSpinDensityMatrices;
	multi_array<double, 5> tempSpinDensityCovarianceMatrices;
	multi_array<double, 3> tempIntensities;
	multi_array<double, 4> tempPhases;
	if(not readFitResultMatrices(sysTree, sysFit, sysMapping, tempSpinDensityMatrices, tempSpinDensityCovarianceMatrices,
	                             tempIntensities, tempPhases)) {
		printErr << "error while reading spin-density matrix from fit result tree in '" << _sysFileNames[idxSystematics-1] << "'." << endl;
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
                                                     vector<Long64_t>& mapping) const
{
	if(not tree or not fit) {
		printErr << "'tree' or 'fit' is not a pointer to a valid object." << endl;
		return false;
	}

	// reset mapping
	mapping.assign(_nrMassBins, numeric_limits<Long64_t>::max());

	// extract data from tree
	const Long64_t nrEntries = tree->GetEntries();

	if(_debug) {
		printDebug << "check that the centers of mass bins of " << nrEntries << " entries in tree are at a known place, "
		           << "and map the " << _nrMassBins << " mass bins to those entries." << endl;
	}

	for(Long64_t idx=0; idx<nrEntries; ++idx) {
		if(tree->GetEntry(idx) == 0) {
			printErr << "error while reading entry " << idx << " from tree." << endl;
			return false;
		}
		//FIXME: this would also be the place to select the best fit in case one file contains more than one fit result per mass bin
		const double mass = fit->massBinCenter();

		if(_debug) {
			printDebug << "entry " << idx << ": center of mass bin at " << mass << " MeV/c^2" << endl;
		}

		bool found = false;
		size_t idxMass=0;
		while(idxMass<_nrMassBins) {
			if(abs(_massBinCenters[idxMass]-mass) < 1000.*numeric_limits<double>::epsilon()) {
				found = true;
				break;
			}
			++idxMass;
		}

		if(not found) {
			printErr << "could not map mass bin centered at " << mass << " MeV/c^2 to a known mass bin." << endl;
			return false;
		}

		if(mapping[idxMass] != numeric_limits<Long64_t>::max()) {
			printErr << "cannot map tree entry " << idx << " to mass bin " << idxMass << " (" << _massBinCenters[idxMass] << " MeV/c^2)  "
			         << "which is already mapped to tree entry " << mapping[idxMass] << "." << endl;
			return false;
		}

		if(_debug) {
			printDebug << "mapping mass bin " << idxMass << " (" << _massBinCenters[idxMass] << " MeV/c^2) to tree entry " << idx << "." << endl;
		}
		mapping[idxMass] = idx;
	} // end loop over entries in tree

	// check that all mass bins are mapped
	for(size_t idx=0; idx<mapping.size(); ++idx) {
		if(mapping[idx] == numeric_limits<Long64_t>::max()) {
			printErr << "mass bin " << idx << " (" << _massBinCenters[idx] << " MeV/c^2) not mapped." << endl;
			return false;
		}
	}

	if(_debug) {
		ostringstream output;
		for(size_t idx=0; idx<mapping.size(); ++idx) {
			output << " " << idx << "->" << mapping[idx];
		}
		printDebug << "etablished mapping:" << output.str() << endl;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readFitResultMassBins(TTree* tree,
                                                    rpwa::fitResult* fit)
{
	if(not tree or not fit) {
		printErr << "'tree' or 'fit' is not a pointer to a valid object." << endl;
		return false;
	}

	// extract data from tree
	const Long64_t nrEntries = tree->GetEntries();
	_massBinCenters.clear();

	if(_debug) {
		printDebug << "getting center of mass bins from " << nrEntries << " entries in tree." << endl;
	}

	for(Long64_t idx=0; idx<nrEntries; ++idx) {
		if(tree->GetEntry(idx) == 0) {
			printErr << "error while reading entry " << idx << " from tree." << endl;
			return false;
		}
		const double newMass = fit->massBinCenter();

		if(_debug) {
			printDebug << "entry " << idx << ": center of mass bin at " << newMass << " MeV/c^2" << endl;
		}

		bool found = false;
		for(size_t idxMass=0; idxMass<_massBinCenters.size(); ++idxMass) {
			if(abs(_massBinCenters[idxMass]-newMass) < 1000.*numeric_limits<double>::epsilon()) {
				found = true;
				if(_debug) {
					printDebug << "this center of mass bin already was encountered before." << endl;
				}
				break;
			}
		}

		if(not found) {
			_massBinCenters.push_back(fit->massBinCenter());
		}
	} // end loop over entries in tree

	// sort mass bins
	sort(_massBinCenters.begin(), _massBinCenters.end());

	_nrMassBins = _massBinCenters.size();

	printInfo << "found " << _nrMassBins << " mass bins, center of first and last mass bins: "
	          << _massBinCenters[0] << " and " << _massBinCenters[_nrMassBins - 1] << " MeV/c^2." << endl;

	_massStep = (_massBinCenters[_nrMassBins - 1] - _massBinCenters[0]) / (_nrMassBins - 1);
	for(size_t idxMass=1; idxMass<_nrMassBins; ++idxMass) {
		if(abs(_massBinCenters[idxMass]-_massBinCenters[idxMass-1] - _massStep) > 1000.*numeric_limits<double>::epsilon()) {
			printErr << "mass distance between bins " << idxMass-1 << " (" << _massBinCenters[idxMass-1] << " MeV/c^2) and "
			         << idxMass << " (" << _massBinCenters[idxMass] << " MeV/c^2) does not agree with nominal distance "
			         << _massStep << " MeV/c^2" << endl;
			return false;
		}
	}
	if(_debug) {
		printDebug << "distance between two mass bins is " << _massStep << " MeV/c^2." << endl;
	}

	_massMin=_massBinCenters[0] - _massStep / 2;
	_massMax=_massBinCenters[_nrMassBins - 1] + _massStep / 2;
	if(_debug) {
		printDebug << "mass bins cover the mass range from " << _massMin << " to " << _massMax << " MeV/c^2." << endl;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readFitResultMatrices(TTree* tree,
                                                    rpwa::fitResult* fit,
                                                    const vector<Long64_t>& mapping,
                                                    multi_array<complex<double>, 3>& spinDensityMatrices,
                                                    multi_array<double, 5>& spinDensityCovarianceMatrices,
                                                    multi_array<double, 3>& intensities,
                                                    multi_array<double, 4>& phases) const
{
	if(not tree or not fit) {
		printErr << "'tree' or 'fit' is not a pointer to a valid object." << endl;
		return false;
	}

	if(_debug) {
		printDebug << "reading spin-density matrices for " << _nrWaves << " waves from fit result." << endl;
	}

	spinDensityMatrices.resize(extents[_nrMassBins][_nrWaves][_nrWaves]);
	spinDensityCovarianceMatrices.resize(extents[_nrMassBins][_nrWaves][_nrWaves][2][2]);

	intensities.resize(extents[_nrMassBins][_nrWaves][2]);
	phases.resize(extents[_nrMassBins][_nrWaves][_nrWaves][2]);

	for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
		if(_debug) {
			printDebug << "reading entry " << mapping[idxMass] << " for mass bin " << idxMass << " (" << _massBinCenters[idxMass] << " MeV/c^2) from tree." << endl;
		}
		// FIXME: in case of reading the fit result for a systematic tree this might happen, so this should be allowed in certain cases
		if(tree->GetEntry(mapping[idxMass]) == 0) {
			printErr << "error while reading entry " << mapping[idxMass] << " from tree." << endl;
			return false;
		}

		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			const int idx = fit->waveIndex(_waveNames[idxWave]);
			if(idx == -1) {
				printErr << "wave '" << _waveNames[idxWave] << "' not in fit result." << endl;
				return false;
			}

			intensities[idxMass][idxWave][0] = fit->intensity(idx);
			intensities[idxMass][idxWave][1] = fit->intensityErr(idx);

			for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
				const int jdx = fit->waveIndex(_waveNames[jdxWave]);
				if(jdx == -1) {
					printErr << "wave '" << _waveNames[jdxWave] << "' not in fit result." << endl;
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

		if(_debug) {
			ostringstream output;
			ostringstream outputCovariance;

			for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
				output << " (";
				outputCovariance << " (";
				for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
					output << " " << spinDensityMatrices[idxMass][idxWave][jdxWave];
					outputCovariance << " (";
					for(size_t idx=0; idx<2; ++idx) {
						outputCovariance << " (";
						for(size_t jdx=0; jdx<2; ++jdx) {
							outputCovariance << " " << spinDensityCovarianceMatrices[idxMass][idxWave][jdxWave][idx][jdx];
						}
						outputCovariance << " )";
					}
					outputCovariance << " )";
				}
				output << " )";
				outputCovariance << " )";
			}

			printDebug << "spin-density matrix: " << output.str() << endl;
			printDebug << "spin-density covariance matrix: " << outputCovariance.str() << endl;
		}
	} // end loop over mass bins

	return true;
}


bool
rpwa::massDepFit::massDepFit::readFitResultIntegrals(TTree* tree,
                                                     rpwa::fitResult* fit,
                                                     const vector<Long64_t>& mapping,
                                                     multi_array<double, 2>& phaseSpaceIntegrals) const
{
	if(not tree or not fit) {
		printErr << "'tree' or 'fit' is not a pointer to a valid object." << endl;
		return false;
	}

	phaseSpaceIntegrals.resize(extents[_nrMassBins][_nrWaves]);

	if(_debug) {
		printDebug << "reading phase-space integrals for " << _nrWaves << " waves from fit result." << endl;
	}

	for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
		if(_debug) {
			printDebug << "reading entry " << mapping[idxMass] << " for mass bin " << idxMass << " (" << _massBinCenters[idxMass] << " MeV/c^2) from tree." << endl;
		}
		if(tree->GetEntry(mapping[idxMass]) == 0) {
			printErr << "error while reading entry " << mapping[idxMass] << " from tree." << endl;
			return false;
		}

		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			const double ps = fit->phaseSpaceIntegral(_waveNames[idxWave]);
			phaseSpaceIntegrals[idxMass][idxWave] = ps;
		}
	}

	if(_debug) {
		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			ostringstream output;
			for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
				output << " " << phaseSpaceIntegrals[idxMass][idxWave];
			}
			printDebug << "phase-space integrals for wave '" << _waveNames[idxWave] << "' (" << idxWave << "):" << output.str() << endl;
		}
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::readPhaseSpaceIntegralMatrices(const vector<string>& overwritePhaseSpace,
                                                             multi_array<double, 2>& phaseSpaceIntegrals) const
{
	phaseSpaceIntegrals.resize(extents[_nrMassBins][_nrWaves]);

	if(_debug) {
		printDebug << "reading phase-space integrals for " << _nrWaves << " waves from integral matrices." << endl;
	}

	for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
		ostringstream sFileName;
		for(size_t idxPart=0; idxPart<overwritePhaseSpace.size(); ++idxPart) {
			sFileName << overwritePhaseSpace[idxPart];
			if(idxPart == 0) {
				sFileName << _massMin + idxMass*_massStep;
			} else if(idxPart == 1) {
				sFileName << _massMin + (idxMass+1)*_massStep;
			}
		}
		const string fileName = sFileName.str();

		if(_debug) {
			printDebug << "reading phase-space integrals for mass bin " << idxMass << " (" << _massBinCenters[idxMass] << " MeV/c^2) from file '" << fileName << "'." << endl;
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
			ostringstream output;
			for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
				output << " " << phaseSpaceIntegrals[idxMass][idxWave];
			}
			printDebug << "phase-space integrals for wave '" << _waveNames[idxWave] << "' (" << idxWave << "):" << output.str() << endl;
		}
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::prepareMassLimits()
{
	if(_debug) {
		printDebug << "determine which mass bins to use in the fit for " << _nrMassBins << " mass bins, center of first and last mass bins: "
		           << _massBinCenters[0] << " and " << _massBinCenters[_nrMassBins - 1] << " MeV/c^2." << endl;
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
			           << "-" << (_waveMassLimits[idxWave].second<0. ? _massMax : _waveMassLimits[idxWave].second) << " MeV/c^2, "
			           << "bin range " << binFirst << "-" << binLast << endl;
		}
		_waveMassBinLimits.push_back(make_pair(binFirst, binLast));
	}

	_wavePairMassBinLimits.resize(extents[_nrWaves][_nrWaves]);
	for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
		for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
			_wavePairMassBinLimits[idxWave][jdxWave] = make_pair(max(_waveMassBinLimits[idxWave].first,  _waveMassBinLimits[jdxWave].first),
			                                                     min(_waveMassBinLimits[idxWave].second, _waveMassBinLimits[jdxWave].second));
		}
	}

	if(_debug) {
		printDebug << "waves and mass limits:" << endl;
		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			ostringstream output;
			for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
				output << _wavePairMassBinLimits[idxWave][jdxWave].first << "-" << _wavePairMassBinLimits[idxWave][jdxWave].second << " ";
			}
			printDebug << _waveNames[idxWave] << " " << _waveMassBinLimits[idxWave].first << "-" << _waveMassBinLimits[idxWave].second
			           << ": " << output.str() << endl;
		}
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::createPlots(const rpwa::massDepFit::model& fitModel,
                                          TFile* outFile,
                                          const bool rangePlotting) const
{
	if(_debug) {
		printDebug << "start creating plots." << endl;
	}

	for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
		TDirectory* outDirectory = NULL;
		if(_nrBins == 1) {
			outDirectory = outFile;
		} else {
			ostringstream name;
			name << "bin" << idxBin;
			outDirectory = outFile->mkdir(name.str().c_str());
		}

		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			if(not createPlotsWave(fitModel, outDirectory, rangePlotting, idxWave, idxBin)) {
				printErr << "error while creating intensity plots for wave '" << _waveNames[idxWave] << "' in bin " << idxBin << "." << endl;
				return false;
			}
		}

		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			for(size_t jdxWave=idxWave+1; jdxWave<_nrWaves; ++jdxWave) {
				if(not createPlotsWavePair(fitModel, outDirectory, rangePlotting, idxWave, jdxWave, idxBin)) {
					printErr << "error while creating intensity plots for wave pair '" << _waveNames[idxWave] << "' and '" << _waveNames[jdxWave] << "' in bin " << idxBin << "." << endl;
					return false;
				}
			}
		}
	}

	if (fitModel.getFsmdFunction() != NULL) {
		TF1* func = fitModel.getFsmdFunction();

		TGraph graph;
		graph.SetName("finalStateMassDependence");
		graph.SetTitle("finalStateMassDependence");
		graph.SetDrawOption("CP");

		Int_t point = -1;
		for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
			++point;
			const double mass = _massBinCenters[idxMass] * MASSSCALE;

			graph.SetPoint(point, mass, pow(func->Eval(_massBinCenters[idxMass]), 2.));
		}

		outFile->cd();
		graph.Write();
	}

	if(_debug) {
		printDebug << "finished creating plots." << endl;
	}

	return true;
}


bool
rpwa::massDepFit::massDepFit::createPlotsWave(const rpwa::massDepFit::model& fitModel,
                                              TDirectory* outDirectory,
                                              const bool rangePlotting,
                                              const size_t idxWave,
                                              const size_t idxBin) const
{
	if(_debug) {
		printDebug << "start creating plots for wave '" << _waveNames[idxWave] << "' in bin " << idxBin << "." << endl;
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

	const vector<pair<size_t, size_t> >& compChannel = fitModel.getComponentChannel(idxWave);
	vector<TGraph*> components;
	for(size_t idxComponents=0; idxComponents<compChannel.size(); ++idxComponents) {
		const size_t idxComponent = compChannel[idxComponents].first;
		TGraph* component = new TGraph;
		component->SetName((_waveNames[idxWave] + "__" + fitModel.getComponent(idxComponent)->getName()).c_str());
		component->SetTitle((_waveNames[idxWave] + "__" + fitModel.getComponent(idxComponent)->getName()).c_str());

		Color_t color = kBlue;
		if(fitModel.getComponent(idxComponent)->getName().find("bkg") != string::npos) {
			color = kMagenta;
		}
		component->SetLineColor(color);
		component->SetMarkerColor(color);

		graphs.Add(component, "CP");
		components.push_back(component);
	}

	double maxP = -numeric_limits<double>::max();
	double maxIE = -numeric_limits<double>::max();
	Int_t point = -1;
	Int_t pointLimit = -1;
	for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
		++point;
		const double mass = _massBinCenters[idxMass] * MASSSCALE;
		const double halfBin = _massStep/2000.;

		data->SetPoint(point, mass, _inIntensities[idxBin][idxMass][idxWave][0]);
		data->SetPointError(point, halfBin, _inIntensities[idxBin][idxMass][idxWave][1]);
		maxIE = max(maxIE, _inIntensities[idxBin][idxMass][idxWave][0]+_inIntensities[idxBin][idxMass][idxWave][1]);

		if(_sysPlotting) {
			double maxSI = -numeric_limits<double>::max();
			double minSI = numeric_limits<double>::max();
			for(size_t idxSystematics=0; idxSystematics<_nrSystematics; ++idxSystematics) {
				maxSI = max(maxSI, _sysIntensities[idxSystematics][idxMass][idxWave][0]);
				minSI = min(minSI, _sysIntensities[idxSystematics][idxMass][idxWave][0]);
			}
			systematics->SetPoint(point, mass, (maxSI+minSI)/2.);
			systematics->SetPointError(point, halfBin, (maxSI-minSI)/2.);
			maxIE = max(maxIE, maxSI);
		}

		const double ps = pow(_inPhaseSpaceIntegrals[idxBin][idxMass][idxWave], 2);
		phaseSpace->SetPoint(point, mass, ps);
		maxP = max(maxP, ps);

		// check that this mass bin should be taken into account for this
		// combination of waves
		if(rangePlotting && (idxMass < _wavePairMassBinLimits[idxWave][idxWave].first || idxMass > _wavePairMassBinLimits[idxWave][idxWave].second)) {
			continue;
		}
		++pointLimit;

		const double intensity = fitModel.intensity(idxWave, idxBin, _massBinCenters[idxMass], idxMass);
		fit->SetPoint(pointLimit, mass, intensity);
		maxIE = max(maxIE, intensity);

		for(size_t idxComponents=0; idxComponents<compChannel.size(); ++idxComponents) {
			const size_t idxComponent = compChannel[idxComponents].first;
			const size_t idxChannel = compChannel[idxComponents].second;

			complex<double> prodAmp = fitModel.getComponent(idxComponent)->val(idxBin, _massBinCenters[idxMass]);
			prodAmp *= fitModel.getComponent(idxComponent)->getChannel(idxChannel).getCouplingPhaseSpace(idxBin, _massBinCenters[idxMass], idxMass);
			prodAmp *= fitModel.calcFsmd(_massBinCenters[idxMass], idxMass);

			components[idxComponents]->SetPoint(pointLimit, mass, norm(prodAmp));
			maxIE = max(maxIE, norm(prodAmp));
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
                                                  TDirectory* outDirectory,
                                                  const bool rangePlotting,
                                                  const size_t idxWave,
                                                  const size_t jdxWave,
                                                  const size_t idxBin) const
{
	if(_debug) {
		printDebug << "start creating plots for wave pair '" << _waveNames[idxWave] << "' and '" << _waveNames[jdxWave] << "' in bin " << idxBin << "." << endl;
	}

	const string phaseName = _waveNames[idxWave] + "__" + _waveNames[jdxWave] + "__phase";
	const string realName = _waveNames[idxWave] + "__" + _waveNames[jdxWave] + "__real";
	const string imagName = _waveNames[idxWave] + "__" + _waveNames[jdxWave] + "__imag";

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
		const double mass = _massBinCenters[idxMass] * MASSSCALE;
		const double halfBin = _massStep/2000.;

		phaseData->SetPoint(point, mass, _inPhases[idxBin][idxMass][idxWave][jdxWave][0]);
		phaseData->SetPointError(point, halfBin, _inPhases[idxBin][idxMass][idxWave][jdxWave][1]);

		realData->SetPoint(point, mass, _inSpinDensityMatrices[idxBin][idxMass][idxWave][jdxWave].real());
		realData->SetPointError(point, halfBin, sqrt(_inSpinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][0]));

		imagData->SetPoint(point, mass, _inSpinDensityMatrices[idxBin][idxMass][idxWave][jdxWave].imag());
		imagData->SetPointError(point, halfBin, sqrt(_inSpinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][1]));

		if(_sysPlotting) {
			const double dataP = _inPhases[idxBin][idxMass][idxWave][jdxWave][0];
			double maxSP = -numeric_limits<double>::max();
			double minSP = numeric_limits<double>::max();

			double maxSR = -numeric_limits<double>::max();
			double minSR = numeric_limits<double>::max();

			double maxSI = -numeric_limits<double>::max();
			double minSI = numeric_limits<double>::max();

			for(size_t idxSystematics=0; idxSystematics<_nrSystematics; ++idxSystematics) {
				double sysP = _sysPhases[idxSystematics][idxMass][idxWave][jdxWave][0];
				if(abs(sysP+360.-dataP) < abs(sysP-dataP)) {
					sysP = sysP+360;
				} else if(abs(sysP-360.-dataP) < abs(sysP-dataP)) {
					sysP = sysP-360;
				}
				maxSP = max(maxSP, sysP);
				minSP = min(minSP, sysP);

				maxSR = max(maxSR, _sysSpinDensityMatrices[idxSystematics][idxMass][idxWave][jdxWave].real());
				minSR = min(minSR, _sysSpinDensityMatrices[idxSystematics][idxMass][idxWave][jdxWave].real());

				maxSI = max(maxSI, _sysSpinDensityMatrices[idxSystematics][idxMass][idxWave][jdxWave].imag());
				minSI = min(minSI, _sysSpinDensityMatrices[idxSystematics][idxMass][idxWave][jdxWave].imag());
			}
			phaseSystematics->SetPoint(point, mass, (maxSP+minSP)/2.);
			phaseSystematics->SetPointError(point, halfBin, (maxSP-minSP)/2.);

			realSystematics->SetPoint(point, mass, (maxSR+minSR)/2.);
			realSystematics->SetPointError(point, halfBin, (maxSR-minSR)/2.);

			imagSystematics->SetPoint(point, mass, (maxSI+minSI)/2.);
			imagSystematics->SetPointError(point, halfBin, (maxSI-minSI)/2.);
		}

		phaseFitAll.SetPoint(point, mass, fitModel.phase(idxWave, jdxWave, idxBin, _massBinCenters[idxMass], idxMass) * TMath::RadToDeg());

		// check that this mass bin should be taken into account for this
		// combination of waves
		if(rangePlotting && (idxMass < _wavePairMassBinLimits[idxWave][jdxWave].first || idxMass > _wavePairMassBinLimits[idxWave][jdxWave].second)) {
			continue;
		}
		++pointLimit;

		const complex<double> element = fitModel.spinDensityMatrix(idxWave, jdxWave, idxBin, _massBinCenters[idxMass], idxMass);
		realFit->SetPoint(pointLimit, mass, element.real());
		imagFit->SetPoint(pointLimit, mass, element.imag());
	}

	// rectify phase graphs
	point = -1;
	pointLimit = -1;
	for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
		++point;
		const double mass = _massBinCenters[idxMass] * MASSSCALE;

		double valueFit=0;
		if(point != 0) {
			int bestOffs = 0;
			double bestDiff = numeric_limits<double>::max();

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
		double bestDiff = numeric_limits<double>::max();

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


void
usage(const string& progName,
      const int     errCode = 0)
{
	cerr << "performs mass-dependent fit" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " [-o outfile -M minimizer -m algorithm -t # -P -R -C -d -q -h] config file" << endl
	     << "    where:" << endl
	     << "        -o file    path to output file (default: 'mDep.result.root')" << endl
	     << "        -M name    minimizer (default: Minuit2)" << endl
	     << "        -m name    minimization algorithm (optional, default: Migrad)" << endl
	     << "                   available minimizers: Minuit:      Migrad, Simplex, Minimize, Migrad_imp" << endl
	     << "                                         Minuit2:     Migrad, Simplex, Combined, Scan, Fumili" << endl
	     << "                                         GSLMultiMin: ConjugateFR, ConjugatePR, BFGS, BFGS2, SteepestDescent" << endl
	     << "                                         GSLMultiFit: -" << endl
	     << "                                         GSLSimAn:    -" << endl
	     << "                                         Linear:      Robust" << endl
	     << "                                         Fumili:      -" << endl
	     << "        -g #       minimizer strategy: 0 = low, 1 = medium, 2 = high effort  (default: 1)" << endl
	     << "        -t #       minimizer tolerance (default: 1e-10)" << endl
	     << "        -P         plotting only - no fit" << endl
	     << "        -R         plot in fit range only" << endl
	     << "        -C         switch OFF covariances between real and imag part" << endl
	     << "        -d         additional debug output (default: false)" << endl
	     << "        -q         run quietly (default: false)" << endl
	     << "        -h         print help" << endl
	     << endl;
	exit(errCode);
}


// changes status of variables (fixed/released)
// fixed values from config remain fixed
// parameters are taken from current status of fitter
// level
// 0 = release only couplings
// 1 = release couplings and masses
// 2 = release couplings, masses and widths
bool
releasePars(Minimizer* minimizer,
            const rpwa::massDepFit::model& compset,
            const string& freeParameters)
{
	// save current state
	const size_t npar = compset.getNrParameters();
	double par[npar];

	if(minimizer->NFree() == 0) {
		compset.getParameters(par);
	} else {
		for(size_t i=0; i<npar; ++i) {
			par[i]=minimizer->X()[i];
		}
	}

	// tokenize freeParameters string (default deparators also include '*')
	char_separator<char> separators(" ,\t\n");
	tokenizer<char_separator<char> > tokenizeFreeParameters(freeParameters, separators);

	// reset minimizer
	minimizer->Clear();

	size_t parcount=0;
	// first add all couplings
	for(size_t idxComponent=0; idxComponent<compset.getNrComponents(); ++idxComponent) {
		const rpwa::massDepFit::component* comp = compset.getComponent(idxComponent);
		for(vector<rpwa::massDepFit::channel>::const_iterator it=comp->getChannels().begin(); it!=comp->getChannels().end(); ++it) {
			for(size_t idxBin=0; idxBin<it->getNrBins(); ++idxBin) {
				ostringstream prefixName;
				prefixName << "coupling__bin"
				           << idxBin
				           << "__"
				           << comp->getName()
				           << "__"
				           << it->getWaveName();

				printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__real") << "') set to " << par[parcount] << endl;
				minimizer->SetVariable(parcount,
				                       prefixName.str() + "__real",
				                       par[parcount],
				                       0.1);
				++parcount;

				if(not it->isAnchor()) {
					printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__imag") << "') set to " << par[parcount] << endl;
					minimizer->SetVariable(parcount,
					                       prefixName.str() + "__imag",
					                       par[parcount],
					                       0.1);
					++parcount;
				}
			}
		} // end loop over channels
	} // end loop over components

	// second add parameters of the components, i.e. mass and width
	for(size_t idxComponent=0; idxComponent<compset.getNrComponents(); ++idxComponent) {
		const rpwa::massDepFit::component* comp = compset.getComponent(idxComponent);
		for(size_t idxParameter=0; idxParameter<comp->getNrParameters(); ++idxParameter) {
			const string name = comp->getName() + "__" + comp->getParameterName(idxParameter);

			bool free = false;
			if(find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "*")!=tokenizeFreeParameters.end()
			   || find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), comp->getName())!=tokenizeFreeParameters.end()
			   || find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), comp->getParameterName(idxParameter))!=tokenizeFreeParameters.end()
			   || find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), comp->getName() + "__" + comp->getParameterName(idxParameter))!=tokenizeFreeParameters.end()) {
				free = true;
			}
			bool fix = not free;
			if(comp->getParameterFixed(idxParameter)) {
				fix = true;
			}

			if(fix) {
				printInfo << "parameter " << parcount << " ('" << name << "') fixed to " << par[parcount] << endl;
				minimizer->SetFixedVariable(parcount,
				                            name,
				                            par[parcount]);
			} else if(comp->getParameterLimitedLower(idxParameter) && comp->getParameterLimitedUpper(idxParameter)) {
				printInfo << "parameter " << parcount << " ('" << name << "') set to " << par[parcount]
				          << " (limited between " << comp->getParameterLimitLower(idxParameter)
				          << " and " << comp->getParameterLimitUpper(idxParameter) << ")" << endl;
				minimizer->SetLimitedVariable(parcount,
				                              name,
				                              par[parcount],
				                              comp->getParameterStep(idxParameter),
				                              comp->getParameterLimitLower(idxParameter),
				                              comp->getParameterLimitUpper(idxParameter));
			} else if(comp->getParameterLimitedLower(idxParameter)) {
				printInfo << "parameter " << parcount << " ('" << name << "') set to " << par[parcount]
				          << " (limited larger than " << comp->getParameterLimitLower(idxParameter) << ")" << endl;
				minimizer->SetLowerLimitedVariable(parcount,
				                                   name,
				                                   par[parcount],
				                                   comp->getParameterStep(idxParameter),
				                                   comp->getParameterLimitLower(idxParameter));
			} else if(comp->getParameterLimitedUpper(idxParameter)) {
				printInfo << "parameter " << parcount << " ('" << name << "') set to " << par[parcount]
				          << " (limited smaller than " << comp->getParameterLimitUpper(idxParameter) << ")" << endl;
				minimizer->SetUpperLimitedVariable(parcount,
				                                   name,
				                                   par[parcount],
				                                   comp->getParameterStep(idxParameter),
				                                   comp->getParameterLimitUpper(idxParameter));
			} else {
				printInfo << "parameter " << parcount << " ('" << name << "') set to " << par[parcount] << endl;
				minimizer->SetVariable(parcount,
				                       name,
				                       par[parcount],
				                       comp->getParameterStep(idxParameter));
			}
			++parcount;
		}
	} // end loop over components

	// set phase space
	const unsigned int nfreeFsmd=compset.getFsmdNrParameters();
	for(unsigned int ifreeFsmd=0; ifreeFsmd<nfreeFsmd; ++ifreeFsmd) {
		const double val = par[parcount];
		double lower,upper;
		compset.getFsmdParameterLimits(ifreeFsmd, lower, upper);
		ostringstream name;
		name << "PSP__" << ifreeFsmd;
		printInfo << "parameter " << parcount << " ('" << name.str() << "') set to " << val << endl;
		minimizer->SetLimitedVariable(parcount,
		                              name.str().c_str(),
		                              val,
		                              0.0001,
		                              lower,
		                              upper);
		++parcount;
	}

	return true;
}

int
main(int    argc,
     char** argv)
{
	printCompilerInfo();
	printLibraryInfo ();
	printGitHash     ();
	cout << endl;

	// --------------------------------------------------------------------------
	// internal parameters
	const string       valTreeName           = "pwa";
	const string       valBranchName         = "fitResult_v2";
	const unsigned int maxNmbOfIterations    = 20000;
	const unsigned int maxNmbOfFunctionCalls = 2000000;
	const bool         runHesse              = true;
	const bool         runMinos              = false;

	// ---------------------------------------------------------------------------
	// parse command line options
	const string progName           = argv[0];
	string       outFileName        = "mDep.result.root";     // output filename
	string       minimizerType[2]   = {"Minuit2", "Migrad"};  // minimizer, minimization algorithm
	int          minimizerStrategy  = 1;                      // minimizer strategy
	double       minimizerTolerance = 1e-10;                  // minimizer tolerance
	bool         onlyPlotting       = false;
	bool         rangePlotting      = false;
	bool         doCov              = true;
	bool         debug              = false;
	bool         quiet              = false;
	extern char* optarg;
	extern int   optind;
	int c;
	while ((c = getopt(argc, argv, "o:M:m:g:t:PRCdqh")) != -1) {
		switch (c) {
		case 'o':
			outFileName = optarg;
			break;
		case 'M':
			minimizerType[0] = optarg;
			break;
		case 'm':
			minimizerType[1] = optarg;
			break;
		case 'g':
			minimizerStrategy = atoi(optarg);
			break;
		case 't':
			minimizerTolerance = atof(optarg);
			break;
		case 'P':
			onlyPlotting=true;
			break;
		case 'R':
			rangePlotting=true;
			break;
		case 'C':
			doCov=false;
			break;
		case 'd':
			debug = true;
			break;
		case 'q':
			quiet = true;
			break;
		case '?':
		case 'h':
			usage(progName, 1);
			break;
		}
	}

	// there must only be one remaining (unhandled) argument which is the
	// configuration file
	if(optind+1 != argc) {
		printErr << "you need to specify exactly one configuration file." << endl;
		usage(progName, 1);
	}
	const string configFileName = argv[optind];

	// report parameters
	printInfo << "running " << progName << " with the following parameters:" << endl
	          << "    path to configuration file ..................... '" << configFileName << "'" << endl
	          << "    path to output file ............................ '" << outFileName << "'" << endl
	          << "    minimizer ...................................... "  << minimizerType[0] << ", " << minimizerType[1] << endl
	          << "    minimizer strategy ............................. "  << minimizerStrategy  << endl
	          << "    minimizer tolerance ............................ "  << minimizerTolerance << endl
	          << "    only plotting .................................. "  << yesNo(onlyPlotting) << endl
	          << "    plot in fit range only ......................... "  << yesNo(rangePlotting) << endl
	          << "    take covariance into account ................... "  << yesNo(doCov) << endl
	          << "    debug .......................................... "  << yesNo(debug) << endl
	          << "    quiet .......................................... "  << yesNo(quiet) << endl;

	rpwa::massDepFit::massDepFit mdepFit;
	mdepFit.setDebug(debug);

	Config configFile;
	if(not parseLibConfigFile(configFileName, configFile, debug)) {
		printErr << "could not read configuration file '" << configFileName << "'." << endl;
		return 1;
	}
	Setting& configRoot = configFile.getRoot();

	// input section
	const Setting* configInput = findLibConfigGroup(configRoot, "input");
	if(not configInput) {
		printErr << "'input' section in configuration file does not exist." << endl;
		return 1;
	}
	if(not mdepFit.readConfigInput(configInput)) {
		printErr << "error while reading 'input' section from configuration file." << endl;
		return 1;
	}

	// extract information from fit results
	if(not mdepFit.readInFiles(valTreeName, valBranchName)) {
		printErr << "error while trying to read fit result." << endl;
		return 1;
	}

	// extract information for systematic errors
	if(not mdepFit.readSystematicsFiles(valTreeName, valBranchName)) {
		printErr << "error while trying to read fit results for systematic errors." << endl;
		return 1;
	}

	// prepare mass limits
	if(not mdepFit.prepareMassLimits()) {
		printErr << "error determine which bins to use in the fit." << endl;
		return 1;
	}

	// set-up fit model (resonances, background, final-state mass dependence
	rpwa::massDepFit::model compset;
	if(not mdepFit.readConfigModel(&configRoot, compset)) {
		printErr << "error while reading fit model from configuration file." << endl;
		return 1;
	}

	if(not compset.init(mdepFit.getWaveNames(), mdepFit.getMassBinCenters(), mdepFit.getAnchorWaveName(), mdepFit.getAnchorComponentName())) {
		printErr << "error while initializing the fit model." << endl;
		return 1;
	}

	if(onlyPlotting) {
		printInfo << "plotting only mode, skipping minimzation." << endl;
	} else {
		// set-up likelihood
		rpwa::massDepFit::likelihood L;
		L.init(&compset,
		       mdepFit.getMassBinCenters(),
		       mdepFit.getInSpinDensityMatrices(),
		       mdepFit.getInSpinDensityCovarianceMatrices(),
		       mdepFit.getWavePairMassBinLimits(),
		       doCov);

		// setup minimizer
		printInfo << "creating and setting up minimizer '" << minimizerType[0] << "' "
		          << "using algorithm '" << minimizerType[1] << "'" << endl;
		Minimizer* minimizer = Factory::CreateMinimizer(minimizerType[0], minimizerType[1]);
		if(not minimizer) {
			printErr << "could not create minimizer. exiting." << endl;
			return 1;
		}
		minimizer->SetFunction        (L);
		minimizer->SetStrategy        (minimizerStrategy);
		minimizer->SetTolerance       (minimizerTolerance);
		minimizer->SetPrintLevel      ((quiet) ? 0 : 3);
		minimizer->SetMaxIterations   (maxNmbOfIterations);
		minimizer->SetMaxFunctionCalls(maxNmbOfFunctionCalls);

		// keep list of parameters to free
		const vector<string> freeParameters = mdepFit.getFreeParameters();
		const size_t nrSteps = freeParameters.size();

		bool success = true;
		for(size_t step=0; step<nrSteps; ++step) {
			// set startvalues
			if(not releasePars(minimizer, compset, freeParameters[step])) {
				printErr << "error while setting start parameters for step " << step << "." << endl;
				return 1;
			}

			// keep track of the time spend in the fit
			TStopwatch stopwatch;

			printInfo << "performing minimization step " << step << ": '" << freeParameters[step] << "' (" << minimizer->NFree() << " free parameters)." << endl;
			stopwatch.Start();
			success &= minimizer->Minimize();
			stopwatch.Stop();

			if(not success) {
				printWarn << "minimization failed." << endl;
			} else {
				printInfo << "minimization successful." << endl;
			}
			printInfo << "minimization took " <<  maxPrecisionAlign(stopwatch.CpuTime()) << " s" << endl;
		}

		if(runHesse) {
			TStopwatch stopwatch;

			printInfo << "calculating Hessian matrix." << endl;
			stopwatch.Start();
			success &= minimizer->Hesse();
			stopwatch.Stop();

			if(not success) {
				printWarn << "calculation of Hessian matrix failed." << endl;
			} else {
				printInfo << "calculation of Hessian matrix successful." << endl;
			}
			printInfo << "calculating Hessian matrix took " <<  maxPrecisionAlign(stopwatch.CpuTime()) << " s" << endl;
		}

		printInfo << "minimizer status summary:" << endl
		          << "    total number of parameters .......................... " << minimizer->NDim()             << endl
		          << "    number of free parameters ........................... " << minimizer->NFree()            << endl
		          << "    maximum allowed number of iterations ................ " << minimizer->MaxIterations()    << endl
		          << "    maximum allowed number of function calls ............ " << minimizer->MaxFunctionCalls() << endl
		          << "    minimizer status .................................... " << minimizer->Status()           << endl
		          << "    minimizer provides error and error matrix ........... " << minimizer->ProvidesError()    << endl
		          << "    minimizer has performed detailed error validation ... " << minimizer->IsValidError()     << endl
		          << "    estimated distance to minimum ....................... " << minimizer->Edm()              << endl
		          << "    statistical scale used for error calculation ........ " << minimizer->ErrorDef()         << endl
		          << "    minimizer strategy .................................. " << minimizer->Strategy()         << endl
		          << "    absolute tolerance .................................. " << minimizer->Tolerance()        << endl;

		// print results
		ostringstream output;
		const unsigned int nmbPar  = L.NDim();
		for (unsigned int i = 0; i<nmbPar; ++i) {
			output << "    parameter [" << setw(3) << i << "] ";
			output << minimizer->VariableName(i) << " " ;
			output << maxPrecisionAlign(minimizer->X()[i]) << " +- " << maxPrecisionAlign(minimizer->Errors()[i]);

			if (runMinos) {  // does not work for all parameters
				double minosErrLow, minosErrUp;
				if (minimizer->GetMinosError(i, minosErrLow, minosErrUp)) {
					output << "    Minos: " << "[" << minosErrLow << ", +" << minosErrUp << "]";
				}
			}
			output << endl;
		}
		printInfo << "minimization result:" << endl
		          << output.str();

		const double* par=minimizer->X();
		double chi2 = 0.;
		if(success) {
			chi2 = L.DoEval(par);
		}
		printInfo << "chi2 = " << maxPrecisionAlign(chi2) << endl;
		compset.setParameters(par);

		const unsigned int nrDataPoints = L.NDataPoints();
		const unsigned int nrFree = minimizer->NFree();
		const unsigned int ndf = nrDataPoints - nrFree;
		printInfo << "ndf = " << nrDataPoints << "-" << nrFree << "=" << ndf << endl;

		double chi2red = chi2/(double)ndf;
		printInfo << "chi2/ndf = " << maxPrecisionAlign(chi2red) << endl;

		if(not mdepFit.updateConfig(&configRoot, compset, minimizer, chi2, ndf, chi2red)) {
			printErr << "error while updating configuration file." << endl;
			return 1;
		}
	}

	string confFileName(outFileName);
	if(extensionFromPath(confFileName) == "root") {
		confFileName = changeFileExtension(confFileName, "conf");
	} else {
		confFileName += ".conf";
	}

	if(debug) {
		printDebug << "name of output configuration file: '" << confFileName << "'." << endl;
	}
	configFile.writeFile(confFileName.c_str());

	string rootFileName(outFileName);
	if(extensionFromPath(rootFileName) != "root") {
		rootFileName += ".root";
	}

	if(debug) {
		printDebug << "name of output ROOT file: '" << confFileName << "'." << endl;
	}
	TFile* outFile = TFile::Open(rootFileName.c_str(), "RECREATE");
	if(not mdepFit.createPlots(compset, outFile, rangePlotting)) {
		printErr << "error while creating plots." << endl;
		return 1;
	}
	outFile->Close();

	return 0;
}
