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
#include <map>
#include <string>

#include <boost/assign/std/vector.hpp>

#include <yaml-cpp/yaml.h>

#include <reportingUtils.hpp>
#include <yamlCppUtils.hpp>

#include "components.h"
#include "data.h"
#include "information.h"
#include "model.h"


namespace {


	// debug flag for functions in this (anonymous) namespace
	bool debug = false;


	std::string
	readAnchorWaveName(const YAML::Node& configRoot)
	{
		if(debug) {
			printDebug << "reading name of wave used as anchor." << std::endl;
		}

		const YAML::Node& configModel = configRoot["model"];
		if(not configModel) {
			printErr << "'model' is not a valid YAML node." << std::endl;
			throw;
		}
		if(not configModel.IsMap()) {
			printErr << "'model' is not a YAML map." << std::endl;
			throw;
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
		                     ("name", rpwa::YamlCppUtils::TypeString);
		if(not checkIfAllVariablesAreThere(configAnchors, mandatoryArguments)) {
			printErr << "'anchorwave' does not contain all required variables." << std::endl;
			throw;
		}

		return configAnchors["name"].as<std::string>();
	}


}


void
rpwa::resonanceFit::read(const YAML::Node& configRoot,
                         rpwa::resonanceFit::informationConstPtr& fitInformation,
                         rpwa::resonanceFit::dataConstPtr& fitData,
                         rpwa::resonanceFit::modelConstPtr& fitModel,
                         rpwa::resonanceFit::parameters& fitParameters,
                         rpwa::resonanceFit::parameters& fitParametersError,
                         std::map<std::string, double>& fitQuality,
                         std::vector<std::string>& freeParameters,
                         const bool useBranchings,
                         const rpwa::resonanceFit::function::useCovarianceMatrix useCovariance,
                         const std::string& valTreeName,
                         const std::string& valBranchName)
{
	// get information of fit quality if a previous fit was stored in the
	// configuration file
	fitQuality = rpwa::resonanceFit::readFitQuality(configRoot);

	// get information for which parameters to release in which order
	freeParameters = rpwa::resonanceFit::readFreeParameters(configRoot);

	// get fit results and waves to use in the resonance fit
	fitInformation = rpwa::resonanceFit::readInformation(configRoot);
	if(not fitInformation) {
		printErr << "error while reading 'input' in configuration file." << std::endl;
		throw;
	}

	fitData = rpwa::resonanceFit::readData(fitInformation,
	                                       readAnchorWaveName(configRoot),
	                                       useCovariance,
	                                       valTreeName,
	                                       valBranchName);

	// set-up fit model (resonances, background, final-state mass-dependence)
	fitModel = rpwa::resonanceFit::readModel(configRoot,
	                                         fitInformation,
	                                         fitData,
	                                         fitParameters,
	                                         fitParametersError,
	                                         useBranchings);

	std::ostringstream output;
	for(size_t idxComponent = 0; idxComponent < fitModel->getNrComponents(); ++idxComponent) {
		output << "    " << fitModel->getComponent(idxComponent)->getName() << std::endl;
	}
	printInfo << "fitting " << fitModel->getNrComponents() << " components to the data:" << std::endl
	          << output.str();
	if(fitModel->getFsmd()) {
		printInfo << "using final-state mass-dependence." << std::endl;
	} else {
		printInfo << "not using final-state mass-dependence." << std::endl;
	}
}


void
rpwa::resonanceFit::read(const YAML::Node& configRoot,
                         const double maxMassBinCenter,
                         rpwa::resonanceFit::informationConstPtr& fitInformation,
                         rpwa::resonanceFit::modelConstPtr& fitModel,
                         rpwa::resonanceFit::parameters& fitParameters,
                         rpwa::resonanceFit::parameters& fitParametersError,
                         std::map<std::string, double>& fitQuality,
                         std::vector<std::string>& freeParameters,
                         const bool useBranchings)
{
	// get information of fit quality if a previous fit was stored in the
	// configuration file
	fitQuality = rpwa::resonanceFit::readFitQuality(configRoot);

	// get information for which parameters to release in which order
	freeParameters = rpwa::resonanceFit::readFreeParameters(configRoot);

	// get fit results and waves to use in the resonance fit
	fitInformation = rpwa::resonanceFit::readInformation(configRoot);
	if(not fitInformation) {
		printErr << "error while reading 'input' in configuration file." << std::endl;
		throw;
	}

	// create a fake data object containing the information required to
	// initialize the model, that is:
	// - wave names used in each bin
	//   this cannot simply be done, for the moment this works only if no
	//   alternative wave names are given
	// - massBinCenters
	// - phaseSpaceIntegrals
	std::vector<size_t> nrMassBins(fitInformation->nrBins());
	boost::multi_array<double, 2> massBinCenters(boost::extents[fitInformation->nrBins()][2]);
	boost::multi_array<std::string, 2> waveNames(boost::extents[fitInformation->nrBins()][fitInformation->nrWaves()]);
	boost::multi_array<double, 3> phaseSpaceIntegrals(boost::extents[fitInformation->nrBins()][2][fitInformation->nrWaves()]);
	for(size_t idxBin = 0; idxBin < fitInformation->nrBins(); ++idxBin) {
		nrMassBins[idxBin] = 2;

		massBinCenters[idxBin][0] = 0.0;
		massBinCenters[idxBin][1] = maxMassBinCenter;

		for(size_t idxWave = 0; idxWave < fitInformation->nrWaves(); ++idxWave) {
			if(fitInformation->getWave(idxWave).waveNameAlternatives().size() != 0) {
				printErr << "reading a configuration file without reading data does only work if no alternative wave names are used. Aborting..." << std::endl;
				throw;
			}
			waveNames[idxBin][idxWave] = fitInformation->getWave(idxWave).waveName();

			phaseSpaceIntegrals[idxBin][0][idxWave] = 0.0;
			phaseSpaceIntegrals[idxBin][1][idxWave] = 1.0;
		}
	}
	const rpwa::resonanceFit::baseDataConstPtr fitData(new baseData(nrMassBins,
	                                                                massBinCenters,
	                                                                waveNames,
	                                                                phaseSpaceIntegrals));

	// set-up fit model (resonances, background, final-state mass-dependence)
	fitModel = rpwa::resonanceFit::readModel(configRoot,
	                                         fitInformation,
	                                         fitData,
	                                         fitParameters,
	                                         fitParametersError,
	                                         useBranchings);
}


void
rpwa::resonanceFit::read(const std::string& configFileName,
                         rpwa::resonanceFit::informationConstPtr& fitInformation,
                         rpwa::resonanceFit::dataConstPtr& fitData,
                         rpwa::resonanceFit::modelConstPtr& fitModel,
                         rpwa::resonanceFit::parameters& fitParameters,
                         rpwa::resonanceFit::parameters& fitParametersError,
                         std::map<std::string, double>& fitQuality,
                         std::vector<std::string>& freeParameters,
                         const bool useBranchings,
                         const rpwa::resonanceFit::function::useCovarianceMatrix useCovariance,
                         const std::string& valTreeName,
                         const std::string& valBranchName)
{
	YAML::Node configRoot;
	if(not rpwa::YamlCppUtils::parseYamlFile(configFileName, configRoot, rpwa::resonanceFit::debug())) {
		printErr << "could not read configuration file '" << configFileName << "'." << std::endl;
		throw;
	}

	read(configRoot,
	     fitInformation,
	     fitData,
	     fitModel,
	     fitParameters,
	     fitParametersError,
	     fitQuality,
	     freeParameters,
	     useBranchings,
	     useCovariance,
	     valTreeName,
	     valBranchName);
}


void
rpwa::resonanceFit::read(const std::string& configFileName,
                         const double maxMassBinCenter,
                         rpwa::resonanceFit::informationConstPtr& fitInformation,
                         rpwa::resonanceFit::modelConstPtr& fitModel,
                         rpwa::resonanceFit::parameters& fitParameters,
                         rpwa::resonanceFit::parameters& fitParametersError,
                         std::map<std::string, double>& fitQuality,
                         std::vector<std::string>& freeParameters,
                         const bool useBranchings)
{
	YAML::Node configRoot;
	if(not rpwa::YamlCppUtils::parseYamlFile(configFileName, configRoot, rpwa::resonanceFit::debug())) {
		printErr << "could not read configuration file '" << configFileName << "'." << std::endl;
		throw;
	}

	read(configRoot,
	     maxMassBinCenter,
	     fitInformation,
	     fitModel,
	     fitParameters,
	     fitParametersError,
	     fitQuality,
	     freeParameters,
	     useBranchings);
}


bool
rpwa::resonanceFit::debug()
{
	return ::debug;
}


void
rpwa::resonanceFit::setDebug(const bool newDebug)
{
	::debug = newDebug;
}
