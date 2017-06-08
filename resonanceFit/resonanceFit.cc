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

#include "data.h"
#include "input.h"
#include "resonanceFitHelper.h"


namespace {


	// debug flag for functions in this (anonymous) namespace
	bool debug = false;


	std::vector<std::string>
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
		if(not configAnchors.IsSequence()) {
			printErr << "'anchorwave' is not a YAML sequence." << std::endl;
			throw;
		}

		const size_t nrBins = configAnchors.size();
		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading anchor waves and components for " << nrBins << " bins from configuration file." << std::endl;
		}

		std::vector<std::string> anchorWaveNames(nrBins);
		for(size_t idxBin = 0; idxBin < nrBins; ++idxBin) {
			const YAML::Node& configAnchor = configAnchors[idxBin];

			if(configAnchor.size() != 2) {
				printErr << "entry " << idxBin << " of 'anchorwave' needs to contain exactly two elements." << std::endl;
				throw;
			}
			if(not checkVariableType(configAnchor[0], rpwa::YamlCppUtils::TypeString) or
			   not checkVariableType(configAnchor[1], rpwa::YamlCppUtils::TypeString)) {
				printErr << "elements of entry " << idxBin << " in 'anchorwave' must be 'string'." << std::endl;
				throw;
			}

			anchorWaveNames[idxBin] = configAnchor[0].as<std::string>();

			if(rpwa::resonanceFit::debug()) {
				printDebug << "read anchor wave '" << anchorWaveNames[idxBin] << "' for entry " << idxBin << "." << std::endl;
			}
		}

		return anchorWaveNames;
	}


}


void
rpwa::resonanceFit::read(const YAML::Node& configRoot,
                         rpwa::resonanceFit::inputConstPtr& fitInput,
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
	fitInput = rpwa::resonanceFit::readInput(configRoot);
	if(not fitInput) {
		printErr << "error while reading 'input' in configuration file." << std::endl;
		throw;
	}

	const std::vector<std::string> anchorWaveNames = readAnchorWaveName(configRoot);
	if(fitInput->nrBins() != anchorWaveNames.size()) {
		printErr << "expected to find " << fitInput->nrBins() << " entries in 'anchorwave', found " << anchorWaveNames.size() << " instead." << std::endl;
		throw;
	}

	fitData = rpwa::resonanceFit::readData(fitInput,
	                                       anchorWaveNames,
	                                       useCovariance,
	                                       valTreeName,
	                                       valBranchName);

	// set-up fit model (resonances, background, final-state mass-dependence)
	fitModel = rpwa::resonanceFit::readModel(configRoot,
	                                         fitInput,
	                                         fitData,
	                                         fitParameters,
	                                         fitParametersError,
	                                         useBranchings);
}


void
rpwa::resonanceFit::read(const YAML::Node& configRoot,
                         const double maxMassBinCenter,
                         rpwa::resonanceFit::inputConstPtr& fitInput,
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
	fitInput = rpwa::resonanceFit::readInput(configRoot);
	if(not fitInput) {
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
	std::vector<size_t> nrWaves(fitInput->nrBins());
	boost::multi_array<std::string, 2> waveNames(boost::extents[fitInput->nrBins()][std::max_element(fitInput->bins().begin(), fitInput->bins().end(), [](const rpwa::resonanceFit::input::bin& binMax, const rpwa::resonanceFit::input::bin& bin){ return binMax.nrWaves() < bin.nrWaves(); })->nrWaves()]);
	std::vector<size_t> nrMassBins(fitInput->nrBins());
	boost::multi_array<double, 2> massBinCenters(boost::extents[fitInput->nrBins()][2]);
	boost::multi_array<double, 3> phaseSpaceIntegrals(boost::extents[fitInput->nrBins()][2][std::max_element(fitInput->bins().begin(), fitInput->bins().end(), [](const rpwa::resonanceFit::input::bin& binMax, const rpwa::resonanceFit::input::bin& bin){ return binMax.nrWaves() < bin.nrWaves(); })->nrWaves()]);
	for(size_t idxBin = 0; idxBin < fitInput->nrBins(); ++idxBin) {
		const rpwa::resonanceFit::input::bin& fitInputBin = fitInput->getBin(idxBin);

		nrWaves[idxBin] = fitInputBin.nrWaves();
		nrMassBins[idxBin] = 2;

		massBinCenters[idxBin][0] = 0.0;
		massBinCenters[idxBin][1] = maxMassBinCenter;

		for(size_t idxWave = 0; idxWave < fitInputBin.nrWaves(); ++idxWave) {
			waveNames[idxBin][idxWave] = fitInputBin.getWave(idxWave).waveName();

			phaseSpaceIntegrals[idxBin][0][idxWave] = 0.0;
			phaseSpaceIntegrals[idxBin][1][idxWave] = 1.0;
		}
	}
	const rpwa::resonanceFit::baseDataConstPtr fitData(new baseData(nrWaves,
	                                                                waveNames,
	                                                                nrMassBins,
	                                                                massBinCenters,
	                                                                phaseSpaceIntegrals));

	// set-up fit model (resonances, background, final-state mass-dependence)
	fitModel = rpwa::resonanceFit::readModel(configRoot,
	                                         fitInput,
	                                         fitData,
	                                         fitParameters,
	                                         fitParametersError,
	                                         useBranchings);
}


void
rpwa::resonanceFit::read(const std::string& configFileName,
                         rpwa::resonanceFit::inputConstPtr& fitInput,
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
	     fitInput,
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
                         rpwa::resonanceFit::inputConstPtr& fitInput,
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
	     fitInput,
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
