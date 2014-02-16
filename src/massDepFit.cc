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
//      minimizes massDepFitLikeli function
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <algorithm>
#include <cassert>
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>
#include <string>

#include <boost/assign/std/vector.hpp>

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
#include <Math/Minimizer.h>
#include <Math/Factory.h>

#include <libconfig.h++>

#include "fitResult.h"
#include "libConfigUtils.hpp"
#include "massDepFit.h"
#include "massDepFitLikeli.h"
#include "massDepFitModel.h"
#include "pwacomponent.h"
#include "reportingUtils.hpp"
#include "reportingUtilsEnvironment.h"


#define MASSSCALE 0.001


using namespace std;
using namespace libconfig;
using namespace ROOT::Math;
using namespace rpwa;
using namespace boost;
using namespace boost::assign;


bool massDepFit::_debug = false;


massDepFit::massDepFit()
	: _sysPlotting(false),
	  _nrMassBins(0),
	  _nrSystematics(0),
	  _nrWaves(0)
{
}


bool
massDepFit::readConfigInput(const Setting* configInput)
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

	return true;
}


bool
massDepFit::readConfigInputFitResults(const Setting* configInputFitResults)
{
	if(not configInputFitResults) {
		printErr << "'configInputFitResults' is not a pointer to a valid object." << endl;
		return false;
	}

	const int nrFitResults = configInputFitResults->getLength();
	if(nrFitResults != 1) {
		printErr << "handling of more than one entry in 'fitresults' not yet supported." << endl;
		return false;
	}

	const Setting* configInputFitResult = &((*configInputFitResults)[0]);

	map<string, Setting::Type> mandatoryArguments;
	insert(mandatoryArguments)
	    ("name", Setting::TypeString);
	if(not checkIfAllVariablesAreThere(configInputFitResult, mandatoryArguments)) {
		printErr << "'fitresults' list in 'input' section in configuration file contains errors." << endl;
		return false;
	}

	configInputFitResult->lookupValue("name", _inFileName);

	if(_debug) {
		printDebug << "read file name of fit results of mass-independent fit: '" << _inFileName << "'." << endl;
	}

	return true;
}


bool
massDepFit::readConfigInputWaves(const Setting* configInputWaves)
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

	return true;
}


bool
massDepFit::readConfigInputSystematics(const Setting* configInputSystematics)
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

	for(int idxSystematic=0; idxSystematic<nrSystematics; ++idxSystematic) {
		const string fileName = (*configInputSystematics)[idxSystematic];
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
massDepFit::readInFile(const string& valTreeName,
                       const string& valBranchName)
{
	if(_debug) {
		printDebug << "reading fit result from file '" << _inFileName << "'." << endl;
	}

	TFile* inFile = TFile::Open(_inFileName.c_str());
	if(not inFile) {
		printErr << "input file '" << _inFileName << "' not found."<< endl;
		return false;
	}
	if(inFile->IsZombie()) {
		printErr << "error while reading input file '" << _inFileName << "'."<< endl;
		delete inFile;
		return false;
	}

	if(_debug) {
		printDebug << "searching for tree '" << valTreeName << "' in file '" << _inFileName << "'." << endl;
	}

	TTree* inTree;
	inFile->GetObject(valTreeName.c_str(), inTree);
	if(not inTree) {
		printErr << "input tree '" << valTreeName << "' not found in input file '" << _inFileName << "'."<< endl;
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
		printErr << "could not extract mass bins from fit result tree in '" << _inFileName << "'." << endl;
		delete inFile;
		return false;
	}

	vector<Long64_t> inMapping;
	if(not checkFitResultMassBins(inTree, inFit, inMapping)) {
		printErr << "error while checking and mapping mass bins from fit result tree in '" << _inFileName << "'." << endl;
		delete inFile;
		return false;
	}

	if(not readFitResultMatrices(inTree, inFit, inMapping, _inSpinDensityMatrices, _inSpinDensityCovarianceMatrices,
	                             _inIntensities, _inPhases)) {
		printErr << "error while reading spin-density matrix from fit result tree in '" << _inFileName << "'." << endl;
		delete inFile;
		return false;
	}
	if(not readFitResultIntegrals(inTree, inFit, inMapping, _inPhaseSpaceIntegrals)) {
		printErr << "error while reading phase-space integrals from fit result tree in '" << _inFileName << "'." << endl;
		delete inFile;
		return false;
	}

	delete inFile;
	return true;
}


bool
massDepFit::readSystematicsFiles(const string& valTreeName,
                                 const string& valBranchName)
{
	if(not _sysPlotting) {
		return true;
	}

	if(_debug) {
		printDebug << "reading fit results for systematic errors from " << _nrSystematics << " files." << endl;
	}

	_sysSpinDensityMatrices.resize(extents[_nrSystematics][_nrMassBins][_nrWaves][_nrWaves]);
	_sysSpinDensityCovarianceMatrices.resize(extents[_nrSystematics][_nrMassBins][_nrWaves][_nrWaves][2][2]);
	_sysIntensities.resize(extents[_nrSystematics][_nrMassBins][_nrWaves][2]);
	_sysPhases.resize(extents[_nrSystematics][_nrMassBins][_nrWaves][_nrWaves][2]);

	_sysSpinDensityMatrices[0] = _inSpinDensityMatrices;
	_sysSpinDensityCovarianceMatrices[0] = _inSpinDensityCovarianceMatrices;
	_sysIntensities[0] = _inIntensities;
	_sysPhases[0] = _inPhases;

	for(size_t idxSystematics=1; idxSystematics<_nrSystematics; ++idxSystematics) {
		readSystematicsFile(idxSystematics, valTreeName, valBranchName);
	}

	return true;
}


bool
massDepFit::readSystematicsFile(const size_t idxSystematics,
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
	boost::multi_array<double, 3> tempIntensities;
	boost::multi_array<double, 4> tempPhases;
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
massDepFit::checkFitResultMassBins(TTree* tree,
                                   rpwa::fitResult* fit,
                                   vector<Long64_t>& mapping) const
{
	if(not tree or not fit) {
		printErr << "'tree' or 'fit' is not a pointer to a valid object." << endl;
		return false;
	}

	// reset mapping
	mapping.assign(_nrMassBins, numeric_limits<size_t>::max());

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

		if(mapping[idxMass] != numeric_limits<size_t>::max()) {
			printErr << "cannat map tree entry " << idx << " to mass bin " << idxMass << " (" << _massBinCenters[idxMass] << " MeV/c^2)  "
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
		if(mapping[idx] == numeric_limits<size_t>::max()) {
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
massDepFit::readFitResultMassBins(TTree* tree,
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

	return true;
}


bool
massDepFit::readFitResultMatrices(TTree* tree,
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
massDepFit::readFitResultIntegrals(TTree* tree,
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
			phaseSpaceIntegrals[idxMass][idxWave] = ps*ps;
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
massDepFit::prepareMassLimits()
{
	if(_debug) {
		printDebug << "determine which mass bins to use in the fit for " << _nrMassBins << " mass bins, center of first and last mass bins: "
		           << _massBinCenters[0] << " and " << _massBinCenters[_nrMassBins - 1] << " MeV/c^2." << endl;
	}

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
void releasePars(Minimizer* minimizer, const pwacompset& compset,
		 const vector<string>& anchorwave_reso,
		 const vector<string>& anchorwave_channel,
		 int level){
  // copy state
  unsigned int npar=minimizer->NDim();
  double par[npar];
  for(unsigned int i=0;i<npar;++i)par[i]=minimizer->X()[i];
  minimizer->Clear();

  unsigned int parcount=0;
  for(unsigned int ic=0;ic<compset.n();++ic){
    const pwacomponent& comp=*compset[ic];
    TString name(comp.name());
    double mmin,mmax,gmin,gmax;
    comp.getLimits(mmin,mmax,gmin,gmax);
    if(comp.fixM() || level==0)minimizer->SetFixedVariable(parcount,
					       (name+"_M").Data() ,
					       par[parcount]);
    else minimizer->SetLimitedVariable(parcount,
				       (name+"_M").Data(),
				       par[parcount],
				       5.0,
				       mmin,mmax);
    if(level==0 && !comp.fixM()) printInfo << minimizer->VariableName(parcount)
			   << " fixed to " << par[parcount] << endl;
    ++parcount;
    if(comp.fixGamma() || level < 2)minimizer->SetFixedVariable(parcount,
						   (name+"_Gamma").Data() ,
						    par[parcount]);
    else minimizer->SetLimitedVariable(parcount,
				       (name+"_Gamma").Data(),
				       par[parcount],
				       5.0,
				       gmin,gmax);
    if(level<2 && !comp.fixGamma()) printInfo << minimizer->VariableName(parcount)
			  << " fixed to " << par[parcount] << endl;
    ++parcount;

    std::map<std::string,pwachannel >::const_iterator it=comp.channels().begin();
    while(it!=comp.channels().end()){
      minimizer->SetVariable(parcount,(name + "_ReC" + it->first).Data() , par[parcount], 10.0);
      ++parcount;
      // fix one phase
      if(find(anchorwave_reso.begin(),anchorwave_reso.end(),name)!=anchorwave_reso.end() && find(anchorwave_channel.begin(),anchorwave_channel.end(),it->first)!=anchorwave_channel.end()){
	minimizer->SetFixedVariable(parcount,(name + "_ImC" + it->first).Data() , 0.0);
      }
      else {minimizer->SetVariable(parcount,(name + "_ImC" + it->first).Data() , par[parcount], 0.10);}
      ++parcount;
      ++it;
    } // end loop over channels
  }// end loop over components
  // set phase space
  unsigned int nfreeFsmd=compset.nFreeFsmdPar();
  for(unsigned int ifreeFsmd=0;ifreeFsmd<nfreeFsmd;++ifreeFsmd){
    double val,lower,upper;
    val=par[parcount];
    compset.getFreeFsmdLimits(ifreeFsmd,lower,upper);
    TString name("PSP_"); name+=+ifreeFsmd;
    minimizer->SetLimitedVariable(parcount, 
				  name.Data(), 
				  val, 0.0001 ,lower,upper);
  }



  const unsigned int nfree=minimizer->NFree();
  printInfo <<  nfree  << " Free Parameters in fit" << endl;


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
	const unsigned int maxNmbOfFunctionCalls = 40000;
	const bool         runHesse              = true;
	const bool         runMinos              = false;

	// ---------------------------------------------------------------------------
	// parse command line options
	const string progName           = argv[0];
	bool         doCov              = true;
	string       outFileName        = "mDep.result.root";     // output filename
	string       minimizerType[2]   = {"Minuit2", "Migrad"};  // minimizer, minimization algorithm
	int          minimizerStrategy  = 1;                      // minimizer strategy
	double       minimizerTolerance = 1e-10;                  // minimizer tolerance
	bool         onlyPlotting       = false;
	bool         rangePlotting      = false;
	bool         debug              = false;
	bool         quiet              = false;
	extern char* optarg;
	extern int   optind;
	int c;
	while ((c = getopt(argc, argv, "o:M:m:g:t:PRCdqh")) != -1)
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

	// there must only be one remaining (unhandled) argument which is the
	// configuration file
	if(optind+1 != argc) {
		printErr << "you need to specify exactly one configuration file." << endl;
		usage(progName, 1);
	}
	const string configFileName = argv[optind];

	massDepFit mdepFit;
	mdepFit.setDebug(debug);

	Config configFile;
	if(not parseLibConfigFile(configFileName, configFile, debug)) {
		printErr << "could not read configuration file '" << configFileName << "'." << endl;
		exit(1);
	}
	const Setting& configRoot = configFile.getRoot();

	// input section
	const Setting* configInput = findLibConfigGroup(configRoot, "input");
	if(not configInput) {
		printErr << "'input' section in configuration file does not exist." << endl;
		exit(1);
	}
	if(not mdepFit.readConfigInput(configInput)) {
		printErr << "error while reading 'input' section from configuration file." << endl;
		exit(1);
	}

	// extract information from fit results
	if(not mdepFit.readInFile(valTreeName, valBranchName)) {
		printErr << "error while trying to read fit result." << endl;
		exit(1);
	}

	// extract information for systematic errors
	if(not mdepFit.readSystematicsFiles(valTreeName, valBranchName)) {
		printErr << "error while trying to read fit results for systematic errors." << endl;
		exit(1);
	}

	// prepare mass limits
	if(not mdepFit.prepareMassLimits()) {
		printErr << "error determine which bins to use in the fit." << endl;
		exit(1);
	}

  printInfo << "creating and setting up likelihood function" << endl;
  printInfo << "doCovariances = " << doCov << endl;

  
  // Setup Component Set (Resonances + Background)
  pwacompset compset;
  bool check=true;

	compset.setWaveList(mdepFit.getWaveNames());

  // overall final-state mass dependence
  TF1* fPS = NULL;
  const Setting* configFsmd = findLibConfigGroup(configRoot, "finalStateMassDependence", false);
  if(configFsmd){
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
      exit(1);
    }
    
    string formula;
    configFsmd->lookupValue("formula", formula);

    fPS = new TF1("fps", formula.c_str(), 900, 3000);
    const unsigned int nrPar = fPS->GetNpar();

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

    for (unsigned int i=0; i<nrPar; ++i) {
      fPS->SetParameter(i, configFsmdValue[i]);
      fPS->SetParError(i, configFsmdError[i]);
      fPS->SetParLimits(i, configFsmdLower[i], configFsmdUpper[i]);

      if (((int)configFsmdFix[i]) != 0) {
        fPS->FixParameter(i, configFsmdValue[i]);
      }
    }
    printInfo << "using final-state mass dependence as defined in the configuration file." << endl;
  } else {
    printInfo << "not using final-state mass dependence." << endl;
  }
  compset.setFuncFsmd(fPS);


  // Resonances
  if(configFile.exists("components.resonances")){
    const Setting &bws = configRoot["components"]["resonances"];
    // loop through breitwigners
    int nbw=bws.getLength();
    printInfo << "found " << nbw << " Resonances in config" << endl;
    for(int ibw = 0; ibw < nbw; ++ibw) {
      const Setting &bw = bws[ibw];
      string jpc;
      string name;
      double mass=-1;double ml,mu;int mfix;
      double width=-1;double wl,wu;int wfix;int wdyn;

      check&=bw.lookupValue("name",     name);
      check&=bw.lookupValue("jpc",       jpc);
      const Setting &massSet = bw["mass"];
      check&=massSet.lookupValue("val",        mass);
      check&=massSet.lookupValue("lower",        ml);
      check&=massSet.lookupValue("upper",        mu);
      check&=massSet.lookupValue("fix",        mfix);
      const Setting &widthSet = bw["width"];
      check&=widthSet.lookupValue("val",       width);
      check&=widthSet.lookupValue("lower",       wl);
      check&=widthSet.lookupValue("upper",       wu);
      check&=widthSet.lookupValue("fix",       wfix);
      bool checkdyn=widthSet.lookupValue("dyn",       wdyn);
      if(!checkdyn)wdyn=0;
      cout << "---------------------------------------------------------------------" << endl;
      cout << name << "    JPC = " << jpc << endl;
      cout << "mass(limits)  = " << mass <<" ("<<ml<<","<<mu<<") MeV/c^2";
      if(mfix==1)cout<<"  -- FIXED";
      cout<< endl;
      cout << "width(limits) = " << width <<" ("<<wl<<","<<wu<<") MeV/c^2";
      if(wfix==1)cout<<"  -- FIXED";
      if(wdyn!=0)cout<<"  -- DYNAMIC WIDTH";
      else cout<<"  -- CONST WIDTH";
      cout<< endl;
      const Setting &channelSet = bw["decaychannels"];
      unsigned int nCh=channelSet.getLength();
      cout << "Decaychannels (coupling):" << endl;
      std::map<std::string,pwachannel > channels;
      for(unsigned int iCh=0;iCh<nCh;++iCh){
	const Setting &ch = channelSet[iCh];
	string amp;
	double cRe=0;
	double cIm=0;
	check&=ch.lookupValue("amp",amp);
	check&=ch.lookupValue("coupling_Re",cRe);
	check&=ch.lookupValue("coupling_Im",cIm);
	complex<double> C(cRe,cIm);
	cout << "   " << amp << "  " << C << endl;
        map<string, size_t>::const_iterator it=mdepFit.getWaveIndices().find(amp);
        const multi_array<double, 2>& inPhaseSpaceIntegrals = mdepFit.getInPhaseSpaceIntegrals();
        multi_array<double, 2>::const_array_view<1>::type view = inPhaseSpaceIntegrals[indices[multi_array<double, 2>::index_range()][it->second]];
	channels[amp]=pwachannel(C,mdepFit.getMassBinCenters(),std::vector<double>(view.begin(), view.end()));
      }// end loop over channels
      if(!check){
	printErr << "Bad config value lookup! Check your config file!" << endl;
	return 1;
      }
      pwacomponent* comp1=new pwacomponent(name,mass,width,channels);
      cerr << "created component" << endl;
      comp1->setLimits(ml,mu,wl,wu);
      comp1->setFixed(mfix,wfix);
      if(wdyn==0)comp1->setConstWidth();
      compset.add(comp1);
      cout << "CHECK val(m0)="<< comp1->val(mass) << endl;
    }// end loop over resonances
  }
  cout << endl;
  // Background components
  if(configFile.exists("components.background")){
    const Setting &bws = configRoot["components"]["background"];
    // loop through breitwigners
    int nbw=bws.getLength();
    printInfo << "found " << nbw << " Background components in config" << endl;
    for(int ibw = 0; ibw < nbw; ++ibw) {
      const Setting &bw = bws[ibw];
      string name;
      double mass=-1;double ml,mu;int mfix;
      double width=-1;double wl,wu;int wfix;

      check&=bw.lookupValue("name",     name);
      const Setting &massSet = bw["m0"];
      check&=massSet.lookupValue("val",        mass);
      check&=massSet.lookupValue("lower",        ml);
      check&=massSet.lookupValue("upper",        mu);
      check&=massSet.lookupValue("fix",        mfix);
      const Setting &widthSet = bw["g"];
      check&=widthSet.lookupValue("val",       width);
      check&=widthSet.lookupValue("lower",       wl);
      check&=widthSet.lookupValue("upper",       wu);
      check&=widthSet.lookupValue("fix",       wfix);
      cout << "---------------------------------------------------------------------" << endl;
      cout << name << endl;
      cout << "mass-offset(limits)  = " << mass <<" ("<<ml<<","<<mu<<") MeV/c^2";
      if(mfix==1)cout<<"  -- FIXED";
      cout<< endl;
      cout << "g(limits)            = " << width <<" ("<<wl<<","<<wu<<") MeV/c^2";
      if(wfix==1)cout<<"  -- FIXED";
      cout<< endl;
      std::map<std::string,pwachannel > channels;
      string amp;
      double cRe=0;
      double cIm=0;
      double mIso1=0;
      double mIso2=0;
      check&=bw.lookupValue("amp",amp);
      check&=bw.lookupValue("coupling_Re",cRe);
      check&=bw.lookupValue("coupling_Im",cIm);
      check&=bw.lookupValue("mIsobar1",mIso1);
      check&=bw.lookupValue("mIsobar2",mIso2);
      complex<double> C(cRe,cIm);
      cout << "Decaychannel (coupling):" << endl;
      cout << "   " << amp << "  " << C << endl;
      cout << "   Isobar masses: " << mIso1<<"  "<< mIso2<< endl;
      map<string, size_t>::const_iterator it=mdepFit.getWaveIndices().find(amp);
      const multi_array<double, 2>& inPhaseSpaceIntegrals = mdepFit.getInPhaseSpaceIntegrals();
      multi_array<double, 2>::const_array_view<1>::type view = inPhaseSpaceIntegrals[indices[multi_array<double, 2>::index_range()][it->second]];
      channels[amp]=pwachannel(C,mdepFit.getMassBinCenters(),std::vector<double>(view.begin(), view.end()));

      if(!check){
	printErr << "Bad config value lookup! Check your config file!" << endl;
	return 1;
      }
      pwabkg* bkg=new pwabkg(name,mass,width,channels);
      bkg->setIsobars(mIso1,mIso2);
      bkg->setLimits(ml,mu,wl,wu);
      bkg->setFixed(mfix,wfix);
      compset.add(bkg);
    }// end loop over background
  }// endif


 cout << "---------------------------------------------------------------------" << endl << endl;

	if(not compset.doMapping()) {
		printErr << "error while mapping the waves to the decay channels and components." << endl;
		exit(1);
	}

 cout << "---------------------------------------------------------------------" << endl << endl;

  // set anchorwave
  vector<string> anchorwave_channel;
  vector<string> anchorwave_reso;
  if(configFile.exists("components.anchorwave")){
    const Setting &anc = configRoot["components"]["anchorwave"];
    // loop through breitwigners
    unsigned int nanc=anc.getLength();
    for(unsigned int ianc=0;ianc<nanc;++ianc){
      string ch,re;
      const Setting &anco = anc[ianc];
      anco.lookupValue("channel",ch);
      anco.lookupValue("resonance",re);
      cout << "Ancorwave: "<< endl;
      cout << "    " << re << endl;
      cout << "    " << ch << endl;
      anchorwave_channel.push_back(ch);
      anchorwave_reso.push_back(re);
    }
  }


    cout << "---------------------------------------------------------------------" << endl << endl;



	massDepFitLikeli L;
	L.init(&compset,
	       mdepFit.getMassBinCenters(),
	       mdepFit.getInSpinDensityMatrices(),
	       mdepFit.getInSpinDensityCovarianceMatrices(),
	       mdepFit.getWavePairMassBinLimits(),
	       doCov);

   const unsigned int nmbPar  = L.NDim();
  // double par[nmbPar];
  // for(unsigned int ip=0;ip<nmbPar;++ip)par[ip]=1.4;


  // TStopwatch watch;
  // L.DoEval(par);
  // watch.Stop();


  //printInfo << "TESTCALL TO LIKELIHOOD takes " <<  maxPrecisionAlign(watch.CpuTime()) << " s" << endl;

  printInfo << nmbPar << " Parameters in fit" << endl;
 
	// ---------------------------------------------------------------------------
	// setup minimizer
	printInfo << "creating and setting up minimizer '" << minimizerType[0] << "' "
	          << "using algorithm '" << minimizerType[1] << "'" << endl;
	Minimizer* minimizer = Factory::CreateMinimizer(minimizerType[0], minimizerType[1]);
	if(not minimizer) { 
		printErr << "could not create minimizer. exiting." << endl;
		throw;
	}
	minimizer->SetFunction        (L);
	minimizer->SetStrategy        (minimizerStrategy);
	minimizer->SetTolerance       (minimizerTolerance);
	minimizer->SetPrintLevel      ((quiet) ? 0 : 3);
	minimizer->SetMaxIterations   (maxNmbOfIterations);
	minimizer->SetMaxFunctionCalls(maxNmbOfFunctionCalls);

  // ---------------------------------------------------------------------------

  // Set startvalues
  unsigned int parcount=0;
  for(unsigned int ic=0;ic<compset.n();++ic){
    const pwacomponent& comp=*compset[ic];
    TString name(comp.name());
    double mmin,mmax,gmin,gmax;
    comp.getLimits(mmin,mmax,gmin,gmax);
    if(comp.fixM())minimizer->SetFixedVariable(parcount++,
					       (name+"_M").Data() ,
					       comp.m0());
    else minimizer->SetLimitedVariable(parcount++,
				       (name+"_M").Data(),
				       comp.m0(),
				       0.10,
				       mmin,mmax);
    if(comp.fixGamma())minimizer->SetFixedVariable(parcount++,
						   (name+"_Gamma").Data() ,
						   comp.gamma());
    else minimizer->SetLimitedVariable(parcount++,
				       (name+"_Gamma").Data(),
				       comp.gamma(),
				       0.01,
				       gmin,gmax);
    std::map<std::string,pwachannel >::const_iterator it=comp.channels().begin();
    while(it!=comp.channels().end()){
      minimizer->SetVariable(parcount++,(name + "_ReC" + it->first).Data() , it->second.C().real(), 0.10);

      // fix one phase
      if(find(anchorwave_reso.begin(),anchorwave_reso.end(),name)!=anchorwave_reso.end() && find(anchorwave_channel.begin(),anchorwave_channel.end(),it->first)!=anchorwave_channel.end()){
	minimizer->SetFixedVariable(parcount++,(name + "_ImC" + it->first).Data() , 0.0);
      }
      else {minimizer->SetVariable(parcount++,(name + "_ImC" + it->first).Data() , it->second.C().imag(), 0.10);}

      ++it;
    } // end loop over channels
  }// end loop over components
  // set phase space
  unsigned int nfreeFsmd=compset.nFreeFsmdPar();
  for(unsigned int ifreeFsmd=0;ifreeFsmd<nfreeFsmd;++ifreeFsmd){
    double val,lower,upper;
    val=compset.getFreeFsmdPar(ifreeFsmd);
    compset.getFreeFsmdLimits(ifreeFsmd,lower,upper);
    TString name("PSP_"); name+=+ifreeFsmd;
    minimizer->SetLimitedVariable(parcount++, 
				  name.Data(), 
				  val, 0.0001 ,lower,upper);
  }



  const unsigned int nfree=minimizer->NFree();
  printInfo <<  nfree  << " Free Parameters in fit" << endl;


  // find minimum of likelihood function
  double chi2=0;
  if(onlyPlotting) printInfo << "Plotting mode, skipping minimzation!" << endl;
  else {
    printInfo << "performing minimization. MASSES AND WIDTHS FIXED" << endl;
    
    // only do couplings
    TStopwatch fitW;
    // releasePars(minimizer,compset,anchorwave_reso,anchorwave_channel,0);
    bool success = minimizer->Minimize();
    if(!success)printWarn << "minimization failed." << endl;
    else printInfo << "minimization successful." << endl;
    printInfo << "Minimization took " <<  maxPrecisionAlign(fitW.CpuTime()) << " s" << endl;
    //release masses
    releasePars(minimizer,compset,anchorwave_reso,anchorwave_channel,1);
    printInfo << "performing minimization. MASSES RELEASED" << endl;
    fitW.Start();
    success &= minimizer->Minimize();
    if(!success)printWarn << "minimization failed." << endl;
    else printInfo << "minimization successful." << endl;
    printInfo << "Minimization took " <<  maxPrecisionAlign(fitW.CpuTime()) << " s" << endl;
    //release widths
    releasePars(minimizer,compset,anchorwave_reso,anchorwave_channel,2);
    printInfo << "performing minimization. ALL RELEASED" << endl;
    fitW.Start();
    success &= minimizer->Minimize();
    printInfo << "Minimization took " <<  maxPrecisionAlign(fitW.CpuTime()) << " s" << endl;

    const double* par=minimizer->X();
    compset.setPar(par);
    cerr << compset << endl;
    if (success){
      printInfo << "minimization finished successfully." << endl;
      chi2=minimizer->MinValue();
    }
    else
      printWarn << "minimization failed." << endl;
    if (runHesse) {
      printInfo << "calculating Hessian matrix." << endl;
      success = minimizer->Hesse();  // comes only with ROOT 5.24+
      if (!success)
	printWarn << "calculation of Hessian matrix failed." << endl;
    }
  }
  printInfo << "minimization stopped after " << minimizer->NCalls() << " function calls. minimizer status summary:" << endl
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


  // ---------------------------------------------------------------------------
  // print results
  //map<TString, double> errormap;
  printInfo << "minimization result:" << endl;
  for (unsigned int i = 0; i< nmbPar; ++i) {
    cout << "    parameter [" << setw(3) << i << "] ";
    cout << minimizer->VariableName(i) << " " ;
      //	 << setw(maxParNameLength); //<< L.parName(i) << " = ";
    //if (parIsFixed[i])
    //  cout << minimizer->X()[i] << " (fixed)" << endl;
    //else {
      cout << setw(12) << maxPrecisionAlign(minimizer->X()[i]) << " +- "
	   << setw(12) << maxPrecisionAlign(minimizer->Errors()[i]);
      //errormap[minimizer]=minimizer->Errors()[i];


      if (runMinos && (i == 156)) {  // does not work for all parameters
	double minosErrLow = 0;
	double minosErrUp  = 0;
	const bool success = minimizer->GetMinosError(i, minosErrLow, minosErrUp);
	if (success)
	  cout << "    Minos: " << "[" << minosErrLow << ", +" << minosErrUp << "]" << endl;
      } else
	cout << endl;
  }

 cout << "---------------------------------------------------------------------" << endl;
 // Reduced chi2

 printInfo << chi2 << " chi2" << endl;
 unsigned int numdata=L.NDataPoints();
 // numDOF
 unsigned int numDOF=numdata-nfree;
 printInfo << numDOF << " degrees of freedom" << endl;
 double redChi2 = chi2/(double)numDOF;
 printInfo << redChi2 << " chi2/nDF" << endl;
 cout << "---------------------------------------------------------------------" << endl;


  // write out results
  // Likelihood and such
 const Setting& fitqualS= configRoot["fitquality"];
 Setting& chi2S=fitqualS["chi2"];
 chi2S=chi2;
 Setting& ndfS=fitqualS["ndf"];
 ndfS=(int)numDOF;
 Setting& redchi2S=fitqualS["redchi2"];
 redchi2S=redChi2;

  if(configFsmd){
    const Setting& configFsmdValue = (*configFsmd)["val"];
    const Setting& configFsmdError = (*configFsmd)["error"];
    const Setting& configFsmdFix = (*configFsmd)["fix"];

    const unsigned int nrPar = fPS->GetNpar();
    unsigned int iPar = 0;
    for(unsigned int i=0; i<nrPar; ++i) {
      if(((int)configFsmdFix[i]) == 0) {
        ostringstream sName;
        sName << "PSP_" << iPar;
        configFsmdValue[i] = compset.getFreeFsmdPar(iPar);
        configFsmdError[i] = minimizer->Errors()[minimizer->VariableIndex(sName.str())];
        ++iPar;
      }
    }
    assert(compset.nFreeFsmdPar() == iPar);
  }
 
  // Setup Component Set (Resonances + Background)
  const Setting& bws= configRoot["components"]["resonances"];
  const Setting& bkgs= configRoot["components"]["background"];
  unsigned int nbws=bws.getLength();
  unsigned int nbkgs=bkgs.getLength();
  // loop over components
  unsigned int nc=compset.n();
  for(unsigned int ic=0;ic<nc;++ic){
    const pwacomponent* comp=compset[ic];
    string name=comp->name();
    // search corresponding setting

    string sname;
    bool found=false;
    for(unsigned int is=0;is<nbws;++is){
      const Setting& bw = bws[is];
      bw.lookupValue("name",     sname);
      if(sname==name){
	found=true;
	// set values to this setting
	Setting& sm = bw["mass"];
	Setting& smval = sm["val"];
	smval = comp->m0();
	Setting& smerr = sm["error"];
	TString merrname=name+"_M";
	smerr=minimizer->Errors()[minimizer->VariableIndex(merrname.Data())];

	Setting& sw = bw["width"];
	Setting& swval = sw["val"];
	swval = comp->gamma();

	Setting& swerr = sw["error"];
	TString werrname=name+"_Gamma";
	swerr=minimizer->Errors()[minimizer->VariableIndex(werrname.Data())];

	cout << name
	     << "   mass="<<double(smval)<<" +- "<<double(smerr)
	     << "   width="<<double(swval)<<" +- "<<double(swerr)<< endl;

	// loop through channel and fix couplings
	const Setting& sChn=bw["decaychannels"];
	unsigned int nCh=sChn.getLength();
	const std::map<std::string,pwachannel >& ch=comp->channels();
	std::map<std::string,pwachannel>::const_iterator it=ch.begin();
	for(;it!=ch.end();++it){
	  std::complex<double> c= it->second.C();
	  string ampname=it->first;
	  // loop through channels in setting
	  for(unsigned int isc=0;isc<nCh;++isc){
	    Setting& sCh=sChn[isc];
	    string amp; sCh.lookupValue("amp",amp);
	    if(amp==ampname){
	      Setting& sRe=sCh["coupling_Re"];
	      sRe=c.real();
	      Setting& sIm=sCh["coupling_Im"];
	       sIm=c.imag();
	      break;
	    } // endif
	  } // end loop through cannels in setting

	} // end loop through channels of component

	break;
      }
    }

    // loop over background settings
    if(!found){
      for(unsigned int is=0;is<nbkgs;++is){
	const Setting& bw = bkgs[is];
	bw.lookupValue("name",     sname);
	if(sname==name){
	  Setting& sm = bw["m0"];
	  Setting& smval = sm["val"];
	  smval = comp->m0();
	  Setting& sw = bw["g"];
	  Setting& swval = sw["val"];
	  swval = comp->gamma();

	  const pwachannel& ch=comp->channels().begin()->second;
	  std::complex<double> c=ch.C();
	  Setting& sRe=bw["coupling_Re"];
	  sRe=c.real();
	  Setting& sIm=bw["coupling_Im"];
	  sIm=c.imag();
	  break;
	}
      }
    }




  }


  
  string outconfig(outFileName);
  outconfig.append(".conf");
  configFile.writeFile(outconfig.c_str());

  cerr << "Fitting finished... Start building graphs ... " << endl;

  int syscolor=kAzure-9;

   std::vector<std::string> wl=compset.getWaveList();
   std::map<std::string, unsigned int> wmap;
   const size_t nrMassBins=mdepFit.getMassBinCenters().size();

   std::vector<TGraphErrors*> datagraphs;
   std::vector<TGraphErrors*> intenssysgraphs;

   std::vector<TMultiGraph*> graphs;

   for(unsigned int iw=0; iw<wl.size();++iw){
     wmap[wl[iw]]=iw;
     graphs.push_back(new TMultiGraph);

     intenssysgraphs.push_back(new TGraphErrors(nrMassBins));
     string name=("sys_");name.append(wl[iw]);
     intenssysgraphs[iw]->SetName(name.c_str());
     intenssysgraphs[iw]->SetTitle(name.c_str());
     intenssysgraphs[iw]->SetLineColor(syscolor);
     intenssysgraphs[iw]->SetFillColor(syscolor);
     intenssysgraphs[iw]->SetDrawOption("2");
     graphs[iw]->Add(intenssysgraphs[iw],"2");



     graphs[iw]->SetName(wl[iw].c_str());
     graphs[iw]->SetTitle(wl[iw].c_str());
     graphs[iw]->SetDrawOption("AP");
     datagraphs.push_back(new TGraphErrors(nrMassBins));
     name="data_";name.append(wl[iw]);
     datagraphs[iw]->SetName(name.c_str());
     datagraphs[iw]->SetTitle(name.c_str());
     datagraphs[iw]->SetDrawOption("AP");
     //datagraphs[iw]->SetLineColor(kRed);
     //datagraphs[iw]->SetMarkerColor(kRed);

     graphs[iw]->Add(datagraphs[iw],"P");
   }




     // build fitgraphs
   unsigned int nbins=nrMassBins;//200;
   //double mmin=1200.;
   //double md=10.;
   std::vector<TGraph*> fitgraphs;
   std::vector<TGraph*> absphasegraphs;


   for(unsigned int iw=0; iw<wl.size();++iw){
     fitgraphs.push_back(new TGraph(nbins));
     string name("fit_");name.append(wl[iw]);
     fitgraphs[iw]->SetName(name.c_str());
     fitgraphs[iw]->SetTitle(name.c_str());
     fitgraphs[iw]->SetLineColor(kRed);
     fitgraphs[iw]->SetLineWidth(2);
     fitgraphs[iw]->SetMarkerColor(kRed);
     fitgraphs[iw]->SetDrawOption("AP");
     //fitgraphs[iw]->SetMarkerStyle(22);
     graphs[iw]->Add(fitgraphs[iw],"cp");
     graphs[iw]->Add(new TGraph(mdepFit.getMassBinCenters().size(), &(mdepFit.getMassBinCenters()[0]), &(mdepFit.getInPhaseSpaceIntegrals()[iw][0])));

     absphasegraphs.push_back(new TGraph(nbins));
     name="absphase_";name.append(wl[iw]);
     absphasegraphs[iw]->SetName(name.c_str());
     absphasegraphs[iw]->SetTitle(name.c_str());
     absphasegraphs[iw]->SetLineColor(kRed);
     absphasegraphs[iw]->SetLineWidth(2);
     absphasegraphs[iw]->SetMarkerColor(kRed);
     absphasegraphs[iw]->SetDrawOption("AP");



   }

   std::vector<TGraph*> compgraphs; // individual components
   // loop over components and build graphs
     for(unsigned int ic=0;ic<compset.n();++ic){
       const pwacomponent* c=compset[ic];
       std::map<std::string,pwachannel >::const_iterator it=c->channels().begin();
       while(it!=c->channels().end()){
	 string name=c->name();name.append("__");
	 name.append(it->first);
	 TGraph* gcomp=new TGraph(nbins);
	 gcomp->SetName(name.c_str());
	 gcomp->SetTitle(name.c_str());
	 unsigned int color=kBlue;
	 if(dynamic_cast<const pwabkg*>(c)!=NULL)color=kMagenta;
	 gcomp->SetLineColor(color);
	 gcomp->SetMarkerColor(color);

	 compgraphs.push_back(gcomp);
	 graphs[wmap[it->first]]->Add(gcomp,"cp");
	 ++it;
       }// end loop over channels

     }// end loop over components


   std::vector<TGraphErrors*> phasedatagraphs;
   std::vector<TGraphErrors*> phasesysgraphs;


   std::vector<TGraphErrors*> realdatagraphs;
   std::vector<TGraphErrors*> realsysgraphs;
   std::vector<TGraphErrors*> imagdatagraphs;
   std::vector<TGraphErrors*> imagsysgraphs;

   std::vector<TGraph*> realfitgraphs;
   std::vector<TGraph*> imagfitgraphs;

    std::vector<TMultiGraph*> phasegraphs;
   std::vector<TMultiGraph*> overlapRegraphs;
   std::vector<TMultiGraph*> overlapImgraphs;

   std::vector<TGraph*> phasefitgraphs;
   unsigned int count=0;
  

  for(unsigned int iw=0; iw<wl.size();++iw){
     for(unsigned int iw2=iw+1; iw2<wl.size();++iw2){



       phasegraphs.push_back(new TMultiGraph);
       overlapImgraphs.push_back(new TMultiGraph);
       overlapRegraphs.push_back(new TMultiGraph);
       string name("dPhi_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       phasegraphs[count]->SetName(name.c_str());
       phasegraphs[count]->SetTitle(name.c_str());
       name="Re_";name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       overlapRegraphs[count]->SetName(name.c_str());
       overlapRegraphs[count]->SetTitle(name.c_str());
       name="Im_";name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       overlapImgraphs[count]->SetName(name.c_str());
       overlapImgraphs[count]->SetTitle(name.c_str());

       phasesysgraphs.push_back(new TGraphErrors(nbins));
       name=("dPhi_sys_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       phasesysgraphs[count]->SetName(name.c_str());
       phasesysgraphs[count]->SetTitle(name.c_str());
       phasesysgraphs[count]->SetLineColor(syscolor);
       phasesysgraphs[count]->SetFillColor(syscolor);
       phasesysgraphs[count]->SetDrawOption("2");
       //phasesysgraphs[count]->SetLineWidth(1);
       //phasesysgraphs[count]->SetFillStyle(1001);
       
       phasegraphs[count]->Add(phasesysgraphs[count],"2");


       phasedatagraphs.push_back(new TGraphErrors(nbins));
       name=("dPhi_data_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       phasedatagraphs[count]->SetName(name.c_str());
       phasedatagraphs[count]->SetTitle(name.c_str());
       phasegraphs[count]->Add(phasedatagraphs[count],"p");


       realsysgraphs.push_back(new TGraphErrors(nbins));
       name=("RE_sys_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       realsysgraphs[count]->SetName(name.c_str());
       realsysgraphs[count]->SetTitle(name.c_str());
       realsysgraphs[count]->SetLineColor(syscolor);
       realsysgraphs[count]->SetFillColor(syscolor);
       realsysgraphs[count]->SetDrawOption("2");
       overlapRegraphs[count]->Add(realsysgraphs[count],"2");

       imagsysgraphs.push_back(new TGraphErrors(nbins));
       name=("IM_sys_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       imagsysgraphs[count]->SetName(name.c_str());
       imagsysgraphs[count]->SetTitle(name.c_str());
       imagsysgraphs[count]->SetLineColor(syscolor);
       imagsysgraphs[count]->SetFillColor(syscolor);
       imagsysgraphs[count]->SetDrawOption("2");
       overlapImgraphs[count]->Add(imagsysgraphs[count],"2");


       realdatagraphs.push_back(new TGraphErrors(nbins));
       name=("RE_data_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       realdatagraphs[count]->SetName(name.c_str());
       realdatagraphs[count]->SetTitle(name.c_str());
       overlapRegraphs[count]->Add(realdatagraphs[count],"p");

       imagdatagraphs.push_back(new TGraphErrors(nbins));
       name=("IM_data_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       imagdatagraphs[count]->SetName(name.c_str());
       imagdatagraphs[count]->SetTitle(name.c_str());
       //imagdatagraphs[count]->SetLineStyle(2);
       overlapImgraphs[count]->Add(imagdatagraphs[count],"p");

       ++count;
     }
   }

   count=0;
   for(unsigned int iw=0; iw<wl.size();++iw){
     for(unsigned int iw2=iw+1; iw2<wl.size();++iw2){
       phasefitgraphs.push_back(new TGraph(nbins));
       string name("dPhi_fit_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       phasefitgraphs[count]->SetName(name.c_str());
       phasefitgraphs[count]->SetTitle(name.c_str());
       phasefitgraphs[count]->SetLineColor(kRed);
       phasefitgraphs[count]->SetMarkerColor(kRed);
       phasefitgraphs[count]->SetDrawOption("P");
       phasefitgraphs[count]->SetMarkerStyle(24);
       phasefitgraphs[count]->SetMarkerSize(0.2);
       phasegraphs[count]->Add(phasefitgraphs[count],"cp");

       realfitgraphs.push_back(new TGraph(nbins));
       name=("Re_fit_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       realfitgraphs[count]->SetName(name.c_str());
       realfitgraphs[count]->SetTitle(name.c_str());
       realfitgraphs[count]->SetLineColor(kRed);
       realfitgraphs[count]->SetMarkerColor(kRed);
       realfitgraphs[count]->SetDrawOption("AP");
       realfitgraphs[count]->SetMarkerStyle(24);
       realfitgraphs[count]->SetMarkerSize(0.2);
       overlapRegraphs[count]->Add(realfitgraphs[count],"cp");

       imagfitgraphs.push_back(new TGraph(nbins));
       name=("Im_fit_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       imagfitgraphs[count]->SetName(name.c_str());
       imagfitgraphs[count]->SetTitle(name.c_str());
       imagfitgraphs[count]->SetLineColor(kRed);
       //imagfitgraphs[count]->SetLineStyle(2);
       imagfitgraphs[count]->SetMarkerColor(kRed);
       imagfitgraphs[count]->SetDrawOption("AP");
       imagfitgraphs[count]->SetMarkerStyle(24);
       imagfitgraphs[count]->SetMarkerSize(0.2);
       overlapImgraphs[count]->Add(imagfitgraphs[count],"cp");

       ++count;
     }
   }

   std::vector<TGraph2D*> phase2d;
   count=0;
   for(unsigned int iw=0; iw<wl.size();++iw){
     for(unsigned int iw2=iw+1; iw2<wl.size();++iw2){
       phase2d.push_back(new TGraph2D(nbins));
       string name("dPhi_2d_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       phase2d[count]->SetName(name.c_str());
       phase2d[count]->SetTitle(name.c_str());
       phase2d[count]->SetLineColor(kRed);
       phase2d[count]->SetMarkerColor(kRed);
       //phasegraphs[count]->Add(phasefitgraphs[count],"cp");
       ++count;
     }
   }






  // get data
   vector<double> prevphase(wl.size());
   double binwidth=MASSSCALE*30; // half binwidth
   //double w=2*30/10;
   for(unsigned int i=0;i<nrMassBins;++i){
     double m=mdepFit.getMassBinCenters()[i];
     double mScaled=m*MASSSCALE;
     //cout << "MASS: "<<m << endl;
     unsigned int count=0;
     for(unsigned int iw=0; iw<wl.size();++iw){
		// check that this mass bin should be taken into account for this
		// combination of waves
		if(rangePlotting && (i < mdepFit.getWavePairMassBinLimits()[iw][iw].first || i > mdepFit.getWavePairMassBinLimits()[iw][iw].second)) {
			continue;
		}

       datagraphs[iw]->SetPoint(i,mScaled,mdepFit.getInIntensities()[i][iw][0]);
       datagraphs[iw]->SetPointError(i,binwidth,mdepFit.getInIntensities()[i][iw][1]);
      fitgraphs[iw]->SetPoint(i,mScaled,compset.intensity(wl[iw],m));      
      double absphase=compset.phase(wl[iw],m)*TMath::RadToDeg();
      if(i>0){
	double absp=absphase+360;
	double absm=absphase-360;
	if(fabs(absphase-prevphase[iw])>fabs(absp-prevphase[iw])){
	  absphase=absp;
	}
	else if(fabs(absphase-prevphase[iw])>fabs(absm-prevphase[iw])){
	  absphase=absm;
	}
      }
      prevphase[iw]=absphase;
      absphasegraphs[iw]->SetPoint(i,mScaled,absphase);      
      if(mdepFit.getSysPlotting()){
	double maxIntens=-10000000;
	double minIntens=10000000;
	for(unsigned int iSys=0;iSys<mdepFit.getNrSystematics();++iSys){
	  // get data
// rename one wave
	  string mywave1=wl[iw];


	  if(mywave1=="1-1++0+pi-_01_rho1700=pi-+_10_pi1300=pi+-_00_sigma.amp" && iSys>0)mywave1="1-1++0+pi-_01_eta11600=pi-+_10_pi1300=pi+-_00_sigma.amp";

	  // check if waves are in fit
	  // FIXME: skip waves that are not in the fit
	  double myI=mdepFit.getSysIntensities()[iSys][i][iw][0];
	  if(maxIntens<myI)maxIntens=myI;
	  if(minIntens>myI)minIntens=myI;
	} // end loop over systematic trees

	intenssysgraphs[iw]->SetPoint(i,mScaled,(maxIntens+minIntens)*0.5);
	intenssysgraphs[iw]->SetPointError(i,binwidth,(maxIntens-minIntens)*0.5);
      }

      //cout << "getting phases" << endl;

       // second loop to get phase differences
       
       for(unsigned int iw2=iw+1; iw2<wl.size();++iw2){
			// check that this mass bin should be taken into account for this
			// combination of waves
			if(rangePlotting && (i < mdepFit.getWavePairMassBinLimits()[iw][iw2].first || i > mdepFit.getWavePairMassBinLimits()[iw][iw2].second)) {
				continue;
			}

	 realdatagraphs[count]->SetPoint(i,
				     mScaled,
				     mdepFit.getInSpinDensityMatrices()[i][iw][iw2].real());
	 realdatagraphs[count]->SetPointError(i,
			    binwidth,
			    sqrt(mdepFit.getInSpinDensityCovarianceMatrices()[i][iw][iw2][0][0]));
	 imagdatagraphs[count]->SetPoint(i,
		       mScaled,
		       mdepFit.getInSpinDensityMatrices()[i][iw][iw2].imag());
     
	 imagdatagraphs[count]->SetPointError(i,
			    binwidth,
			    sqrt(mdepFit.getInSpinDensityCovarianceMatrices()[i][iw][iw2][1][1]));


	 double dataphi=mdepFit.getInPhases()[i][iw][iw2][0];

	 phasedatagraphs[count]->SetPoint(i,mScaled,dataphi);
	   
	 TVector2 v;v.SetMagPhi(1,dataphi/TMath::RadToDeg());
	 phase2d[count]->SetPoint(i,v.X(),v.Y(),mScaled);

	 phasedatagraphs[count]->SetPointError(i,binwidth,mdepFit.getInPhases()[i][iw][iw2][1]);
	 double fitphase=compset.phase(wl[iw],wl[iw2],m)*TMath::RadToDeg();

	 if(mdepFit.getSysPlotting()){
	   //cerr << "start sysplotting" << endl;
	   // loop over systematics files
	   double maxPhase=-10000;
	   double minPhase=10000;
	   double maxRe=-10000000;
	   double minRe=10000000;
	   double maxIm=-10000000;
	   double minIm=10000000;
	  
	   for(unsigned int iSys=0;iSys<mdepFit.getNrSystematics();++iSys){
	     //cerr << iSys;
	   // get data
	     // rename one wave
	     string mywave1=wl[iw];
	     string mywave2=wl[iw2];

if(mywave1=="1-1++0+pi-_01_rho1700=pi-+_10_pi1300=pi+-_00_sigma.amp" && iSys>0)mywave1="1-1++0+pi-_01_eta11600=pi-+_10_pi1300=pi+-_00_sigma.amp";
if(mywave2=="1-1++0+pi-_01_rho1700=pi-+_10_pi1300=pi+-_00_sigma.amp" && iSys>0)mywave2="1-1++0+pi-_01_eta11600=pi-+_10_pi1300=pi+-_00_sigma.amp";

	     // check if waves are in fit
	     // FIXME: skip pairs where at least one wave is not in the fit

	     double myphi=mdepFit.getSysPhases()[iSys][i][iw][iw2][0];
	     double myphiplus=myphi+360;
	     double myphiminus=myphi-360;
	     // translate by 2pi to get closest solution to fit
	     if(fabs(myphiplus-dataphi)<fabs(myphi-dataphi)){
		 myphi=myphiplus;
		 //cout << "myphiminus" << endl;
	     }
	     if(fabs(myphiminus-dataphi)<fabs(myphi-dataphi)){
	       myphi=myphiminus;
	       //cout << "myphiplus" << endl;
	     }

	     if(myphi>maxPhase)maxPhase=myphi;
	     if(myphi<minPhase)minPhase=myphi;

	     // real and imag part:
	     complex<double> r=mdepFit.getSysSpinDensityMatrices()[iSys][i][iw][iw2];
	     if(maxRe<r.real())maxRe=r.real();
	     if(minRe>r.real())minRe=r.real();
	     if(maxIm<r.imag())maxIm=r.imag();
	     if(minIm>r.imag())minIm=r.imag();
	   }// end loop over sys trees
	   // cerr << "loop over systrees finished" << endl;

	   phasesysgraphs[count]->SetPoint(i,mScaled,(maxPhase+minPhase)*0.5);
	   phasesysgraphs[count]->SetPointError(i,binwidth,(maxPhase-minPhase)*0.5);
	   
	   realsysgraphs[count]->SetPoint(i,mScaled,(maxRe+minRe)*0.5);
	   realsysgraphs[count]->SetPointError(i,binwidth,(maxRe-minRe)*0.5);
	   imagsysgraphs[count]->SetPoint(i,mScaled,(maxIm+minIm)*0.5);
	   imagsysgraphs[count]->SetPointError(i,binwidth,(maxIm-minIm)*0.5);
	

	 }// end if sysplotting
	 //cerr << "sysplotting finished" << endl;

	 phasefitgraphs[count]->SetPoint(i,mScaled,fitphase);

	 complex<double> fitval=compset.overlap(wl[iw],wl[iw2],m);

	 realfitgraphs[count]->SetPoint(i,mScaled,fitval.real());
	 imagfitgraphs[count]->SetPoint(i,mScaled,fitval.imag());


	 count++;
       }// end inner loop over waves

     } // end loop over waves
     //cerr << "outer loop over waves finished" << endl;
     // loop over components to fill individual graphs
     unsigned int compcount=0;
       for(unsigned int ic=0;ic<compset.n();++ic){
	 const pwacomponent* c=compset[ic];
			for(std::map<std::string, pwachannel>::const_iterator itChan=c->channels().begin(); itChan!=c->channels().end(); ++itChan, ++compcount) {
				// check that this mass bin should be taken into account for this
				// combination of waves
				const size_t iw = wmap[itChan->first];
				if(rangePlotting && (i < mdepFit.getWavePairMassBinLimits()[iw][iw].first || i > mdepFit.getWavePairMassBinLimits()[iw][iw].second)) {
					continue;
				}

				const double I=norm(c->val(m)*itChan->second.C()*sqrt(itChan->second.ps(m))*compset.calcFsmd(m));
				compgraphs[compcount]->SetPoint(i,mScaled,I);
			} // end loop over channels
       }// end loop over components

       cerr << "Finished plotting mass-bin " << m << endl;
       //mprev=m;
   }
   cerr << "Finished Loop Over DataBins" << endl;



//    for(unsigned int im=0;im<nbins;++im){ // fine loop in masses -> fits

//      double m=mmin+im*md;
//      for(unsigned int iw=0; iw<wl.size();++iw){
//        fitgraphs[iw]->SetPoint(im,m,compset.intensity(wl[iw],m));
//      } // end loop over waves


//    }// end loop over mass bins





   TFile* outfile=TFile::Open(outFileName.c_str(),"RECREATE");
   for(unsigned int iw=0; iw<wl.size();++iw){
     TGraph* g=(TGraph*)graphs[iw]->GetListOfGraphs()->At(3);
      for(unsigned int ib=0;ib<nbins;++ib){
        double m,ps;g->GetPoint(ib,m,ps);
        g->SetPoint(ib,m*MASSSCALE,ps*500.);
       }

     graphs[iw]->Write();
     //absphasegraphs[iw]->Write();
   }


 for(unsigned int iw=0; iw<phasegraphs.size();++iw){

/// rectivfy phase graphs

   unsigned int refbin=6;
   double m;
   double predph;
   phasedatagraphs[iw]->GetPoint(refbin,m,predph);
   double prefph;
   phasefitgraphs[iw]->GetPoint(refbin,m,prefph);
   for(unsigned int ib=refbin+1;ib<nbins;++ib){
     double dph; phasedatagraphs[iw]->GetPoint(ib,m,dph);
     double fph; phasefitgraphs[iw]->GetPoint(ib,m,fph);
     double dp,dm;dp=dph+360;dm=dph-360;
     double fp,fm;fp=fph+360;fm=fph-360;
     if(1){
       if(fabs(fp-prefph)<fabs(fph-prefph) && fabs(fp-prefph)<fabs(fm-prefph))
         phasefitgraphs[iw]->SetPoint(ib,m,fp);
       else if(fabs(fm-prefph)<fabs(fph-prefph) && fabs(fm-prefph)<fabs(fp-prefph))
         phasefitgraphs[iw]->SetPoint(ib,m,fm);
       phasefitgraphs[iw]->GetPoint(ib,m,prefph);

       if(fabs(dp-prefph)<fabs(dph-prefph) && fabs(dp-prefph)<fabs(dm-prefph))
         phasedatagraphs[iw]->SetPoint(ib,m,dp);
       else if(fabs(dm-prefph)<fabs(dph-prefph) && fabs(dm-prefph)<fabs(dp-prefph))
         phasedatagraphs[iw]->SetPoint(ib,m,dm);



       phasedatagraphs[iw]->GetPoint(ib,m,predph);



   // put systematic error closest to fit/data
       double sph; phasesysgraphs[iw]->GetPoint(ib,m,sph);
       double sp,sm;sp=sph+360;sm=sph-360;
       if(fabs(sp-prefph)<fabs(sph-prefph) && fabs(sp-prefph)<fabs(sm-prefph))
	 phasesysgraphs[iw]->SetPoint(ib,m,sp);
       else if(fabs(sm-prefph)<fabs(sph-prefph) && fabs(sm-prefph)<fabs(sp-prefph))
	 phasesysgraphs[iw]->SetPoint(ib,m,sm);


     }
   }
   // backward:
   phasedatagraphs[iw]->GetPoint(refbin,m,predph);
   phasefitgraphs[iw]->GetPoint(refbin,m,prefph);
   for(unsigned int i=0;i<refbin;++i){
       unsigned int ib=refbin-i-1;
       double dph; phasedatagraphs[iw]->GetPoint(ib,m,dph);
     double fph; phasefitgraphs[iw]->GetPoint(ib,m,fph);
     double dp,dm;dp=dph+360;dm=dph-360;
     double fp,fm;fp=fph+360;fm=fph-360;
     if(1){

       if(fabs(fp-prefph)<fabs(fph-prefph) && fabs(fp-prefph)<fabs(fm-prefph)){
         phasefitgraphs[iw]->SetPoint(ib,m,fp);

       }
       else if(fabs(fm-prefph)<fabs(fph-prefph) && fabs(fm-prefph)<fabs(fp-prefph)){
         phasefitgraphs[iw]->SetPoint(ib,m,fm);

       }

       phasefitgraphs[iw]->GetPoint(ib,m,prefph);



       if(fabs(dp-prefph)<fabs(dph-prefph) && fabs(dp-prefph)<fabs(dm-prefph)){
         phasedatagraphs[iw]->SetPoint(ib,m,dp);
       }

       else if(fabs(dm-prefph)<fabs(dph-prefph) && fabs(dm-prefph)<fabs(dp-prefph)){
         phasedatagraphs[iw]->SetPoint(ib,m,dm);
       }

       phasedatagraphs[iw]->GetPoint(ib,m,predph);


       // put systematic error closest to fit
       double sph; phasesysgraphs[iw]->GetPoint(ib,m,sph);
       double sp,sm;sp=sph+360;sm=sph-360;
       if(fabs(sp-predph)<fabs(sph-predph) && fabs(sp-predph)<fabs(sm-predph))
	 phasesysgraphs[iw]->SetPoint(ib,m,sp);
       else if(fabs(sm-predph)<fabs(sph-predph) && fabs(sm-predph)<fabs(sp-predph))
	 phasesysgraphs[iw]->SetPoint(ib,m,sm);


     }
   }


     phasegraphs[iw]->Write();
     phasesysgraphs[iw]->Write();
     overlapRegraphs[iw]->Write();
     overlapImgraphs[iw]->Write();
     //phase2d[iw]->Write();
 } // end loop over waves

 if (fPS != NULL) {
   fPS->Write("funcFinalStateMassDependence");
 }

outfile->Close();

   return 0;

}
