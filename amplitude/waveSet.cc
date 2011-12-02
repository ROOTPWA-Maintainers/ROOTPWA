///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
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
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      handles wave sets for fitting
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "libconfig.h++"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TList.h"

#include "libConfigUtils.hpp"
#include "waveDescription.h"
#include "waveSet.h"

  
using namespace std;
using namespace libconfig;
using namespace rpwa;


#ifdef USE_STD_COMPLEX_TREE_LEAFS
ClassImp(waveSet);
#endif

bool waveSet::_debug = false;


waveSet::waveSet()
	: TObject            (),
	  _decayAmpTreeNames (),
	  _decayAmpMassRanges(),
	  _decayAmpFileNames (),
	  _decayAmpTrees     (),
	  _waveDescs         (),
	  _decayAmps         ()
{
	//waveSet::Class()->IgnoreTObjectStreamer();  // don't store TObject's fBits and fUniqueID
}


waveSet::~waveSet()
{ }


void
waveSet::clear()
{
	_decayAmpTreeNames.clear ();
	_decayAmpMassRanges.clear();
	_decayAmpFileNames.clear ();
	_decayAmpTrees.clear     ();
	_waveDescs.clear         ();
	_decayAmps.clear         ();
}


waveSet&
waveSet::operator =(const waveSet& set)
{
	if (this != &set) {
		TObject::operator   =(set);
		_decayAmpTreeNames  = set._decayAmpTreeNames;
		_decayAmpMassRanges = set._decayAmpMassRanges;
		_decayAmpFileNames  = set._decayAmpFileNames;
		_decayAmpTrees      = set._decayAmpTrees;
		_waveDescs          = set._waveDescs;
		_decayAmps          = set._decayAmps;
	}
	return *this;
}


void
waveSet::setDecayAmpFileNames(const vector<string>& ampFileNames)
{
	_decayAmpFileNames.clear();
	_decayAmpFileNames = ampFileNames;
}


bool
waveSet::buildWaveSet(const string& waveSetFileName)
{
	if (waveSetFileName == "") {
		printWarn << "empty wave set file name. cannot build wave set." << endl;
		return false;
	}

	libconfig::Config config;
	if (not parseLibConfigFile(waveSetFileName, config, _debug)) {
		printWarn << "problems reading wave set file '" << waveSetFileName << "'. "
		          << "cannot construct wave set." << endl;
		return false;
	}

	const Setting& configRoot = config.getRoot();
	printInfo << "constructing wave set from file '" << waveSetFileName << "'" << endl;

	// find waveSet group
	const Setting* configWaveSet = findLibConfigGroup(configRoot, "waveSet", true);
	if (not configWaveSet) {
		printWarn << "cannot construct wave set." << endl;
		return false;
	}

	// find list of decay amplitude specifications
	const Setting* configDecayAmps = findLibConfigList(*configWaveSet, "decayAmplitudes", true);
	if (not configDecayAmps) {
		printWarn << "cannot construct wave set." << endl;
		return false;
	}

	// read list of decay amplitude tree names
	const int nmbEntries = configDecayAmps->getLength();
	_decayAmpTreeNames.clear ();
	_decayAmpMassRanges.clear();
	_decayAmpTreeNames.resize (nmbEntries, "");
	_decayAmpMassRanges.resize(nmbEntries, make_pair(0, numeric_limits<double>::infinity()));
	bool success = true;
	for (int i = 0; i < nmbEntries; ++i) {
		// get tree name
		if (not (*configDecayAmps)[i].lookupValue("treeName", _decayAmpTreeNames[i])) {
			printWarn << "entry [" << i << "] in decay amplitude list does not specify a tree name";
			success = false;
		}
		// get mass range
		const Setting* configMassRange = findLibConfigArray((*configDecayAmps)[i], "massRange", false);
		if (configMassRange) {
			try {
				if (configMassRange->getLength() == 1)
					_decayAmpMassRanges[i].first  = (*configMassRange)[0];
				else if (configMassRange->getLength() == 2) {
					_decayAmpMassRanges[i].first  = (*configMassRange)[0];
					_decayAmpMassRanges[i].second = (*configMassRange)[1];
				} else {
					printWarn << "cannot read mass range from entry [" << i << "] in decay amplitude list. "
					          << "array length is neither 1 nor 2." << endl;
					success = false;
				}
			} catch (const SettingTypeException&) {
				try {  // accept also integer values
					if (configMassRange->getLength() == 1)
						_decayAmpMassRanges[i].first  = (unsigned int)(*configMassRange)[0];
					else if (configMassRange->getLength() == 2) {
						_decayAmpMassRanges[i].first  = (unsigned int)(*configMassRange)[0];
						_decayAmpMassRanges[i].second = (unsigned int)(*configMassRange)[1];
					}
				} catch (const SettingTypeException& settingEx) {
					printWarn << "mass range value at '" << settingEx.getPath() << "' is not of type double. "
					          << "using default mass range." << endl;
					success = false;
				}
			}
		}
	}

	if (success)
		printSucc << "constructed wave set from file '" << waveSetFileName << "'" << endl;
	else
		printWarn << "problems constructing wave set from file '" << waveSetFileName << "'" << endl;
	return success;
}


bool
waveSet::getDecayAmplitudeTrees()
{
	_decayAmpTrees.clear();
	_waveDescs.clear    ();
	if (_decayAmpFileNames.size() == 0) {
		printWarn << "array with decay amplitude file names is empty. cannot get amplitude trees." << endl;
		return false;
	}

#ifdef USE_STD_COMPLEX_TREE_LEAFS

	bool success = true;
	// force loading predefined std::complex dictionary
	// see http://root.cern.ch/phpBB3/viewtopic.php?f=5&t=9618&p=50164
	gROOT->ProcessLine("#include <complex>");

	// get trees from .root files
	_decayAmpTrees.resize(nmbWaves(), 0);
	_waveDescs.resize    (nmbWaves(), 0);
	unsigned int countOpenFails = 0;
	unsigned int countDuplTrees = 0;
	for (unsigned int ampFileIndex = 0; ampFileIndex < _decayAmpFileNames.size(); ++ampFileIndex) {
		// open .root file
		if (_debug)
			printDebug << "opening .root decay amplitude file '" << _decayAmpFileNames[ampFileIndex] << "'"
			           << endl;
		TFile* ampFile = TFile::Open(_decayAmpFileNames[ampFileIndex].c_str(), "READ");
		if (not ampFile or ampFile->IsZombie()) {
			printWarn << "cannot open decay amplitude file '" << _decayAmpFileNames[ampFileIndex] << "'. "
			          << "skipping file." << endl;
			++countOpenFails;
			continue;			
		}
		// get all amplitude trees from file
		vector<pair<TTree*, unsigned int> > trees;
		for (unsigned int i = 0; i < nmbWaves(); ++i) {
			if (_debug)
				printDebug << "looking for tree '" << _decayAmpTreeNames[i] << "' "
				           << "in file '" << _decayAmpFileNames[ampFileIndex] << "'" << endl;
			TTree* tree = 0;
			ampFile->GetObject(_decayAmpTreeNames[i].c_str(), tree);
			if (tree) {
				trees.push_back(make_pair(tree, i));
				if (_debug)
					printDebug << "found decay amplitude tree '" << _decayAmpTreeNames[i] << "' "
					           << "in file '" << _decayAmpFileNames[ampFileIndex] << "'" << endl;
			} else if (_debug)
				printDebug << "did not find decay amplitude tree '" << _decayAmpTreeNames[i] << "' "
				           << "in file '" << _decayAmpFileNames[ampFileIndex] << "'" << endl;
		}
		// store tree pointers and make sure trees do not already exist
		for (unsigned int i = 0; i < trees.size(); ++i) {
			if (_decayAmpTrees[trees[i].second]) {
				printWarn << "tree '" << trees[i].first->GetName() << "' exists in file '"
				          << _decayAmpTrees[trees[i].second]->GetDirectory()->GetName() << "' and in file '"
					// GetCurrentFile()
				          << trees[i].first->GetDirectory()->GetName() << "'. skipping tree from latter file."
				          << endl;
				++countDuplTrees;
				continue;
			}
			_decayAmpTrees[trees[i].second] = trees[i].first;
		}		
	}

	// get corresponding wave descriptions from trees
	unsigned int countMissingTrees     = 0;
	unsigned int countMissingWaveDescs = 0;
	for (unsigned int i = 0; i < nmbWaves(); ++i) {
		if (not _decayAmpTrees[i]) {
			printWarn << "did not find decay amplitude tree '" << _decayAmpTreeNames[i] << "' "
			          << "in any of the given decay amplitude files" << endl;
			++countMissingTrees;
			continue;
		}
		TList* userInfo = _decayAmpTrees[i]->GetUserInfo();
		if (not userInfo) {
			printWarn << "tree '" << _decayAmpTreeNames[i] << "' does not have user info." << endl;
			++countMissingWaveDescs;
			continue;
		}
		waveDescription* waveDesc
			= static_cast<waveDescription*>(userInfo->FindObject("rpwa::waveDescription"));
		if (not waveDesc) {
			printInfo << "user info of tree '" << _decayAmpTreeNames[i] << "' "
			          << "does not contain wave description." << endl;
			++countMissingWaveDescs;
			continue;
		}
		_waveDescs[i] = waveDesc;
	}

	// report
	if (countOpenFails > 0) {
		success = false;
		printWarn << "problems opening " << countOpenFails << " out of " << _decayAmpFileNames.size()
		          << " decay amplitude files." << endl;
	}
	if (countDuplTrees > 0) {
		success = false;
		printWarn << "found " << countDuplTrees << " duplicate trees. "
		          << "decay amplitude data were ignored." << endl;
	}
	if (countMissingTrees > 0) {
		success = false;
		printWarn << countMissingTrees << " out of " << nmbWaves() << " trees are missing. "
		          << "decay amplitude data incomplete." << endl;
	}
	if (countMissingWaveDescs > 0) {
		success = false;
		printWarn << countMissingWaveDescs << " out of " << nmbWaves() << " wave descriptions "
		          << "are missing. wave set data incomplete." << endl;
	}
	if (    (countOpenFails    == 0) and (countDuplTrees        == 0)
	    and (countMissingTrees == 0) and (countMissingWaveDescs == 0))
		printSucc << "opened all " << _decayAmpFileNames.size() << " decay amplitude files "
		          << "and found all " << nmbWaves() << " decay amplitude trees "
		          << "and wave descriptions" << endl;

	return success;

#else

	printWarn << "cannot read trees, because your ROOT version does not support "
	          << "std::complex tree leafs. please consider updating your ROOT installation.";
	return false;

#endif
}


bool
waveSet::constructDecayAmps()
{
	_decayAmps.clear();
	if (_waveDescs.size() != nmbWaves()) {
		printWarn << "size of wave description array (= " << _waveDescs.size() << ") "
		          << "!= number of waves (= " << nmbWaves() << "). "
		          << "cannot construct isobar decay amplitudes" << endl;
		return false;
	}
	unsigned int countFail = 0;
	_decayAmps.resize(nmbWaves(), isobarAmplitudePtr());
	for (unsigned int i = 0; i < nmbWaves(); ++i) {
		if (_waveDescs[i])
			if (not _waveDescs[i]->constructAmplitude(_decayAmps[i])) {
//!!! improve reporting
				printWarn << "problems constructing decay amplitude for wave description[" << i << "]" << endl;
				++countFail;
			}
	}
	if (countFail == 0) {
		printSucc << "constructed decay amplitude objects for all " << nmbWaves() << " waves" << endl;
		return true;
	}
	printWarn << "failed to construct " << countFail << " out of " << nmbWaves()
	          << " decay amplitude objects" << endl;
	return false;
}


bool
waveSet::buildWaveNames()
{
	return false;
}


ostream&
waveSet::print(ostream& out) const
{
	out << "wave set" << endl;
	for (unsigned int i = 0; i < nmbWaves(); ++i)
		out << "    [" << i << "]: tree name = '" << _decayAmpTreeNames[i] << "', mass range = "
		    << "[" << _decayAmpMassRanges[i].first << ", " << _decayAmpMassRanges[i].second << "] MeV/c^2"
		    << endl;
	return out;
}
