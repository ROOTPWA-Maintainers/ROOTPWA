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
	: TObject           (),
	  decayAmpTreeNames (),
	  decayAmpMassRanges(),
	  decayAmpTrees     (),
	  decayAmpWaveDescs ()
{
	//waveSet::Class()->IgnoreTObjectStreamer();  // don't store TObject's fBits and fUniqueID
}


waveSet::~waveSet()
{ }


void
waveSet::clear()
{
	decayAmpTreeNames.clear ();
	decayAmpMassRanges.clear();
	decayAmpTrees.clear     ();
	decayAmpWaveDescs.clear ();
}


waveSet&
waveSet::operator =(const waveSet& set)
{
	if (this != &set) {
		TObject::operator  =(set);
		decayAmpTreeNames  = set.decayAmpTreeNames;
		decayAmpMassRanges = set.decayAmpMassRanges;
		decayAmpTrees      = set.decayAmpTrees;
		decayAmpWaveDescs  = set.decayAmpWaveDescs;
	}
	return *this;
}


bool
waveSet::parseWaveSetFile(const string& waveSetFileName)
{
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

	// read decay amplitude list
	const int nmbEntries = configDecayAmps->getLength();
	decayAmpTreeNames.clear ();
	decayAmpMassRanges.clear();
	decayAmpTreeNames.resize (nmbEntries);
	decayAmpMassRanges.resize(nmbEntries);
	bool success = true;
	for (int i = 0; i < nmbEntries; ++i) {
		// get tree name
		string treeName = "";
		if (not (*configDecayAmps)[i].lookupValue("treeName", treeName)) {
			printWarn << "entry [" << i << "] in decay amplitude list does not specify a tree name.";
			success = false;
		}
		decayAmpTreeNames[i] = treeName;
		// get mass range
		pair<double, double> massRange(0, numeric_limits<double>::infinity());
		const Setting*       configMassRange = findLibConfigArray((*configDecayAmps)[i], "massRange", false);
		if (configMassRange) {
			if (configMassRange->getLength() == 1)
				massRange.first  = (*configMassRange)[0];
			else if (configMassRange->getLength() == 2) {
				massRange.first  = (*configMassRange)[0];
				massRange.second = (*configMassRange)[1];
			} else {
				printWarn << "cannot read mass range from entry [" << i << "] in decay amplitude list. "
				          << "array length is neither 1 nor 2." << endl;
				success = false;
			}
		}
		decayAmpMassRanges[i] = massRange;
	}	

	if (success)
		printSucc << "constructed wave set from file '" << waveSetFileName << "'" << endl;
	else
		printSucc << "problems constructing wave set from file '" << waveSetFileName << "'" << endl;
	return success;
}


bool
waveSet::getDecayAmplitudeTrees(const vector<string>& ampFileNames)
{
	decayAmpTrees.clear    ();
	decayAmpWaveDescs.clear();
	bool success = true;

#ifdef USE_STD_COMPLEX_TREE_LEAFS
	// force loading predefined std::complex dictionary
	// see http://root.cern.ch/phpBB3/viewtopic.php?f=5&t=9618&p=50164
	gROOT->ProcessLine("#include <complex>");

	// get trees from .root files
	decayAmpTrees.resize    (nmbDecayAmps(), 0);
	decayAmpWaveDescs.resize(nmbDecayAmps(), 0);
	unsigned int countOpenFails = 0;
	unsigned int countDuplTrees = 0;
	for (unsigned int ampFileIndex = 0; ampFileIndex < ampFileNames.size(); ++ampFileIndex) {
		// open .root file
		if (_debug)
			printDebug << "opening .root decay amplitude file '" << ampFileNames[ampFileIndex] << "'"
			           << endl;
		TFile* ampFile = TFile::Open(ampFileNames[ampFileIndex].c_str(), "READ");
		if (not ampFile or ampFile->IsZombie()) {
			printWarn << "cannot open decay amplitude file '" << ampFileNames[ampFileIndex] << "'. "
			          << "skipping file." << endl;
			++countOpenFails;
			continue;			
		}
		// get all amplitude trees from file
		vector<pair<TTree*, unsigned int> > trees;
		for (unsigned int i = 0; i < nmbDecayAmps(); ++i) {
			if (_debug)
				printDebug << "looking for tree '" << decayAmpTreeNames[i] << "' "
				           << "in file '" << ampFileNames[ampFileIndex] << "'" << endl;
			TTree* tree = 0;
			ampFile->GetObject(decayAmpTreeNames[i].c_str(), tree);
			if (tree) {
				trees.push_back(make_pair(tree, i));
				if (_debug)
					printDebug << "found decay amplitude tree '" << decayAmpTreeNames[i] << "' "
					           << "in file '" << ampFileNames[ampFileIndex] << "'" << endl;
			} else if (_debug)
				printDebug << "did not find decay amplitude tree '" << decayAmpTreeNames[i] << "' "
				           << "in file '" << ampFileNames[ampFileIndex] << "'" << endl;
		}
		// store tree pointers and make sure trees do not already exist
		for (unsigned int i = 0; i < trees.size(); ++i) {
			if (decayAmpTrees[trees[i].second]) {
				printWarn << "tree '" << trees[i].first->GetName() << "' exists in file '"
				          << decayAmpTrees[trees[i].second]->GetDirectory()->GetName() << "' and in file '"
					// GetCurrentFile()
				          << trees[i].first->GetDirectory()->GetName() << "'. skipping tree from latter file."
				          << endl;
				++countDuplTrees;
				continue;
			}
			decayAmpTrees[trees[i].second] = trees[i].first;
		}		
	}

	// get corresponding wave descriptions from trees
	unsigned int countMissingTrees     = 0;
	unsigned int countMissingWaveDescs = 0;
	for (unsigned int i = 0; i < nmbDecayAmps(); ++i) {
		if (not decayAmpTrees[i]) {
			printWarn << "did not find decay amplitude tree '" << decayAmpTreeNames[i] << "' "
			          << "in any of the given decay amplitude files" << endl;
			++countMissingTrees;
			continue;
		}
		TList* userInfo = decayAmpTrees[i]->GetUserInfo();
		if (not userInfo) {
			printWarn << "tree '" << decayAmpTreeNames[i] << "' does not have user info." << endl;
			++countMissingWaveDescs;
			continue;
		}
		waveDescription* waveDesc
			= static_cast<waveDescription*>(userInfo->FindObject("rpwa::waveDescription"));
		if (not waveDesc) {
			printInfo << "user info of tree '" << decayAmpTreeNames[i] << "' "
			          << "does not contain wave description." << endl;
			++countMissingWaveDescs;
			continue;
		}
		decayAmpWaveDescs[i] = waveDesc;
	}

	// report
	if (countOpenFails > 0) {
		success = false;
		printWarn << "problems opening " << countOpenFails << " out of " << ampFileNames.size()
		          << " decay amplitude files." << endl;
	}
	if (countDuplTrees > 0) {
		success = false;
		printWarn << "found " << countDuplTrees << " duplicate trees. "
		          << "decay amplitude data were ignored." << endl;
	}
	if (countMissingTrees > 0) {
		success = false;
		printWarn << countMissingTrees << " out of " << nmbDecayAmps() << " trees are missing. "
		          << "decay amplitude data incomplete." << endl;
	}
	if (countMissingWaveDescs > 0) {
		success = false;
		printWarn << countMissingWaveDescs << " out of " << nmbDecayAmps() << " wave descriptions "
		          << "are missing. wave set data incomplete." << endl;
	}
	if (    (countOpenFails    == 0) and (countDuplTrees        == 0)
	    and (countMissingTrees == 0) and (countMissingWaveDescs == 0))
		printSucc << "could open all " << ampFileNames.size() << " decay amplitude files "
		          << "and found all " << nmbDecayAmps() << " decay amplitude trees "
		          << "and wave descriptions" << endl;

#endif

	return success;
}


ostream&
waveSet::print(ostream& out) const
{
	out << "wave set" << endl;
	for (unsigned int i = 0; i < decayAmpTreeNames.size(); ++i)
		out << "    [" << i << "]: tree name = '" << decayAmpTreeNames[i] << "', mass range = "
		    << "[" << decayAmpMassRanges[i].first << ", " << decayAmpMassRanges[i].second << "] MeV/c^2"
		    << endl;
	return out;
}
