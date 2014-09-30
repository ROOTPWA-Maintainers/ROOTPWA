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
//
// Description:
//      creates tree with absolute or relative differences of two trees
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <string>
#include <cassert>
#include <algorithm>

#include <boost/progress.hpp>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TVector3.h"

#include "reportingUtils.hpp"
#include "reportingUtilsRoot.hpp"
#include "evtTreeHelper.h"


using namespace std;
using namespace boost;
using namespace rpwa;


void
fillDiffLeafs(const TClonesArray& partNames,
              TClonesArray*       partMomenta[3],
              const unsigned int  index,
              const bool          absoluteDiff,
              const bool          debug = false)
{
	const TVector3* momenta[2] = {((TVector3*)(*(partMomenta[0]))[index]),
	                              ((TVector3*)(*(partMomenta[1]))[index])};
	TVector3 diff;
	if (absoluteDiff)
		diff = *(momenta[0]) - *(momenta[1]);
	else {
		diff.SetX((momenta[0]->X() - momenta[1]->X()) / momenta[0]->X());
		diff.SetY((momenta[0]->Y() - momenta[1]->Y()) / momenta[0]->Y());
		diff.SetZ((momenta[0]->Z() - momenta[1]->Z()) / momenta[0]->Z());
	}
	new((*(partMomenta[2]))[index]) TVector3(diff);

	if (debug) {
		const TObjString name = *((TObjString*)partNames[index]);
		printDebug << "particle[" << index << "]: '" << name.GetString().Data() << "', "
		           << *(momenta[0]) << " vs. " << *(momenta[1]) << ", "
		           << ((absoluteDiff) ? "absolute" : "relative" ) << " difference " << diff << endl;
	}
}


bool
createDiffTree(const string&  inFileNamePatternA       = "testEvents.root",
               const string&  inFileNamePatternB       = "testEvents2.root",
               const string&  outFileName              = "testDiffTree.root",
               const bool     absoluteDiff             = true,
               const long int maxNmbEvents             = -1,
               const string&  inTreeName               = "rootPwaEvtTree",
               const string&  prodKinPartNamesObjName  = "prodKinParticles",
               const string&  prodKinMomentaLeafName   = "prodKinMomenta",
               const string&  decayKinPartNamesObjName = "decayKinParticles",
               const string&  decayKinMomentaLeafName  = "decayKinMomenta",
               const bool     debug                    = false,
               const long int treeCacheSize            = 25000000)  // 25 MByte ROOT tree read cache
{
	// open input files
	const string  inFileNamePatterns[2] = {inFileNamePatternA, inFileNamePatternB};
	TTree*        inTrees           [2] = {0, 0};
	TClonesArray* prodKinPartNames  [2] = {0, 0};
	TClonesArray* decayKinPartNames [2] = {0, 0};
	for (unsigned int i = 0; i < 2; ++i) {
		vector<string> rootFileNames;
		rootFileNames.push_back(inFileNamePatterns[i]);
		vector<string> evtFileNames;
		vector<TTree*> trees;
		if (not openRootEvtFiles(trees, prodKinPartNames[i], decayKinPartNames[i],
		                         rootFileNames, evtFileNames,
		                         inTreeName, prodKinPartNamesObjName, prodKinMomentaLeafName,
		                         decayKinPartNamesObjName, decayKinMomentaLeafName, debug)) {
			printErr << "problems opening input file(s). exiting." << endl;
			return false;
		}
		inTrees[i] = trees[0];
	}

	// check compatibility of input data
	const long int nmbEventsInTrees[2] = {inTrees[0]->GetEntries(), inTrees[1]->GetEntries()};
	if (nmbEventsInTrees[0] != nmbEventsInTrees[1])
		printWarn << "trees have different number of events: " << nmbEventsInTrees[0] << " vs. "
		          << nmbEventsInTrees[1] << endl;
	const int nmbProdKinPart [2] = {prodKinPartNames [0]->GetEntriesFast(),
	                                prodKinPartNames [1]->GetEntriesFast()};
	const int nmbDecayKinPart[2] = {decayKinPartNames[0]->GetEntriesFast(),
	                                decayKinPartNames[1]->GetEntriesFast()};
	assert(nmbProdKinPart [0] == nmbProdKinPart [1]);
	assert(nmbDecayKinPart[0] == nmbDecayKinPart[1]);
	for (int i = 0; i < nmbProdKinPart[0]; i++)
		assert(   ((TObjString*)(*(prodKinPartNames[0]))[i])->GetString()
		       == ((TObjString*)(*(prodKinPartNames[1]))[i])->GetString());
	for (int i = 0; i < nmbDecayKinPart[0]; i++)
		assert(   ((TObjString*)(*(decayKinPartNames[0]))[i])->GetString()
		       == ((TObjString*)(*(decayKinPartNames[1]))[i])->GetString());

	// create output file
	printInfo << "creating output file '" << outFileName << "'" << endl;
	TFile* outFile = TFile::Open(outFileName.c_str(), "RECREATE");
	if (not outFile) {
		printWarn << "cannot open output file '" << outFileName << "'" << endl;
		return false;
	}

	// copy particle names to output file
	prodKinPartNames [0]->Write(prodKinPartNamesObjName.c_str (), TObject::kSingleKey);
	decayKinPartNames[0]->Write(decayKinPartNamesObjName.c_str(), TObject::kSingleKey);

	// create output tree
	TTree* tree = new TTree(inTreeName.c_str(), inTreeName.c_str());
	if (not tree) {
		printErr << "problems creating tree '" << inTreeName << "' "
		         << "in file '" << outFileName << "'" << endl;
		return false;
	}

	// create leaf variables
	TClonesArray* prodKinMomenta [3] = {0, 0, new TClonesArray("TVector3", 0)};
	TClonesArray* decayKinMomenta[3] = {0, 0, new TClonesArray("TVector3", 0)};

	// connect leaf variables to input tree branches
	for (unsigned int i = 0; i < 2; ++i) {
		inTrees[i]->SetBranchAddress(prodKinMomentaLeafName.c_str(),  &prodKinMomenta [i]);
		inTrees[i]->SetBranchAddress(decayKinMomentaLeafName.c_str(), &decayKinMomenta[i]);
	}

	// connect leaf variables to output tree branches
	const int splitLevel = 99;
	const int bufSize    = 256000;
	tree->Branch(prodKinMomentaLeafName.c_str(),  "TClonesArray", &prodKinMomenta [2], bufSize, splitLevel);
	tree->Branch(decayKinMomentaLeafName.c_str(), "TClonesArray", &decayKinMomenta[2], bufSize, splitLevel);

	// loop over events
	long int nmbEvents = min(nmbEventsInTrees[0], nmbEventsInTrees[1]);
	if (maxNmbEvents > 0)
		nmbEvents = min(maxNmbEvents, nmbEvents);
	progress_display* progressIndicator = (not debug) ? new progress_display(nmbEvents, cout, "") : 0;
	for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
		if (progressIndicator)
			++(*progressIndicator);

		if ((inTrees[0]->LoadTree(eventIndex) < 0) or (inTrees[1]->LoadTree(eventIndex) < 0))
			break;
		inTrees[0]->GetEntry(eventIndex);
		inTrees[1]->GetEntry(eventIndex);

		assert(    (prodKinMomenta[0]->GetEntriesFast() == nmbProdKinPart[0])
		       and (prodKinMomenta[1]->GetEntriesFast() == nmbProdKinPart[0]));
		if (debug)
			printDebug << "event[" << eventIndex << "]: " << nmbProdKinPart[0]
			           << " production kinematics particles:" << endl;
		prodKinMomenta[2]->Clear();
		for (int i = 0; i < nmbProdKinPart[0]; ++i)
			fillDiffLeafs(*(prodKinPartNames[0]), prodKinMomenta, i, absoluteDiff, debug);

		assert(    (decayKinMomenta[0]->GetEntriesFast() == nmbDecayKinPart[0])
		       and (decayKinMomenta[1]->GetEntriesFast() == nmbDecayKinPart[0]));
		if (debug)
			printDebug << "event[" << eventIndex << "]: " << nmbDecayKinPart[0]
			           << " decay kinematics particles:" << endl;
		decayKinMomenta[2]->Clear();
		for (int i = 0; i < nmbDecayKinPart[0]; ++i)
			fillDiffLeafs(*(decayKinPartNames[0]), decayKinMomenta, i, absoluteDiff, debug);

		tree->Fill();
		if (debug)
			cout << endl;
	}

	tree->OptimizeBaskets(treeCacheSize, 1, "d");
	tree->Write();
	outFile->Close();
	printSucc << "wrote " << nmbEvents << " events to file '" << outFileName << "'" << endl;
	return true;
}
