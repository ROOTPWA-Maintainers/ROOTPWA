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

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TVector3.h"

#include "utilities.h"


using namespace std;


void
createDiffLeafs(TClonesArray*      partNames  [3],
		TClonesArray*      partMomenta[3],
		const unsigned int index,
		const bool         absoluteDiff,
		const bool         debug = false)
{
  const string   names  [2] = {((TObjString*)(*(partNames[0]))[index])->GetString().Data(),
			       ((TObjString*)(*(partNames[1]))[index])->GetString().Data()};
  const TVector3 momenta[2] = {*((TVector3*)(*(partMomenta[0]))[index]),
			       *((TVector3*)(*(partMomenta[1]))[index])};
  assert(names[0] == names[1]);
  new((*(partNames[2]))[index]) TObjString(names[0].c_str());
  
  TVector3 diff;
  if (absoluteDiff)
    diff = momenta[0] - momenta[1];
  else {
    diff.SetX((momenta[0].X() - momenta[1].X()) / momenta[0].X());
    diff.SetY((momenta[0].Y() - momenta[1].Y()) / momenta[0].Y());
    diff.SetZ((momenta[0].Z() - momenta[1].Z()) / momenta[0].Z());
  }
  new((*(partMomenta[2]))[index]) TVector3(diff);
  
  if (debug) {
	cout << "        particle[" << index << "]: " << names[0] << " vs. " << names[1] << ", "
	     << momenta[0] << " vs. " << momenta[1] << endl;
  }
}


bool
createDiffTree(const string&  inFileNamePatternA    = "testEvents.root",
	       const string&  inFileNamePatternB    = "testEvents2.root",
	       const string&  outFileName           = "testDiffTree.root",
	       const bool     absoluteDiff          = true,
	       const long int maxNmbEvents          = -1,
	       const string&  inTreeName            = "rootPwaEvtTree",
	       const string&  leafNameIsPartNames   = "initialStateNames",
	       const string&  leafNameIsPartMomenta = "initialStateMomenta",
	       const string&  leafNameFsPartNames   = "finalStateNames",
	       const string&  leafNameFsPartMomenta = "finalStateMomenta",
	       const bool     debug                 = false)
{
  
  // open input files
  TChain chain[2];
  printInfo << "opening input file(s) '" << inFileNamePatternA << "'" << endl;
  chain[0].SetName(inTreeName.c_str());
  if (chain[0].Add(inFileNamePatternA.c_str()) < 1) {
    printWarn << "no events in input file(s) '" << inFileNamePatternA << "'" << endl;
    return false;
  }
  printInfo << "opening input file(s) '" << inFileNamePatternB << "'" << endl;
  chain[1].SetName(inTreeName.c_str());
  if (chain[1].Add(inFileNamePatternB.c_str()) < 1) {
    printWarn << "no events in input file(s) '" << inFileNamePatternB << "'" << endl;
    return false;
  }
  const long int nmbEventsChain[2] = {chain[0].GetEntries(), chain[1].GetEntries()};
  if (nmbEventsChain[0] != nmbEventsChain[1])
    printWarn << "trees have different number of events: " << nmbEventsChain[0] << " vs. "
	      << nmbEventsChain[1] << endl;
  chain[0].GetListOfFiles()->ls();
  chain[1].GetListOfFiles()->ls();

  // create output file
  printInfo << "creating output file '" << outFileName << "'" << endl;
  TFile* outFile = TFile::Open(outFileName.c_str(), "RECREATE");
  if (!outFile) {
    printWarn << "cannot open output file '" << outFileName << "'" << endl;
    return false;
  }

  // create output tree
  TTree* tree = new TTree(inTreeName.c_str(), inTreeName.c_str());
  if (!tree) {
    printErr << "problems creating tree '" << inTreeName << "' "
	     << "in file '" << outFileName << "'" << endl;
    return false;
  }

  // create leaf variables
  TClonesArray* initialStateNames  [3] = {0, 0, new TClonesArray("TObjString", 0)};
  TClonesArray* initialStateMomenta[3] = {0, 0, new TClonesArray("TVector3",   0)};
  TClonesArray* finalStateNames    [3] = {0, 0, new TClonesArray("TObjString", 0)};
  TClonesArray* finalStateMomenta  [3] = {0, 0, new TClonesArray("TVector3",   0)};

  // connect leaf variables to input tree branches
  for (unsigned int i = 0; i < 2; ++i) {
    chain[i].SetBranchAddress(leafNameIsPartNames.c_str(),   &initialStateNames  [i]);
    chain[i].SetBranchAddress(leafNameIsPartMomenta.c_str(), &initialStateMomenta[i]);
    chain[i].SetBranchAddress(leafNameFsPartNames.c_str(),   &finalStateNames    [i]);
    chain[i].SetBranchAddress(leafNameFsPartMomenta.c_str(), &finalStateMomenta  [i]);
  }

  // connect leaf variables to output tree branches
  const int split   = 0;
  const int bufSize = 256000;
  tree->Branch(leafNameIsPartNames.c_str(),   "TClonesArray", &initialStateNames  [2], bufSize, split);
  tree->Branch(leafNameIsPartMomenta.c_str(), "TClonesArray", &initialStateMomenta[2], bufSize, split);
  tree->Branch(leafNameFsPartNames.c_str(),   "TClonesArray", &finalStateNames    [2], bufSize, split);
  tree->Branch(leafNameFsPartMomenta.c_str(), "TClonesArray", &finalStateMomenta  [2], bufSize, split);

  // loop over events
  long int nmbEvents = min(nmbEventsChain[0], nmbEventsChain[1]);
  if (maxNmbEvents > 0)
    nmbEvents = min(maxNmbEvents, nmbEvents);
  for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
    if (!debug)
      progressIndicator(eventIndex, nmbEvents);

    if ((chain[0].LoadTree(eventIndex) < 0) || (chain[1].LoadTree(eventIndex) < 0))
      break;
    chain[0].GetEntry(eventIndex);
    chain[1].GetEntry(eventIndex);

    const unsigned int nmbIsPart[2] = {initialStateNames[0]->GetEntriesFast(),
				       initialStateNames[1]->GetEntriesFast()};
    const unsigned int nmbFsPart[2] = {finalStateNames[0]->GetEntriesFast(),
				       finalStateNames[1]->GetEntriesFast()};
    assert(nmbIsPart[0] == nmbIsPart[1]);
    assert(nmbFsPart[0] == nmbFsPart[1]);

    if (debug)
      printInfo << "event[" << eventIndex << "]: " << nmbIsPart[0]
		<< " initial state particles:" << endl;
    initialStateNames  [2]->Clear();
    initialStateMomenta[2]->Clear();
    for (unsigned int i = 0; i < nmbIsPart[0]; ++i)
      createDiffLeafs(initialStateNames, initialStateMomenta, i, absoluteDiff, debug);
    
    if (debug)
      printInfo << "event[" << eventIndex << "]: " << nmbFsPart[0]
		<< " final state particles:" << endl;
    finalStateNames  [2]->Clear();
    finalStateMomenta[2]->Clear();
    for (unsigned int i = 0; i < nmbFsPart[0]; ++i)
      createDiffLeafs(finalStateNames, finalStateMomenta, i, absoluteDiff, debug);

    tree->Fill();
    if (debug)
      cout << endl;
  }

  tree->Write();
  tree->OptimizeBaskets(10000000, 1.1, "d");
  outFile->Close();
  printInfo << "wrote " << nmbEvents << " events to file '" << outFileName << "'" << endl;
  return true;
}
