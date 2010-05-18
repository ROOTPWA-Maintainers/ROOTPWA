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
//      program that converts standrad ASCII PWA2000 .evt files into
//      new ROOT tree format
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <fstream>
#include <sstream>
#include <string>
#include <cassert>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TVector3.h"

#include "pputil.h"
#include "particleProperties.h"
#include "utilities.h"


using namespace std;
using namespace rpwa;


string
getParticleName(const int id,
		const int charge)
{
  assert((charge == -1) || (charge == 0) || (charge == +1));
  string name = id2name((Geant_ID)id);
  name = particleProperties::stripChargeFromName(name);
  stringstream n;
  n << name << sign(charge);
  return n.str();
}


bool
convertEvtToTree(//const string&  evtFileName           = "1500.1540.3pi.evt",
		 const string&  evtFileName           = "1720.1840.5pi.evt",
		 const string&  outFileName           = "testEvents.root",
		 // const string&  evtFileName           = "testTree.evt",
		 // const string&  outFileName           = "testEvents2.root",
		 const long int maxNmbEvents          = -1,
		 const string&  outTreeName           = "rootPwaEvtTree",
		 const string&  leafNameIsPartNames   = "initialStateNames",
		 const string&  leafNameIsPartMomenta = "initialStateMomenta",
		 const string&  leafNameFsPartNames   = "finalStateNames",
		 const string&  leafNameFsPartMomenta = "finalStateMomenta",
		 const bool     debug                 = false)
{
  // open input file
  printInfo << "opening input file '" << evtFileName << "'" << endl;
  ifstream evtFile(evtFileName.c_str());
  if (!evtFile || !evtFile.good()) {
    printWarn << "cannot open input file '" << evtFileName << "'" << endl;
    return false;
  }

  // create output file
  printInfo << "creating output file '" << outFileName << "'" << endl;
  TFile* outFile = TFile::Open(outFileName.c_str(), "RECREATE");
  if (!outFile) {
    printErr << "cannot open output file '" << outFileName << "'" << endl;
    return false;
  }

  // create tree
  TTree* tree = new TTree(outTreeName.c_str(), outTreeName.c_str());
  if (!tree) {
    printErr << "problems creating tree '" << outTreeName << "' "
	     << "in file '" << outFileName << "'" << endl;
    return false;
  }

  // create leaf variables
  TClonesArray* initialStateNames   = new TClonesArray("TObjString", 1);
  TClonesArray* initialStateMomenta = new TClonesArray("TVector3",   1);
  TClonesArray* finalStateNames     = new TClonesArray("TObjString", 0);
  TClonesArray* finalStateMomenta   = new TClonesArray("TVector3",   0);

  // connect leaf variables to tree branches
  const int split   = 0;
  const int bufSize = 256000;
  tree->Branch(leafNameIsPartNames.c_str(),   "TClonesArray", &initialStateNames,   bufSize, split);
  tree->Branch(leafNameIsPartMomenta.c_str(), "TClonesArray", &initialStateMomenta, bufSize, split);
  tree->Branch(leafNameFsPartNames.c_str(),   "TClonesArray", &finalStateNames,     bufSize, split);
  tree->Branch(leafNameFsPartMomenta.c_str(), "TClonesArray", &finalStateMomenta,   bufSize, split);

  // loop over events and fill tree
  bool     success     = true;
  long int countEvents = 0;
  long int countLines  = 0;
  while (evtFile.good()) {

    // read event
    string line;

    // read number of particles
    int nmbParticles = 0;
    if (getline(evtFile, line)) {
      ++countLines;
      stringstream lineStream(line);
      int n;
      if (lineStream >> n)
    	nmbParticles = n;
      else {
    	printWarn << "event " << countEvents + 1 << ": error reading number of particles "
		  << "from line " << countLines << ": " << line << endl;
    	success = false;
      }
    } else
      break;
    assert(nmbParticles > 0);
    if (debug)
      printInfo << "# of particles = " << nmbParticles << endl;
    
    // read initial state particle data (beam only)
    initialStateNames->Clear  ();
    initialStateMomenta->Clear();
    if (getline(evtFile, line)) {
      ++countLines;
      stringstream lineStream(line);
      int    id = 0, charge = 0;
      double momX = 0, momY = 0, momZ = 0, E = 0;
      if (lineStream >> id >> charge >> momX >> momY >> momZ >> E) {
      	const string name = getParticleName(id , charge);
      	new((*initialStateNames  )[0]) TObjString(name.c_str());
      	new((*initialStateMomenta)[0]) TVector3  (momX, momY, momZ);
      } else {
      	printWarn << "event " << countEvents + 1 << ": error reading beam data "
      		  << "from line " << countLines << ": " << line << endl;
      	success = false;
      }
    } else
      break;
    assert(initialStateNames->GetEntriesFast() == initialStateMomenta->GetEntriesFast());
    const unsigned int nmbIsPart = initialStateNames->GetEntriesFast();
    if (debug) {
      printInfo << nmbIsPart << " initial state particles:" << endl;
      for (unsigned int i = 0; i < nmbIsPart; ++i)
      	cout << "        particle[" << i << "]: "
      	     << ((TObjString*)(*initialStateNames)[i])->GetString() << "; "
      	     << *((TVector3*)(*initialStateMomenta)[i]) << endl;
    }

    // read final state particle data
    finalStateNames->Clear  ();
    finalStateMomenta->Clear();
    for (int i = 0; i < nmbParticles - 1; ++i) {
      if (getline(evtFile, line)) {
    	++countLines;
    	stringstream lineStream(line);
    	int    id, charge;
    	double momX, momY, momZ, E;
    	if (lineStream >> id >> charge >> momX >> momY >> momZ >> E) {
	  const string name = getParticleName(id , charge);
    	  new((*finalStateNames  )[i]) TObjString(name.c_str());
    	  new((*finalStateMomenta)[i]) TVector3  (momX, momY, momZ);
    	} else {
    	  printWarn << "event " << countEvents + 1 << ": error reading final state "
		    << "particle[" << i << "] data from line " << countLines << ": " << line << endl;
    	  success = false;
    	}
      } else
	break;
    }
    assert(finalStateNames->GetEntriesFast() == finalStateMomenta->GetEntriesFast());
    const unsigned int nmbFsPart = finalStateNames->GetEntriesFast();
    if (debug) {
      printInfo << nmbFsPart << " final state particles:" << endl;
      for (unsigned int i = 0; i < nmbFsPart; ++i)
      	cout << "        particle[" << i << "]: "
      	     << ((TObjString*)(*finalStateNames)[i])->GetString() << "; "
      	     << *((TVector3*)(*finalStateMomenta)[i]) << endl;
    }

    tree->Fill();
    ++countEvents;
    if ((maxNmbEvents > 0) && (countEvents >= maxNmbEvents))
      break;
    if (debug)
      cout << endl;
  }

  tree->Write();
  tree->OptimizeBaskets(10000000, 1.1, "d");
  outFile->Close();
  printInfo << "read " << countLines << " lines from file '" << evtFileName << "' "
	    << "and wrote " << countEvents << " events to file '" << outFileName << "'" << endl;
  return success;
}
