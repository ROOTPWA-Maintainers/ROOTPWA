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
//      checks validity of amplitude specified by key file by
//      performing a number of consistency checks
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <iostream>
#include <fstream>
#include <unistd.h>
#include <vector>
#include <complex>

#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TStopwatch.h"

#include "svnVersion.h"
#include "utilities.h"
#include "particleDataTable.h"
#include "evtTreeHelper.h"
#include "keyFileParser.h"


using namespace std;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{
  cerr << "usage:" << endl
       << progName
       << " -d test data [-p PDG file -t tree name -l leaf names -v -h] key file(s)" << endl
       << "    where:" << endl
       << "        -d file    path to file with test data (.evt or ROOT format)" << endl
       << "        -p file    path to particle data table file (default: ./particleDataTable.txt)" << endl
       << "        -t name    name of tree in ROOT data files (default: rootPwaEvtTree)"
       << "        -l names   semicolon separated tree leaf names (default: 'prodKinParticles;prodKinMomenta;decayKinParticles;decayKinMomenta')"
       << "        -v         verbose; print debug output (default: false)" << endl
       << "        -h         print help" << endl
       << endl;
  exit(errCode);
}


int
main(int    argc,
     char** argv)
{
  printCompilerInfo();
  printSvnVersion();

  // parse command line options
  const string progName     = argv[0];
  string       dataFileName = "";
  string       pdgFileName  = "./particleDataTable.txt";
  string       inTreeName   = "rootPwaEvtTree";
  string       leafNames    = "prodKinParticles;prodKinMomenta;"
                              "decayKinParticles;decayKinMomenta";
  bool         debug        = false;
  extern char* optarg;
  extern int   optind;
  int          c;
  while ((c = getopt(argc, argv, "d:p:t:l:vh")) != -1)
    switch (c) {
    case 'd':
      dataFileName = optarg;
      break;
    case 'p':
      pdgFileName = optarg;
      break;
    case 't':
      inTreeName = optarg;
      break;
    case 'l':
      leafNames = optarg;
      break;
    case 'v':
      debug = true;
      break;
    case 'h':
    default:
      usage(progName);
    }

  // get key file names
  if (optind >= argc) {
    printErr << "you need to specify at least one key file to process. aborting." << endl;;
    usage(progName, 1);
  }
  vector<string> keyFileNames;
  while (optind < argc) {
    const string fileName = argv[optind++];
    keyFileNames.push_back(fileName);
  }

  // get leaf names
  const vector<string> leafNameTokens            = tokenizeString(leafNames, ";");
  const string         prodKinParticlesLeafName  = leafNameTokens[0];
  const string         prodKinMomentaLeafName    = leafNameTokens[1];
  const string         decayKinParticlesLeafName = leafNameTokens[2];
  const string         decayKinMomentaLeafName   = leafNameTokens[3];
  // printInfo << "using the following leaf names:" << endl
  // 	    << "        production kinematics: particle names = '" << prodKinParticlesLeafName << "', "
  // 	    << "momenta = '" << prodKinMomentaLeafName << "'" << endl
  // 	    << "        decay kinematics:      particle names = '" << decayKinParticlesLeafName << "', "
  // 	    << "momenta = '" << decayKinMomentaLeafName << "'" << endl;

  // determine test data format
  bool rootDataFormat = false;
  if (dataFileName.substr(dataFileName.length() - 5) == ".root")
    rootDataFormat = true;

  TTree* tree = 0;
  if (rootDataFormat) {
    // open root file and build chain
    TChain* chain = new TChain(inTreeName.c_str());
    printInfo << "opening ROOT input file '" << dataFileName << "'" << endl;
    if (chain->Add(dataFileName.c_str()) < 1)
      printWarn << "no events in ROOT input file '" << dataFileName << "'" << endl;
    chain->GetListOfFiles()->ls();
    tree = chain;
  } else {
    // convert .evt file to root tree
    printInfo << "opening .evt input file '" << dataFileName << "'" << endl;
    ifstream evtFile(dataFileName.c_str());
    if (!evtFile || !evtFile.good()) {
      printErr << "cannot open .evt input file '" << dataFileName << "'. aborting." << endl;
      exit(1);
    }
    // create tree
    tree = new TTree(inTreeName.c_str(), inTreeName.c_str());
    if (!tree) {
      printErr << "problems creating tree '" << inTreeName << "'. aborting." << endl;
      exit(1);
    }
    if (not fillTreeFromEvt(evtFile, *tree, -1,
			    prodKinParticlesLeafName,  prodKinMomentaLeafName,
			    decayKinParticlesLeafName, decayKinMomentaLeafName, debug)) {
      printErr << "problems creating tree from .evt input file '" << dataFileName << "' "
	       << "aborting." << endl;
      exit(1);
    }
  }

  // initialize particle data table
  particleDataTable& pdt = particleDataTable::instance();
  pdt.readFile(pdgFileName);

  // loop over key files
  map<string, vector<string> > keyFileErrors;  // maps error description to key files
  for (unsigned int i = 0; i < keyFileNames.size(); ++i) {
    printInfo << "checking key file '" << keyFileNames[i] << "'" << endl;
    // parse key file and create decay topology and amplitude instances
    keyFileParser&         parser = keyFileParser::instance();
    isobarDecayTopologyPtr decayTopo;
    if (not parser.parse(keyFileNames[i], decayTopo)) {
      printWarn << "problems constructing decay topology from key file '" << keyFileNames[i] << "'. "
		<< "skipping." << endl;
      keyFileErrors["parsing errors"].push_back(keyFileNames[i]);
      continue;
    }
    // check topology
    if (not decayTopo->checkTopology())
      keyFileErrors["problematic topology"].push_back(keyFileNames[i]);
    if (not decayTopo->checkConsistency())
      keyFileErrors["inconsistent decay"].push_back(keyFileNames[i]);
    isobarHelicityAmplitude amplitude(decayTopo);
    parser.setAmplitudeOptions(amplitude);
    
    // read data from tree, calculate amplitudes, and write them to output file
    vector<complex<double> > ampValues;
    if (not processTree(*tree, *decayTopo, amplitude, ampValues,
			prodKinParticlesLeafName,  prodKinMomentaLeafName,
			decayKinParticlesLeafName, decayKinMomentaLeafName, debug))
      printWarn << "problems reading tree" << endl;
  }

  // count key files with errors
  unsigned int countKeyFileErr = 0;
  for (unsigned int i = 0; i < keyFileNames.size(); ++i)
    for (map<string, vector<string> >::const_iterator entry = keyFileErrors.begin();
	 entry != keyFileErrors.end(); ++entry) {
      bool foundErr = false;
      for (unsigned int j = 0; j < entry->second.size(); ++j)
	if (entry->second[j] == keyFileNames[i]) {
	  ++countKeyFileErr;
	  foundErr = true;
	  break;
	}
      if (foundErr)
	break;
    }
  
  printInfo << keyFileNames.size() - countKeyFileErr << " of " << keyFileNames.size()
	    << " keyfile(s) passed all tests" << endl;
  if (countKeyFileErr > 0) {
    printInfo << countKeyFileErr << " problematic keyfile(s):" << endl;
    for (map<string, vector<string> >::const_iterator entry = keyFileErrors.begin();
	 entry != keyFileErrors.end(); ++entry) {
      cout << "        " << entry->first << ":" << endl;
      for (unsigned int j = 0; j < entry->second.size(); ++j)
	cout << "            " << entry->second[j] << endl;
    }
  }

  // clean up
  delete tree;
  return 0;
}
