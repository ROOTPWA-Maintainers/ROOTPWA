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
//      reads in data files in .evt or tree format, calculates
//      amplitudes for each event based on given key file, and writes
//      out amplitudes in PWA2000 binary or ascii format
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
       << " -k key file [-p PDG file -t tree name -o output file -a -l leaf names -v -h] input data file(s)" << endl
       << "    where:" << endl
       << "        -k file    path to key file" << endl
       << "        -p file    path to particle data table file (default: ./particleDataTable.txt)" << endl
       << "        -t name    name of tree in ROOT data files (default: rootPwaEvtTree)"
       << "        -o file    path to amplitude file (default: ./out.amp)" << endl
       << "        -a         write amplitudes in ASCII format (default: binary)" << endl
       << "        -l names   semicolon separated tree leaf names (default: 'initialStateNames;initialStateMomenta;finalStateNames;finalStateMomenta')"
       << "        -v         verbose; print debug output (default: false)" << endl
       << "        -h         print help" << endl
       << endl;
  exit(errCode);
}


int
main(int    argc,
     char** argv)
{
  // parse command line options
  const string progName     = argv[0];
  string       keyFileName  = "";
  string       pdgFileName  = "./particleDataTable.txt";
  string       inTreeName   = "rootPwaEvtTree";
  string       ampFileName  = "./out.amp";
  bool         asciiOutput  = false;
  string       leafNames    = "initialStateNames;initialStateMomenta;"
                              "finalStateNames;finalStateMomenta";
  bool         debug        = false;
  extern char* optarg;
  extern int   optind;
  int          c;
  while ((c = getopt(argc, argv, "k:p:t:o:al:vh")) != -1)
    switch (c) {
    case 'k':
      keyFileName = optarg;
      break;
    case 'p':
      pdgFileName = optarg;
      break;
    case 't':
      inTreeName = optarg;
      break;
    case 'o':
      ampFileName = optarg;
      break;
    case 'a':
      asciiOutput = true;
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

  // get input file names
  if (optind >= argc) {
    printErr << "you need to specify at least one data file to process. exiting." << endl;;
    usage(progName, 1);
  }
  vector<string> rootFileNames;
  vector<string> evtFileNames;
  while (optind < argc) {
    const string fileName = argv[optind++];
    if (fileName.substr(fileName.length() - 5) == ".root")
      rootFileNames.push_back(fileName);
    else if (fileName.substr(fileName.length() - 4) == ".evt")
      evtFileNames.push_back(fileName);
    else
      printWarn << "input file '" << fileName << "' is neither a .root nor a .evt file. "
		<< "skipping." << endl;
  }

  if ((rootFileNames.size() == 0) and (evtFileNames.size() == 0)) {
    printErr << "specified input files are neither .root nor .evt files. exiting.";
    usage(progName, 1);
  }

  // get leaf names
  const vector<string> leafNameTokens        = tokenizeString(leafNames, ";");
  const string         leafNameIsPartNames   = leafNameTokens[0];
  const string         leafNameIsPartMomenta = leafNameTokens[1];
  const string         leafNameFsPartNames   = leafNameTokens[2];
  const string         leafNameFsPartMomenta = leafNameTokens[3];
  printInfo << "using the following leaf names:" << endl
	    << "        production kinematics: particle names = '" << leafNameIsPartNames << "', "
	    << "momenta = '" << leafNameIsPartMomenta << "'" << endl
	    << "        decay final state:     particle names = '" << leafNameFsPartNames << "', "
	    << "momenta = '" << leafNameFsPartMomenta << "'" << endl;

  // open root files and build chain
  TChain* chain = 0;
  if (rootFileNames.size() > 0) {
    chain = new TChain(inTreeName.c_str());
    for (unsigned int i = 0; i < rootFileNames.size(); ++i) {
      printInfo << "opening ROOT input file '" << rootFileNames[i] << "'" << endl;
      if (chain->Add(rootFileNames[i].c_str()) < 1)
	printWarn << "no events in ROOT input file '" << rootFileNames[i] << "'" << endl;
    }
    chain->GetListOfFiles()->ls();
  }

  // convert .evt files to root trees
  vector<TTree*> trees;
  if (chain)
    trees.push_back(chain);
  for (unsigned int i = 0; i < evtFileNames.size(); ++i) {
    printInfo << "opening .evt input file '" << evtFileNames[i] << "'" << endl;
    ifstream evtFile(evtFileNames[i].c_str());
    if (!evtFile || !evtFile.good()) {
      printWarn << "cannot open .evt input file '" << evtFileNames[i] << "'. skipping." << endl;
      continue;
    }
    printInfo << "converting .evt input file '" << evtFileNames[i] << "' "
	      << "into memory resident tree. this might reduce performance. "
	      << "ROOT input format is recommended." << endl;
    // create tree
    TTree* tree = new TTree(inTreeName.c_str(), inTreeName.c_str());
    if (!tree) {
      printErr << "problems creating tree '" << inTreeName << "'. skipping." << endl;
      continue;
    }
    if (fillTreeFromEvt(evtFile, *tree, -1,
			leafNameIsPartNames, leafNameIsPartMomenta,
			leafNameFsPartNames, leafNameFsPartMomenta, debug))
      trees.push_back(tree);
    else {
      printWarn << "problems creating tree from .evt input file '" << evtFileNames[i] << "' "
		<< "skipping." << endl;
    }
  }

  // initialize particle data table
  particleDataTable& pdt = particleDataTable::instance();
  pdt.readFile(pdgFileName);

  // parse key file and create decay topology and amplitude instances
  if (keyFileName == "") {
    printErr << "no key file specified. exiting." << endl;
    usage(progName, 1);
  }
  keyFileParser&         parser = keyFileParser::instance();
  isobarDecayTopologyPtr decayTopo;
  if (not parser.parse(keyFileName, decayTopo)) {
    printErr << "problems constructing decay topology from key file '" << keyFileName << "'. "
	     << "exiting." << endl;
    exit(1);
  }
  printInfo << *decayTopo;
  decayTopo->checkTopology();
  decayTopo->checkConsistency();
  isobarHelicityAmplitude amplitude(*decayTopo);
  parser.setAmplitudeOptions(amplitude);
  
  // create output file for amplitudes
  printInfo << "creating amplitude file '" << ampFileName << "'" << endl;
  ofstream ampFile(ampFileName.c_str());
  if (!ampFile) {
    printErr << "cannot create amplitude file '" << ampFileName << "'. exiting." << endl;
    exit(1);
  }

  // read data from tree(s) and calculate amplitudes
  TStopwatch timer;
  timer.Reset();
  timer.Start();
  long int countEvents = 0;
  for (unsigned int i = 0; i < trees.size(); ++i) {
    // create branch pointers and leaf variables
    TBranch*      initialStateNamesBr   = 0;
    TBranch*      initialStateMomentaBr = 0;
    TBranch*      finalStateNamesBr     = 0;
    TBranch*      finalStateMomentaBr   = 0;
    TClonesArray* initialStateNames     = 0;
    TClonesArray* initialStateMomenta   = 0;
    TClonesArray* finalStateNames       = 0;
    TClonesArray* finalStateMomenta     = 0;
	
    // connect leaf variables to tree branches
    trees[i]->SetBranchAddress(leafNameIsPartNames.c_str  (), &initialStateNames,   &initialStateNamesBr  );
    trees[i]->SetBranchAddress(leafNameIsPartMomenta.c_str(), &initialStateMomenta, &initialStateMomentaBr);
    trees[i]->SetBranchAddress(leafNameFsPartNames.c_str  (), &finalStateNames,     &finalStateNamesBr    );
    trees[i]->SetBranchAddress(leafNameFsPartMomenta.c_str(), &finalStateMomenta,   &finalStateMomentaBr  );

    // loop over events
    const long int nmbEvents = trees[i]->GetEntries();
    printInfo << "processing ";
    if (chain and (i == 0)) 
      cout << "chain of .root files";
    else
      cout << ".evt tree[" << ((chain) ? i : i + 1) << "]";
    cout << endl;
    for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
      progressIndicator(eventIndex, nmbEvents);
      
      if (trees[i]->LoadTree(eventIndex) < 0)
	break;
      // read only required branches
      initialStateNamesBr->GetEntry  (eventIndex);
      initialStateMomentaBr->GetEntry(eventIndex);
      finalStateNamesBr->GetEntry    (eventIndex);
      finalStateMomentaBr->GetEntry  (eventIndex);

      if (   not initialStateNames or not initialStateMomenta
	  or not finalStateNames   or not finalStateMomenta) {
	printWarn << "at least one data array is null pointer: "
		  << "initialStateNames = "   << initialStateNames   << ", "
		  << "initialStateMomenta = " << initialStateMomenta << ", "
		  << "finalStateNames = "     << finalStateNames     << ", "
		  << "finalStateMomenta = "   << finalStateMomenta   << ". "
		  << "skipping event." << endl;
	continue;
      }

      if (decayTopo->readData(*initialStateNames, *initialStateMomenta,
			      *finalStateNames,   *finalStateMomenta)) {
      	const complex<double> amp = amplitude();
	if (asciiOutput)
	  ampFile  << setprecision(numeric_limits<double>::digits10 + 1) << amp << endl;
	else
	  ampFile.write((char*)(&amp), sizeof(complex<double>));
	++countEvents;
      }
    }
  }
  
  timer.Stop();
  printInfo << "successfully calculated amplitudes for " << countEvents << " events and "
	    << "wrote them to '" << ampFileName << "'" << endl;
  printInfo << "this job consumed: ";
  timer.Print();
 
  // clean up
  for (unsigned int i = 0; i < trees.size(); ++i)
    delete trees[i];
  trees.clear();

  return 0;
}
