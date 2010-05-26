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
       << "        -l names   semicolon separated tree leaf names (default: 'prodKinParticles;prodKinMomenta;decayKinParticles;decayKinMomenta')"
       << "        -v         verbose; print debug output (default: false)" << endl
       << "        -h         print help" << endl
       << endl;
  exit(errCode);
}


long int
processTree(TTree&                         tree,
	    const isobarDecayTopologyPtr&  decayTopo,
	    const isobarHelicityAmplitude& amplitude,
	    ostream&                       ampFile,
	    const bool                     asciiOutput,
	    const string&                  prodKinParticlesLeafName  = "prodKinParticles",
	    const string&                  prodKinMomentaLeafName    = "prodKinMomenta",
	    const string&                  decayKinParticlesLeafName = "decayKinParticles",
	    const string&                  decayKinMomentaLeafName   = "decayKinMomenta",
	    const bool                     debug                     = false)
{
  // create branch pointers and leaf variables
  TBranch*      prodKinParticlesBr  = 0;
  TBranch*      prodKinMomentaBr    = 0;
  TBranch*      decayKinParticlesBr = 0;
  TBranch*      decayKinMomentaBr   = 0;
  TClonesArray* prodKinParticles    = 0;
  TClonesArray* prodKinMomenta      = 0;
  TClonesArray* decayKinParticles   = 0;
  TClonesArray* decayKinMomenta     = 0;
	
  // connect leaf variables to tree branches
  tree.SetBranchAddress(prodKinParticlesLeafName.c_str(),  &prodKinParticles,  &prodKinParticlesBr );
  tree.SetBranchAddress(prodKinMomentaLeafName.c_str(),    &prodKinMomenta,    &prodKinMomentaBr   );
  tree.SetBranchAddress(decayKinParticlesLeafName.c_str(), &decayKinParticles, &decayKinParticlesBr);
  tree.SetBranchAddress(decayKinMomentaLeafName.c_str(),   &decayKinMomenta,   &decayKinMomentaBr  );

  // loop over events
  const long int nmbEvents   = tree.GetEntries();
  long int       countEvents = 0;
  for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
    progressIndicator(eventIndex, nmbEvents);
      
    if (tree.LoadTree(eventIndex) < 0)
      break;
    // read only required branches
    prodKinParticlesBr->GetEntry (eventIndex);
    prodKinMomentaBr->GetEntry   (eventIndex);
    decayKinParticlesBr->GetEntry(eventIndex);
    decayKinMomentaBr->GetEntry  (eventIndex);

    if (   not prodKinParticles  or not prodKinMomenta
	or not decayKinParticles or not decayKinMomenta) {
      printWarn << "at least one of the input data arrays is a null pointer: "
		<< "        production kinematics: particle names = " << prodKinParticles << ", "
		<< "momenta = " << prodKinMomenta << endl
		<< "        decay kinematics:      particle names = " << decayKinParticles << ", "
		<< "momenta = " << decayKinMomenta << endl
		<< "skipping event." << endl;
      continue;
    }

    if (decayTopo->readData(*prodKinParticles,  *prodKinMomenta,
			    *decayKinParticles, *decayKinMomenta)) {
      const complex<double> amp = amplitude();
      if (asciiOutput)
	ampFile  << setprecision(numeric_limits<double>::digits10 + 1) << amp << endl;
      else
	ampFile.write((char*)(&amp), sizeof(complex<double>));
      ++countEvents;
    }
  }
  return countEvents;
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
  string       leafNames    = "prodKinParticles;prodKinMomenta;"
                              "decayKinParticles;decayKinMomenta";
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
  const vector<string> leafNameTokens            = tokenizeString(leafNames, ";");
  const string         prodKinParticlesLeafName  = leafNameTokens[0];
  const string         prodKinMomentaLeafName    = leafNameTokens[1];
  const string         decayKinParticlesLeafName = leafNameTokens[2];
  const string         decayKinMomentaLeafName   = leafNameTokens[3];
  printInfo << "using the following leaf names:" << endl
	    << "        production kinematics: particle names = '" << prodKinParticlesLeafName << "', "
	    << "momenta = '" << prodKinMomentaLeafName << "'" << endl
	    << "        decay kinematics:      particle names = '" << decayKinParticlesLeafName << "', "
	    << "momenta = '" << decayKinMomentaLeafName << "'" << endl;

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
			prodKinParticlesLeafName,  prodKinMomentaLeafName,
			decayKinParticlesLeafName, decayKinMomentaLeafName, debug))
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
  isobarHelicityAmplitude amplitude(decayTopo);
  parser.setAmplitudeOptions(amplitude);
  
  // create output file for amplitudes
  printInfo << "creating amplitude file '" << ampFileName << "'; "
	    << ((asciiOutput) ? "ASCII" : "binary") << " mode" << endl;
  ofstream ampFile(ampFileName.c_str());
  if (!ampFile) {
    printErr << "cannot create amplitude file '" << ampFileName << "'. exiting." << endl;
    exit(1);
  }

  // read data from tree(s), calculate amplitudes, and write them to output file
  TStopwatch timer;
  timer.Reset();
  timer.Start();
  long int countEvents = 0;
  for (unsigned int i = 0; i < trees.size(); ++i) {
    printInfo << "processing ";
    if (chain and (i == 0)) 
      cout << "chain of .root files";
    else
      cout << ".evt tree[" << ((chain) ? i : i + 1) << "]";
    cout << endl;
    countEvents += processTree(*trees[i], decayTopo, amplitude, ampFile, asciiOutput,
			       prodKinParticlesLeafName,  prodKinMomentaLeafName,
			       decayKinParticlesLeafName, decayKinMomentaLeafName, debug);
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
