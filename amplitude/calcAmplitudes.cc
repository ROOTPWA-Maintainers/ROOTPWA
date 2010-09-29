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

#include "svnVersion.h"
#include "utilities.h"
#include "particleDataTable.h"
#include "evtTreeHelper.h"
#include "keyFileParser.h"
#include "isobarHelicityAmplitude.h"


using namespace std;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{
	cerr << "calculates amplitudes for events in input data files and writes them to file" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " -k key file [-n max. # of events -p PDG file -o output file -a -t tree name -l leaf names -r target particle name -v -h] input data file(s) (.evt or .root format)" << endl
	     << "    where:" << endl
	     << "        -k file    path to key file" << endl
	     << "        -n #       maximum number of events to read (default: all)" << endl
	     << "        -p file    path to particle data table file (default: ./particleDataTable.txt)" << endl
	     << "        -o file    path to amplitude file (default: ./out.amp)" << endl
	     << "        -a         write amplitudes in ASCII format (default: binary)" << endl
	     << "        -t name    name of tree in ROOT data files (default: rootPwaEvtTree)" << endl
	     << "        -l names   semicolon separated tree leaf names (default: 'prodKinParticles;prodKinMomenta;decayKinParticles;decayKinMomenta')" << endl
	     << "        -r         target particle name for .evt files (default: 'p+')" << endl
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
	const string progName           = argv[0];
	string       keyFileName        = "";
	long int     maxNmbEvents       = -1;
	string       pdgFileName        = "./particleDataTable.txt";
	string       ampFileName        = "./out.amp";
	bool         asciiOutput        = false;
	string       inTreeName         = "rootPwaEvtTree";
	string       leafNames          = "prodKinParticles;prodKinMomenta;"
		                                "decayKinParticles;decayKinMomenta";
	string       targetParticleName = "p+";
	bool         debug              = false;
	extern char* optarg;
	extern int   optind;
	int          c;
	while ((c = getopt(argc, argv, "k:n:p:o:at:l:r:vh")) != -1)
		switch (c) {
		case 'k':
			keyFileName = optarg;
			break;
		case 'n':
			maxNmbEvents = atol(optarg);
			break;
		case 'p':
			pdgFileName = optarg;
			break;
		case 'o':
			ampFileName = optarg;
			break;
		case 'a':
			asciiOutput = true;
			break;
		case 't':
			inTreeName = optarg;
			break;
		case 'l':
			leafNames = optarg;
			break;
		case 'r':
			targetParticleName = optarg;
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
		printErr << "you need to specify at least one data file to process. aborting." << endl;;
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
		printErr << "none of the specified input files is a .root or .evt file. aborting.";
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
		if (not evtFile or not evtFile.good()) {
			printWarn << "cannot open .evt input file '" << evtFileNames[i] << "'. skipping." << endl;
			continue;
		}
		printInfo << "converting .evt input file '" << evtFileNames[i] << "' "
		          << "into memory resident tree. this might reduce performance. "
		          << "ROOT input format is recommended." << endl;
		// create tree
		TTree* tree = new TTree(inTreeName.c_str(), inTreeName.c_str());
		if (not tree) {
			printErr << "problems creating tree '" << inTreeName << "'. skipping." << endl;
			continue;
		}
		if (fillTreeFromEvt(evtFile, *tree, -1,
		                    prodKinParticlesLeafName,  prodKinMomentaLeafName,
		                    decayKinParticlesLeafName, decayKinMomentaLeafName,
		                    targetParticleName, debug))
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
		printErr << "no key file specified. aborting." << endl;
		usage(progName, 1);
	}
	keyFileParser&         parser = keyFileParser::instance();
	isobarDecayTopologyPtr decayTopo;
	if (not parser.parse(keyFileName) or not parser.constructDecayTopology(decayTopo)) {
		printErr << "problems constructing decay topology from key file '" << keyFileName << "'. "
		         << "aborting." << endl;
		exit(1);
	}
	decayTopo->checkTopology();
	decayTopo->checkConsistency();
	isobarHelicityAmplitude amplitude(decayTopo);
	parser.setAmplitudeOptions(amplitude);
	printInfo << amplitude;
  
	// create output file for amplitudes
	printInfo << "creating amplitude file '" << ampFileName << "'; "
	          << ((asciiOutput) ? "ASCII" : "binary") << " mode" << endl;
	ofstream ampFile(ampFileName.c_str());
	if (not ampFile) {
		printErr << "cannot create amplitude file '" << ampFileName << "'. aborting." << endl;
		exit(1);
	}
  
	// read data from tree(s), calculate amplitudes
	TStopwatch timer;
	timer.Reset();
	timer.Start();
	vector<complex<double> > ampValues;
	for (unsigned int i = 0; i < trees.size(); ++i) {
		printInfo << "processing ";
		if (chain and (i == 0)) 
			cout << "chain of .root files";
		else
			cout << ".evt tree[" << ((chain) ? i : i + 1) << "]";
		cout << endl;
		if (not processTree(*trees[i], *decayTopo, amplitude, ampValues, maxNmbEvents - ampValues.size(),
		                    prodKinParticlesLeafName,  prodKinMomentaLeafName,
		                    decayKinParticlesLeafName, decayKinMomentaLeafName))
			printWarn << "problems reading tree" << endl;
	}
	printInfo << "successfully calculated amplitudes for " << ampValues.size() << " events" << endl;
	timer.Stop();
	printInfo << "this job consumed: ";
	timer.Print();
  
	// write amplitudes to output file
	for (unsigned int i = 0; i < ampValues.size(); ++i)
		if (asciiOutput)
			ampFile  << setprecision(numeric_limits<double>::digits10 + 1) << ampValues[i] << endl;
		else
			ampFile.write((char*)(&ampValues[i]), sizeof(complex<double>));
	printInfo << "successfully wrote " << ampValues.size() << " amplitude values to "
	          << "'" << ampFileName << "' in " << ((asciiOutput) ? "ASCII" : "binary")
	          << " mode" << endl;
  
	// clean up
	for (unsigned int i = 0; i < trees.size(); ++i)
		delete trees[i];
	trees.clear();

	return 0;
}
