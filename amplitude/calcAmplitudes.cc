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

#include <boost/tokenizer.hpp>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TStopwatch.h"

#include "fileUtils.hpp"
#include "particleDataTable.h"
#include "evtTreeHelper.h"
#include "waveDescription.h"
#include "amplitudeTreeLeaf.h"


using namespace std;
using namespace boost;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{
	cerr << "calculates amplitudes for events in input data files and writes them to file" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " -k key file [-n max. # of events -p PDG file -o output file -m amplitude leaf name "
	     << "-a -t tree name -l leaf names -r target particle name -v -h] "
	     << "input data file(s) (.evt or .root format)" << endl
	     << "    where:" << endl
	     << "        -k file    path to key file" << endl
	     << "        -n #       maximum number of events to read (default: all)" << endl
	     << "        -p file    path to particle data table file (default: ./particleDataTable.txt)" << endl
	     << "        -o file    path to amplitude file (.amp or .root format; default: ./out.root)" << endl
	     << "        -m         amplitude leaf name (default: 'amplitude')" << endl
	     << "        -a         write .amp files in ASCII format (default: binary)" << endl
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

#if AMPLITUDETREELEAF_ENABLED
	// force loading predefined std::complex dictionary
	// see http://root.cern.ch/phpBB3/viewtopic.php?f=5&t=9618&p=50164
	gROOT->ProcessLine("#include <complex>");
#endif
	
	// parse command line options
	const string progName           = argv[0];
	string       keyFileName        = "";
	long int     maxNmbEvents       = -1;
	string       pdgFileName        = "./particleDataTable.txt";
	string       ampFileName        = "./out.root";
	string       ampLeafName        = "amplitude";
	bool         asciiOutput        = false;
	string       inTreeName         = "rootPwaEvtTree";
	string       leafNames          = "prodKinParticles;prodKinMomenta;"
		                                "decayKinParticles;decayKinMomenta";
	string       targetParticleName = "p+";
	bool         debug              = false;
	extern char* optarg;
	extern int   optind;
	int          c;
	while ((c = getopt(argc, argv, "k:n:p:o:m:at:l:r:vh")) != -1)
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
		case 'm':
			ampLeafName = optarg;
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
		const string fileExt  = extensionFromPath(fileName);
		if (fileExt == "root")
			rootFileNames.push_back(fileName);
		else if (fileExt == "evt")
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
	typedef tokenizer<char_separator<char> > tokenizer;
	char_separator<char> separator(";");
	tokenizer            leafNameTokens(leafNames, separator);
	tokenizer::iterator  leafNameToken             = leafNameTokens.begin();
	const string         prodKinParticlesLeafName  = *leafNameToken;
	const string         prodKinMomentaLeafName    = *(++leafNameToken);
	const string         decayKinParticlesLeafName = *(++leafNameToken);
	const string         decayKinMomentaLeafName   = *(++leafNameToken);
	printInfo << "using the following leaf names:" << endl
	          << "        production kinematics: "
	          << "particle names = '" << prodKinParticlesLeafName << "', "
	          << "momenta = '" << prodKinMomentaLeafName << "'" << endl
	          << "        decay kinematics:      "
	          << "particle names = '" << decayKinParticlesLeafName << "', "
	          << "momenta = '" << decayKinMomentaLeafName << "'" << endl;
  
	// open root files and build chain
	TChain* inChain = 0;
	if (rootFileNames.size() > 0) {
		inChain = new TChain(inTreeName.c_str());
		for (unsigned int i = 0; i < rootFileNames.size(); ++i) {
			printInfo << "opening ROOT input file '" << rootFileNames[i] << "'" << endl;
			if (inChain->Add(rootFileNames[i].c_str()) < 1)
				printWarn << "no events in ROOT input file '" << rootFileNames[i] << "'" << endl;
		}
		inChain->GetListOfFiles()->ls();
	}

	// convert .evt files to root trees
	vector<TTree*> inTrees;
	if (inChain)
		inTrees.push_back(inChain);
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
			inTrees.push_back(tree);
		else {
			printWarn << "problems creating tree from .evt input file '" << evtFileNames[i] << "' "
			          << "skipping." << endl;
		}
	}

	// initialize particle data table
	particleDataTable::readFile(pdgFileName);
  
	// parse key file and create amplitude instance
	if (keyFileName == "") {
		printErr << "no key file specified. aborting." << endl;
		usage(progName, 1);
	}
	waveDescription    waveDesc;
	isobarAmplitudePtr amplitude;
	if (   not waveDesc.parseKeyFile(keyFileName)
	    or not waveDesc.constructAmplitude(amplitude)) {
		printErr << "problems constructing decay topology from key file '" << keyFileName << "'. "
		         << "aborting." << endl;
		exit(1);
	}
	printInfo << *amplitude;

	// create output file for amplitudes
	bool writeRootFormat = false;
#if AMPLITUDETREELEAF_ENABLED
	const string ampFileExt = extensionFromPath(ampFileName);
	if (ampFileExt == "root")
		writeRootFormat = true;
	else if (ampFileExt == "amp")
		writeRootFormat = false;
	else {
		printErr << "specified amplitude file '" << ampFileName << "' is neither a .root "
		         << "nor a .amp file. aborting." << endl;
		usage(progName);
	}
#endif
	printInfo << "creating amplitude file '" << ampFileName << "'";
	ofstream*          ampFilePlain = 0;
	TFile*             ampFileRoot  = 0;
#if AMPLITUDETREELEAF_ENABLED
	amplitudeTreeLeaf* ampTreeLeaf  = 0;
	TTree*             ampTree      = 0;
	if (writeRootFormat) {
		cout << endl;
		ampFileRoot = TFile::Open(ampFileName.c_str(), "RECREATE");
		// write wave description
		const string waveName = waveDesc.waveNameFromTopologyOld(*(amplitude->decayTopology()));
		waveDesc.Write((waveName + ".key").c_str());
		// create output tree
		ampTreeLeaf = new amplitudeTreeLeaf();
		ampTreeLeaf->setNmbIncohSubAmps(1);
		const string ampTreeName = waveName + ".amp";
		ampTree = new TTree(ampTreeName.c_str(), ampTreeName.c_str());
		ampTree->Branch(ampLeafName.c_str(), &ampTreeLeaf);
	} else
#endif
	{
		cout << "; " << ((asciiOutput) ? "ASCII" : "binary") << " mode" << endl;
		ampFilePlain = new ofstream(ampFileName.c_str());
	}
	if ((writeRootFormat and not ampFileRoot)
	    or (not writeRootFormat and not ampFilePlain and not *ampFilePlain)) {
		printErr << "cannot create amplitude file '" << ampFileName << "'. aborting." << endl;
		exit(1);
	}
  
	// read data from tree(s), calculate amplitudes
	TStopwatch timer;
	timer.Reset();
	timer.Start();
	vector<complex<double> > ampValues;
	for (unsigned int i = 0; i < inTrees.size(); ++i) {
		printInfo << "processing ";
		if (inChain and (i == 0)) 
			cout << "chain of .root files";
		else
			cout << ".evt tree[" << ((inChain) ? i : i + 1) << "]";
		cout << endl;
		if (not processTree(*inTrees[i], amplitude, ampValues, maxNmbEvents - ampValues.size(),
		                    prodKinParticlesLeafName,  prodKinMomentaLeafName,
		                    decayKinParticlesLeafName, decayKinMomentaLeafName))
			printWarn << "problems reading tree" << endl;
	}
	printInfo << "successfully calculated amplitudes for " << ampValues.size() << " events" << endl;
	timer.Stop();
	printInfo << "this job consumed: ";
	timer.Print();
  
	// write amplitudes to output file
	for (unsigned int i = 0; i < ampValues.size(); ++i) {
#if AMPLITUDETREELEAF_ENABLED
		if (ampFileRoot) {
			ampTreeLeaf->setIncohSubAmp(ampValues[i]);
			ampTree->Fill();
		}
#endif
		if (ampFilePlain) {
			if (asciiOutput)
				*ampFilePlain << setprecision(numeric_limits<double>::digits10 + 1) << ampValues[i] << endl;
			else
				ampFilePlain->write((char*)(&ampValues[i]), sizeof(complex<double>));
		}
	}
#if AMPLITUDETREELEAF_ENABLED
	if (ampFileRoot) {
		ampTree->Write();
		ampFileRoot->Close();
		delete ampFileRoot;
	}
#endif
	if (ampFilePlain) {
		ampFilePlain->close();
		delete ampFilePlain;
	}
	printInfo << "successfully wrote " << ampValues.size() << " amplitude values to "
	          << "'" << ampFileName << "'";
	if (not writeRootFormat)
		cout << " in " << ((asciiOutput) ? "ASCII" : "binary")
		     << " mode";
	cout << endl;
  
	// clean up
	for (unsigned int i = 0; i < inTrees.size(); ++i)
		delete inTrees[i];
	inTrees.clear();

	return 0;
}
