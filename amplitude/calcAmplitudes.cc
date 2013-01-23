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
//      reads in data files in .evt or ROOT tree format, calculates
//      amplitudes for each event based on given key file, and writes
//      out amplitudes in ROOT tree or PWA2000 binary or ascii format
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
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{

	cerr << "calculates decay amplitudes for given wave for events in input data files and" << endl
	     << "writes amplitudes to file" << endl
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
	     << "        -l names   semicolon separated object/leaf names in input data" << endl
	     << "                   (default: 'prodKinParticles;prodKinMomenta;decayKinParticles;decayKinMomenta')" << endl
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
	printLibraryInfo ();
	printGitHash     ();
	cout << endl;

#ifdef USE_STD_COMPLEX_TREE_LEAFS
	// force loading predefined std::complex dictionary
	// see http://root.cern.ch/phpBB3/viewtopic.php?f=5&t=9618&p=50164
	gROOT->ProcessLine("#include <complex>");
#endif
	
	// parse command line options
	const string   progName                 = argv[0];
	string         keyFileName              = "";
	long int       maxNmbEvents             = -1;
	string         pdgFileName              = "./particleDataTable.txt";
	string         ampFileName              = "./out.root";
	string         ampLeafName              = "amplitude";
	bool           asciiOutput              = false;
	string         inTreeName               = "rootPwaEvtTree";
	string         leafNames                = "prodKinParticles;prodKinMomenta;decayKinParticles;decayKinMomenta";
	bool           newKeyFileNameConvention = false;
	bool           debug                    = false;
	const long int treeCacheSize            = 1000000;  // 1 MByte ROOT tree read cache size
	extern char*   optarg;
	extern int     optind;
	int            c;
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

	// get object and leaf names for event data
	string prodKinPartNamesObjName,  prodKinMomentaLeafName;
	string decayKinPartNamesObjName, decayKinMomentaLeafName;
	parseLeafAndObjNames(leafNames, prodKinPartNamesObjName, prodKinMomentaLeafName,
	                     decayKinPartNamesObjName, decayKinMomentaLeafName);

	// open .root and .evt files for event data
	vector<TTree*> inTrees;
	TClonesArray*  prodKinPartNames  = 0;
	TClonesArray*  decayKinPartNames = 0;
	if (not openRootEvtFiles(inTrees, prodKinPartNames, decayKinPartNames,
	                         rootFileNames, evtFileNames,
	                         inTreeName, prodKinPartNamesObjName, prodKinMomentaLeafName,
	                         decayKinPartNamesObjName, decayKinMomentaLeafName, debug)) {
		printErr << "problems opening input file(s). aborting." << endl;
		exit(1);
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
#ifdef USE_STD_COMPLEX_TREE_LEAFS
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
#ifdef USE_STD_COMPLEX_TREE_LEAFS
	amplitudeTreeLeaf* ampTreeLeaf  = 0;
	TTree*             ampTree      = 0;
	if (writeRootFormat) {
		cout << endl;
		ampFileRoot = TFile::Open(ampFileName.c_str(), "RECREATE");
		// write wave description
		const string waveName = waveDesc.waveNameFromTopology(*(amplitude->decayTopology()),
		                                                      newKeyFileNameConvention);
		waveDesc.Write((waveName + ".key").c_str());
		// create output tree
		ampTreeLeaf = new amplitudeTreeLeaf();
		const string ampTreeName = waveName + ".amp";
		ampTree = new TTree(ampTreeName.c_str(), ampTreeName.c_str());
		const int splitLevel = 99;
		const int bufSize    = 256000;
		ampTree->Branch(ampLeafName.c_str(), &ampTreeLeaf, bufSize, splitLevel);
	} else
#endif
		{
			cout << "; " << ((asciiOutput) ? "ASCII" : "binary") << " mode" << endl;
			ampFilePlain = new ofstream(ampFileName.c_str());
		}
	if (   (writeRootFormat and not ampFileRoot)
	    or (not writeRootFormat and not ampFilePlain and not *ampFilePlain)) {
		printErr << "cannot create amplitude file '" << ampFileName << "'. aborting." << endl;
		exit(1);
	}
  
	// read data from tree(s) and calculate decay amplitudes
	TStopwatch timer;
	timer.Reset();
	timer.Start();
	vector<complex<double> > ampValues;
	for (unsigned int i = 0; i < inTrees.size(); ++i) {
		printInfo << "processing ";
		if ((rootFileNames.size() > 0) and (i == 0)) 
			cout << "chain of .root files";
		else
			cout << ".evt tree[" << ((rootFileNames.size() > 0) ? i : i + 1) << "]";
		cout << endl;
		if (not processTree(*inTrees[i], *prodKinPartNames, *decayKinPartNames,
		                    amplitude, ampValues, maxNmbEvents - ampValues.size(),
		                    prodKinMomentaLeafName, decayKinMomentaLeafName))
			printWarn << "problems reading tree" << endl;
	}
	printSucc << "calculated amplitudes for " << ampValues.size() << " events" << endl;
	timer.Stop();
	printInfo << "this job consumed: ";
	timer.Print();
  
	// write amplitudes to output file
	for (unsigned int i = 0; i < ampValues.size(); ++i) {
#ifdef USE_STD_COMPLEX_TREE_LEAFS
		if (ampFileRoot) {
			ampTreeLeaf->setAmp(ampValues[i]);
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
#ifdef USE_STD_COMPLEX_TREE_LEAFS
	if (ampFileRoot) {
		printInfo << "optimizing tree" << endl;
		//ampTree->Print();
		ampTree->OptimizeBaskets(treeCacheSize, 1, "d");
		ampTree->Write();
		//ampTree->Print();
		ampFileRoot->Close();
		delete ampFileRoot;
	}
#endif
	if (ampFilePlain) {
		ampFilePlain->close();
		delete ampFilePlain;
	}
	printSucc << "wrote " << ampValues.size() << " amplitude values to "
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
