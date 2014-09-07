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

#include "factorial.hpp"
#include "dFunction.hpp"


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
	     << "        -b #       maximum number of events in a single parallelization block  (default: 50000)" << endl
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
	vector<string> keyFileNames;
	long int       maxNmbEvents             = -1;
	long int       numParallelEvents        = 50000;
	string         pdgFileName              = "./particleDataTable.txt";
	vector<string> ampFileName; //              = "./out.root";
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
	while ((c = getopt(argc, argv, "k:n:p:o:m:at:l:r:vhb:")) != -1)
		switch (c) {
		case 'k':
			keyFileNames.push_back(optarg);
			break;
		case 'n':
			maxNmbEvents = atol(optarg);
			break;
		case 'b':
			numParallelEvents = atol(optarg);
			break;
		case 'p':
			pdgFileName = optarg;
			break;
		case 'o':
			ampFileName.push_back(optarg);
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

	// caches are initialized here to allow multi-threaded access
	cout << "Init factiorial cache ..." << std::endl;
	initFactorial<double>();
	cout << "Init dFunction cache ..." << std::endl;
	initDFunction<double>();

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
	if (keyFileNames.size() == 0) {
		printErr << "no key file specified. aborting." << endl;
		usage(progName, 1);
	}
	if (keyFileNames.size() != ampFileName.size()) {
		printErr << "numbers of key files and output files do not match. aborting." << endl;
		usage(progName, 1);
	}
	vector<waveDescription> waveDesc(keyFileNames.size());
	vector<isobarAmplitudePtr> amplitude(keyFileNames.size());
	for(unsigned int i = 0; i < keyFileNames.size(); ++i) {
		if (not waveDesc[i].parseKeyFile(keyFileNames[i]) or not waveDesc[i].constructAmplitude(amplitude[i])) {
			printErr << "problems constructing decay topology from key file '" << keyFileNames[i] << "'. "
					 << "aborting." << endl;
			exit(1);
		}
		printInfo << *(amplitude[i]);
	}

	// create output file for amplitudes
	vector<bool>      writeRootFormat(keyFileNames.size(), false);
	vector<ofstream*> ampFilePlain   (keyFileNames.size(), NULL);
	vector<TFile*>    ampFileRoot    (keyFileNames.size(), NULL);
#ifdef USE_STD_COMPLEX_TREE_LEAFS
	vector<amplitudeTreeLeaf*> ampTreeLeaf(keyFileNames.size(), NULL);
	vector<TTree*>             ampTree    (keyFileNames.size(), NULL);
#endif
	for(unsigned int i = 0; i < keyFileNames.size(); ++i) {
#ifdef USE_STD_COMPLEX_TREE_LEAFS
		const string ampFileExt = extensionFromPath(ampFileName[i]);
		if (ampFileExt == "root")
			writeRootFormat[i] = true;
		else if (ampFileExt == "amp")
			writeRootFormat[i] = false;
		else {
			printErr << "specified amplitude file '" << ampFileName[i] << "' is neither a .root "
					 << "nor a .amp file. aborting." << endl;
			usage(progName);
		}
#endif
		printInfo << "creating amplitude file '" << ampFileName[i] << "'";
#ifdef USE_STD_COMPLEX_TREE_LEAFS
		if (writeRootFormat[i]) {
			cout << endl;
			ampFileRoot[i] = TFile::Open(ampFileName[i].c_str(), "RECREATE");
			// write wave description
			const string waveName = waveDesc[i].waveNameFromTopology(*(amplitude[i]->decayTopology()),
																     newKeyFileNameConvention);
			waveDesc[i].Write((waveName + ".key").c_str());
			// create output tree
			ampTreeLeaf[i] = new amplitudeTreeLeaf();
			const string ampTreeName = waveName + ".amp";
			ampTree[i] = new TTree(ampTreeName.c_str(), ampTreeName.c_str());
			const int splitLevel = 99;
			const int bufSize    = 256000;
			ampTree[i]->Branch(ampLeafName.c_str(), &(ampTreeLeaf[i]), bufSize, splitLevel);
		} else
#endif
		{
			cout << "; " << ((asciiOutput) ? "ASCII" : "binary") << " mode" << endl;
			ampFilePlain[i] = new ofstream(ampFileName[i].c_str());
		}
		if (   (writeRootFormat[i] and not ampFileRoot[i])
			or (not writeRootFormat[i] and not ampFilePlain[i] and not *(ampFilePlain[i]))) {
			printErr << "cannot create amplitude file '" << ampFileName[i] << "'. aborting." << endl;
			exit(1);
		}
	}

	// read data from tree(s) and calculate decay amplitudes
	TStopwatch timer;
	timer.Reset();
	timer.Start();
	vector<vector<complex<double> > > ampValues(keyFileNames.size()); // [keyfile][event]
	for (unsigned int i = 0; i < inTrees.size(); ++i) {
		printInfo << "processing ";
		if ((rootFileNames.size() > 0) and (i == 0))
			cout << "chain of .root files";
		else
			cout << ".evt tree[" << ((rootFileNames.size() > 0) ? i : i + 1) << "]";
		cout << endl;
		if (not processTree(*inTrees[i], *prodKinPartNames, *decayKinPartNames,
		                    amplitude, ampValues, maxNmbEvents - ampValues[0].size(),
		                    numParallelEvents, prodKinMomentaLeafName, decayKinMomentaLeafName))
			printWarn << "problems reading tree" << endl;
	}
	printSucc << "calculated amplitudes for " << ampValues.size() << " topologies and "
			  << ampValues[0].size() << " events" << endl;
	timer.Stop();
	printInfo << "this job consumed: ";
	timer.Print();

	// write amplitudes to output file
	for (unsigned int k = 0; k < ampValues.size(); ++k) {

		for (unsigned int i = 0; i < ampValues[k].size(); ++i) {
#ifdef USE_STD_COMPLEX_TREE_LEAFS
			if (ampFileRoot[k]) {
				ampTreeLeaf[k]->setAmp(ampValues[k][i]);
				ampTree[k]->Fill();
			}
#endif
			if (ampFilePlain[k]) {
				if (asciiOutput)
					*(ampFilePlain[k]) << setprecision(numeric_limits<double>::digits10 + 1) << ampValues[k][i] << endl;
				else
					ampFilePlain[k]->write((char*)(&ampValues[k][i]), sizeof(complex<double>));
			}
		}

#ifdef USE_STD_COMPLEX_TREE_LEAFS
		if (ampFileRoot[k]) {
			printInfo << "optimizing tree" << endl;
			//ampTree->Print();
			ampTree[k]->OptimizeBaskets(treeCacheSize, 1, "d");
			ampTree[k]->Write();
			//ampTree->Print();
			ampFileRoot[k]->Close();
			delete ampFileRoot[k];
		}
#endif

		if (ampFilePlain[k]) {
			ampFilePlain[k]->close();
			delete ampFilePlain[k];
		}
		printSucc << "wrote " << ampValues[k].size() << " amplitude values to "
				  << "'" << ampFileName[k] << "'";
		if (not writeRootFormat[k])
			cout << " in " << ((asciiOutput) ? "ASCII" : "binary")
				 << " mode";
		cout << endl;

	}
  
	// clean up
	for (unsigned int i = 0; i < inTrees.size(); ++i)
		delete inTrees[i];
	inTrees.clear();

	return 0;
}
