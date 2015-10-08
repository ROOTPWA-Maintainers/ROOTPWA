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
//      out amplitudes in ROOT tree format
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

#include "amplitudeFileWriter.h"
#include "calcAmplitude.h"
#include "evtTreeHelper.h"
#include "fileUtils.hpp"
#include "particleDataTable.h"
#include "reportingUtilsEnvironment.h"
#include "waveDescription.h"


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
	     << " -k key file [-n max. # of events -p PDG file -o output file -h] "
	     << "input data file(s) (.root format)" << endl
	     << "    where:" << endl
	     << "        -k file    path to key file" << endl
	     << "        -n #       maximum number of events to read (default: all)" << endl
	     << "        -p file    path to particle data table file (default: ./particleDataTable.txt)" << endl
	     << "        -o file    path to amplitude file (.root format; default: constructed from topology)" << endl
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

#if ROOT_VERSION_CODE < ROOT_VERSION(6, 0, 0)
	// force loading predefined std::complex dictionary
	// see http://root.cern.ch/phpBB3/viewtopic.php?f=5&t=9618&p=50164
	gROOT->ProcessLine("#include <complex>");
#endif

	// parse command line options
	const string   progName                 = argv[0];
	string         keyFileName              = "";
	long int       maxNmbEvents             = -1;
	string         pdgFileName              = "./particleDataTable.txt";
	string         ampFileName              = "";
	extern char*   optarg;
	extern int     optind;
	int            c;
	while ((c = getopt(argc, argv, "k:n:p:o:h")) != -1)
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
		case 'h':
		default:
			usage(progName);
		}

	// get input file names
	if (optind >= argc) {
		printErr << "you need to specify at least one data file to process. Aborting..." << endl;
		usage(progName, 1);
	}
	vector<string> eventFileNames;
	while (optind < argc) {
		const string fileName = argv[optind++];
		const string fileExt  = extensionFromPath(fileName);
		if (fileExt != "root")
			printWarn << "input file '" << fileName << "' is not a .root file. skipping." << endl;
		eventFileNames.push_back(fileName);
	}
	if (eventFileNames.size() == 0) {
		printErr << "none of the specified input files is a .root file. Aborting..." << endl;
		usage(progName, 1);
	}

	// open files containing event data
	vector<const eventMetadata*> eventMetas;
	for (size_t i=0; i<eventFileNames.size(); ++i) {
		const string eventFileName = eventFileNames[i];
		TFile* file = TFile::Open(eventFileName.c_str());
		if (file == NULL || file->IsZombie()) {
			printErr << "event file '" << eventFileName << "' cannot be opened. Aborting..." << endl;
			exit(1);
		}

		const eventMetadata* eventMeta = eventMetadata::readEventFile(file);
		if (eventMeta == NULL) {
			printErr << "cannot read event data from event file '" << eventFileName << "'. Aborting..." << endl;
			exit(1);
		}

		eventMetas.push_back(eventMeta);
	}
	if (eventMetas.size() != 1) {
		printErr << "cannot handle more than one event file. Aborting..." << endl;
		exit(1);
	}

	// check file name of amplitude file
	if (ampFileName != "") {
		const string ampFileExt = extensionFromPath(ampFileName);
		if (ampFileExt != "root") {
			printErr << "specified amplitude file '" << ampFileName << "' is not a .root file. Aborting..." << endl;
			usage(progName);
		}
	}

	// initialize particle data table
	particleDataTable::readFile(pdgFileName);

	// parse key file and create amplitude instance
	if (keyFileName == "") {
		printErr << "no key file specified. Aborting..." << endl;
		usage(progName, 1);
	}
	vector<waveDescriptionPtr> waveDescs = waveDescription::parseKeyFile(keyFileName);
	if (waveDescs.size() != 1) {
		printErr << "problems reading wave description from key file '" << keyFileName << "'. "
		         << "Aborting..." << endl;
		exit(1);
	}
	waveDescriptionPtr waveDesc = waveDescs[0];
	isobarAmplitudePtr amplitude;
	if (not waveDesc->constructAmplitude(amplitude)) {
		printErr << "problems constructing decay topology from key file '" << keyFileName << "'. "
		         << "Aborting..." << endl;
		exit(1);
	}
	printInfo << *amplitude;

	// read data from tree(s) and calculate decay amplitudes
	TStopwatch timer;
	timer.Reset();
	timer.Start();
	const vector<complex<double> > ampValues = hli::calcAmplitude(*eventMetas[0], amplitude, maxNmbEvents);
	printSucc << "calculated amplitudes for " << ampValues.size() << " events" << endl;
	timer.Stop();
	printInfo << "this job consumed: ";
	timer.Print();

	// create output file for amplitudes
	if (ampFileName == "") {
		ampFileName = waveDescription::waveNameFromTopology(*(amplitude->decayTopology())) + ".root";
	}
	printInfo << "creating amplitude file '" << ampFileName << "'" << endl;
	TFile* ampFile = TFile::Open(ampFileName.c_str(), "RECREATE");
	if (not ampFile) {
		printErr << "cannot create amplitude file '" << ampFileName << "'. Aborting..." << endl;
		exit(1);
	}

	// write amplitudes to output file
	const string waveName = fileNameNoExtFromPath(ampFileName);

	amplitudeFileWriter ampFileWriter;
	ampFileWriter.initialize(*ampFile, eventMetas, waveDesc->keyFileContent(), waveName);
	ampFileWriter.addAmplitudes(ampValues);
	ampFileWriter.finalize();

	ampFile->Close();
	delete ampFile;

	printSucc << "wrote " << ampValues.size() << " amplitude values to "
	          << "'" << ampFileName << "'" << endl;

	return 0;
}
