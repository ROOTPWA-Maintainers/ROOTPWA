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
//      calculates integral matrix for set of amplitudes
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <stdlib.h>
#include <unistd.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TIterator.h"
#include "TKey.h"

#include "reportingUtils.hpp"
#include "reportingUtilsEnvironment.h"
#include "fileUtils.hpp"
#include "ampIntegralMatrix.h"
#include "amplitudeTreeLeaf.h"


using namespace std;


void
usage(const string& progName,
      const int     errCode = 0)
{
	cerr << "computes integral matrix for given amplitude files" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " [-i binfilename -j rootfilename -n max. # of events -v -h] "
	     << "amplitude files" << endl
	     << "    where:" << endl
	     << "        -i name    input bin file name"  << endl
	     << "        -j name    input root file name"  << endl
	     << "        -n #       maximum number of events to process (default: all)"               << endl
	     << "        -v         verbose; print debug output (default: false)"                     << endl
	     << "        -h         print help" << endl
	     << endl;
	exit(errCode);
}


int
main(int    argc,
     char** argv)
{
	rpwa::printCompilerInfo();
	rpwa::printLibraryInfo ();
	rpwa::printGitHash     ();
	cout << endl;

	// force loading predefined std::complex dictionary
	// see http://root.cern.ch/phpBB3/viewtopic.php?f=5&t=9618&p=50164
	gROOT->ProcessLine("#include <complex>");

	// parse command line options
	const string progName          = argv[0];
	string       inputBinFilePath  = "";
	string       inputRootFilePath = "";
	unsigned int maxNmbEvents      = 0;
	bool         debug             = false;
	extern char* optarg;
	int          c;
	while ((c = getopt(argc, argv, "i:j:n:vh")) != -1) {
		switch (c) {
		case 'i':
			inputBinFilePath = optarg;
			break;
		case 'j':
			inputRootFilePath = optarg;
			break;
		case 'n':
			maxNmbEvents = atoi(optarg);
			break;
		case 'v':
			debug = true;
			break;
		case 'h':
		default:
			usage(progName);
		}
	}

	printInfo << "opening bin file '" << inputBinFilePath << "'." << endl;
	ifstream* binFile = new ifstream();
	binFile->open(inputBinFilePath.c_str());
	vector<complex<double> > binAmps(0);
	unsigned int evt = 0;

	while(not binFile->eof()) {
		complex<double> amp = 0;
		binFile->read((char*)&(amp), sizeof(complex<double>));
		binAmps.push_back(amp);
		evt++;
	}
	binFile->close();
	printInfo << "done with bin file..." << endl;
	printInfo << "opening root file '" << inputRootFilePath << "'." << endl;
	TFile* rootFile = TFile::Open(inputRootFilePath.c_str(), "READ");
	if (not rootFile or rootFile->IsZombie()) {
		printErr << "could open amplitude file '" << inputRootFilePath << "'. skipping file." << endl;
		return 1;
	}
	// find all amplitude trees with name matching *.amp and corresponding .key wave description
	TIterator* keys        = rootFile->GetListOfKeys()->MakeIterator();
	bool       foundAmpKey = false;
	TTree*     tree     = 0;
	while (TKey* k = static_cast<TKey*>(keys->Next())) {
		if (!k) {
			printWarn << "NULL pointer to TKey in file '" << inputRootFilePath << "'. skipping TKey." << endl;
			continue;
		}
		const string     keyName  = k->GetName();
		if (rpwa::extensionFromPath(keyName) == "amp") {
			foundAmpKey = true;
			rootFile->GetObject(keyName.c_str(), tree);
			printInfo << "found TTree with keyname '" << keyName << "'." << endl;
		}
	}
	if (not foundAmpKey)
		printWarn << "no TKey in file '" << inputRootFilePath << "' matches '*.amp'. skipping file." << endl;
	rpwa::amplitudeTreeLeaf* leaf = 0;
	tree->SetBranchAddress("amplitude", &leaf);
	const int nmbEvt = tree->GetEntries();
	vector<complex<double> > rootAmps(nmbEvt, 0);
	for (int evt = 0; evt < nmbEvt; evt++) {
		tree->GetEvent(evt);
		const unsigned int nmbSubAmps = leaf->nmbIncohSubAmps();
		if (nmbSubAmps < 1) {
			printErr << "amplitude object for wave "
			         << "does not contain any amplitude values "
			         << "at event " << evt << " of total " << nmbEvt << ". aborting." << endl;
			return 1;
		}
		// get all incoherent subamps
		rootAmps[evt] = leaf->incohSubAmp(0);
	}
	if (debug) printInfo << "debug: " << rpwa::yesNo(debug) << endl;
	printInfo << "comparing amplitudes: binAmps / rootAmps" << endl;
	const int maxEvt = maxNmbEvents;
	unsigned int nmbEvents = (maxEvt) ? min(nmbEvt, maxEvt) : nmbEvt;
	cout << "looping over " << nmbEvents << " events." << endl;
	bool success = true;
	for (unsigned int evt = 0; evt < nmbEvents; evt++) {
		if (debug) {
			cout << "event #" << evt << endl
				 << rpwa::maxPrecisionAlign(binAmps[evt].real())  << "  "
				 << rpwa::maxPrecisionAlign(binAmps[evt].imag())  << endl
				 << rpwa::maxPrecisionAlign(rootAmps[evt].real()) << "  "
				 << rpwa::maxPrecisionAlign(rootAmps[evt].imag()) << endl;
		}
		if(binAmps[evt] != rootAmps[evt]) {
			cout << "event #" << evt << endl
				 << rpwa::maxPrecisionAlign(binAmps[evt].real())  << "  "
				 << rpwa::maxPrecisionAlign(binAmps[evt].imag())  << endl
				 << rpwa::maxPrecisionAlign(rootAmps[evt].real()) << "  "
				 << rpwa::maxPrecisionAlign(rootAmps[evt].imag()) << endl;
			cout << "amplitudes in event #" << evt << " not equal!" << endl;
			success = false;
		}
	}
	if (success) {
		cout << "Amplitude files are equal! Success!" << endl;
		return 0;
	} else {
		cout << "Amplitude files are different! :(" << endl;
		return 1;
	}
}
