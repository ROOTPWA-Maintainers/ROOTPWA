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
//      compares two amplitude files and calculates differences
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <vector>
#include <complex>

#include "TFile.h"
#include "TTree.h"
#include "TObjString.h"
#include "TSystem.h"

#include "reportingUtils.hpp"


using namespace std;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{
	cerr << "calculates difference of two amplitude files an writes result into tree" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " [-o out file -v -h] .amp file 1 .amp file 2" << endl
	     << "    where:" << endl
	     << "        -o file    path to output ROOT file (default: none)" << endl
	     << "        -v         verbose; print debug output (default: false)" << endl
	     << "        -h         print help" << endl
	     << endl;
	exit(errCode);
}


// tries to extract mass bin boundaries from amplitude file path
// assumes standard directory layout:
// .../[massBinMin].[massBinMax]/{PSP,}AMPS/*.amp
bool
massBinFromPath(const string& path,
                unsigned int& massBinMin,
                unsigned int& massBinMax)
{
	// strip .amp file
	size_t pos;;
	if ((pos = path.find_last_of("/")) == string::npos)
		return false;
	string massBinDir = path.substr(0, pos);
	// strip {PSP,}AMP
	if ((pos = massBinDir.find_last_of("/")) == string::npos)
		return false;
	massBinDir = massBinDir.substr(0, pos);
	// extract [massBinMin].[massBinMax]
	if ((pos = massBinDir.find_last_of("/")) == string::npos)
		return false;
	massBinDir = massBinDir.substr(pos + 1);
	// extract mass bins
	if ((pos = massBinDir.find_last_of(".")) == string::npos)
		return false;
	massBinMin = atoi(massBinDir.substr(0, pos ).c_str());
	massBinMax = atoi(massBinDir.substr(pos + 1).c_str());
	if ((massBinMin == 0) or (massBinMax == 0))
		return false;
	return true;
}


int
main(int    argc,
     char** argv)
{
	printCompilerInfo();
	printLibraryInfo ();
	printGitHash     ();
	cout << endl;

	// parse command line options
	const string progName        = argv[0];
	string       outFileName     = "";
	bool         debug           = false;
	string       ampFileNames[2] = {"", ""};
	extern char* optarg;
	extern int   optind;
	int          c;
	while ((c = getopt(argc, argv, "o:vh")) != -1)
		switch (c) {
		case 'o':
			outFileName = optarg;
			break;
		case 'v':
			debug = true;
			break;
		case 'h':
		default:
			usage(progName);
		}

	// get amplitude file names
	if (argc - optind != 2) {
		printErr << "you have to specify two amplitude files. aborting." << endl;;
		usage(progName, 1);
	}
	ampFileNames[0] = argv[optind++];
	ampFileNames[1] = argv[optind++];

	// open amplitude files
	printInfo << "comparing amplitude files '" << ampFileNames[0] << "' and "
	          << "'" << ampFileNames[1] << endl;
	ifstream ampFiles[2];
	for (unsigned int i = 0; i < 2; ++i) {
		ampFiles[i].open(ampFileNames[i].c_str());
		if (not ampFiles[i]) {
			printErr << "cannot open amplitude file '" << ampFileNames[i] << "'. aborting." << endl;
			exit(1);
		}
	}

	// read amplitudes into memory
	vector<complex<double> > amps[2];
	for (unsigned int i = 0; i < 2; ++i) {
		complex<double> amp;
		while (ampFiles[i].read((char*) &amp, sizeof(complex<double>)))
			amps[i].push_back(amp);
	}
	if (amps[0].size() != amps[1].size())
		printWarn << "the two amplitude files have different number of amplitudes "
		          << "(" << amps[0].size() << " vs. " << amps[1].size() << ")." << endl;
	const unsigned int nmbAmps = min(amps[0].size(), amps[1].size());

	// open output file
	TFile* outFile = 0;
	if (outFileName != "") {
		printInfo << "writing difference tree to '" << outFileName << "'" << endl;
		outFile = TFile::Open(outFileName.c_str(), "RECREATE");
	}

	// create tree
	TTree* tree = new TTree("ampDiffTree", "ampDiffTree");
	if (not tree) {
		printErr << "problems creating tree 'ampDiffTree' " << "in file '" << outFileName << "'" << endl;
		return false;
	}

	// create leaf variables
	const string cmd     = "basename '" + ampFileNames[0] + "'";
	TObjString*  ampName = new TObjString(gSystem->GetFromPipe(cmd.c_str()));
	UInt_t       eventNmb;
	UInt_t       massBinMin = 0, massBinMax = 0;  // [MeV/c^2]
	Double_t     valReal[2], valImag[2];
	Double_t     absDiffReal, absDiffImag;
	Double_t     relDiffReal, relDiffImag;

	// get mass bin boundaries from file paths
	if (not     massBinFromPath(ampFileNames[0], massBinMin, massBinMax)
	    and not massBinFromPath(ampFileNames[1], massBinMin, massBinMax))
		printWarn << "cannot determine mass bin boundaries from file paths" << endl;
	else
		printInfo << "extracted mass bin boundaries [" << massBinMin << ", " << massBinMax << "] from file paths" << endl;

	// connect leaf variables to tree branches
	const int split   = 0;
	const int bufSize = 256000;
	tree->Branch("ampName",     "TObjString", &ampName,  bufSize, split);
	tree->Branch("eventNmb",    &eventNmb,    "eventNmb/i");
	tree->Branch("massBinMin",  &massBinMin,  "massBinMin/i");
	tree->Branch("massBinMax",  &massBinMax,  "massBinMax/i");
	tree->Branch("valReal",     valReal,      "valReal[2]/D");
	tree->Branch("valImag",     valImag,      "valImag[2]/D");
	tree->Branch("absDiffReal", &absDiffReal, "absDiffReal/D");
	tree->Branch("absDiffImag", &absDiffImag, "absDiffImag/D");
	tree->Branch("relDiffReal", &relDiffReal, "relDiffReal/D");
	tree->Branch("relDiffImag", &relDiffImag, "relDiffImag/D");

	// compare amplitudes
	double   maxAbsDiff      = 0;
	long int maxAbsDiffIndex = -1;
	double   maxRelDiff      = 0;
	long int maxRelDiffIndex = -1;
	for (unsigned int i = 0; i < nmbAmps; ++i) {
		const complex<double> absDiff = amps[0][i] - amps[1][i];
		const complex<double> relDiff = complex<double>(absDiff.real() / amps[0][i].real(),
		                                                absDiff.imag() / amps[0][i].imag());
		// fill tree
		if (outFile) {
			eventNmb    = i;
			valReal[0]  = amps[0][i].real();
			valImag[0]  = amps[0][i].imag();
			valReal[1]  = amps[1][i].real();
			valImag[1]  = amps[1][i].imag();
			absDiffReal = absDiff.real();
			absDiffImag = absDiff.imag();
			relDiffReal = relDiff.real();
			relDiffImag = relDiff.imag();
			tree->Fill();
		}
		// print amplitudes
		if (debug) {
			const unsigned int nmbDigits = numeric_limits<double>::digits10 + 1;
			ostringstream s;
			s.precision(nmbDigits);
			s.setf(ios_base::scientific, ios_base::floatfield);
			s << "    " << setw(49) << amps[0][i] << " - " << setw(49) << amps[1][i]
			  << " = " << setw(49) << absDiff	<< ", relative "
			  << "(" << setw(23) << relDiff.real() << ", " << setw(23) << relDiff.imag() << " )";
			cout << s.str() << endl;
		}
		// compute maximum deviations
		const double absMax = max(fabs(absDiffReal), fabs(absDiffImag));
		const double relMax = max(fabs(relDiffReal), fabs(relDiffImag));
		if (absMax > maxAbsDiff) {
			maxAbsDiff      = absMax;
			maxAbsDiffIndex = i;
		}
		if (relMax > maxRelDiff) {
			maxRelDiff      = relMax;
			maxRelDiffIndex = i;
		}
	}
	printInfo << "maximum observed deviations: absolute = " << maxPrecision(maxAbsDiff) << " "
	          << "(event " << maxAbsDiffIndex << "), relative = " << maxPrecision(maxRelDiff) << " "
	          << "(event " << maxRelDiffIndex << ")" << endl;

	if (outFile) {
		outFile->Write();
		outFile->Close();
		delete outFile;
		printInfo << "wrote difference tree to '" << outFileName << "'" << endl;
	}
	return 0;
}
