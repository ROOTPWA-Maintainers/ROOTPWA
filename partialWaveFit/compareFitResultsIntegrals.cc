///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
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
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------
//
// Description:
//      fitting program for rootpwa
//      minimizes pwaLikelihood function
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <iostream>

#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>

#include "ampIntegralMatrix.h"
#include "fileUtils.hpp"
#include "fitResult.h"
#include "reportingUtilsEnvironment.h"

using namespace std;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{

	cerr << "compares normalization integrals from two fit result files" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " -l inputFileLeft -r inputFileRight [-h]" << endl
	     << "    where:" << endl
	     << "        -l         input fit result file left" << endl
	     << "        -r         output fit result file" << endl
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

	const string treeName = "pwa";
	const string branchName = "fitResult_v2";

	const string progName            = argv[0];
	string inputFileNameLeft = "";
	string inputFileNameRight = "";
	extern char* optarg;
	// extern int optind;
	int c;
	while ((c = getopt(argc, argv, "l:r:h")) != -1)
		switch (c) {
		case 'l':
			inputFileNameLeft = optarg;
			break;
		case 'r':
			inputFileNameRight = optarg;
			break;
		case 'h':
			usage(progName);
			break;
		}

	TFile* inputFileLeft = TFile::Open(inputFileNameLeft.c_str(), "READ");
	if(not inputFileLeft || inputFileLeft->IsZombie()) {
		printErr << "could not open input file '" << inputFileNameLeft << "'. Aborting..." << endl;
		return 1;
	}

	TFile* inputFileRight = TFile::Open(inputFileNameRight.c_str(), "READ");
	if(not inputFileRight || inputFileRight->IsZombie()) {
		printErr << "could not open input file '" << inputFileNameRight << "'. Aborting..." << endl;
		return 1;
	}

	TTree* inResultTreeLeft = 0;
	inputFileLeft->GetObject(treeName.c_str(), inResultTreeLeft);
	if(not inResultTreeLeft) {
		printErr << "could not find input tree with name '" << treeName << "' in input file '" << inputFileLeft << "'. Aborting..." << endl;
		return 1;
	}

	TTree* inResultTreeRight = 0;
	inputFileRight->GetObject(treeName.c_str(), inResultTreeRight);
	if(not inResultTreeRight) {
		printErr << "could not find input tree with name '" << treeName << "' in input file '" << inputFileRight << "'. Aborting..." << endl;
		return 1;
	}

	fitResult* inResultLeft = 0;
	inResultTreeLeft->SetBranchAddress(branchName.c_str(), &inResultLeft);

	fitResult* inResultRight = 0;
	inResultTreeRight->SetBranchAddress(branchName.c_str(), &inResultRight);

	if(inResultTreeLeft->GetEntries() != inResultTreeRight->GetEntries()) {
		printErr << "not the same number of entries in input TTrees. Aborting..." << endl;
		return 1;
	}

	if(inResultTreeLeft->GetEntries() > 1) {
		printWarn << "more than one fit result in input TTrees, hoping that they are ordered correctly..." << endl;
	}

	for(long i = 0; i < inResultTreeLeft->GetEntries(); ++i) {
		inResultTreeLeft->GetEntry(i);
		inResultTreeRight->GetEntry(i);
		printInfo << "comparing phase space integral matrix..." << endl << endl;
		if(inResultLeft->normIntegralMatrix().equalToPrecision(inResultRight->normIntegralMatrix(), 0., true)) {
			printSucc << "identical." << endl << endl;
		} else {
			printErr << "difference found." << endl << endl;
		}
		printInfo << "comparing accepted integral matrix..." << endl << endl;
		if (inResultLeft->acceptedNormIntegralMatrix().equalToPrecision(inResultRight->acceptedNormIntegralMatrix(), 0., true)) {
			printSucc << "identical." << endl << endl;
		} else {
			printErr << "difference found." << endl << endl;
		}
		printInfo << "comparing phaseSpaceIntegral vector..." << endl << endl;
		const vector<double>& phaseSpaceIntegralLeft = inResultLeft->phaseSpaceIntegralVector();
		const vector<double>& phaseSpaceIntegralRight = inResultRight->phaseSpaceIntegralVector();
		if(phaseSpaceIntegralLeft.size() != phaseSpaceIntegralRight.size()) {
			printErr << "size mismatch." << endl << endl;
		} else {
			bool success = true;
			for(unsigned int j = 0; j < phaseSpaceIntegralLeft.size(); ++j) {
				if(phaseSpaceIntegralLeft[j] != phaseSpaceIntegralRight[j]) {
					printErr << "difference found in element " << j << " (" << phaseSpaceIntegralLeft[j] << " != " << phaseSpaceIntegralRight[j] << ")." << endl;
					success = false;
				}
			}
			if(success) {
				printSucc << "identical." << endl;
			}
		}
	}

	return 0;
}
