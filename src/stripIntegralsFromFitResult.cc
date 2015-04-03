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

#include "fitResult.h"
#include "reportingUtilsEnvironment.h"

using namespace std;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{

	cerr << "removes normalization integrals from fit result file" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " -i inputFile -o outputFile [-c -h]" << endl
	     << "    where:" << endl
	     << "        -i         input fit result file" << endl
	     << "        -o         output fit result file" << endl
	     << "        -c         also strip covariance matrix (default: false)" << endl
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

	// force loading predefined std::complex dictionary
	// see http://root.cern.ch/phpBB3/viewtopic.php?f=5&t=9618&p=50164
	gROOT->ProcessLine("#include <complex>");

	const string treeName = "pwa";
	const string branchName = "fitResult_v2";

	const string progName            = argv[0];
	string inputFileName = "";
	string outputFileName = "";
	bool stripCovarianceMatrix = false;
	extern char* optarg;
	// extern int optind;
	int c;
	while ((c = getopt(argc, argv, "i:o:ch")) != -1)
		switch (c) {
		case 'i':
			inputFileName = optarg;
			break;
		case 'o':
			outputFileName = optarg;
			break;
		case 'c':
			stripCovarianceMatrix = true;
			break;
		case 'h':
			usage(progName);
			break;
		}

	TFile* inputFile = TFile::Open(inputFileName.c_str(), "READ");
	if(not inputFile || inputFile->IsZombie()) {
		printErr << "could not open input file '" << inputFileName << "'. aborting." << endl;
		return 1;
	}

	TFile* outputFile = TFile::Open(outputFileName.c_str(), "NEW");
	if(not outputFile || outputFile->IsZombie()) {
		printErr << "could not open output file '" << outputFileName << "'. aborting." << endl;
		return 1;
	}

	TTree* inResultTree = 0;
	inputFile->GetObject(treeName.c_str(), inResultTree);
	if(not inResultTree) {
		printErr << "could not find input tree with name '" << treeName << "' in input file '" << inputFileName << "'. aborting." << endl;
		return 1;
	}

	fitResult* inResult = 0;
	inResultTree->SetBranchAddress(branchName.c_str(), &inResult);

	outputFile->cd();
	TTree* outResultTree = new TTree(treeName.c_str(), treeName.c_str());
	fitResult* outResult = 0;
	outResultTree->Branch(branchName.c_str(), &outResult);

	for(long i = 0; i < inResultTree->GetEntries(); ++i) {
		inResultTree->GetEntry(i);
		const unsigned int             nmbEvents              = inResult->nmbEvents();
		const unsigned int             normNmbEvents          = inResult->normNmbEvents();
		const double                   massBinCenter          = inResult->massBinCenter();
		const double                   logLikelihood          = inResult->logLikelihood();
		const int                      rank                   = inResult->rank();

		const vector<TComplex>&        prodAmpsTComplex       = inResult->prodAmps();
		const unsigned int             nmbProdAmps            = prodAmpsTComplex.size();
		vector<complex<double> >       prodAmps(nmbProdAmps);
		for(unsigned int i = 0; i < nmbProdAmps; ++i) {
			prodAmps[i] = complex<double>(prodAmpsTComplex[i].Re(), prodAmpsTComplex[i].Im());
		}

		const vector<string>&          prodAmpNames           = inResult->prodAmpNames();

		const vector<pair<int, int> >& fitParCovMatrixIndices = inResult->fitParCovIndices();

		const complexMatrix            normIntegral(0, 0);
		const complexMatrix            accIntegral(0, 0);
		const vector<double>           phaseSpaceIntegral;

		const bool                     converged              = inResult->converged();
		const bool                     hasHessian             = inResult->hasHessian();

		outResult->reset();
		outResult->fill(nmbEvents,
		                normNmbEvents,
		                massBinCenter,
		                logLikelihood,
		                rank,
		                prodAmps,
		                prodAmpNames,
		                stripCovarianceMatrix ? TMatrixT<double>(0, 0) : inResult->fitParCovMatrix(),
		                fitParCovMatrixIndices,
		                normIntegral,
		                accIntegral,
		                phaseSpaceIntegral,
		                converged,
		                hasHessian);

		outResultTree->Fill();

	}

	outputFile->cd();
	outResultTree->Write();
	outputFile->Close();

	if(outResult) {
		delete outResult;
	}

	return 0;
}
