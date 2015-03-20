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
#include <iomanip>
#include <vector>
#include <string>
#include <complex>
#include <cassert>
#include <time.h>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TComplex.h"
#include "TRandom3.h"
#include "TMatrixT.h"
#include "TVectorT.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TStopwatch.h"
#include "Minuit2/Minuit2Minimizer.h"

#include "reportingUtilsEnvironment.h"
#include "conversionUtils.hpp"
#include "pwaLikelihood.h"
#include "fitResult.h"

#include "amplitudeTreeLeaf.h"


using namespace std;
using namespace ROOT::Math;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{

	cerr << "add analytical covariance matrix to a fit result file." << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " -l # -u # -w wavelist -i fitResultFile [-d amplitude directory -R -N -n normfile"
	     << " [-a normfile] -r rank [-t # -T treename -B branchname -f -q -h]" << endl
	     << "    where:" << endl
	     << "        -l #       lower edge of mass bin [MeV/c^2]" << endl
	     << "        -u #       upper edge of mass bin [MeV/c^2]" << endl
	     << "        -w file    path to wavelist file" << endl
	     << "        -i file    path to fit result input file with parameters" << endl
	     << "        -o file    path to fit result output file" << endl
	     << "        -C         use half-Cauchy priors" << endl
	     << "        -d dir     path to directory with decay amplitude files (default: '.')" << endl
#ifdef USE_STD_COMPLEX_TREE_LEAFS
	     << "        -R         use .root amplitude files (default: false)" << endl
#else
	     << "        -R         use .root amplitude files [not supported; ROOT version too low]" << endl
#endif
	     << "        -N         use normalization of decay amplitudes (default: false)" << endl
	     << "        -n file    path to normalization integral file (default: 'norm.int')" << endl
	     << "        -a file    path to acceptance integral file (default: 'norm.int')" << endl
	     << "        -A #       number of input events to normalize acceptance to" << endl
	     << "        -r #       rank of spin density matrix (default: 1)" << endl
	     << "        -t #       minimizer tolerance (default: 0.001)" << endl
	     << "        -T name    tree name (default: 'pwa')" << endl
	     << "        -B name    branch name (default: 'fitResult_v2')" << endl
	     << "        -q         run quietly (default: false)" << endl
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

	// ---------------------------------------------------------------------------
	// parse command line options
	const string progName            = argv[0];
	double       massBinMin          = 0;                      // [MeV/c^2]
	double       massBinMax          = 0;                      // [MeV/c^2]
	string       waveListFileName    = "";                     // wavelist filename
	string       fitResultFileName   = "";                     // fit result filename
	string       outputFileName      = "fitResult.root";       // output filename
	string       treeName            = "pwa";                  // fit result tree name
	string       branchName          = "fitResult_v2";         // fit result tree branch name
	bool         cauchyPriors        = false;
	string       ampDirName          = ".";                    // decay amplitude directory name
	bool         useRootAmps         = false;                  // if true .root amplitude files are read
	//bool         useRootAmps         = true;                 // if true .root amplitude files are read
	bool         useNormalizedAmps   = false;                  // if true normalized amplitudes are used
	string       normIntFileName     = "";                     // file with normalization integrals
	string       accIntFileName      = "";                     // file with acceptance integrals
	unsigned int numbAccEvents       = 0;                      // number of events used for acceptance integrals
	unsigned int rank                = 1;                      // rank of fit
	double       minimizerTolerance  = 0.001;                  // minimizer tolerance
	bool         quiet               = false;
	extern char* optarg;
	// extern int optind;
	int c;
	while ((c = getopt(argc, argv, "l:u:w:i:o:Cd:RNn:a:A:r:t:T:B:qh")) != -1)
		switch (c) {
		case 'l':
			massBinMin = atof(optarg);
			break;
		case 'u':
			massBinMax = atof(optarg);
			break;
		case 'w':
			waveListFileName = optarg;
			break;
		case 'i':
			fitResultFileName = optarg;
			break;
		case 'o':
			outputFileName = optarg;
			break;
		case 'C':
			cauchyPriors = true;
			break;
		case 'd':
			ampDirName = optarg;
			break;
		case 'R':
#ifdef USE_STD_COMPLEX_TREE_LEAFS
			useRootAmps = true;
#endif
			break;
		case 'N':
			useNormalizedAmps = true;
			break;
		case 'n':
			normIntFileName = optarg;
			break;
		case 'a':
			accIntFileName = optarg;
			break;
		case 'A':
			numbAccEvents = atoi(optarg);
			break;
		case 'r':
			rank = atoi(optarg);
			break;
		case 't':
			minimizerTolerance = atof(optarg);
			break;
		case 'T':
			treeName = optarg;
			break;
		case 'B':
			branchName = optarg;
			break;
		case 'q':
			quiet = true;
			break;
		case 'h':
			usage(progName);
			break;
		}
	if (normIntFileName.length() <= 1) {
		normIntFileName = "norm.int";
		printWarn << "using default normalization integral file '" << normIntFileName << "'" << endl;
	}
	if (accIntFileName.length() <= 1) {
		accIntFileName = "norm.int";
		printWarn << "using default acceptance normalization integral file "
		          << "'" << accIntFileName << "'" << endl;
	}
	if (waveListFileName.length() <= 1) {
		printErr << "no wavelist file specified. aborting." << endl;
		usage(progName, 1);
	}
	// report parameters
	printInfo << "running " << progName << " with the following parameters:" << endl;
	cout << "    mass bin [" <<  massBinMin << ", " <<  massBinMax << "] MeV/c^2" << endl
	     << "    path to wave list file ......................... '" << waveListFileName << "'" << endl
	     << "    path to amplitude directory .................... '" << ampDirName       << "'" << endl
	     << "    use .root amplitude files ...................... "  << yesNo(useRootAmps)      << endl
	     << "    use normalization .............................. "  << yesNo(useNormalizedAmps) << endl
	     << "    path to file with normalization integral ....... '" << normIntFileName  << "'" << endl
	     << "    path to file with acceptance integral .......... '" << accIntFileName   << "'" << endl
	     << "    number of acceptance norm. events .............. "  << numbAccEvents    << endl
	     << "    rank of spin density matrix .................... "  << rank                    << endl
	     << "    minimizer tolerance ............................ "  << minimizerTolerance << endl
	     << "    quiet .......................................... "  << yesNo(quiet) << endl;

	// ---------------------------------------------------------------------------
	// setup likelihood function
	printInfo << "creating and setting up likelihood function" << endl;
	pwaLikelihood<complex<double> > L;
	if (quiet)
		L.setQuiet();
	L.useNormalizedAmps(useNormalizedAmps);
	L.init(rank, waveListFileName, normIntFileName, accIntFileName,
	       ampDirName, numbAccEvents, useRootAmps);
	if (not quiet)
		cout << L << endl;
	const unsigned int nmbPar  = L.NDim();

	if(cauchyPriors) {
		L.setPriorType(L.HALF_CAUCHY);
	}

	printInfo << "using prior: ";
	switch(L.priorType())
	{
		case pwaLikelihood<complex<double> >::FLAT:
			cout << "flat" << endl;
			break;
		case pwaLikelihood<complex<double> >::HALF_CAUCHY:
			cout << "half-cauchy" << endl;
			break;
	}

	TFile* fitResultFile = TFile::Open(fitResultFileName.c_str(), "READ");
	if(not fitResultFile) {
		printErr << "could not open fit result file '" << fitResultFileName << "'. Aborting..." << endl;
		return 1;
	}

	TTree* fitResultTree = 0;
	fitResultFile->GetObject(treeName.c_str(), fitResultTree);
	fitResultTree = fitResultTree->CloneTree();
	if(not fitResultTree) {
		printErr << "could not find result tree in fit result file '" << fitResultFileName << "'. Aborting..." << endl;
		return 1;
	}

	fitResult* result = 0;
	fitResultTree->SetBranchAddress(branchName.c_str(), &result);

	if(fitResultTree->GetEntries() != 1) {
		printErr << "result tree has more than one entry, NOT IMPLEMENTED" << endl;
		return 1;
	}
	fitResultTree->GetEntry(0);

	std::vector<double> pars(nmbPar);
	for(unsigned int i = 0; i < nmbPar; ++i) {
		const string parName = L.parName(i);
		pars[i] = result->fitParameter(parName);
		printInfo << "setting parameter[" << i << "] (name: '" << parName << "') = " << pars[i] << endl;
	}

	TMatrixT<double> hessian = L.HessianAnalytically(&pars[0]);
	TVectorT<double> eigenvalues;
	hessian.EigenVectors(eigenvalues);
	if (not quiet) {
		printInfo << "analytical Hessian eigenvalues:" << endl;
	}
	for(int i=0; i<eigenvalues.GetNrows(); i++) {
		if (not quiet) {
			cout << "	" << eigenvalues[i] << endl;
		}
		if (eigenvalues[i] <= 0.) {
			printWarn << "eigenvalue " << i << " of Hessian is non-positive (" << eigenvalues[i] << ")." << endl;
		}
	}
	const TMatrixT<double> covMatrix = L.CovarianceMatrixAnalytically(hessian);
	if (not quiet) {
		printInfo << "analytical covariance matrix:" << endl;
		covMatrix.Print();
	}
	TMatrixT<double> oldCovMatrix = result->fitParCovMatrix();

	printInfo << "writing result to '" << outputFileName << "'" << endl;
	if(oldCovMatrix.GetNcols() > 0 and oldCovMatrix.GetNrows() > 0) {
		printWarn << "fit result file already has a valid covariance matrix." << endl;
	}
	// open output file and create tree for writing
	TFile* outFile = new TFile(outputFileName.c_str(), "UPDATE");
	if ((not outFile) or outFile->IsZombie())
		printWarn << "cannot open output file '" << outputFileName << "'. "
				  << "no results will be written." << endl;
	else {
		TTree* tree = new TTree(treeName.c_str(), treeName.c_str());
		fitResult newResult;
		const unsigned int& nmbEvents = result->nmbEvents();
		const unsigned int& normNmbEvents = result->normNmbEvents();
		const double& massBinCenter = result->massBinCenter();
		const double& logLikelihood = result->logLikelihood();
		const int& rank = result->rank();
		const vector<TComplex> prodAmpsTComplex = result->prodAmps();
		const unsigned int nmbProdAmps = prodAmpsTComplex.size();
		vector<complex<double> > prodAmps(nmbProdAmps);
		for(unsigned int i = 0; i < nmbProdAmps; ++i) {
			prodAmps[i] = complex<double>(prodAmpsTComplex[i].Re(), prodAmpsTComplex[i].Im());
		}
		const vector<string>& prodAmpNames = result->prodAmpNames();
		const vector<pair<int, int> >& fitParCovMatrixIndices = result->fitParCovIndices();
		complexMatrix normIntegral = result->normIntegralMatrix();
		complexMatrix accIntegral = result->acceptedNormIntegralMatrix();
		vector<double> phaseSpaceIntegral = result->phaseSpaceIntegralVector();
		const bool& converged = result->converged();
		const bool& hasHesse = result->hasHessian();
		newResult.fill(
			nmbEvents,
			normNmbEvents,
			massBinCenter,
			logLikelihood,
			rank,
			prodAmps,
			prodAmpNames,
			covMatrix,
			fitParCovMatrixIndices,
			normIntegral,
			accIntegral,
			phaseSpaceIntegral,
			converged,
			hasHesse
		);
		tree->Branch(branchName.c_str(), &newResult);
		tree->Fill();
		tree->Write();
		outFile->Close();
		printSucc << "covariance matrix successfully added to output file '" << outputFileName <<"'." << endl;
	}

	fitResultFile->Close();

	return 0;
}
