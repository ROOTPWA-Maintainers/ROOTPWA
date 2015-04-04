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
#include "partialWaveFitHelper.h"
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
	     << " [-d amplitude directory -R] -i infile [-o outfile -N -n normfile"
	     << " -a normfile -A # normalisation events -C -q -h]" << endl
	     << "    where:" << endl
	     << "        -d dir     path to directory with decay amplitude files (default: '.')" << endl
	     << "        -R         use .root amplitude files (default: false)" << endl
	     << "        -i file    path to input file" << endl
	     << "        -o file    path to output file (default: 'fitresult.root')" << endl
	     << "        -N         use normalization of decay amplitudes (default: false)" << endl
	     << "        -n file    path to normalization integral file (default: 'norm.int')" << endl
	     << "        -a file    path to acceptance integral file (default: 'norm.int')" << endl
	     << "        -A #       number of input events to normalize acceptance to (default: use number of events from acceptance integral file)" << endl
	     << "        -C         use half-Cauchy priors" << endl
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

#if ROOT_VERSION_CODE < ROOT_VERSION(6, 0, 0)
	// force loading predefined std::complex dictionary
	// see http://root.cern.ch/phpBB3/viewtopic.php?f=5&t=9618&p=50164
	gROOT->ProcessLine("#include <complex>");
#endif

	// ---------------------------------------------------------------------------
	// internal parameters
	const string       valTreeName           = "pwa";
	const string       valBranchName         = "fitResult_v2";

	// ---------------------------------------------------------------------------
	// parse command line options
	const string progName            = argv[0];
	string       ampDirName          = ".";                    // decay amplitude directory name
	bool         useRootAmps         = false;                  // if true .root amplitude files are read
	//bool         useRootAmps         = true;                   // if true .root amplitude files are read
	string       inFileName          = "";                     // input filename
	string       outFileName         = "fitresult.root";       // output filename
	bool         useNormalizedAmps   = false;                  // if true normalized amplitudes are used
	string       normIntFileName     = "";                     // file with normalization integrals
	string       accIntFileName      = "";                     // file with acceptance integrals
	unsigned int numbAccEvents       = 0;                      // number of events used for acceptance integrals
	bool         cauchyPriors        = false;
	bool         quiet               = false;
	extern char* optarg;
	// extern int optind;
	int c;
	while ((c = getopt(argc, argv, "d:Ri:o:Nn:a:A:Cqh")) != -1)
		switch (c) {
		case 'd':
			ampDirName = optarg;
			break;
		case 'R':
			useRootAmps = true;
			break;
		case 'i':
			inFileName = optarg;
			break;
		case 'o':
			outFileName = optarg;
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
		case 'C':
			cauchyPriors = true;
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
	if (inFileName.length() <= 1) {
		printErr << "no input file name specified. aborting." << endl;
		usage(progName, 1);
	}
	// report parameters
	printInfo << "running " << progName << " with the following parameters:" << endl;
	cout << "    path to amplitude directory .................... '" << ampDirName               << "'" << endl
	     << "    use .root amplitude files ...................... "  << yesNo(useRootAmps)       << endl
	     << "    path to input file ............................. '" << inFileName               << "'" << endl
	     << "    path to output file ............................ '" << outFileName              << "'" << endl
	     << "    use normalization .............................. "  << yesNo(useNormalizedAmps) << endl
	     << "        path to file with normalization integral ... '" << normIntFileName          << "'" << endl
	     << "        path to file with acceptance integral ...... '" << accIntFileName           << "'" << endl
	     << "        number of acceptance norm. events .......... "  << numbAccEvents            << endl
	     << "    use half-Cauchy priors ......................... '" << yesNo(cauchyPriors)      << "'" << endl
	     << "    quiet .......................................... '" << yesNo(quiet)             << "'" << endl;

	TFile* inFile = TFile::Open(inFileName.c_str(), "READ");
	if(not inFile || inFile->IsZombie()) {
		printErr << "could not open input file '" << inFileName << "'. aborting." << endl;
		return 1;
	}

	TTree* inTree = 0;
	inFile->GetObject(valTreeName.c_str(), inTree);
	if(not inTree) {
		printErr << "could not find result tree '" << valTreeName << "' in input file '" << inFileName << "'. aborting." << endl;
		return 1;
	}

	fitResult* result = 0;
	inTree->SetBranchAddress(valBranchName.c_str(), &result);

	if(inTree->GetEntries() != 1) {
		printErr << "result tree '" << valTreeName << "' has more than one entry, NOT IMPLEMENTED." << endl;
		return 1;
	}
	inTree->GetEntry(0);

	// create temporary file to store the wavelist
	char tempFileName[] = "XXXXXX";
	close(mkstemp(tempFileName));
	const string waveListFileName(tempFileName);
	ofstream waveListFile(waveListFileName.c_str());
	rpwa::partialWaveFitHelper::extractWaveList(*result, waveListFile);
	waveListFile.close();

	printInfo << "parameters extracted from input fit result" << endl;
	cout << "    mass bin centered at ........................... "  << result->massBinCenter() << " MeV/c^2" << endl
	     << "    path to temporary wave list file ............... '" << waveListFileName        << "'" << endl
	     << "    rank of spin density matrix .................... "  << result->rank()          << endl;

	// ---------------------------------------------------------------------------
	// setup likelihood function
	printInfo << "creating and setting up likelihood function" << endl;
	pwaLikelihood<complex<double> > L;
	if (quiet)
		L.setQuiet();
	L.useNormalizedAmps(useNormalizedAmps);
	if (cauchyPriors)
		L.setPriorType(L.HALF_CAUCHY);
	L.init(result->rank(), waveListFileName, normIntFileName, accIntFileName,
	       ampDirName, numbAccEvents, useRootAmps);
	remove(waveListFileName.c_str());
	if (not quiet)
		cout << L << endl;
	const unsigned int nmbPar = L.NDim();

	unsigned int maxParNameLength = 0;  // maximum length of parameter names
	for(unsigned int i = 0; i < nmbPar; ++i) {
		const string& parName = L.parName(i);
		if (parName.length() > maxParNameLength)
			maxParNameLength = parName.length();
	}
	std::vector<double> pars(nmbPar);
	for(unsigned int i = 0; i < nmbPar; ++i) {
		const string& parName = L.parName(i);
		pars[i] = result->fitParameter(parName);
		printInfo << "setting parameter [" << setw(3) << i << "] " << setw(maxParNameLength) << parName << " = " << maxPrecisionAlign(pars[i]) << endl;
	}

	// analytically calculate Hessian
	const TMatrixT<double> hessian = L.Hessian(pars.data());
	// create and check Hessian eigenvalues
	TVectorT<double> eigenvalues;
	hessian.EigenVectors(eigenvalues);
	if (not quiet) {
		printInfo << "eigenvalues of (analytic) Hessian:" << endl;
	}
	for(int i=0; i<eigenvalues.GetNrows(); ++i) {
		if (not quiet) {
			cout << "    " << maxPrecision(eigenvalues[i]) << endl;
		}
		if (eigenvalues[i] <= 0.) {
			printWarn << "eigenvalue " << i << " of Hessian is not positive (" << maxPrecisionAlign(eigenvalues[i]) << ")." << endl;
		}
	}
	const TMatrixT<double> covMatrix = L.CovarianceMatrix(hessian);
	if (not quiet) {
		printInfo << "(analytic) covariance matrix:" << endl;
		covMatrix.Print();
	}

	const TMatrixT<double> oldCovMatrix = result->fitParCovMatrix();
	if(oldCovMatrix.GetNcols() > 0 and oldCovMatrix.GetNrows() > 0) {
		printWarn << "fit result from input already has a covariance matrix. it will be overwritten." << endl;
	}

	printInfo << "writing result to '" << outFileName << "'" << endl;
	// open output file and create tree for writing
	TFile* outFile = new TFile(outFileName.c_str(), "NEW");
	if ((not outFile) or outFile->IsZombie()) {
		printErr << "cannot open output file '" << outFileName << "'. "
		         << "no results will be written." << endl;
		return 1;
	} else {
		TTree* tree = new TTree(valTreeName.c_str(), valTreeName.c_str());
		fitResult* newResult = new fitResult();
		tree->Branch(valBranchName.c_str(), &newResult);

		const unsigned int             nmbEvents              = result->nmbEvents();
		const unsigned int             normNmbEvents          = result->normNmbEvents();
		const double                   massBinCenter          = result->massBinCenter();
		const double                   logLikelihood          = result->logLikelihood();
		const int                      rank                   = result->rank();

		const vector<TComplex>&        prodAmpsTComplex       = result->prodAmps();
		const unsigned int             nmbProdAmps            = prodAmpsTComplex.size();
		vector<complex<double> >       prodAmps(nmbProdAmps);
		for(unsigned int i = 0; i < nmbProdAmps; ++i) {
			prodAmps[i] = complex<double>(prodAmpsTComplex[i].Re(), prodAmpsTComplex[i].Im());
		}

		const vector<string>&          prodAmpNames           = result->prodAmpNames();

		const vector<pair<int, int> >& fitParCovMatrixIndices = result->fitParCovIndices();
		const complexMatrix&           normIntegral           = result->normIntegralMatrix();
		const complexMatrix&           accIntegral            = result->acceptedNormIntegralMatrix();
		const vector<double>&          phaseSpaceIntegral     = result->phaseSpaceIntegralVector();
		const bool                     converged              = result->converged();
		const bool                     hasHessian             = result->hasHessian();
		newResult->fill(nmbEvents,
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
		                hasHessian);

		// write result to file
		tree->Fill();
		tree->Write();
		outFile->Close();

		printSucc << "covariance matrix successfully added to output file '" << outFileName << "'." << endl;
	}

	inFile->Close();

	return 0;
}
