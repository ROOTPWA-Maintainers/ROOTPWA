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
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TStopwatch.h"
#include "Minuit2/Minuit2Minimizer.h"

#include <nlopt.hpp>

#include "reportingUtilsEnvironment.h"
#include "conversionUtils.hpp"
#include "pwaLikelihood.h"
#include "fitResult.h"

#include "amplitudeTreeLeaf.h"


using namespace std;
using namespace ROOT::Math;
using namespace rpwa;


double rpwaNloptFunc(unsigned n, const double* x, double* gradient, void* func_data)
{
	pwaLikelihood<complex<double> >* L = (pwaLikelihood<complex<double> >*)func_data;
	if(n != L->nmbPars()) {
		printErr << "parameter mismatch between NLopt and pwaLikelihood. Aborting..." << endl;
		throw;
	}
	double likeli;
	if(gradient) {
		L->FdF(x, likeli, gradient);
	} else {
		likeli = L->DoEval(x);
	}
	return likeli;
}


void
usage(const string& progName,
      const int     errCode = 0)
{

	cerr << "performs PWA fit for given mass bin and list of waves" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " -l # -u # -w wavelist [-d amplitude directory -R -o outfile -N -n normfile"
	     << " [-a normfile] -r rank [-t # -m # -q -h]" << endl
	     << "    where:" << endl
	     << "        -l #       lower edge of mass bin [MeV/c^2]" << endl
	     << "        -u #       upper edge of mass bin [MeV/c^2]" << endl
	     << "        -w file    path to wavelist file" << endl
	     << "        -d dir     path to directory with decay amplitude files (default: '.')" << endl
#ifdef USE_STD_COMPLEX_TREE_LEAFS
	     << "        -R         use .root amplitude files (default: false)" << endl
#else
	     << "        -R         use .root amplitude files [not supported; ROOT version too low]" << endl
#endif
	     << "        -o file    path to output file (default: 'fitresult.root')" << endl
	     << "        -s #       seed for random start values (default: 1234567)" << endl
	     << "        -x #       use fixed instead of random start values (default: 0.01)" << endl

	     << "        -N         use normalization of decay amplitudes (default: false)" << endl
	     << "        -n file    path to normalization integral file (default: 'norm.int')" << endl
	     << "        -a file    path to acceptance integral file (default: 'norm.int')" << endl
	     << "        -A #       number of input events to normalize acceptance to" << endl
	     << "        -r #       rank of spin density matrix (default: 1)" << endl
	     << "        -t #       relative parameter tolerance (default: 0.001)" << endl
	     << "        -m #       absolute likelihood tolerance (default: 0.001)" << endl
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
	// internal parameters
	const string       valTreeName           = "pwa";
	const string       valBranchName         = "fitResult_v2";
	double             defaultStartValue     = 0.01;
	bool               useFixedStartValues   = false;
	const unsigned int maxNmbOfIterations    = 50000;
	int                startValSeed          = 1234567;

	// ---------------------------------------------------------------------------
	// parse command line options
	const string progName            = argv[0];
	double       massBinMin          = 0;                      // [MeV/c^2]
	double       massBinMax          = 0;                      // [MeV/c^2]
	string       waveListFileName    = "";                     // wavelist filename
	string       ampDirName          = ".";                    // decay amplitude directory name
	bool         useRootAmps         = false;                  // if true .root amplitude files are read
	//bool         useRootAmps         = true;                   // if true .root amplitude files are read
	string       outFileName         = "fitresult.root";       // output filename
	bool         useNormalizedAmps   = false;                  // if true normalized amplitudes are used
	string       normIntFileName     = "";                     // file with normalization integrals
	string       accIntFileName      = "";                     // file with acceptance integrals
	unsigned int numbAccEvents       = 0;                      // number of events used for acceptance integrals
	unsigned int rank                = 1;                      // rank of fit
	double       minimizerTolerance  = 1e-4;                  // minimizer tolerance
	double       likelihoodTolerance = 1e-6;                   // tolerance of likelihood function
	bool         quiet               = false;
	extern char* optarg;
	// extern int optind;
	int c;
	while ((c = getopt(argc, argv, "l:u:w:d:Ro:s:x::Nn:a:A:r:t:m:qh")) != -1)
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
		case 'd':
			ampDirName = optarg;
			break;
		case 'R':
#ifdef USE_STD_COMPLEX_TREE_LEAFS
			useRootAmps = true;
#endif
			break;
		case 'o':
			outFileName = optarg;
			break;
		case 's':
			startValSeed = atoi(optarg);
			break;
		case 'x':
			if (optarg)
				defaultStartValue = atof(optarg);
			useFixedStartValues = true;
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
		case 'm':
			likelihoodTolerance = atof(optarg);
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
	     << "    path to output file ............................ '" << outFileName      << "'" << endl;
	if (useFixedStartValues)
		cout << "    using fixed instead of random start values ..... " << defaultStartValue << endl;
	cout << "    use normalization .............................. "  << yesNo(useNormalizedAmps) << endl
	     << "        path to file with normalization integral ... '" << normIntFileName  << "'" << endl
	     << "        path to file with acceptance integral ...... '" << accIntFileName   << "'" << endl
	     << "        number of acceptance norm. events .......... "  << numbAccEvents    << endl
	     << "    rank of spin density matrix .................... "  << rank                    << endl
	     << "    relative parameter tolerance.................... "  << minimizerTolerance << endl
	     << "    absolute likelihood tolerance................... "  << likelihoodTolerance << endl
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
	const unsigned int nmbEvts = L.nmbEvents();
	const double sqrtNmbEvts = sqrt((double)nmbEvts);
	const double massBinCenter  = (massBinMin + massBinMax) / 2;

	L.setPriorType(L.HALF_CAUCHY);

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

	nlopt::opt optimizer(nlopt::LD_LBFGS, nmbPar);
	optimizer.set_min_objective(&rpwaNloptFunc, &L);
	optimizer.set_lower_bounds(-2.*sqrtNmbEvts);
	optimizer.set_upper_bounds(2.*sqrtNmbEvts);

	vector<double>params(nmbPar, defaultStartValue);
	if(not useFixedStartValues) {
		TRandom3 random(startValSeed);
		for(unsigned int i = 0; i < params.size(); ++i)
		{
			params[i] = random.Uniform(defaultStartValue, sqrtNmbEvts);
			cout << "    setting parameter [" << setw(3) << i << "] = "
			     << maxPrecisionAlign(params[i]) << endl;
		}
	}

	optimizer.set_xtol_rel(minimizerTolerance);
	optimizer.set_maxeval(maxNmbOfIterations);
	optimizer.set_ftol_abs(likelihoodTolerance);

	{
		printInfo << "optimizer parameters:" << endl;
		cout << "    absolute likelihood tolerance.............. " << optimizer.get_ftol_abs() << endl;
		cout << "    relative likelihood tolerance.............. " << optimizer.get_ftol_rel() << endl;
		cout << "    maximum number of likelihood evaluations... " << optimizer.get_maxeval() << endl;
		cout << "    maximum time for optimization.............. " << optimizer.get_maxtime() << endl;
		cout << "    stop when this likelihood is found......... " << optimizer.get_stopval() << endl;
		cout << "    relative parameter tolerance............... " << optimizer.get_xtol_rel() << endl;
	}

	double likeli;
	printInfo << "starting minimization." << endl;
	nlopt::result result = optimizer.optimize(params, likeli);
	bool converged = true;
	if(result < 0) {
		converged = false;
	}
	if(converged) {
		printSucc << "minimization succeeded." << endl;
	} else {
		printWarn << "minimization failed." << endl;
	}

	// ---------------------------------------------------------------------------
	// print results
	printInfo << "minimization result:" << endl;
	vector<unsigned int> parIndices = L.orderedParIndices();
	for (unsigned int i = 0; i< parIndices.size(); ++i) {
		const unsigned int parIndex = parIndices[i];
		cout << "    parameter [" << setw(3) << i << "] "
		     << setw(12) << L.parName(parIndex) << " = ";
		cout << setw(12) << maxPrecisionAlign(params[parIndex]) << " +- "
		     << setw(12) << "[not available]";
		cout << endl;
	}
	printInfo << "function call summary:" << endl;
	L.printFuncInfo(cout);

	printInfo << "optimizer status:" << endl;
	cout << "    ";
	switch(result) {
		case nlopt::SUCCESS:
			cout << "success!" << endl;
			break;
		case nlopt::STOPVAL_REACHED:
			cout << "likelihood threshold reached." << endl;
			break;
		case nlopt::FTOL_REACHED:
			cout << "function tolerance reached." << endl;
			break;
		case nlopt::XTOL_REACHED:
			cout << "parameter tolerance reached." << endl;
			break;
		case nlopt::MAXEVAL_REACHED:
			cout << "maximum number of evaluations reached." << endl;
			break;
		case nlopt::MAXTIME_REACHED:
			cout << "time limit reached." << endl;
			break;

		case nlopt::FAILURE:
			cout << "generic error." << endl;
			break;
		case nlopt::INVALID_ARGS:
			cout << "invalid arguments." << endl;
			break;
		case nlopt::OUT_OF_MEMORY:
			cout << "out of memory." << endl;
			break;
		case nlopt::ROUNDOFF_LIMITED:
			cout << "roundoff errors limited progress." << endl;
			break;
		case nlopt::FORCED_STOP:
			cout << "forced stop." << endl;
			break;
	}

	// ---------------------------------------------------------------------------
	// write out result
	printInfo << "writing result to '" << outFileName << "'" << endl;
	{
		// open output file and create tree for writing
		TFile* outFile = new TFile(outFileName.c_str(), "UPDATE");
		if ((not outFile) or outFile->IsZombie())
			printWarn << "cannot open output file '" << outFileName << "'. "
			          << "no results will be written." << endl;
		else {
			// check whether output tree already exists
			TTree*     tree;
			fitResult* result       = new fitResult();
			outFile->GetObject(valTreeName.c_str(), tree);
			if (not tree) {
				printInfo << "file '" << outFileName << "' is empty. "
				          << "creating new tree '" << valTreeName << "' for PWA result." << endl;
				tree = new TTree(valTreeName.c_str(), valTreeName.c_str());
				tree->Branch(valBranchName.c_str(), &result);
			} else {
				tree->SetBranchAddress(valBranchName.c_str(), &result);
			}

			{
				// get data structures to construct fitResult
				vector<std::complex<double> > prodAmps;                // production amplitudes
				vector<string>                prodAmpNames;            // names of production amplitudes used in fit
				vector<pair<int,int> >        fitParCovMatrixIndices;  // indices of fit parameters for real and imaginary part in covariance matrix matrix
				L.buildProdAmpArrays(&params[0], prodAmps, fitParCovMatrixIndices, prodAmpNames, true);
				TMatrixT<double> fitParCovMatrix(0, 0);          // nlopt does not calculate the covariance matrix, so do not save it
				const unsigned int nmbWaves = L.nmbWaves() + 1;  // flat wave is not included in L.nmbWaves()
				complexMatrix normIntegral(nmbWaves, nmbWaves);  // normalization integral over full phase space without acceptance
				complexMatrix accIntegral (nmbWaves, nmbWaves);  // normalization integral over full phase space with acceptance
				vector<double> phaseSpaceIntegral;
				L.getIntegralMatrices(normIntegral, accIntegral, phaseSpaceIntegral);
				const int normNmbEvents = (useNormalizedAmps) ? 1 : L.nmbEvents();  // number of events to normalize to

				cout << "filling fitResult:" << endl
				     << "    number of fit parameters ............... " << nmbPar                        << endl
				     << "    number of production amplitudes ........ " << prodAmps.size()               << endl
				     << "    number of production amplitude names ... " << prodAmpNames.size()           << endl
				     << "    number of wave names ................... " << nmbWaves                      << endl
				     << "    number of cov. matrix indices .......... " << fitParCovMatrixIndices.size() << endl
				     << "    dimension of covariance matrix ......... " << fitParCovMatrix.GetNrows() << " x " << fitParCovMatrix.GetNcols() << endl
				     << "    dimension of normalization matrix ...... " << normIntegral.nRows()       << " x " << normIntegral.nCols()       << endl
				     << "    dimension of acceptance matrix ......... " << accIntegral.nRows()        << " x " << accIntegral.nCols()        << endl;
				result->fill(L.nmbEvents(),
				             normNmbEvents,
				             massBinCenter,
				             likeli,
				             rank,
				             prodAmps,
				             prodAmpNames,
				             fitParCovMatrix,
				             fitParCovMatrixIndices,
				             normIntegral,  // contains the sqrt of the integral matrix diagonal elements!!!
				             accIntegral,
				             phaseSpaceIntegral,
				             converged,
				             false);
				//printDebug << *result;
			}

			// write result to file
			tree->Fill();
			tree->Write("", TObject::kOverwrite);
			outFile->Close();
		}
	}

	return 0;
}
