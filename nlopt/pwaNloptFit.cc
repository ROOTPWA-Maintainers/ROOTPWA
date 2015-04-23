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
#include "TVectorT.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TStopwatch.h"
#include "Minuit2/Minuit2Minimizer.h"

#include <nlopt.hpp>

#include "reportingUtilsEnvironment.h"
#include "conversionUtils.hpp"
#include "pwaLikelihood.h"
#include "fitResult.h"
#include "partialWaveFitHelper.h"
#include "amplitudeTreeLeaf.h"


using namespace std;
using namespace ROOT::Math;
using namespace rpwa;


double rpwaNloptFunc(unsigned n, const double* x, double* gradient, void* func_data)
{
	pwaLikelihood<complex<double> >* L = (pwaLikelihood<complex<double> >*)func_data;
	if(n != L->nmbPars()) {
		printErr << "parameter mismatch between NLopt and pwaLikelihood. Aborting..." << endl;
		return 1;
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
	     << " -l # -u # -w wavelist [-d amplitude directory -R -o outfile -s seed -x [startvalue] -N -n normfile"
	     << " -a normfile -A # normalisation events -r rank -t # -m # -C -q -z -h]" << endl
	     << "    where:" << endl
	     << "        -l #       lower edge of mass bin [MeV/c^2]" << endl
	     << "        -u #       upper edge of mass bin [MeV/c^2]" << endl
	     << "        -w file    path to wavelist file" << endl
	     << "        -d dir     path to directory with decay amplitude files (default: '.')" << endl
	     << "        -R         use .root amplitude files (default: false)" << endl
	     << "        -o file    path to output file (default: 'fitresult.root')" << endl
	     << "        -s #       seed for random start values (default: 1234567)" << endl
	     << "        -x #       use fixed instead of random start values (default: 0.01)" << endl

	     << "        -N         use normalization of decay amplitudes (default: false)" << endl
	     << "        -n file    path to normalization integral file (default: 'norm.int')" << endl
	     << "        -a file    path to acceptance integral file (default: 'norm.int')" << endl
	     << "        -A #       number of input events to normalize acceptance to (default: use number of events from acceptance integral file)" << endl
	     << "        -r #       rank of spin density matrix (default: 1)" << endl
	     << "        -t #       relative parameter tolerance (default: 0.0001)" << endl
	     << "        -m #       absolute likelihood tolerance (default: 0.000001)" << endl
	     << "        -C         use half-Cauchy priors (default: false)" << endl
	     << "        -q         run quietly (default: false)" << endl
	     << "        -z         save space by not saving integral and covariance matrices (default: false)" << endl
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
	double       minimizerTolerance  = 1e-4;                   // minimizer tolerance
	double       likelihoodTolerance = 1e-6;                   // tolerance of likelihood function
	bool         cauchy              = false;
	bool         quiet               = false;
	bool         saveSpace           = false;
	extern char* optarg;
	// extern int optind;
	int c;
	while ((c = getopt(argc, argv, "l:u:w:d:Ro:s:x::Nn:a:A:r:t:m:Cqzh")) != -1)
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
			useRootAmps = true;
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
		case 'C':
			cauchy = true;
			break;
		case 'q':
			quiet = true;
			break;
		case 'z':
			saveSpace = true;
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
		printErr << "no wavelist file specified. Aborting..." << endl;
		usage(progName, 1);
	}
	// report parameters
	printInfo << "running " << progName << " with the following parameters:" << endl;
	cout << "    mass bin [" <<  massBinMin << ", " <<  massBinMax << "] MeV/c^2" << endl
	     << "    path to wave list file ......................... '" << waveListFileName << "'" << endl
	     << "    path to amplitude directory .................... '" << ampDirName       << "'" << endl
	     << "    use .root amplitude files ...................... "  << yesNo(useRootAmps)      << endl
	     << "    path to output file ............................ '" << outFileName      << "'" << endl
	     << "    seed for random start values ................... "  << startValSeed            << endl;
	if (useFixedStartValues)
		cout << "    using fixed instead of random start values ..... " << defaultStartValue << endl;
	cout << "    use normalization .............................. "  << yesNo(useNormalizedAmps) << endl
	     << "        path to file with normalization integral ... '" << normIntFileName  << "'" << endl
	     << "        path to file with acceptance integral ...... '" << accIntFileName   << "'" << endl
	     << "        number of acceptance norm. events .......... "  << numbAccEvents    << endl
	     << "    rank of spin density matrix .................... "  << rank                    << endl
	     << "    relative parameter tolerance.................... "  << minimizerTolerance << endl
	     << "    absolute likelihood tolerance................... "  << likelihoodTolerance << endl
	     << "    using half-Cauchy priors........................ "  << yesNo(cauchy) << endl
	     << "    saving integral and covariance matrices......... "  << yesNo(not saveSpace) << endl
	     << "    quiet .......................................... "  << yesNo(quiet) << endl;

	// ---------------------------------------------------------------------------
	// setup likelihood function
	const double massBinCenter  = (massBinMin + massBinMax) / 2;
	printInfo << "creating and setting up likelihood function" << endl;
	pwaLikelihood<complex<double> > L;
	if (quiet)
		L.setQuiet();
	L.useNormalizedAmps(useNormalizedAmps);
	L.init(rank, massBinCenter, waveListFileName, normIntFileName, accIntFileName,
	       ampDirName, numbAccEvents, useRootAmps);
	if (not quiet)
		cout << L << endl;
	const unsigned int nmbPar  = L.NDim();
	const unsigned int nmbEvts = L.nmbEvents();
	const double sqrtNmbEvts = sqrt((double)nmbEvts);

	if (cauchy)
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

	unsigned int maxParNameLength = 0;
	vector<double>params(nmbPar, defaultStartValue);
	if(not useFixedStartValues) {
		TRandom3 random(startValSeed);
		for(unsigned int i = 0; i < params.size(); ++i)
		{
			if(L.parFixed(i)) {
				printErr << "thresholds are not implemented for fits with NLopt. Aborting..." << endl;
				return 1;
			}
			const string& parName = L.parName(i);
			if(parName.length() > maxParNameLength) {
				maxParNameLength = parName.length();
			}
			params[i] = random.Uniform(defaultStartValue, sqrtNmbEvts);
			if(random.Rndm() > 0.5) {
				params[i] *= -1.;
			}
			cout << "    setting parameter [" << setw(3) << i << "] "
			     << setw(maxParNameLength) << parName << " = " << maxPrecisionAlign(params[i]) << endl;
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
	std::vector<double> correctParams = L.CorrectParamSigns(params.data());
	double newLikelihood = L.DoEval(correctParams.data());
	if(likeli != newLikelihood) {
		printErr << "Flipping signs according to sign conventions changed the likelihood (from " << likeli << " to " << newLikelihood << ")." << endl;
		return 1;
	} else {
		printInfo << "Likelihood unchanged at " << newLikelihood << " by flipping signs according to conventions." << endl;
	}
	TMatrixT<double> fitParCovMatrix(0, 0);
	bool hasHessian = false;
	if(not saveSpace) {
		TMatrixT<double> hessian = L.Hessian(correctParams.data());
		fitParCovMatrix.ResizeTo(nmbPar, nmbPar);
		fitParCovMatrix = L.CovarianceMatrix(hessian);
		TVectorT<double> eigenvalues;
		TMatrixT<double> eigenvectors;
		rpwa::partialWaveFitHelper::getEigenvectors(L, hessian, eigenvectors, eigenvalues);
		if (not quiet) {
			printInfo << "analytical Hessian eigenvalues:" << endl;
		}
		for(int i=0; i<eigenvalues.GetNrows(); ++i) {
			if (not quiet) {
				cout << "    " << maxPrecisionAlign(eigenvalues[i]) << endl;
			}
			if (eigenvalues[i] <= 0.) {
				printWarn << "eigenvalue " << i << " of (analytic) Hessian is not positive (" << maxPrecisionAlign(eigenvalues[i]) << ")." << endl;
				converged = false;
			}
		}
		if(converged) hasHessian = true;
	}
	if(converged) {
		printSucc << "minimization succeeded." << endl;
	} else {
		printWarn << "minimization failed." << endl;
	}

	// ---------------------------------------------------------------------------
	// print results
	printInfo << "minimization result:" << endl;
	for (unsigned int i = 0; i < nmbPar; ++i) {
		cout << "    parameter [" << setw(3) << i << "] "
		     << setw(maxParNameLength) << L.parName(i) << " = "
		     << setw(12) << maxPrecisionAlign(correctParams[i]) << " +- ";
		if(not saveSpace) {
			cout << setw(12) << maxPrecisionAlign(sqrt(fitParCovMatrix(i, i)));
		} else {
			cout << setw(12) << "[not available]";
		}
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
				L.buildProdAmpArrays(correctParams.data(), prodAmps, fitParCovMatrixIndices, prodAmpNames, true);
				complexMatrix normIntegral(0, 0);  // normalization integral over full phase space without acceptance
				complexMatrix accIntegral (0, 0);  // normalization integral over full phase space with acceptance
				const unsigned int nmbWaves = L.nmbWaves() + 1;  // flat wave is not included in L.nmbWaves()
				vector<double> phaseSpaceIntegral;
				if(not saveSpace) {
					normIntegral.resizeTo(nmbWaves, nmbWaves);  // normalization integral over full phase space without acceptance
					accIntegral.resizeTo(nmbWaves, nmbWaves);  // normalization integral over full phase space with acceptance
					L.getIntegralMatrices(normIntegral, accIntegral, phaseSpaceIntegral);
				}
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
				             normIntegral,
				             accIntegral,
				             phaseSpaceIntegral,  // contains the sqrt of the integral matrix diagonal elements!!!
				             converged,
				             hasHessian);
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
