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

#include "reportingUtilsEnvironment.h"
#include "conversionUtils.hpp"
#include "pwaLikelihood.h"
#include "fitResult.h"
#ifdef USE_CUDA
#include "complex.cuh"
#include "likelihoodInterface.cuh"
#endif
#include "amplitudeTreeLeaf.h"


using namespace std;
using namespace ROOT::Math;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{

	cerr << "performs PWA fit for given mass bin and list of waves" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " -l # -u # -w wavelist [-d amplitude directory -R -o outfile -S start value file -s seed -x [start value] -N -n normfile"
	     << " -a normfile -A # normalisation events -r rank -M minimizer -m algorithm -g strategy -t # -e -c -H -q -h]" << endl
	     << "    where:" << endl
	     << "        -l #       lower edge of mass bin [MeV/c^2]" << endl
	     << "        -u #       upper edge of mass bin [MeV/c^2]" << endl
	     << "        -w file    path to wavelist file" << endl
	     << "        -d dir     path to directory with decay amplitude files (default: '.')" << endl
	     << "        -R         use .root amplitude files (default: false)" << endl
	     << "        -o file    path to output file (default: 'fitresult.root')" << endl
	     << "        -S file    path to file with start values (default: none; highest priority)" << endl
	     << "        -s #       seed for random start values (default: 1234567)" << endl
	     << "        -x #       use fixed instead of random start values (default: 0.01)" << endl

	     << "        -N         use normalization of decay amplitudes (default: false)" << endl
	     << "        -n file    path to normalization integral file (default: 'norm.int')" << endl
	     << "        -a file    path to acceptance integral file (default: 'norm.int')" << endl
	     << "        -A #       number of input events to normalize acceptance to (default: use number of events from acceptance integral file)" << endl
	     << "        -r #       rank of spin density matrix (default: 1)" << endl
	     << "        -M name    minimizer (default: Minuit2)" << endl
	     << "        -m name    minimization algorithm (optional, default: Migrad)" << endl
	     << "                   available minimizers: Minuit:      Migrad, Simplex, Minimize, Migrad_imp" << endl
	     << "                                         Minuit2:     Migrad, Simplex, Combined, Scan, Fumili" << endl
	     << "                                         GSLMultiMin: ConjugateFR, ConjugatePR, BFGS, BFGS2, SteepestDescent" << endl
	     << "                                         GSLMultiFit: -" << endl
	     << "                                         GSLSimAn:    -" << endl
	     << "                                         Linear:      Robust" << endl
	     << "                                         Fumili:      -" << endl
	     << "        -g #       minimizer strategy: 0 = low, 1 = medium, 2 = high effort  (default: 1)" << endl
	     << "        -t #       minimizer tolerance (default: 1e-10)" << endl
#if ROOT_VERSION_CODE >= ROOT_VERSION(5, 34, 19)
	     << "        -e         set minimizer storage level to 1 (only available for Minuit2, default: 0)" << endl
#else
	     << "        -e         set minimizer storage level to 1 [not supported and not required; ROOT version too low to switch off minimizer storage]" << endl
#endif
#ifdef USE_CUDA
	     << "        -c         enable CUDA acceleration (default: off)" << endl
#else
	     << "        -c         enable CUDA acceleration [not supported by your platform]" << endl
#endif
	     << "        -H         check analytical Hessian eigenvalues (default: false)" << endl
	     << "        -q         run quietly (default: false)" << endl
	     << "        -h         print help" << endl
	     << endl;
	exit(errCode);
}


ostream&
printMinimizerStatus(ostream&   out,
                     Minimizer& minimizer)
{
	out << "minimization stopped after " << minimizer.NCalls() << " function calls. "
	    << "minimizer status summary:" << endl
	    << "    total number of parameters .......................... " << minimizer.NDim()                 << endl
	    << "    number of free parameters ........................... " << minimizer.NFree()                << endl
	    << "    maximum allowed number of iterations ................ " << minimizer.MaxIterations()        << endl
	    << "    maximum allowed number of function calls ............ " << minimizer.MaxFunctionCalls()     << endl
	    << "    minimizer status .................................... " << minimizer.Status()               << endl
	    << "    minimizer provides error and error matrix ........... " << yesNo(minimizer.ProvidesError()) << endl
	    << "    minimizer has performed detailed error validation ... " << yesNo(minimizer.IsValidError())  << endl
	    << "    estimated distance to minimum ....................... " << minimizer.Edm()                  << endl
	    << "    statistical scale used for error calculation ........ " << minimizer.ErrorDef()             << endl
	    << "    strategy ............................................ " << minimizer.Strategy()             << endl
	    << "    absolute tolerance .................................. " << minimizer.Tolerance()            << endl
	    << "    precision ........................................... " << minimizer.Precision()            << endl;
	return out;
}


ostream&
operator <<(ostream&   out,
            Minimizer& minimizer)
{
	return printMinimizerStatus(out, minimizer);
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
	double             startValStep          = 0.0005;
	const unsigned int maxNmbOfIterations    = 20000;
	const unsigned int maxNmbOfFunctionCalls = 40000;
	const bool         runHesse              = true;
	const bool         runMinos              = false;
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
	string       startValFileName    = "";                     // file with start values
	bool         useNormalizedAmps   = false;                  // if true normalized amplitudes are used
	string       normIntFileName     = "";                     // file with normalization integrals
	string       accIntFileName      = "";                     // file with acceptance integrals
	unsigned int numbAccEvents       = 0;                      // number of events used for acceptance integrals
	unsigned int rank                = 1;                      // rank of fit
	string       minimizerType[2]    = {"Minuit2", "Migrad"};  // minimizer, minimization algorithm
	int          minimizerStrategy   = 1;                      // minimizer strategy
	double       minimizerTolerance  = 1e-10;                  // minimizer tolerance
#if ROOT_VERSION_CODE >= ROOT_VERSION(5, 34, 19)
	bool         saveMinimizerMemory = true;
#endif
	bool         cudaEnabled         = false;                  // if true CUDA kernels are activated
	bool         checkHessian        = false;                  // if true checks analytical Hessian eigenvalues
	bool         quiet               = false;
	extern char* optarg;
	// extern int optind;
	int c;
	while ((c = getopt(argc, argv, "l:u:w:d:Ro:S:s:x::Nn:a:A:r:M:m:g:t:ecHqh")) != -1)
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
		case 'S':
			startValFileName = optarg;
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
		case 'M':
			minimizerType[0] = optarg;
			break;
		case 'm':
			minimizerType[1] = optarg;
			break;
		case 'g':
			minimizerStrategy = atoi(optarg);
			break;
		case 't':
			minimizerTolerance = atof(optarg);
			break;
		case 'e':
#if ROOT_VERSION_CODE >= ROOT_VERSION(5, 34, 19)
			saveMinimizerMemory = false;
#endif
			break;
		case 'c':
#ifdef USE_CUDA
			cudaEnabled = true;
#endif
			break;
		case 'H':
			checkHessian = true;
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
	     << "    path to file with start values ................. '" << startValFileName << "'" << endl
	     << "    seed for random start values ................... "  << startValSeed            << endl;
	if (useFixedStartValues)
		cout << "    using fixed instead of random start values ..... " << defaultStartValue << endl;
	cout << "    use normalization .............................. "  << yesNo(useNormalizedAmps) << endl
	     << "        path to file with normalization integral ... '" << normIntFileName  << "'" << endl
	     << "        path to file with acceptance integral ...... '" << accIntFileName   << "'" << endl
	     << "        number of acceptance norm. events .......... "  << numbAccEvents    << endl
	     << "    rank of spin density matrix .................... "  << rank                    << endl
	     << "    minimizer ...................................... "  << minimizerType[0] << ", " << minimizerType[1] << endl
	     << "    minimizer strategy ............................. "  << minimizerStrategy  << endl
	     << "    minimizer tolerance ............................ "  << minimizerTolerance << endl
	     << "    CUDA acceleration .............................. "  << enDisabled(cudaEnabled) << endl
	     << "    check analytical Hessian eigenvalues............ "  << yesNo(checkHessian) << endl
	     << "    quiet .......................................... "  << yesNo(quiet) << endl;

	// ---------------------------------------------------------------------------
	// setup likelihood function
	const double massBinCenter  = (massBinMin + massBinMax) / 2;
	printInfo << "creating and setting up likelihood function" << endl;
	pwaLikelihood<complex<double> > L;
	if (quiet)
		L.setQuiet();
	L.useNormalizedAmps(useNormalizedAmps);
#ifdef USE_CUDA
	L.enableCuda(cudaEnabled);
#endif
	L.init(rank, massBinCenter, waveListFileName, normIntFileName, accIntFileName,
	       ampDirName, numbAccEvents, useRootAmps);
	if (not quiet)
		cout << L << endl;
	const unsigned int nmbPar  = L.NDim();
	const unsigned int nmbEvts = L.nmbEvents();

	// ---------------------------------------------------------------------------
	// setup minimizer
	printInfo << "creating and setting up minimizer '" << minimizerType[0] << "' "
	          << "using algorithm '" << minimizerType[1] << "'" << endl;
	Minimizer* minimizer = Factory::CreateMinimizer(minimizerType[0], minimizerType[1]);
	if (not minimizer) {
		printErr << "could not create minimizer. exiting." << endl;
		return 1;
	}

	// special for Minuit2
#if ROOT_VERSION_CODE >= ROOT_VERSION(5, 34, 19)
	if(saveMinimizerMemory and dynamic_cast<ROOT::Minuit2::Minuit2Minimizer*>(minimizer)) {
			((ROOT::Minuit2::Minuit2Minimizer*)minimizer)->SetStorageLevel(0);
			printInfo << "Minuit2 storage level set to 0." << endl;
	}
#endif

	minimizer->SetFunction        (L);
	minimizer->SetStrategy        (minimizerStrategy);
	minimizer->SetTolerance       (minimizerTolerance);

	// setting the ErrorDef to 1 since the ROOT interface does not
	// Propagate the value. Will do the error rescaling by hand below.
	minimizer->SetErrorDef(1);
	minimizer->SetPrintLevel      ((quiet) ? 0 : 3);
	minimizer->SetMaxIterations   (maxNmbOfIterations);
	minimizer->SetMaxFunctionCalls(maxNmbOfFunctionCalls);

	// ---------------------------------------------------------------------------
	// read in fitResult with start values
	printInfo << "reading start values from '" << startValFileName << "'" << endl;
	fitResult*   startFitResult = NULL;
	bool         startValValid  = false;
	TFile*       startValFile   = NULL;
	if (startValFileName.length() <= 2)
		printWarn << "start value file name '" << startValFileName << "' is invalid. "
		          << "using default start values." << endl;
	else {
		// open root file
		startValFile = TFile::Open(startValFileName.c_str(), "READ");
		if (not startValFile or startValFile->IsZombie())
			printWarn << "cannot open start value file '" << startValFileName << "'. "
			          << "using default start values." << endl;
		else {
			// get tree with start values
			TTree* tree;
			startValFile->GetObject(valTreeName.c_str(), tree);
			if (not tree)
				printWarn << "cannot find start value tree '"<< valTreeName << "' in file "
				          << "'" << startValFileName << "'" << endl;
			else {
				startFitResult = new fitResult();
				tree->SetBranchAddress(valBranchName.c_str(), &startFitResult);
				// find tree entry which is closest to mass bin center
				unsigned int bestIndex = 0;
				double       bestMass  = 0;
				for (unsigned int i = 0; i < tree->GetEntriesFast(); ++i) {
					tree->GetEntry(i);
					if (fabs(massBinCenter - startFitResult->massBinCenter()) <= fabs(massBinCenter - bestMass)) {
						bestIndex = i;
						bestMass  = startFitResult->massBinCenter();
					}
				}
				tree->GetEntry(bestIndex);
				startValValid = true;
			}
		}
	}

	// ---------------------------------------------------------------------------
	// set start parameter values
	printInfo << "setting start values for " << nmbPar << " parameters" << endl
	          << "    parameter naming scheme is: V[rank index]_[IGJPCMR][isobar spec]" << endl;
	unsigned int maxParNameLength = 0;       // maximum length of parameter names
	{
		for (unsigned int i = 0; i < nmbPar; ++i) {
			const string parName = L.parName(i);
			if (parName.length() > maxParNameLength)
				maxParNameLength = parName.length();
		}
		// use local instance of random number generator so that other
		// code has no chance of tampering with gRandom and thus cannot
		// affect the reproducability of the start values
		TRandom3     random(startValSeed);
		const double sqrtNmbEvts = sqrt((double)nmbEvts);
		bool         success     = true;
		for (unsigned int i = 0; i < nmbPar; ++i) {
			const string parName = L.parName(i);

			double startVal;
			if (startValValid) {
				// get parameter value from fitResult
				assert(startFitResult);
				startVal = startFitResult->fitParameter(parName);
			} else {
				startVal = (useFixedStartValues) ? defaultStartValue : random.Uniform(defaultStartValue, sqrtNmbEvts);
				if(random.Rndm() > 0.5) {
					startVal *= -1.;
				}
			}

			// check if parameter needs to be fixed
			if (not L.parFixed(i)) {
				if (startVal == 0) {
					cout << "    read start value 0 for parameter " << parName << ". "
					     << "using default start value." << endl;
					startVal = (useFixedStartValues) ? defaultStartValue : random.Uniform(defaultStartValue, sqrtNmbEvts);
					if(random.Rndm() > 0.5) {
						startVal *= -1.;
					}
				}
				cout << "    setting parameter [" << setw(3) << i << "] "
				     << setw(maxParNameLength) << parName << " = " << maxPrecisionAlign(startVal) << endl;
				if (not minimizer->SetVariable(i, parName, startVal, startValStep))
					success = false;
			} else {
				cout << "    fixing parameter  [" << setw(3) << i << "] "
				     << setw(maxParNameLength) << parName << " = 0" << endl;
				if (not minimizer->SetFixedVariable(i, parName, 0.))  // fix this parameter to 0
					success = false;
			}
		}
		if (not success) {
			printErr << "something went wrong when setting log likelihood parameters. Aborting..." << endl;
			return 1;
		}
		// cleanup
		if(startValFile) {
			startValFile->Close();
			delete startValFile;
			startValFile = NULL;
		}
	}

	// ---------------------------------------------------------------------------
	// find minimum of likelihood function
	bool converged  = false;
	bool hasHessian = false;
	std::vector<double> correctParams;
	printInfo << "performing minimization" << endl;
	{
		TStopwatch timer;
		timer.Start();
		converged = minimizer->Minimize();
		timer.Stop();

		correctParams = L.CorrectParamSigns(minimizer->X());
		const double newLikelihood = L.DoEval(correctParams.data());
		if(minimizer->MinValue() != newLikelihood) {
			printErr << "Flipping signs according to sign conventions changed the likelihood (from " << maxPrecisionAlign(minimizer->MinValue()) << " to " << maxPrecisionAlign(newLikelihood) << ")." << endl;
			return 1;
		} else {
			printInfo << "Likelihood unchanged at " << maxPrecisionAlign(newLikelihood) << " by flipping signs according to conventions." << endl;
		}

		if (checkHessian) {
			// analytically calculate Hessian
			TMatrixT<double> hessian = L.Hessian(correctParams.data());
			// create reduced hessian without fixed parameters
			TMatrixT<double> reducedHessian(nmbPar - L.nmbParsFixed(), nmbPar - L.nmbParsFixed());
			unsigned int iReduced = 0;
			for(unsigned int i = 0; i < nmbPar; ++i) {
				unsigned int jReduced = 0;
				if (not L.parFixed(i)) {
					for(unsigned int j = 0; j < nmbPar; ++j) {
						if (not L.parFixed(j)) {
							reducedHessian[iReduced][jReduced] = hessian[i][j];
							jReduced++;
						}
					}
					iReduced++;
				}
			}
			// create and check Hessian eigenvalues
			TVectorT<double> eigenvalues;
			reducedHessian.EigenVectors(eigenvalues);
			if (not quiet) {
				printInfo << "eigenvalues of (analytic) Hessian:" << endl;
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
		}
		if (converged)
			printInfo << "minimization finished successfully. " << flush;
		else
			printWarn << "minimization failed. " << flush;
		cout << "used " << flush;
		timer.Print();
		printInfo << *minimizer;
		if (runHesse) {
			printInfo << "calculating Hessian matrix" << endl;
			timer.Start();
			hasHessian = minimizer->Hesse();
			timer.Stop();
			if (hasHessian) {
				printInfo << "successfully calculated Hessian matrix. " << flush;
			} else {
				printWarn << "calculation of Hessian matrix failed. " << flush;
				converged = false;
			}
			cout << "used " << flush;
			timer.Print();
			printInfo << *minimizer;
		}
	}

	// ---------------------------------------------------------------------------
	// print results
	printInfo << "minimization result:" << endl;
	const double inverseOfSqrtTwo = 1. / sqrt(2.);
	for (unsigned int i = 0; i < nmbPar; ++i) {
		cout << "    parameter [" << setw(3) << i << "] "
		     << setw(maxParNameLength) << L.parName(i) << " = ";
		if (L.parFixed(i))
			cout << correctParams[i] << " (fixed)" << endl;
		else {
			cout << setw(12) << maxPrecisionAlign(correctParams[i]) << " +- "
			     << setw(12) << maxPrecisionAlign(inverseOfSqrtTwo * minimizer->Errors()[i]);
			if (runMinos) {
				double minosErrLow = 0;
				double minosErrUp  = 0;
				const bool success = minimizer->GetMinosError(i, minosErrLow, minosErrUp);
				if (success)
					cout << "    Minos: " << "[" << minosErrLow << ", +" << minosErrUp << "]" << endl;
			} else
				cout << endl;
		}
	}
	printInfo << "function call summary:" << endl;
	L.printFuncInfo(cout);
#ifdef USE_CUDA
	printInfo << "total CUDA kernel time: "
	          << cuda::likelihoodInterface<cuda::complex<double> >::kernelTime() << " sec" << endl;
#endif

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
				TMatrixT<double> fitParCovMatrix(nmbPar, nmbPar);  // covariance matrix of fit parameters
				for(unsigned int i = 0; i < nmbPar; ++i)
					for(unsigned int j = 0; j < nmbPar; ++j)
					  // The factor 0.5 is needed because
					  // MINUIT by default assumes a Chi2
					  // function and not a loglikeli
					  // (see Minuit manual!)
					  // Note: SetErrorDef in ROOT does not work
						fitParCovMatrix[i][j] = 0.5* minimizer->CovMatrix(i, j);
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
				             minimizer->MinValue(),
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

	if (minimizer)
		delete minimizer;
	return 0;
}
