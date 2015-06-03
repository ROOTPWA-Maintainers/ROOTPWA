#include "pwaFit.h"

#include <complex>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TComplex.h"
#include "TRandom3.h"
#include "TVectorT.h"
#include <boost/progress.hpp>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "partialWaveFitHelper.h"

#include "reportingUtilsEnvironment.h"
#include "conversionUtils.hpp"
#include "pwaLikelihood.h"

#include <reportingUtils.hpp>


using namespace std;
using namespace ROOT::Math;
using namespace rpwa;

const std::string valTreeName   = "pwa";
const std::string valBranchName = "fitResult_v2";

rpwa::fitResultPtr
rpwa::hli::pwaFit(const rpwa::pwaLikelihood<std::complex<double> >& L,
                  const unsigned int              seed=0,
                  const double                    massBinMin=0,
                  const double                    massBinMax=0,
                  const std::string               startValFileName="",
                  const bool                      checkHessian=false,
                  const bool                      verbose=false)
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
	std::string  minimizerType[2]         = {"Minuit2", "Migrad"};  // minimizer, minimization algorithm
	int          minimizerStrategy        = 1;                      // minimizer strategy
	double       minimizerTolerance       = 1e-10;                  // minimizer tolerance
#if ROOT_VERSION_CODE >= ROOT_VERSION(5, 34, 19)
	bool         saveMinimizerMemory      = true;
#endif
//	bool         cudaEnabled              = false;                  // if true CUDA kernels are activated
	bool         quiet               = not verbose;

	// report parameters
	printInfo << "running pwaFit with the following parameters:" << endl;
	cout << "    mass bin [" <<  massBinMin << ", " <<  massBinMax << "] MeV/c^2" << endl
	     << "    path to file with start values ................. '" << startValFileName << "'" << endl
	     << "    seed for random start values ................... "  << seed            << endl;
	if (useFixedStartValues)
		cout << "    using fixed instead of random start values ..... " << defaultStartValue << endl;
	cout << "    rank of spin density matrix .................... "  << L.rank()                    << endl
	     << "    minimizer ...................................... "  << minimizerType[0] << ", " << minimizerType[1] << endl
	     << "    minimizer strategy ............................. "  << minimizerStrategy  << endl
	     << "    minimizer tolerance ............................ "  << minimizerTolerance << endl
//	     << "    CUDA acceleration .............................. "  << enDisabled(cudaEnabled) << endl
	     << "    check analytical Hessian eigenvalues............ "  << yesNo(checkHessian) << endl
	     << "    quiet .......................................... "  << yesNo(quiet) << endl;

	// ---------------------------------------------------------------------------
	// setup likelihood function
	const double massBinCenter  = (massBinMin + massBinMax) / 2;
	printInfo << "creating and setting up likelihood function" << endl;
//	if (quiet)
//		L.setQuiet();
//	L.useNormalizedAmps(true);
//#ifdef USE_CUDA
//	L.enableCuda(cudaEnabled);
//#endif
//	L.init(rank, massBinCenter, waveNames, waveThresholds, normMatrix, accMatrix, ampTrees, accNormEvents);
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
		printErr << "could not create minimizer. exiting." << std::endl;
		throw;
	}

	// special for Minuit2
#if ROOT_VERSION_CODE >= ROOT_VERSION(5, 34, 19)
	if(saveMinimizerMemory and dynamic_cast<ROOT::Minuit2::Minuit2Minimizer*>(minimizer)) {
			((ROOT::Minuit2::Minuit2Minimizer*)minimizer)->SetStorageLevel(0);
			printInfo << "Minuit2 storage level set to 0." << std::endl;
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
	printInfo << "reading start values from '" << startValFileName << "'" << std::endl;
	fitResult*   startFitResult = NULL;
	bool         startValValid  = false;
	TFile*       startValFile   = NULL;
	if (startValFileName.length() <= 2)
		printWarn << "start value file name '" << startValFileName << "' is invalid. "
		          << "using default start values." << std::endl;
	else {
		// open root file
		startValFile = TFile::Open(startValFileName.c_str(), "READ");
		if (not startValFile or startValFile->IsZombie())
			printWarn << "cannot open start value file '" << startValFileName << "'. "
			          << "using default start values." << std::endl;
		else {
			// get tree with start values
			TTree* tree;
			startValFile->GetObject(valTreeName.c_str(), tree);
			if (not tree)
				printWarn << "cannot find start value tree '"<< valTreeName << "' in file "
				          << "'" << startValFileName << "'" << std::endl;
			else {
				startFitResult = new rpwa::fitResult();
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
	printInfo << "setting start values for " << nmbPar << " parameters" << std::endl
	          << "    parameter naming scheme is: V[rank index]_[IGJPCMR][isobar spec]" << std::endl;
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
		TRandom3     random(seed);
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
			throw;
		}
		// cleanup
//		if(startValFile) {
//			startValFile->Close();
//			delete startValFile;
//			startValFile = NULL;
//		}
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
			throw;
		} else {
			printInfo << "Likelihood unchanged at " << maxPrecisionAlign(newLikelihood) << " by flipping signs according to conventions." << endl;
		}

		if (checkHessian) {
			// analytically calculate Hessian
			TMatrixT<double> hessian = L.Hessian(correctParams.data());
			// create and check Hessian eigenvalues
			TVectorT<double> eigenvalues;
			TMatrixT<double> eigenvectors;
			rpwa::partialWaveFitHelper::getEigenvectors(L, hessian, eigenvectors, eigenvalues);
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
//		printInfo << *minimizer;
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
//			printInfo << *minimizer;
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
	// get data structures to construct fitResult
	fitResult* result       = new fitResult();
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
	const int normNmbEvents = 1;  // number of events to normalize to

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
	             L.rank(),
	             prodAmps,
	             prodAmpNames,
	             fitParCovMatrix,
	             fitParCovMatrixIndices,
	             normIntegral,
	             accIntegral,
	             phaseSpaceIntegral,  // contains the sqrt of the integral matrix diagonal elements!!!
	             converged,
	             hasHessian);
	if (minimizer)
		delete minimizer;
	return rpwa::fitResultPtr(result);
}
