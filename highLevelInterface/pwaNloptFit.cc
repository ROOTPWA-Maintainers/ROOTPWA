#include "pwaNloptFit.h"

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

#include <nlopt.hpp>

#include "reportingUtilsEnvironment.h"
#include "conversionUtils.hpp"
#include "pwaLikelihood.h"

#include <reportingUtils.hpp>

using namespace std;
using namespace ROOT::Math;
using namespace rpwa;

const std::string valTreeName   = "pwa";
const std::string valBranchName = "fitResult_v2";

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

rpwa::fitResultPtr
rpwa::hli::pwaNloptFit(std::map<std::string, TTree*>&  ampTrees,
                       const rpwa::ampIntegralMatrix&  normMatrix,
                       rpwa::ampIntegralMatrix&        accMatrix,
                       const std::vector<std::string>& waveNames,
                       const std::vector<double>&      waveThresholds,
                       const double                    massBinMin=0,
                       const double                    massBinMax=0,
                       const unsigned int              seed=0,
                       const bool                      cauchy=false,
                       const std::string               startValFileName="",
                       const unsigned int              accEventsOverride=0,
                       const unsigned int              rank=1,
                       const bool                      checkHessian=false,
                       const bool                      saveSpace=false,
                       const bool                      verbose=false)
{

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
	const double             minimizerTolerance    = 1e-4;                   // minimizer tolerance
	const double             likelihoodTolerance   = 1e-6;                   // tolerance of likelihood function
#if ROOT_VERSION_CODE >= ROOT_VERSION(5, 34, 19)
//	bool         saveMinimizerMemory      = true;
#endif
//	bool         cudaEnabled              = false;                  // if true CUDA kernels are activated
	bool         quiet               = not verbose;

	// report parameters
	printInfo << "running pwaNloptFit with the following parameters:" << endl;
	cout << "    mass bin [" <<  massBinMin << ", " <<  massBinMax << "] MeV/c^2" << endl
	     << "    path to file with start values ................. '" << startValFileName << "'" << endl
	     << "    seed for random start values ................... "  << seed            << endl;
	if (useFixedStartValues)
		cout << "    using fixed instead of random start values ..... " << defaultStartValue << endl;
	cout << "    rank of spin density matrix .................... "  << rank                    << endl
		     << "    minimizer tolerance ............................ "  << minimizerTolerance << endl
		     << "    likelihood tolerance ........................... "  << likelihoodTolerance << endl
//	     << "    CUDA acceleration .............................. "  << enDisabled(cudaEnabled) << endl
	     << "    check analytical Hessian eigenvalues............ "  << yesNo(checkHessian) << endl
	     << "    quiet .......................................... "  << yesNo(quiet) << endl;

	// ---------------------------------------------------------------------------
	// setup likelihood function
	const double massBinCenter  = (massBinMin + massBinMax) / 2;
	unsigned int accNormEvents = 0;
	if (accEventsOverride == 0) {
		accNormEvents = normMatrix.nmbEvents();
	}
	else {
		accNormEvents = accEventsOverride;
	}
	printInfo << "creating and setting up likelihood function" << endl;
	pwaLikelihood<complex<double> > L;
	if (quiet)
		L.setQuiet();
	L.useNormalizedAmps(true);
#ifdef USE_CUDA
	L.enableCuda(cudaEnabled);
#endif
	L.init(rank, massBinCenter, waveNames, waveThresholds, normMatrix, accMatrix, ampTrees, accNormEvents);
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
		TRandom3 random(seed);
		for(unsigned int i = 0; i < params.size(); ++i)
		{
			if(L.parFixed(i)) {
				printErr << "thresholds are not implemented for fits with NLopt. Aborting..." << endl;
				throw;
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
		throw;
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
		hessian.EigenVectors(eigenvalues);
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

	rpwa::fitResult* fitResult = new rpwa::fitResult();
	vector<std::complex<double> > prodAmps;                // production amplitudes
	vector<string>                prodAmpNames;            // names of production amplitudes used in fit
	vector<pair<int,int> >        fitParCovMatrixIndices;  // indices of fit parameters for real and imaginary part in covariance matrix matrix
	L.buildProdAmpArrays(correctParams.data(), prodAmps, fitParCovMatrixIndices, prodAmpNames, true);
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
	fitResult->fill(L.nmbEvents(),
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
	return rpwa::fitResultPtr(fitResult);
}
