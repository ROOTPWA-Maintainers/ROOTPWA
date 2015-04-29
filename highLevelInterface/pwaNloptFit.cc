#include "pwaNloptFit.h"

#include <complex>

#include <TFile.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <TVectorT.h>

#include <nlopt.hpp>

#include <conversionUtils.hpp>
#include <partialWaveFitHelper.h>
#include <pwaLikelihood.h>
#include <reportingUtils.hpp>


using namespace std;
using namespace rpwa;


double
rpwaNloptFunc(unsigned int n, const double* x, double* gradient, void* func_data)
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


fitResultPtr
rpwa::hli::pwaNloptFit(map<string, TTree*>&     ampTrees,
                       const ampIntegralMatrix& normMatrix,
                       ampIntegralMatrix&       accMatrix,
                       const vector<string>&    waveNames,
                       const vector<double>&    waveThresholds,
                       const double             massBinMin = 0.,
                       const double             massBinMax = 0.,
                       const unsigned int       seed = 0,
                       const string&            startValFileName = "",
                       const unsigned int       accEventsOverride = 0,
                       const unsigned int       rank = 1,
                       const bool               cauchy = false,
                       const double             cauchyWidth = 0.5,
                       const bool               checkHessian = false,
                       const bool               saveSpace = false,
                       const bool               verbose = false)
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
	const double       defaultStartValue     = 0.01;
	const bool         useFixedStartValues   = false;
	const unsigned int maxNmbOfIterations    = 50000;
	const double       minimizerTolerance    = 1e-4;                   // minimizer tolerance
	const double       likelihoodTolerance   = 1e-6;                   // tolerance of likelihood function
	const bool         cudaEnabled           = false;                  // if true CUDA kernels are activated
	const bool         quiet                 = not verbose;

	// report parameters
	printInfo << "running pwaNloptFit with the following parameters:" << endl;
	cout << "    mass bin [" << massBinMin << ", " << massBinMax << "] MeV/c^2" << endl
	     << "    seed for random start values ................... "  << seed                    << endl;
	if (useFixedStartValues)
		cout << "    using fixed instead of random start values ..... " << defaultStartValue << endl;
	cout << "    rank of spin density matrix .................... "  << rank                    << endl
	     << "    using half-Cauchy priors........................ "  << yesNo(cauchy) << endl;
	if(cauchy) {
		cout << "    width of cauchy priors.......................... "  << cauchyWidth << endl;
	}
	cout << "    check analytical Hessian eigenvalues............ "  << yesNo(checkHessian)     << endl
	     << "    minimizer tolerance ............................ "  << minimizerTolerance      << endl
	     << "    likelihood tolerance ........................... "  << likelihoodTolerance     << endl
	     << "    saving integral and covariance matrices......... "  << yesNo(not saveSpace)    << endl
	     << "    CUDA acceleration .............................. "  << enDisabled(cudaEnabled) << endl
	     << "    quiet .......................................... "  << yesNo(quiet)            << endl;

	// ---------------------------------------------------------------------------
	// setup likelihood function
	const double massBinCenter = (massBinMin + massBinMax) / 2;
	unsigned int accNormEvents = 0;
	if (accEventsOverride == 0) {
		accNormEvents = normMatrix.nmbEvents();
	} else {
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

	if (cauchy) {
		L.setPriorType(L.HALF_CAUCHY);
		L.setCauchyWidth(cauchyWidth);
	}

	printInfo << "using prior: ";
	switch(L.priorType())
	{
		case pwaLikelihood<complex<double> >::FLAT:
			cout << "flat" << endl;
			break;
		case pwaLikelihood<complex<double> >::HALF_CAUCHY:
			cout      << "half-cauchy" << endl;
			printInfo << "cauchy width: " << L.cauchyWidth() << endl;
			break;
	}

	// ---------------------------------------------------------------------------
	// setup minimizer
	nlopt::opt optimizer(nlopt::LD_LBFGS, nmbPar);
	optimizer.set_min_objective(&rpwaNloptFunc, &L);

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

	// ---------------------------------------------------------------------------
	// set start parameter values
	printInfo << "setting start values for " << nmbPar << " parameters" << endl
	          << "    parameter naming scheme is: V[rank index]_[wave name]" << endl;
	unsigned int maxParNameLength = 0;       // maximum length of parameter names
	vector<double>params(nmbPar);
	{
		for (unsigned int i = 0; i < nmbPar; ++i) {
			const string& parName = L.parName(i);
			if (parName.length() > maxParNameLength)
				maxParNameLength = parName.length();
		}
		// use local instance of random number generator so that other
		// code has no chance of tampering with gRandom and thus cannot
		// affect the reproducability of the start values
		TRandom3     random(seed);
		const double sqrtNmbEvts = sqrt((double)nmbEvts);
		for (unsigned int i = 0; i < nmbPar; ++i) {
			const string& parName = L.parName(i);

			double startVal;
			{
				startVal = (useFixedStartValues) ? defaultStartValue : random.Uniform(defaultStartValue, sqrtNmbEvts);
				if(random.Rndm() > 0.5) {
					startVal *= -1.;
				}
			}

			// check if parameter needs to be fixed
			if (not L.parFixed(i)) {
				cout << "    setting parameter [" << setw(3) << i << "] "
				     << setw(maxParNameLength) << parName << " = " << maxPrecisionAlign(startVal) << endl;
				params[i] = startVal;
			} else {
				printErr << "thresholds are not implemented for fits with NLopt. Aborting..." << endl;
				throw;
			}
		}
	}

	// ---------------------------------------------------------------------------
	// find minimum of likelihood function
	bool converged  = false;
	bool hasHessian = false;
	double likeli;
	vector<double> correctParams;
	TMatrixT<double> fitParCovMatrix(0, 0);
	printInfo << "performing minimization." << endl;
	{
		TStopwatch timer;
		timer.Start();
		nlopt::result result = optimizer.optimize(params, likeli);
		if(result < 0) {
			converged = false;
		}
		timer.Stop();

		correctParams = L.CorrectParamSigns(params.data());
		const double newLikelihood = L.DoEval(correctParams.data());
		if(likeli != newLikelihood) {
			printErr << "Flipping signs according to sign conventions changed the likelihood (from " << maxPrecisionAlign(likeli) << " to " << maxPrecisionAlign(newLikelihood) << ")." << endl;
			throw;
		} else {
			printInfo << "Likelihood unchanged at " << maxPrecisionAlign(newLikelihood) << " by flipping signs according to conventions." << endl;
		}
		if (checkHessian || not saveSpace) {
			TMatrixT<double> hessian = L.Hessian(correctParams.data());
			// create and check Hessian eigenvalues
			TVectorT<double> eigenvalues;
			TMatrixT<double> eigenvectors;
			partialWaveFitHelper::getEigenvectors(L, hessian, eigenvectors, eigenvalues);
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
			if (not saveSpace) {
				fitParCovMatrix.ResizeTo(nmbPar, nmbPar);
				fitParCovMatrix = L.CovarianceMatrix(hessian);
				if(converged) hasHessian = true;
			}
		}
		if (converged) {
			printSucc << "minimization finished successfully. " << flush;
		} else {
			printWarn << "minimization failed. " << flush;
		}
		cout << "used " << flush;
		timer.Print();
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
#ifdef USE_CUDA
	printInfo << "total CUDA kernel time: "
	          << cuda::likelihoodInterface<cuda::complex<double> >::kernelTime() << " sec" << endl;
#endif

	// get data structures to construct fitResult
	vector<complex<double> > prodAmps;                // production amplitudes
	vector<string>           prodAmpNames;            // names of production amplitudes used in fit
	vector<pair<int,int> >   fitParCovMatrixIndices;  // indices of fit parameters for real and imaginary part in covariance matrix matrix
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

	fitResult* result = new fitResult();
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

	return fitResultPtr(result);
}
