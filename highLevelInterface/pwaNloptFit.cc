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


namespace {
	/***
	 * Add fixed parameters to the given list of fitter parameters
	 * @param pars vector of fitter parameters to which the fixed parameters will be added
	 * @param L pwaLikelihood used to determine the parameter mapping
	 * @return parameters including the fixes ones
	 */
	std::vector<double>
	addFixedParams(const double* pars, const rpwa::pwaLikelihood<std::complex<double>>& L);

	/***
	 * Add fixed parameters to the given list of fitter parameters
	 * @param pars vector of fitter parameters to which the fixed parameters will be added
	 * @param L pwaLikelihood used to determine the parameter mapping
	 * @param parsIncludingFixed output list of parameters including the fixed ones
	 */
	void
	addFixedParams(const double* pars, const rpwa::pwaLikelihood<std::complex<double>>& L, double* parsIncludingFixed);
}

using namespace std;
using namespace rpwa;


double
rpwaNloptFunc(unsigned int nmbParsNlopt, const double* parsNlopt, double* gradientNlopt, void* func_data)
{
	const pwaLikelihood<complex<double> >* L = (const pwaLikelihood<complex<double> >*)func_data;
	if(nmbParsNlopt != L->nmbPars() - L->nmbParsFixed()) {
		printErr << "parameter mismatch between NLopt and pwaLikelihood:"
		         << "Likelihood needs "  << L->nmbPars() << " parameters including " << L->nmbParsFixed() << " fixed parameters. Got " << nmbParsNlopt << " parameters."
		         << "Aborting..." << endl;
		throw;
	}
	if (L->nmbParsFixed() != 0) {
		static double* parsIncludingFixed = new double[L->nmbPars()];
		addFixedParams(parsNlopt, *L, parsIncludingFixed);
		parsNlopt = parsIncludingFixed;
	}
	double likeli;
	if (gradientNlopt) {
		static double* gradientLocal = new double[L->nmbPars()]; // local gradient buffer used in case of fixed parameters
		// the likelihood gradient is the nlopt gradient if there are no fixed parameters
		// use the local gradient buffer otherwise
		double* gradientLikelihood = (L->nmbParsFixed() == 0) ? gradientNlopt : gradientLocal;
		L->FdF(parsNlopt, likeli, gradientLikelihood);
		if (L->nmbParsFixed() != 0) { // copy gradient form local gradient buffer to nlopt gradient in case of fixed parameters
			unsigned int iNlopt = 0;
			for (unsigned int iRootPWA = 0; iRootPWA < L->nmbPars(); ++iRootPWA) {
				if (not L->parameter(iRootPWA).fixed()) {
					gradientNlopt[iNlopt] = gradientLikelihood[iRootPWA];
					++iNlopt;
				}
			}
		}

	} else {
		likeli = L->DoEval(parsNlopt);
	}
	return likeli;
}


fitResultPtr
rpwa::hli::pwaNloptFit(const pwaLikelihood<complex<double> >& L,
                       const multibinBoundariesType&          multibinBoundaries,
                       const unsigned int                     seed,
                       const string&                          startValFileName,
                       const bool                             checkHessian,
                       const bool                             saveSpace,
                       const bool                             verbose)
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
	const bool         quiet                 = not verbose;

	// report parameters
	printInfo << "running pwaNloptFit with the following parameters:" << endl;
	for (const auto& bin: multibinBoundaries) {
		char prevFill = std::cout.fill('.');
		cout << "    " << bin.first << " bin " << std::setw((bin.first.length() < 45) ? (45 - bin.first.length()) : 0) << " ["
		     << bin.second.first << ", " << bin.second.second << "]" << endl;
		std::cout.fill(prevFill);
	}
	cout << "    seed for random start values ................... "  << seed                    << endl
	     << "    path to file with start values ................. '" << startValFileName << "'" << endl;
	if (useFixedStartValues)
		cout << "    using fixed instead of random start values ..... " << defaultStartValue << endl;
	cout << "    check analytical Hessian eigenvalues............ "  << yesNo(checkHessian)     << endl
	     << "    minimizer tolerance ............................ "  << minimizerTolerance      << endl
	     << "    likelihood tolerance ........................... "  << likelihoodTolerance     << endl
	     << "    saving integral and covariance matrices......... "  << yesNo(not saveSpace)    << endl
	     << "    quiet .......................................... "  << yesNo(quiet)            << endl;

	// ---------------------------------------------------------------------------
	// setup likelihood function
	if (not quiet) {
		printInfo << "likelihood initialized with the following parameters:" << endl;
		cout << L << endl;

		printInfo << "using prior: ";
		switch (L.priorType()) {
			case pwaLikelihood<complex<double> >::FLAT:
				cout << "flat" << endl;
				break;
			case pwaLikelihood<complex<double> >::HALF_CAUCHY:
				cout      << "half-cauchy" << endl;
				printInfo << "cauchy width: " << L.cauchyWidth() << endl;
				break;
		}
	}

	const unsigned int nmbPars      = L.NDim();
	const unsigned int nmbParsFixed = L.nmbParsFixed();
	const unsigned int nmbEvts      = L.nmbEvents();

	// ---------------------------------------------------------------------------
	// setup minimizer
	nlopt::opt optimizer(nlopt::LD_LBFGS, nmbPars - nmbParsFixed);
	optimizer.set_min_objective(&rpwaNloptFunc, &const_cast<pwaLikelihood<complex<double> >&>(L));

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
	// read in fitResult with start values
	if (startValFileName.length() > 0) {
		printErr << "cannot use start values from fitResult with pwaNloptFit. Aborting..." << std::endl;
		throw;
	}

	// ---------------------------------------------------------------------------
	// set start parameter values
	printInfo << "setting start values for " << nmbPars << " parameters" << endl
	          << "    parameter naming scheme is: V[rank index]_[wave name]" << endl;
	unsigned int maxParNameLength = 0;       // maximum length of parameter names
	vector<double> params;
	params.reserve(nmbPars - nmbParsFixed);
	{
		for (unsigned int i = 0; i < nmbPars; ++i) {
			const string parName = L.parameter(i).parName();
			if (parName.length() > maxParNameLength)
				maxParNameLength = parName.length();
		}
		// use local instance of random number generator so that other
		// code has no chance of tampering with gRandom and thus cannot
		// affect the reproducability of the start values
		TRandom3     random(seed);
		const double sqrtNmbEvts = sqrt((double)nmbEvts);
		for (unsigned int i = 0; i < nmbPars; ++i) {
			const string parName = L.parameter(i).parName();

			double startVal;
			{
				startVal = (useFixedStartValues) ? defaultStartValue : random.Uniform(defaultStartValue, sqrtNmbEvts);
				if (random.Rndm() > 0.5) {
					startVal *= -1.;
				}
			}

			// check if parameter needs to be fixed
			if (not L.parameter(i).fixed()) {
				if (not quiet)
					cout << "    setting parameter [" << setw(3) << i << "] "
					     << setw(maxParNameLength) << parName << " = " << maxPrecisionAlign(startVal) << endl;
				params.push_back(startVal);
			} else if (not quiet){
				cout << "    fixing parameter  [" << setw(3) << i << "] "
				     << setw(maxParNameLength) << parName << " = 0" << endl;
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
		if (result > 0) {
			converged = true;
		}
		timer.Stop();

		// add fixed parameters to params
		params = addFixedParams(params.data(), L);

		correctParams = L.CorrectParamSigns(params.data());
		const double newLikelihood = L.DoEval(correctParams.data());
		if (likeli != newLikelihood) {
			printErr << "Flipping signs according to sign conventions changed the likelihood (from " << maxPrecisionAlign(likeli) << " to " << maxPrecisionAlign(newLikelihood) << ")." << endl;
			throw;
		} else {
			printInfo << "Likelihood unchanged at " << maxPrecisionAlign(newLikelihood) << " by flipping signs according to conventions." << endl;
		}

		if (checkHessian or not saveSpace) {
			// analytically calculate Hessian
			TMatrixT<double> hessian = L.Hessian(correctParams.data());
			// calculate and check Hessian eigenvalues
			vector<pair<TVectorT<double>, double> > eigenVectors = L.HessianEigenVectors(hessian);
			if (not quiet) {
				printInfo << "eigenvalues of (analytic) Hessian:" << endl;
			}
			for(size_t i=0; i<eigenVectors.size(); ++i) {
				if (not quiet) {
					cout << "    " << maxPrecisionAlign(eigenVectors[i].second) << endl;
				}
				if (eigenVectors[i].second <= 0.) {
					printWarn << "eigenvalue " << i << " of (analytic) Hessian is not positive (" << maxPrecisionAlign(eigenVectors[i].second) << ")." << endl;
					converged = false;
				}
			}
			if (not saveSpace) {
				fitParCovMatrix.ResizeTo(nmbPars, nmbPars);
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
		switch (result) {
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
	if (not quiet) {
		printInfo << "minimization result:" << endl;
		for (unsigned int i = 0; i < nmbPars; ++i) {
			cout << "    parameter [" << setw(3) << i << "] "
			     << setw(maxParNameLength) << L.parameter(i).parName() << " = "
			     << setw(12) << maxPrecisionAlign(correctParams[i]) << " +- ";
			if(not saveSpace) {
				cout << setw(12) << maxPrecisionAlign(sqrt(fitParCovMatrix(i, i)));
			} else {
				cout << setw(12) << "[not available]";
			}
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
	const unsigned int nmbWaves = L.nmbWaves() + 1;   // flat wave is not included in L.nmbWaves()
	vector<complex<double> > prodAmps;                // production amplitudes
	vector<string>           prodAmpNames;            // names of production amplitudes used in fit
	vector<pair<int,int> >   fitParCovMatrixIndices;  // indices of fit parameters for real and imaginary part in covariance matrix matrix
	L.buildProdAmpArrays(correctParams.data(), prodAmps, fitParCovMatrixIndices, prodAmpNames, true);
	complexMatrix normIntegral(0, 0);                 // normalization integral over full phase space without acceptance
	complexMatrix accIntegral (0, 0);                 // normalization integral over full phase space with acceptance
	vector<double> phaseSpaceIntegral;
	if (not saveSpace) {
		L.getIntegralMatrices(normIntegral, accIntegral, phaseSpaceIntegral, true);
	}
	const int normNmbEvents = (L.normalizedAmpsUsed()) ? 1 : L.nmbEvents();  // number of events to normalize to

	cout << "filling fitResult:" << endl
	     << "    number of fit parameters ............... " << nmbPars                       << endl
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
	             multibinBoundaries,
	             likeli,
	             L.rank(),
	             prodAmps,
	             prodAmpNames,
	             (saveSpace) ? nullptr : &fitParCovMatrix,
	             fitParCovMatrixIndices,
	             (saveSpace) ? nullptr : &normIntegral,
	             (saveSpace) ? nullptr : &accIntegral,
	             (saveSpace) ? nullptr : &phaseSpaceIntegral,  // contains the sqrt of the integral matrix diagonal elements!!!
	             converged,
	             hasHessian);
	return fitResultPtr(result);
}


namespace {
std::vector<double>
addFixedParams(const double* pars, const rpwa::pwaLikelihood<std::complex<double>>& L) {
	vector<double> paramsIncludingFixed(L.nmbPars());
	addFixedParams(pars, L, paramsIncludingFixed.data());
	return paramsIncludingFixed;
}

void
addFixedParams(const double* pars, const rpwa::pwaLikelihood<std::complex<double>>& L, double* parsIncludingFixed) {
	unsigned int iNlopt = 0;
	for (unsigned int iRootPWA = 0; iRootPWA < L.nmbPars(); ++iRootPWA) {
		if (L.parameter(iRootPWA).fixed()) {
			parsIncludingFixed[iRootPWA] = 0.0;
		} else {
			parsIncludingFixed[iRootPWA] = pars[iNlopt];
			++iNlopt;
		}
	}
}
}
