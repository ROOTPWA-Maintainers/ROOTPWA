#include "pwaFit.h"

#include <complex>

#include <Math/Factory.h>
#include <Math/Minimizer.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <TVectorT.h>

#include <conversionUtils.hpp>
#include <partialWaveFitHelper.h>
#include <pwaLikelihood.h>
#include <reportingUtils.hpp>


using namespace std;
using namespace ROOT::Math;
using namespace rpwa;


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


namespace ROOT {
	namespace Math {
		ostream&
		operator<<(ostream& out,
		           Minimizer& minimizer)
		{
			return printMinimizerStatus(out, minimizer);
		}
	}
}


fitResultPtr
rpwa::hli::pwaFit(const pwaLikelihood<complex<double> >& L,
                  const double                           massBinMin,
                  const double                           massBinMax,
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
	const double       startValStep          = 0.0005;
	const unsigned int maxNmbOfIterations    = 20000;
	const unsigned int maxNmbOfFunctionCalls = 40000;
	const bool         runHesse              = true;
	const bool         runMinos              = false;
	const string       minimizerType[2]      = {"Minuit2", "Migrad"};  // minimizer, minimization algorithm
	const int          minimizerStrategy     = 1;                      // minimizer strategy
	const double       minimizerTolerance    = 1e-10;                  // minimizer tolerance
#if ROOT_VERSION_CODE >= ROOT_VERSION(5, 34, 19)
	const bool         saveMinimizerMemory   = true;
#endif
	const bool         quiet                 = not verbose;

	// report parameters
	printInfo << "running pwaFit with the following parameters:" << endl;
	cout << "    mass bin [" << massBinMin << ", " << massBinMax << "] GeV/c^2" << endl
	     << "    seed for random start values ................... "  << seed                    << endl
	     << "    path to file with start values ................. '" << startValFileName << "'" << endl;
	if (useFixedStartValues)
		cout << "    using fixed instead of random start values ..... " << defaultStartValue << endl;
	cout << "    check analytical Hessian eigenvalues............ "  << yesNo(checkHessian)     << endl
	     << "    minimizer ...................................... '" << minimizerType[0] << ", " << minimizerType[1] << "'" << endl
	     << "    minimizer strategy ............................. "  << minimizerStrategy       << endl
	     << "    minimizer tolerance ............................ "  << minimizerTolerance      << endl
	     << "    saving integral and covariance matrices......... "  << yesNo(not saveSpace)    << endl
	     << "    quiet .......................................... "  << yesNo(quiet)            << endl;

	// ---------------------------------------------------------------------------
	// setup likelihood function
	if (not quiet) {
		printInfo << "likelihood initialized with the following parameters:" << endl;
		cout << L << endl;

		printInfo << "using prior: ";
		switch(L.priorType()) {
			case pwaLikelihood<complex<double> >::FLAT:
				cout << "flat" << endl;
				break;
			case pwaLikelihood<complex<double> >::HALF_CAUCHY:
				cout      << "half-cauchy" << endl;
				printInfo << "cauchy width: " << L.cauchyWidth() << endl;
				break;
		}
	}

	const double massBinCenter = (massBinMin + massBinMax) / 2;
	const unsigned int nmbPar  = L.NDim();
	const unsigned int nmbEvts = L.nmbEvents();

	// ---------------------------------------------------------------------------
	// setup minimizer
	printInfo << "creating and setting up minimizer '" << minimizerType[0] << "' "
	          << "using algorithm '" << minimizerType[1] << "'" << endl;
	Minimizer* minimizer = Factory::CreateMinimizer(minimizerType[0], minimizerType[1]);
	if (not minimizer) {
		printErr << "could not create minimizer. exiting." << endl;
		throw;
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
	          << "    parameter naming scheme is: V[rank index]_[wave name]" << endl;
	unsigned int maxParNameLength = 0;       // maximum length of parameter names
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
		bool         success     = true;
		for (unsigned int i = 0; i < nmbPar; ++i) {
			const string& parName = L.parName(i);

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
	vector<double> correctParams;
	TMatrixT<double> fitParCovMatrix(0, 0);
	printInfo << "performing minimization." << endl;
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
		}
		if (converged) {
			printSucc << "minimization finished successfully. " << flush;
		} else {
			printWarn << "minimization failed. " << flush;
		}
		cout << "used " << flush;
		timer.Print();
		printInfo << *minimizer;
		if (runHesse && not saveSpace) {
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

			fitParCovMatrix.ResizeTo(nmbPar, nmbPar);
			for(unsigned int i = 0; i < nmbPar; ++i)
				for(unsigned int j = 0; j < nmbPar; ++j)
					// The factor 0.5 is needed because
					// MINUIT by default assumes a Chi2
					// function and not a loglikeli
					// (see Minuit manual!)
					// Note: SetErrorDef in ROOT does not work
					fitParCovMatrix[i][j] = 0.5 * minimizer->CovMatrix(i, j);
		}
	}

	// ---------------------------------------------------------------------------
	// print results
	printInfo << "minimization result:" << endl;
	for (unsigned int i = 0; i < nmbPar; ++i) {
		cout << "    parameter [" << setw(3) << i << "] "
		     << setw(maxParNameLength) << L.parName(i) << " = ";
		if (L.parFixed(i))
			cout << correctParams[i] << " (fixed)";
		else {
			cout << setw(12) << maxPrecisionAlign(correctParams[i]) << " +- ";
			if(runHesse && not saveSpace) {
				cout << setw(12) << maxPrecisionAlign(sqrt(fitParCovMatrix(i, i)));
			} else {
				cout << setw(12) << "[not available]";
			}
			if (runMinos) {
				double minosErrLow = 0;
				double minosErrUp  = 0;
				const bool success = minimizer->GetMinosError(i, minosErrLow, minosErrUp);
				if (success)
					cout << "    Minos: " << "[" << minosErrLow << ", +" << minosErrUp << "]";
			}
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

	return fitResultPtr(result);
}
