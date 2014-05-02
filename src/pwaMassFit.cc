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
//      fitting program for massdependent fit rootpwa
//      minimizes rpwa::massDepFit::likelihood function
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <boost/tokenizer.hpp>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <TFile.h>
#include <TStopwatch.h>

#include "conversionUtils.hpp"
#include "fileUtils.hpp"
#include "libConfigUtils.hpp"
#include "massDepFit.h"
#include "massDepFitComponents.h"
#include "massDepFitLikeli.h"
#include "massDepFitModel.h"
#include "reportingUtils.hpp"
#include "reportingUtilsEnvironment.h"


using namespace std;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{
	cerr << "performs mass-dependent fit" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " [-o outfile -M minimizer -m algorithm -t # -P -R -C -d -q -h] config file" << endl
	     << "    where:" << endl
	     << "        -o file    path to output file (default: 'mDep.result.root')" << endl
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
	     << "        -P         plotting only - no fit" << endl
	     << "        -R         plot in fit range only" << endl
	     << "        -A         fit to the production amplitudes (default: spin-density matrix)" << endl
	     << "        -C         fit to spin-density matrix:   switch OFF covariances between real and imag part" << endl
	     << "                   fit to production amplitudes: switch OFF covariances between amplitudes" << endl
	     << "        -d         additional debug output (default: false)" << endl
	     << "        -q         run quietly (default: false)" << endl
	     << "        -h         print help" << endl
	     << endl;
	exit(errCode);
}


// changes status of variables (fixed/released)
// * fixed values from config remain
// * additional parameters can be freed with freeParameters
// * couplings are always free
bool
releasePars(ROOT::Math::Minimizer* minimizer,
            const rpwa::massDepFit::model& compset,
            const string& freeParameters)
{
	// save current state
	const size_t npar = compset.getNrParameters();
	double par[npar];

	if(minimizer->NFree() == 0) {
		compset.getParameters(par);
	} else {
		for(size_t i=0; i<npar; ++i) {
			par[i]=minimizer->X()[i];
		}
	}

	// tokenize freeParameters string (default deparators also include '*')
	boost::char_separator<char> separators(" ,\t\n");
	boost::tokenizer<boost::char_separator<char> > tokenizeFreeParameters(freeParameters, separators);

	// reset minimizer
	minimizer->Clear();

	size_t parcount=0;
	// first add all couplings
	for(size_t idxComponent=0; idxComponent<compset.getNrComponents(); ++idxComponent) {
		const rpwa::massDepFit::component* comp = compset.getComponent(idxComponent);
		for(size_t idxChannel=0; idxChannel<comp->getNrChannels(); ++idxChannel) {
			// if branchings are used, not every channel has its own coupling
			if(idxChannel != comp->getChannelIdxCoupling(idxChannel)) {
				continue;
			}

			const rpwa::massDepFit::channel& channel = comp->getChannel(idxChannel);
			for(size_t idxBin=0; idxBin<channel.getNrBins(); ++idxBin) {
				ostringstream prefixName;
				prefixName << "coupling__bin"
				           << idxBin
				           << "__"
				           << comp->getName()
				           << "__"
				           << channel.getWaveName();

				printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__real") << "') set to " << par[parcount] << endl;
				minimizer->SetVariable(parcount,
				                       prefixName.str() + "__real",
				                       par[parcount],
				                       0.1);
				++parcount;

				if(not channel.isAnchor()) {
					printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__imag") << "') set to " << par[parcount] << endl;
					minimizer->SetVariable(parcount,
					                       prefixName.str() + "__imag",
					                       par[parcount],
					                       0.1);
					++parcount;
				}
			}
		} // end loop over channels
	} // end loop over components

	// second add parameters of the components, i.e. mass and width
	for(size_t idxComponent=0; idxComponent<compset.getNrComponents(); ++idxComponent) {
		const rpwa::massDepFit::component* comp = compset.getComponent(idxComponent);
		for(size_t idxParameter=0; idxParameter<comp->getNrParameters(); ++idxParameter) {
			const string name = comp->getName() + "__" + comp->getParameterName(idxParameter);

			bool free = false;
			if(find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "*")!=tokenizeFreeParameters.end()
			   || find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), comp->getName())!=tokenizeFreeParameters.end()
			   || find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), comp->getParameterName(idxParameter))!=tokenizeFreeParameters.end()
			   || find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), comp->getName() + "__" + comp->getParameterName(idxParameter))!=tokenizeFreeParameters.end()) {
				free = true;
			}
			bool fix = not free;
			if(comp->getParameterFixed(idxParameter)) {
				fix = true;
			}

			if(fix) {
				printInfo << "parameter " << parcount << " ('" << name << "') fixed to " << par[parcount] << endl;
				minimizer->SetFixedVariable(parcount,
				                            name,
				                            par[parcount]);
			} else if(comp->getParameterLimitedLower(idxParameter) && comp->getParameterLimitedUpper(idxParameter)) {
				printInfo << "parameter " << parcount << " ('" << name << "') set to " << par[parcount]
				          << " (limited between " << comp->getParameterLimitLower(idxParameter)
				          << " and " << comp->getParameterLimitUpper(idxParameter) << ")" << endl;
				minimizer->SetLimitedVariable(parcount,
				                              name,
				                              par[parcount],
				                              comp->getParameterStep(idxParameter),
				                              comp->getParameterLimitLower(idxParameter),
				                              comp->getParameterLimitUpper(idxParameter));
			} else if(comp->getParameterLimitedLower(idxParameter)) {
				printInfo << "parameter " << parcount << " ('" << name << "') set to " << par[parcount]
				          << " (limited larger than " << comp->getParameterLimitLower(idxParameter) << ")" << endl;
				minimizer->SetLowerLimitedVariable(parcount,
				                                   name,
				                                   par[parcount],
				                                   comp->getParameterStep(idxParameter),
				                                   comp->getParameterLimitLower(idxParameter));
			} else if(comp->getParameterLimitedUpper(idxParameter)) {
				printInfo << "parameter " << parcount << " ('" << name << "') set to " << par[parcount]
				          << " (limited smaller than " << comp->getParameterLimitUpper(idxParameter) << ")" << endl;
				minimizer->SetUpperLimitedVariable(parcount,
				                                   name,
				                                   par[parcount],
				                                   comp->getParameterStep(idxParameter),
				                                   comp->getParameterLimitUpper(idxParameter));
			} else {
				printInfo << "parameter " << parcount << " ('" << name << "') set to " << par[parcount] << endl;
				minimizer->SetVariable(parcount,
				                       name,
				                       par[parcount],
				                       comp->getParameterStep(idxParameter));
			}
			++parcount;
		}
	} // end loop over components

	// set phase space
	const unsigned int nfreeFsmd=compset.getFsmdNrParameters();
	for(unsigned int ifreeFsmd=0; ifreeFsmd<nfreeFsmd; ++ifreeFsmd) {
		const double val = par[parcount];
		double lower,upper;
		compset.getFsmdParameterLimits(ifreeFsmd, lower, upper);
		ostringstream name;
		name << "PSP__" << ifreeFsmd;
		printInfo << "parameter " << parcount << " ('" << name.str() << "') set to " << val << endl;
		minimizer->SetLimitedVariable(parcount,
		                              name.str().c_str(),
		                              val,
		                              0.0001,
		                              lower,
		                              upper);
		++parcount;
	}

	return true;
}

int
main(int    argc,
     char** argv)
{
	printCompilerInfo();
	printLibraryInfo ();
	printGitHash     ();
	cout << endl;

	// --------------------------------------------------------------------------
	// internal parameters
	const string       valTreeName           = "pwa";
	const string       valBranchName         = "fitResult_v2";
	const unsigned int maxNmbOfIterations    = 20000;
	const unsigned int maxNmbOfFunctionCalls = 2000000;
	const bool         runHesse              = true;
	const bool         runMinos              = false;

	// ---------------------------------------------------------------------------
	// parse command line options
	const string progName           = argv[0];
	string       outFileName        = "mDep.result.root";     // output filename
	string       minimizerType[2]   = {"Minuit2", "Migrad"};  // minimizer, minimization algorithm
	int          minimizerStrategy  = 1;                      // minimizer strategy
	double       minimizerTolerance = 1e-10;                  // minimizer tolerance
	bool         onlyPlotting       = false;
	bool         rangePlotting      = false;
	bool         doProdAmp          = false;
	bool         doCov              = true;
	bool         debug              = false;
	bool         quiet              = false;
	extern char* optarg;
	extern int   optind;
	int c;
	while ((c = getopt(argc, argv, "o:M:m:g:t:PRACdqh")) != -1) {
		switch (c) {
		case 'o':
			outFileName = optarg;
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
		case 'P':
			onlyPlotting = true;
			break;
		case 'R':
			rangePlotting = true;
			break;
		case 'A':
			doProdAmp = true;
			break;
		case 'C':
			doCov = false;
			break;
		case 'd':
			debug = true;
			break;
		case 'q':
			quiet = true;
			break;
		case '?':
		case 'h':
			usage(progName, 1);
			break;
		}
	}

	// there must only be one remaining (unhandled) argument which is the
	// configuration file
	if(optind+1 != argc) {
		printErr << "you need to specify exactly one configuration file." << endl;
		usage(progName, 1);
	}
	const string configFileName = argv[optind];

	// report parameters
	printInfo << "running " << progName << " with the following parameters:" << endl
	          << "    path to configuration file ..................... '" << configFileName << "'" << endl
	          << "    path to output file ............................ '" << outFileName << "'" << endl
	          << "    minimizer ...................................... "  << minimizerType[0] << ", " << minimizerType[1] << endl
	          << "    minimizer strategy ............................. "  << minimizerStrategy  << endl
	          << "    minimizer tolerance ............................ "  << minimizerTolerance << endl
	          << "    only plotting .................................. "  << yesNo(onlyPlotting) << endl
	          << "    plot in fit range only ......................... "  << yesNo(rangePlotting) << endl
	          << "    fit to production amplitudes ................... "  << yesNo(doProdAmp) << endl
	          << "    take covariance into account ................... "  << yesNo(doCov) << endl
	          << "    debug .......................................... "  << yesNo(debug) << endl
	          << "    quiet .......................................... "  << yesNo(quiet) << endl;

	rpwa::massDepFit::massDepFit mdepFit;
	mdepFit.setDebug(debug);

	libconfig::Config configFile;
	if(not parseLibConfigFile(configFileName, configFile, debug)) {
		printErr << "could not read configuration file '" << configFileName << "'." << endl;
		return 1;
	}
	libconfig::Setting& configRoot = configFile.getRoot();

	// input section
	const libconfig::Setting* configInput = findLibConfigGroup(configRoot, "input");
	if(not configInput) {
		printErr << "'input' section in configuration file does not exist." << endl;
		return 1;
	}
	if(not mdepFit.readConfigInput(configInput)) {
		printErr << "error while reading 'input' section from configuration file." << endl;
		return 1;
	}

	// extract information from fit results
	if(not mdepFit.readInFiles(valTreeName, valBranchName)) {
		printErr << "error while trying to read fit result." << endl;
		return 1;
	}

	// extract information for systematic errors
	if(not mdepFit.readSystematicsFiles(valTreeName, valBranchName)) {
		printErr << "error while trying to read fit results for systematic errors." << endl;
		return 1;
	}

	// prepare mass limits
	if(not mdepFit.prepareMassLimits()) {
		printErr << "error determine which bins to use in the fit." << endl;
		return 1;
	}

	// set-up fit model (resonances, background, final-state mass dependence
	rpwa::massDepFit::model compset;
	if(not mdepFit.readConfigModel(&configRoot, compset)) {
		printErr << "error while reading fit model from configuration file." << endl;
		return 1;
	}

	if(not compset.init(mdepFit.getWaveNames(), mdepFit.getMassBinCenters(), mdepFit.getAnchorWaveName(), mdepFit.getAnchorComponentName())) {
		printErr << "error while initializing the fit model." << endl;
		return 1;
	}

	// set-up likelihood
	rpwa::massDepFit::likelihood L;
	if(not L.init(&compset,
	              mdepFit.getMassBinCenters(),
	              mdepFit.getInProductionAmplitudes(),
	              mdepFit.getInProductionAmplitudesCovariance(),
	              mdepFit.getInSpinDensityMatrices(),
	              mdepFit.getInSpinDensityCovarianceMatrices(),
	              mdepFit.getWavePairMassBinLimits(),
	              doProdAmp,
	              doCov)) {
		printErr << "error while initializing the likelihood calculator." << endl;
		return 1;
	}

	if(onlyPlotting) {
		printInfo << "plotting only mode, skipping minimzation." << endl;

		// still calculate a chi2 from the values in the configuration file
		const size_t npar = compset.getNrParameters();
		double par[npar];
		compset.getParameters(par);

		printInfo << "chi2 (valid only if fit was successful) = " << maxPrecisionAlign(L.DoEval(par)) << endl;
	} else {
		// setup minimizer
		printInfo << "creating and setting up minimizer '" << minimizerType[0] << "' "
		          << "using algorithm '" << minimizerType[1] << "'" << endl;
		std::auto_ptr<ROOT::Math::Minimizer> minimizer(ROOT::Math::Factory::CreateMinimizer(minimizerType[0], minimizerType[1]));
		if(minimizer.get() == NULL) {
			printErr << "could not create minimizer. exiting." << endl;
			return 1;
		}
		minimizer->SetFunction        (L);
		minimizer->SetStrategy        (minimizerStrategy);
		minimizer->SetTolerance       (minimizerTolerance);
		minimizer->SetPrintLevel      ((quiet) ? 0 : 3);
		minimizer->SetMaxIterations   (maxNmbOfIterations);
		minimizer->SetMaxFunctionCalls(maxNmbOfFunctionCalls);

		// keep list of parameters to free
		const vector<string> freeParameters = mdepFit.getFreeParameters();
		const size_t nrSteps = freeParameters.size();

		bool success = true;
		for(size_t step=0; step<nrSteps; ++step) {
			// set startvalues
			if(not releasePars(minimizer.get(), compset, freeParameters[step])) {
				printErr << "error while setting start parameters for step " << step << "." << endl;
				return 1;
			}

			// keep track of the time spend in the fit
			TStopwatch stopwatch;

			printInfo << "performing minimization step " << step << ": '" << freeParameters[step] << "' (" << minimizer->NFree() << " free parameters)." << endl;
			stopwatch.Start();
			success &= minimizer->Minimize();
			stopwatch.Stop();

			if(not success) {
				printWarn << "minimization failed." << endl;
			} else {
				printInfo << "minimization successful." << endl;
			}
			printInfo << "minimization took " <<  maxPrecisionAlign(stopwatch.CpuTime()) << " s" << endl;
		}

		if(runHesse) {
			TStopwatch stopwatch;

			printInfo << "calculating Hessian matrix." << endl;
			stopwatch.Start();
			success &= minimizer->Hesse();
			stopwatch.Stop();

			if(not success) {
				printWarn << "calculation of Hessian matrix failed." << endl;
			} else {
				printInfo << "calculation of Hessian matrix successful." << endl;
			}
			printInfo << "calculating Hessian matrix took " <<  maxPrecisionAlign(stopwatch.CpuTime()) << " s" << endl;
		}

		printInfo << "minimizer status summary:" << endl
		          << "    total number of parameters .......................... " << minimizer->NDim()             << endl
		          << "    number of free parameters ........................... " << minimizer->NFree()            << endl
		          << "    maximum allowed number of iterations ................ " << minimizer->MaxIterations()    << endl
		          << "    maximum allowed number of function calls ............ " << minimizer->MaxFunctionCalls() << endl
		          << "    minimizer status .................................... " << minimizer->Status()           << endl
		          << "    minimizer provides error and error matrix ........... " << minimizer->ProvidesError()    << endl
		          << "    minimizer has performed detailed error validation ... " << minimizer->IsValidError()     << endl
		          << "    estimated distance to minimum ....................... " << minimizer->Edm()              << endl
		          << "    statistical scale used for error calculation ........ " << minimizer->ErrorDef()         << endl
		          << "    minimizer strategy .................................. " << minimizer->Strategy()         << endl
		          << "    absolute tolerance .................................. " << minimizer->Tolerance()        << endl;

		// print results
		ostringstream output;
		const unsigned int nmbPar  = L.NDim();
		for(unsigned int i = 0; i<nmbPar; ++i) {
			output << "    parameter [" << setw(3) << i << "] ";
			output << minimizer->VariableName(i) << " " ;
			output << maxPrecisionAlign(minimizer->X()[i]) << " +- " << maxPrecisionAlign(minimizer->Errors()[i]);

			if(runMinos) {  // does not work for all parameters
				double minosErrLow, minosErrUp;
				if(minimizer->GetMinosError(i, minosErrLow, minosErrUp)) {
					output << "    Minos: " << "[" << minosErrLow << ", +" << minosErrUp << "]";
				}
			}
			output << endl;
		}
		printInfo << "minimization result:" << endl
		          << output.str();

		const double* par=minimizer->X();
		double chi2 = 0.;
		if(success) {
			chi2 = L.DoEval(par);
		} else {
			printInfo << "chi2 (if fit were successful) = " << maxPrecisionAlign(L.DoEval(par)) << endl;
		}
		printInfo << "chi2 = " << maxPrecisionAlign(chi2) << endl;
		compset.setParameters(par);

		const unsigned int nrDataPoints = L.NDataPoints();
		const unsigned int nrFree = minimizer->NFree();
		const unsigned int ndf = nrDataPoints - nrFree;
		printInfo << "ndf = " << nrDataPoints << "-" << nrFree << "=" << ndf << endl;

		double chi2red = chi2/(double)ndf;
		printInfo << "chi2/ndf = " << maxPrecisionAlign(chi2red) << endl;

		if(not mdepFit.updateConfig(&configRoot, compset, minimizer.get(), chi2, ndf, chi2red)) {
			printErr << "error while updating configuration file." << endl;
			return 1;
		}
	}

	string confFileName(outFileName);
	if(extensionFromPath(confFileName) == "root") {
		confFileName = changeFileExtension(confFileName, "conf");
	} else {
		confFileName += ".conf";
	}

	if(debug) {
		printDebug << "name of output configuration file: '" << confFileName << "'." << endl;
	}
	configFile.writeFile(confFileName.c_str());

	string rootFileName(outFileName);
	if(extensionFromPath(rootFileName) != "root") {
		rootFileName += ".root";
	}

	if(debug) {
		printDebug << "name of output ROOT file: '" << confFileName << "'." << endl;
	}
	std::auto_ptr<TFile> outFile(TFile::Open(rootFileName.c_str(), "RECREATE"));
	if(outFile.get() == NULL || outFile->IsZombie()) {
		printErr << "error while creating ROOT file '" << rootFileName << "' for plots of fit result."<< endl;
		return 1;
	}
	if(not mdepFit.createPlots(compset, outFile.get(), rangePlotting)) {
		printErr << "error while creating plots." << endl;
		return 1;
	}
	outFile->Close();

	return 0;
}
