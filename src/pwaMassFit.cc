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
#include <Minuit2/Minuit2Minimizer.h>
#include <TFile.h>
#include <TROOT.h>
#include <TStopwatch.h>

#include "conversionUtils.hpp"
#include "fileUtils.hpp"
#include "libConfigUtils.hpp"
#include "massDepFit.h"
#include "massDepFitComponents.h"
#include "massDepFitFsmd.h"
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
	     << " [-o outfile -M minimizer -m algorithm -g # -t # -P -R -A -B -C -d -q -h] config file" << endl
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
	     << "        -B         use branchings (reducing number of couplings)" << endl
	     << "        -C         fit to spin-density matrix:   switch OFF covariances between real and imag part" << endl
	     << "                   fit to production amplitudes: switch OFF covariances between amplitudes" << endl
	     << "        -d         additional debug output (default: false)" << endl
	     << "        -q         run quietly (default: false)" << endl
	     << "        -h         print help" << endl
	     << endl;
	exit(errCode);
}


// changes status of variables (fixed/released)
// * couplings are always free
// * additional parameters can be freed with freeParameters
// * branchings also have to be freed explicitely
// * fixed values from config remain fixed
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
				           << "__";
				if(compset.useBranchings() && comp->getNrChannels() > 1) {
					const std::string waveQN = channel.getWaveName().substr(0, 7);
					prefixName << waveQN;
				} else {
					prefixName << channel.getWaveName();
				}

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

	// second eventually add all branchings
	if(compset.useBranchings()) {
		for(size_t idxComponent=0; idxComponent<compset.getNrComponents(); ++idxComponent) {
			const rpwa::massDepFit::component* comp = compset.getComponent(idxComponent);
			for(size_t idxChannel=0; idxChannel<comp->getNrChannels(); ++idxChannel) {
				// if branchings are used, not every channel has its own coupling
				if(idxChannel != comp->getChannelIdxBranching(idxChannel) || comp->getNrChannels() <= 1) {
					continue;
				}

				const rpwa::massDepFit::channel& channel = comp->getChannel(idxChannel);
				const std::string waveDecay = channel.getWaveName().substr(7);
				ostringstream prefixName;
				prefixName << "branching__"
				           << comp->getName()
				           << "__"
				           << waveDecay;

				bool free = false;
				if(find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "*")!=tokenizeFreeParameters.end()
				   || find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "branching")!=tokenizeFreeParameters.end() ) {
					free = true;
				}
				bool fix = not free;

				if(idxChannel == 0) {
					printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__real") << "') fixed to " << par[parcount] << endl;
					minimizer->SetFixedVariable(parcount,
					                            prefixName.str() + "__real",
					                            par[parcount]);
					++parcount;
				} else if (fix) {
					printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__real") << "') fixed to " << par[parcount] << endl;
					minimizer->SetFixedVariable(parcount,
					                            prefixName.str() + "__real",
					                            par[parcount]);
					++parcount;

					printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__imag") << "') fixed to " << par[parcount] << endl;
					minimizer->SetFixedVariable(parcount,
					                            prefixName.str() + "__imag",
					                            par[parcount]);
					++parcount;
				} else {
					printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__real") << "') set to " << par[parcount] << endl;
					minimizer->SetVariable(parcount,
					                       prefixName.str() + "__real",
					                       par[parcount],
					                       0.1);
					++parcount;

					printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__imag") << "') set to " << par[parcount] << endl;
					minimizer->SetVariable(parcount,
					                       prefixName.str() + "__imag",
					                       par[parcount],
					                       0.1);
					++parcount;
				}
			} // end loop over channels
		} // end loop over components
	}

	// third add parameters of the components, i.e. mass and width
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

	// set parameters for final-state mass-dependence
	if(compset.getFsmd() != NULL) {
		const rpwa::massDepFit::fsmd* fsmd = compset.getFsmd();
		for(size_t idxParameter=0; idxParameter<fsmd->getNrParameters(); ++idxParameter) {
			std::ostringstream name;
			name << "PSP__" << idxParameter;

			const bool fix = fsmd->getParameterFixed(idxParameter);

			if(fix) {
				printInfo << "parameter " << parcount << " ('" << name.str() << "') fixed to " << par[parcount] << endl;
				minimizer->SetFixedVariable(parcount,
				                            name.str(),
				                            par[parcount]);
			} else if(fsmd->getParameterLimitedLower(idxParameter) && fsmd->getParameterLimitedUpper(idxParameter)) {
				printInfo << "parameter " << parcount << " ('" << name.str() << "') set to " << par[parcount]
				          << " (limited between " << fsmd->getParameterLimitLower(idxParameter)
				          << " and " << fsmd->getParameterLimitUpper(idxParameter) << ")" << endl;
				minimizer->SetLimitedVariable(parcount,
				                              name.str(),
				                              par[parcount],
				                              fsmd->getParameterStep(idxParameter),
				                              fsmd->getParameterLimitLower(idxParameter),
				                              fsmd->getParameterLimitUpper(idxParameter));
			} else if(fsmd->getParameterLimitedLower(idxParameter)) {
				printInfo << "parameter " << parcount << " ('" << name.str() << "') set to " << par[parcount]
				          << " (limited larger than " << fsmd->getParameterLimitLower(idxParameter) << ")" << endl;
				minimizer->SetLowerLimitedVariable(parcount,
				                                   name.str(),
				                                   par[parcount],
				                                   fsmd->getParameterStep(idxParameter),
				                                   fsmd->getParameterLimitLower(idxParameter));
			} else if(fsmd->getParameterLimitedUpper(idxParameter)) {
				printInfo << "parameter " << parcount << " ('" << name.str() << "') set to " << par[parcount]
				          << " (limited smaller than " << fsmd->getParameterLimitUpper(idxParameter) << ")" << endl;
				minimizer->SetUpperLimitedVariable(parcount,
				                                   name.str(),
				                                   par[parcount],
				                                   fsmd->getParameterStep(idxParameter),
				                                   fsmd->getParameterLimitUpper(idxParameter));
			} else {
				printInfo << "parameter " << parcount << " ('" << name.str() << "') set to " << par[parcount] << endl;
				minimizer->SetVariable(parcount,
				                       name.str(),
				                       par[parcount],
				                       fsmd->getParameterStep(idxParameter));
			}
			++parcount;
		}
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

	// force loading predefined std::complex dictionary
	// see http://root.cern.ch/phpBB3/viewtopic.php?f=5&t=9618&p=50164
	gROOT->ProcessLine("#include <complex>");

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
	bool         doBranching        = false;
	bool         doCov              = true;
	bool         debug              = false;
	bool         quiet              = false;
	extern char* optarg;
	extern int   optind;
	int c;
	while ((c = getopt(argc, argv, "o:M:m:g:t:PRABCdqh")) != -1) {
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
		case 'B':
			doBranching = true;
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
	          << "    use branchings ................................. "  << yesNo(doBranching) << endl
	          << "    take covariance into account ................... "  << yesNo(doCov) << endl
	          << "    debug .......................................... "  << yesNo(debug) << endl
	          << "    quiet .......................................... "  << yesNo(quiet) << endl;

	rpwa::massDepFit::massDepFit mdepFit;
	mdepFit.setDebug(debug);

	rpwa::massDepFit::model compset;
	compset.useBranchings(doBranching);

	rpwa::massDepFit::likelihood L;
	L.fitProductionAmplitudes(doProdAmp);
	L.useCovariance(doCov);

	libconfig::Config configFile;
	if(not parseLibConfigFile(configFileName, configFile, debug)) {
		printErr << "could not read configuration file '" << configFileName << "'." << endl;
		return 1;
	}
	libconfig::Setting& configRoot = configFile.getRoot();

	// read configuration file
	if(not mdepFit.readConfig(&configRoot, compset, valTreeName, valBranchName)) {
		printErr << "error while reading configuration file '" << configFileName << "'." << endl;
		return 1;
	}

	// set-up fit model and likelihood
	if(not mdepFit.init(compset, L)) {
		printErr << "error while reading configuration file '" << configFileName << "'." << endl;
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

		// special for Minuit2
		if(dynamic_cast<ROOT::Minuit2::Minuit2Minimizer*>(minimizer.get())) {
#if ROOT_VERSION_CODE >= ROOT_VERSION(5, 34, 19)
			((ROOT::Minuit2::Minuit2Minimizer*)minimizer.get())->SetStorageLevel(0);
#endif
		}

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
