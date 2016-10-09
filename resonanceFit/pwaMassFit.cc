///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010-2012 Sebastian Neubert (TUM)
//    Copyright 2014,2015 Sebastian Uhl (TUM)
//
//    This file is part of ROOTPWA
//
//    ROOTPWA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ROOTPWA is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ROOTPWA.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      program to extract resonance parameters from a partial-wave fit
//
//-------------------------------------------------------------------------


#include <fstream>

#include <yaml-cpp/yaml.h>

#include <TFile.h>
#include <TROOT.h>
#include <TStopwatch.h>

#include "conversionUtils.hpp"
#include "fileUtils.hpp"
#include "massDepFit.h"
#include "massDepFitCache.h"
#include "massDepFitComponents.h"
#include "massDepFitFsmd.h"
#include "massDepFitFunction.h"
#include "massDepFitMinimizerRoot.h"
#include "massDepFitModel.h"
#include "massDepFitParameters.h"
#include "reportingUtils.hpp"
#include "reportingUtilsEnvironment.h"
#include "yamlCppUtils.hpp"


void
usage(const std::string& progName,
      const int     errCode = 0)
{
	std::cerr << "performs mass-dependent fit" << std::endl
	          << std::endl
	          << "usage:" << std::endl
	          << progName
	          << " [-o outfile -c # -M minimizer -m algorithm -g # -t # -P -R -F # -A -B -C # -d -q -h] config file" << std::endl
	          << "    where:" << std::endl
	          << "        -o file    path to output file (default: 'mDep.result.root')" << std::endl
	          << "        -c #       maximal number of function calls (default: depends on number of parameters)" << std::endl
	          << "        -M name    minimizer (default: Minuit2)" << std::endl
	          << "        -m name    minimization algorithm (optional, default: Migrad)" << std::endl
	          << "                   available minimizers: Minuit:      Migrad, Simplex, Minimize, Migrad_imp" << std::endl
	          << "                                         Minuit2:     Migrad, Simplex, Combined, Scan, Fumili" << std::endl
	          << "                                         GSLMultiMin: ConjugateFR, ConjugatePR, BFGS, BFGS2, SteepestDescent" << std::endl
	          << "                                         GSLMultiFit: -" << std::endl
	          << "                                         GSLSimAn:    -" << std::endl
	          << "                                         Linear:      Robust" << std::endl
	          << "                                         Fumili:      -" << std::endl
	          << "        -g #       minimizer strategy: 0 = low, 1 = medium, 2 = high effort  (default: 1)" << std::endl
	          << "        -t #       minimizer tolerance (default: 1e-10)" << std::endl
	          << "        -P         plotting only - no fit" << std::endl
	          << "        -R         plot in fit range only" << std::endl
	          << "        -F #       finer binning for plotting (default: 1)" << std::endl
	          << "        -A         fit to the production amplitudes (default: spin-density matrix)" << std::endl
	          << "        -B         use branchings (reducing number of couplings)" << std::endl
	          << "        -C #       part of the covariance matrix to use:" << std::endl
	          << "                       1 = only diagonal elements" << std::endl
	          << "                       2 = take covariance between real and imaginary part of the same complex number into account" << std::endl
	          << "                       3 = full covariance matrix (not available while fitting to the spin-density matrix)" << std::endl
	          << "        -d         additional debug output (default: false)" << std::endl
	          << "        -q         run quietly (default: false)" << std::endl
	          << "        -h         print help" << std::endl
	          << std::endl;
	exit(errCode);
}


int
main(int    argc,
     char** argv)
{
	rpwa::printCompilerInfo();
	rpwa::printLibraryInfo ();
	rpwa::printGitHash     ();
	std::cout << std::endl;

#if ROOT_VERSION_CODE < ROOT_VERSION(6, 0, 0)
	// force loading predefined std::complex dictionary
	// see http://root.cern.ch/phpBB3/viewtopic.php?f=5&t=9618&p=50164
	gROOT->ProcessLine("#include <complex>");
#endif

	// --------------------------------------------------------------------------
	// internal parameters
	const std::string  valTreeName           = "pwa";
	const std::string  valBranchName         = "fitResult_v2";

	// ---------------------------------------------------------------------------
	// parse command line options
	const std::string progName                 = argv[0];
	std::string       outFileName              = "mDep.result.root";     // output filename
	unsigned int      maxNumberOfFunctionCalls = 0;
	std::string       minimizerType[2]         = {"Minuit2", "Migrad"};  // minimizer, minimization algorithm
	int               minimizerStrategy        = 1;                      // minimizer strategy
	double            minimizerTolerance       = 1e-10;                  // minimizer tolerance
	bool              onlyPlotting             = false;
	bool              rangePlotting            = false;
	size_t            extraBinning             = 1;
	bool              doProdAmp                = false;
	bool              doBranching              = false;
	bool              debug                    = false;
	bool              quiet                    = false;

	rpwa::massDepFit::function::useCovarianceMatrix doCov = rpwa::massDepFit::function::useCovarianceMatrixDefault;

	extern char* optarg;
	extern int   optind;
	int c;
	while ((c = getopt(argc, argv, "o:c:M:m:g:t:PRF:ABC:dqh")) != -1) {
		switch (c) {
		case 'o':
			outFileName = optarg;
			break;
		case 'c':
			maxNumberOfFunctionCalls = atoi(optarg);
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
		case 'F':
			extraBinning = atoi(optarg);
			break;
		case 'A':
			doProdAmp = true;
			break;
		case 'B':
			doBranching = true;
			break;
		case 'C':
			{
				int cov = atoi(optarg);
				if      (cov == 1) { doCov = rpwa::massDepFit::function::useDiagnalElementsOnly;        }
				else if (cov == 2) { doCov = rpwa::massDepFit::function::useComplexDiagnalElementsOnly; }
				else if (cov == 3) { doCov = rpwa::massDepFit::function::useFullCovarianceMatrix;       }
				else               { usage(progName, 1); }
			}
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
		printErr << "you need to specify exactly one configuration file." << std::endl;
		usage(progName, 1);
	}
	const std::string configFileName = argv[optind];

	// report parameters
	printInfo << "running " << progName << " with the following parameters:" << std::endl
	          << "    path to configuration file ..................... '" << configFileName << "'" << std::endl
	          << "    path to output file ............................ '" << outFileName << "'" << std::endl
	          << "    number of calls to fit function ................ "  << maxNumberOfFunctionCalls << std::endl
	          << "    minimizer ...................................... "  << minimizerType[0] << ", " << minimizerType[1] << std::endl
	          << "    minimizer strategy ............................. "  << minimizerStrategy  << std::endl
	          << "    minimizer tolerance ............................ "  << minimizerTolerance << std::endl
	          << "    only plotting .................................. "  << rpwa::yesNo(onlyPlotting) << std::endl
	          << "    plot in fit range only ......................... "  << rpwa::yesNo(rangePlotting) << std::endl
	          << "    fit to production amplitudes ................... "  << rpwa::yesNo(doProdAmp) << std::endl
	          << "    use branchings ................................. "  << rpwa::yesNo(doBranching) << std::endl
	          << "    debug .......................................... "  << rpwa::yesNo(debug) << std::endl
	          << "    quiet .......................................... "  << rpwa::yesNo(quiet) << std::endl;

	rpwa::massDepFit::massDepFit mdepFit;
	mdepFit.setDebug(debug);

	rpwa::massDepFit::modelPtr fitModel(new rpwa::massDepFit::model);
	rpwa::massDepFit::functionPtr fitFunction(new rpwa::massDepFit::function(doProdAmp, doCov));

	YAML::Node configRoot;
	if(not rpwa::YamlCppUtils::parseYamlFile(configFileName, configRoot, debug)) {
		printErr << "could not read configuration file '" << configFileName << "'." << std::endl;
		return 1;
	}

	// read configuration file
	rpwa::massDepFit::parameters fitParameters;
	rpwa::massDepFit::parameters fitParametersError;
	std::map<std::string, double> fitQuality;
	if(not mdepFit.readConfig(configRoot, fitModel, fitParameters, fitParametersError, fitQuality, doBranching, valTreeName, valBranchName)) {
		printErr << "error while reading configuration file '" << configFileName << "'." << std::endl;
		return 1;
	}

	// set-up fit model and fit function
	if(not mdepFit.init(fitModel, fitFunction)) {
		printErr << "error while reading configuration file '" << configFileName << "'." << std::endl;
		return 1;
	}

	rpwa::massDepFit::cache cache(mdepFit.getNrWaves(),
	                              fitModel->getNrComponents()+1,           // nr components + final-state mass-dependence
	                              fitModel->getMaxChannelsInComponent(),
	                              mdepFit.getNrBins(),
	                              mdepFit.getMaxMassBins());

	if(onlyPlotting) {
		printInfo << "plotting only mode, skipping minimzation." << std::endl;

		printInfo << "chi2 (valid only if fit was successful) = " << rpwa::maxPrecisionAlign(fitFunction->chiSquare(fitParameters, cache)) << std::endl;
	} else {
		rpwa::massDepFit::minimizerRoot minimizer(fitModel,
		                                          fitFunction,
		                                          mdepFit.getFreeParameters(),
		                                          maxNumberOfFunctionCalls,
		                                          minimizerType,
		                                          minimizerStrategy,
		                                          minimizerTolerance,
		                                          quiet);

		TStopwatch stopwatch;

		stopwatch.Start();
		fitQuality = minimizer.minimize(fitParameters, fitParametersError, cache);
		stopwatch.Stop();

		printInfo << "minimization took " << rpwa::maxPrecisionAlign(stopwatch.CpuTime()) << " s" << std::endl;
	}

	std::string confFileName(outFileName);
	if(rpwa::extensionFromPath(confFileName) == "root") {
		confFileName = rpwa::changeFileExtension(confFileName, "yaml");
	} else {
		confFileName += ".yaml";
	}

	if(debug) {
		printDebug << "name of output configuration file: '" << confFileName << "'." << std::endl;
	}
	std::ofstream configFile(confFileName.c_str());
	if(not mdepFit.writeConfig(configFile, fitModel, fitParameters, fitParametersError, fitQuality)) {
		printErr << "error while writing result to configuration file." << std::endl;
		return 1;
	}
	configFile.close();

	std::string rootFileName(outFileName);
	if(rpwa::extensionFromPath(rootFileName) != "root") {
		rootFileName += ".root";
	}

	if(debug) {
		printDebug << "name of output ROOT file: '" << rootFileName << "'." << std::endl;
	}
	std::unique_ptr<TFile> outFile(TFile::Open(rootFileName.c_str(), "RECREATE"));
	if(not outFile or outFile->IsZombie()) {
		printErr << "error while creating ROOT file '" << rootFileName << "' for plots of fit result."<< std::endl;
		return 1;
	}
	if(not mdepFit.createPlots(fitModel, fitParameters, cache, outFile.get(), rangePlotting, extraBinning)) {
		printErr << "error while creating plots." << std::endl;
		return 1;
	}
	outFile->Close();

	return 0;
}
