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
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TStopwatch.h"
#include "Minuit2/Minuit2Minimizer.h"

#include "reportingUtilsEnvironment.h"
#include "conversionUtils.hpp"
#include "pwaLikelihood.h"
#include "fitResult.h"

#include "amplitudeTreeLeaf.h"


using namespace std;
using namespace ROOT::Math;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{

	cerr << "calculate likelihood at a point in parameter space given by a fit result file" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " -f fitResultFile [-d amplitude directory -R -N -n normfile"
	     << " [-a normfile] -r rank [-q -h]" << endl
	     << "    where:" << endl
	     << "        -f file    path to fit result file with parameters" << endl
	     << "        -C         use half-Cauchy priors" << endl
	     << "        -d dir     path to directory with decay amplitude files (default: '.')" << endl
#ifdef USE_STD_COMPLEX_TREE_LEAFS
	     << "        -R         use .root amplitude files (default: false)" << endl
#else
	     << "        -R         use .root amplitude files [not supported; ROOT version too low]" << endl
#endif
	     << "        -N         use normalization of decay amplitudes (default: false)" << endl
	     << "        -n file    path to normalization integral file (default: 'norm.int')" << endl
	     << "        -a file    path to acceptance integral file (default: 'norm.int')" << endl
	     << "        -A #       number of input events to normalize acceptance to" << endl
	     << "        -r #       rank of spin density matrix (default: 1)" << endl
	     << "        -q         run quietly (default: false)" << endl
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

	// ---------------------------------------------------------------------------
	// parse command line options
	const string progName            = argv[0];
	string       fitResultFileName   = "";
	bool         cauchyPriors        = false;
	string       ampDirName          = ".";                    // decay amplitude directory name
	bool         useRootAmps         = false;                  // if true .root amplitude files are read
	//bool         useRootAmps         = true;                   // if true .root amplitude files are read
	bool         useNormalizedAmps   = false;                  // if true normalized amplitudes are used
	string       normIntFileName     = "";                     // file with normalization integrals
	string       accIntFileName      = "";                     // file with acceptance integrals
	unsigned int numbAccEvents       = 0;                      // number of events used for acceptance integrals
	unsigned int rank                = 1;                      // rank of fit
	bool         quiet               = false;
	extern char* optarg;
	// extern int optind;
	int c;
	while ((c = getopt(argc, argv, "f:Cd:RNn:a:A:r:qh")) != -1)
		switch (c) {
		case 'f':
			fitResultFileName = optarg;
			break;
		case 'C':
			cauchyPriors = true;
			break;
		case 'd':
			ampDirName = optarg;
			break;
		case 'R':
#ifdef USE_STD_COMPLEX_TREE_LEAFS
			useRootAmps = true;
#endif
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

	// report parameters
	printInfo << "running " << progName << " with the following parameters:" << endl;
	cout << "    path to amplitude directory .................... '" << ampDirName       << "'" << endl
	     << "    use .root amplitude files ...................... "  << yesNo(useRootAmps)      << endl
	     << "    use normalization .............................. "  << yesNo(useNormalizedAmps) << endl
	     << "    path to file with normalization integral ... '" << normIntFileName  << "'" << endl
	     << "    path to file with acceptance integral ...... '" << accIntFileName   << "'" << endl
	     << "    number of acceptance norm. events .......... "  << numbAccEvents    << endl
	     << "    rank of spin density matrix .................... "  << rank                    << endl
	     << "    quiet .......................................... "  << yesNo(quiet) << endl;


	TFile* fitResultFile = TFile::Open(fitResultFileName.c_str(), "READ");
	if(not fitResultFile) {
		printErr << "could not open fit result file '" << fitResultFileName << "'. Aborting..." << endl;
		return 1;
	}

	TTree* fitResultTree = 0;
	fitResultFile->GetObject("pwa", fitResultTree);
	if(not fitResultTree) {
		printErr << "could not find result tree in fit result file '" << fitResultFileName << "'. Aborting..." << endl;
		return 1;
	}

	fitResult* result = 0;
	fitResultTree->SetBranchAddress("fitResult_v2", &result);

	if(fitResultTree->GetEntries() != 1) {
		printErr << "result tree has more than one entry, NOT IMPLEMENTED" << endl;
		return 1;
	}
	fitResultTree->GetEntry(0);

	char tempFileName[] = "XXXXXX";
	close(mkstemp(tempFileName));
	ofstream waveListFile(tempFileName);
	string waveListFileName(tempFileName);

	for(unsigned int i = 0; i < result->nmbWaves(); ++i) {
		const string& waveName(result->waveName(i).Data());
		if(waveName != "flat") {
			waveListFile << waveName << "\n";
		}
	}
	waveListFile.close();

	printInfo << "production amplitudes in fit result file: " << endl;
	for(unsigned int i = 0; i < result->nmbProdAmps(); ++i) {
		cout << "    " << result->prodAmpName(i) << endl;
	}

	// ---------------------------------------------------------------------------
	// setup likelihood function
	printInfo << "creating and setting up likelihood function" << endl;
	pwaLikelihood<complex<double> > L;
	if (quiet)
		L.setQuiet();
	L.useNormalizedAmps(useNormalizedAmps);
	L.init(rank, waveListFileName, normIntFileName, accIntFileName,
	       ampDirName, numbAccEvents, useRootAmps);
	remove(waveListFileName.c_str());
	waveListFileName = "";
	if (not quiet)
		cout << L << endl;
	const unsigned int nmbPar  = L.NDim();

	if(cauchyPriors) {
		L.setPriorType(L.HALF_CAUCHY);
	}

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

	double* pars = new double[nmbPar];
	for(unsigned int i = 0; i < nmbPar; ++i) {
		const string parName = L.parName(i);
		pars[i] = result->fitParameter(parName);
		printInfo << "setting parameter[" << i << "] (name: '" << parName << "') = " << pars[i] << endl;
	}

	cout << endl;
	printSucc << "Likelihood is " << setprecision(10) << L.DoEval(pars) << endl;
	cout << endl;

	return 0;
}
