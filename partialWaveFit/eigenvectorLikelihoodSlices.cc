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
#include "TMatrixT.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TVectorT.h"
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
	     << " -l # -u # -w wavelist -f fitResultFile -o outputFile [-d amplitude directory -R -N -n normfile"
	     << " [-a normfile] -r rank [-t # -q -h]" << endl
	     << "    where:" << endl
	     << "        -l #       lower edge of mass bin [MeV/c^2]" << endl
	     << "        -u #       upper edge of mass bin [MeV/c^2]" << endl
	     << "        -w file    path to wavelist file" << endl
	     << "        -f file    path to fit result file with parameters" << endl
	     << "        -o file    path to output file" << endl
	     << "        -C         use half-Cauchy priors" << endl
	     << "        -d dir     path to directory with decay amplitude files (default: '.')" << endl
	     << "        -R         use .root amplitude files (default: false)" << endl
	     << "        -N         use normalization of decay amplitudes (default: false)" << endl
	     << "        -n file    path to normalization integral file (default: 'norm.int')" << endl
	     << "        -a file    path to acceptance integral file (default: 'norm.int')" << endl
	     << "        -A #       number of input events to normalize acceptance to" << endl
	     << "        -r #       rank of spin density matrix (default: 1)" << endl
	     << "        -t #       minimizer tolerance (default: 0.001)" << endl
	     << "        -q         run quietly (default: false)" << endl
	     << "        -h         print help" << endl
	     << endl;
	exit(errCode);
}

std::vector<double> getEigenvector(const TMatrixT<double>& eigenvectors, const unsigned int& col) {
	std::vector<double> eigenvector(eigenvectors.GetNrows());
	for(int r = 0; r < eigenvectors.GetNrows(); r++) {
		eigenvector[r] = eigenvectors[r][col];
	}
	return eigenvector;
}

std::vector<double> getParsForPoint(const std::vector<double>& eigenvector, const std::vector<double>& maximum, const double& eigenvalue, const int& p) {
	std::vector<double> newPoint(maximum.size(), 0.);
	for(unsigned int r = 0; r < maximum.size(); r++) {
		newPoint[r] = maximum[r] + eigenvector[r] * p / 50 * std::sqrt(eigenvalue);
	}
	return newPoint;
}

bool checkNorm(const std::vector<double>& eigenvector) {
	double norm = 0.;
	for(unsigned int r = 0; r < eigenvector.size(); r++) {
		norm += eigenvector[r] * eigenvector[r];
	}
	if(std::fabs(norm - 1.) < 1e-14) {
		return true;
	} else {
		return false;
	}
}

void normalizeVector(std::vector<double>& eigenvector) {
	double norm = 0.;
	for (unsigned int r = 0; r < eigenvector.size(); r++) {
		norm += eigenvector[r] * eigenvector[r];
	}
	for (unsigned int r = 0; r < eigenvector.size(); r++) {
		eigenvector[r] /= std::sqrt(norm);
	}
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
	double       massBinMin          = 0;                      // [MeV/c^2]
	double       massBinMax          = 0;                      // [MeV/c^2]
	string       waveListFileName    = "";                     // wavelist filename
	string       fitResultFileName   = "";                     // fit result filename
	string       outputFileName      = "eigenSlices.root";     // output filename
	bool         cauchyPriors        = false;
	string       ampDirName          = ".";                    // decay amplitude directory name
	bool         useRootAmps         = false;                  // if true .root amplitude files are read
	//bool         useRootAmps         = true;                   // if true .root amplitude files are read
	bool         useNormalizedAmps   = false;                  // if true normalized amplitudes are used
	string       normIntFileName     = "";                     // file with normalization integrals
	string       accIntFileName      = "";                     // file with acceptance integrals
	unsigned int numbAccEvents       = 0;                      // number of events used for acceptance integrals
	unsigned int rank                = 1;                      // rank of fit
	double       minimizerTolerance  = 0.001;                   // minimizer tolerance
	bool         quiet               = false;
	extern char* optarg;
	// extern int optind;
	int c;
	while ((c = getopt(argc, argv, "l:u:w:f:o:Cd:RNn:a:A:r:t:qh")) != -1)
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
		case 'f':
			fitResultFileName = optarg;
			break;
		case 'o':
			outputFileName = optarg;
			break;
		case 'C':
			cauchyPriors = true;
			break;
		case 'd':
			ampDirName = optarg;
			break;
		case 'R':
			useRootAmps = true;
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
		case 't':
			minimizerTolerance = atof(optarg);
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
		printErr << "no wavelist file specified. aborting." << endl;
		usage(progName, 1);
	}
	// report parameters
	printInfo << "running " << progName << " with the following parameters:" << endl;
	cout << "    mass bin [" <<  massBinMin << ", " <<  massBinMax << "] MeV/c^2" << endl
	     << "    path to wave list file ......................... '" << waveListFileName << "'" << endl
	     << "    path to amplitude directory .................... '" << ampDirName       << "'" << endl
	     << "    use .root amplitude files ...................... "  << yesNo(useRootAmps)      << endl
	     << "    use normalization .............................. "  << yesNo(useNormalizedAmps) << endl
	     << "    path to file with normalization integral ... '" << normIntFileName  << "'" << endl
	     << "    path to file with acceptance integral ...... '" << accIntFileName   << "'" << endl
	     << "    number of acceptance norm. events .......... "  << numbAccEvents    << endl
	     << "    rank of spin density matrix .................... "  << rank                    << endl
	     << "    minimizer tolerance ............................ "  << minimizerTolerance << endl
	     << "    quiet .......................................... "  << yesNo(quiet) << endl;

	// ---------------------------------------------------------------------------
	// setup likelihood function
	printInfo << "creating and setting up likelihood function" << endl;
	pwaLikelihood<complex<double> > L;
	if (quiet)
		L.setQuiet();
	L.useNormalizedAmps(useNormalizedAmps);
	L.init(rank, waveListFileName, normIntFileName, accIntFileName,
	       ampDirName, numbAccEvents, useRootAmps);
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

	std::vector<double> pars(nmbPar);
	for(unsigned int i = 0; i < nmbPar; ++i) {
		const string parName = L.parName(i);
		pars[i] = result->fitParameter(parName);
		printInfo << "setting parameter[" << i << "] (name: '" << parName << "') = " << pars[i] << endl;
	}

	TMatrixT<double> covMatrixMinuit = result->fitParCovMatrix();
	TMatrixT<double> covMatrixAna = L.CovarianceMatrixAnalytically(pars.data());
	TVectorT<double> eigenvaluesMinuit;
	TMatrixT<double> eigenvectorsMinuit = covMatrixMinuit.EigenVectors(eigenvaluesMinuit);
	TVectorT<double> eigenvaluesAna;
	TMatrixT<double> eigenvectorsAna = covMatrixAna.EigenVectors(eigenvaluesAna);
	double maximum = L.DoEval(pars.data());

	const unsigned int nmbGraphs = nmbPar;
	TFile* outFile = new TFile(outputFileName.c_str(), "RECREATE");
	for(unsigned int par = 0; par < nmbGraphs; par++) {
		std::stringstream canvasName;
		canvasName << "eigenvectorSlice" << par;
		TCanvas cnvs(canvasName.str().c_str());
		cnvs.cd();
		TGraph graphLikeli;
		std::vector<double> eigenvector = getEigenvector(eigenvectorsMinuit, par);
		if (not checkNorm(eigenvector)) {
			printWarn << "eigenvector " << par << " is not normalized. Normalizing..." << endl;
			normalizeVector(eigenvector);
			// check again
			if (not checkNorm(eigenvector)) {
				printErr << "could not normalize eigenvector! Aborting..." << endl;
				return 1;
			}
		}
		printInfo << "Eigenvector [" << par << "]: " << eigenvector << endl;
		printInfo << "Eigenvalue  [" << par << "]:" << eigenvaluesMinuit[par] << endl;

		for(int p = -50; p < 50; p++) {
			std::vector<double> newPars = getParsForPoint(eigenvector, pars, eigenvaluesMinuit[par], p);
			const double likeli = L.DoEval(newPars.data()) - maximum;
			graphLikeli.SetPoint(p + 50, p / 50. * std::sqrt(eigenvaluesMinuit[par]), likeli);
		}

		std::stringstream minuitGaussFormula;
		minuitGaussFormula << "-1*TMath::Gaus(x, 0, TMath::Sqrt(" << (eigenvaluesMinuit[par]) << "))+" << (0+1);
		TF1 gaussMinuit("minuitGauss", minuitGaussFormula.str().c_str(), -100, 100);
		cout << minuitGaussFormula.str() << endl;

		std::stringstream anaGaussFormula;
		anaGaussFormula << "-1*TMath::Gaus(x, 0, TMath::Sqrt(" << (eigenvaluesAna[par]) << "))+" << (0+1);
		TF1 gaussAna("analyticGauss", anaGaussFormula.str().c_str(), -100, 100);
		cout << anaGaussFormula.str() << endl;

		graphLikeli.Draw("A*");
		gaussMinuit.Draw("Lsame");
		gaussMinuit.SetLineColor(kBlue);
		gaussAna.Draw("Lsame");
		cnvs.Write();
	}
	printInfo << "Setting Minuit line color to blue, analytical solution to red." << endl;
	printInfo << "Likelihood at maximum is " << setprecision(10) << L.DoEval(pars.data()) << endl;
	outFile->Close();
	printSucc << "Slices successfully written to file '" << outputFileName <<"'." << endl;

	return 0;
}
