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
#include "partialWaveFitHelper.h"
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
	     << " [-d amplitude directory -R] -f fitResultFile [-o outfile -N -n normfile"
	     << " -a normfile -A # normalisation events -C -P width -q -h]" << endl
	     << "    where:" << endl
	     << "        -d dir     path to directory with decay amplitude files (default: '.')" << endl
	     << "        -R         use .root amplitude files (default: false)" << endl
	     << "        -f file    path to fit result file with parameters" << endl
	     << "        -o file    path to output file (default: 'eigenSlices.root')" << endl
	     << "        -N         use normalization of decay amplitudes (default: false)" << endl
	     << "        -n file    path to normalization integral file (default: 'norm.int')" << endl
	     << "        -a file    path to acceptance integral file (default: 'norm.int')" << endl
	     << "        -A #       number of input events to normalize acceptance to (default: use number of events from acceptance integral file)" << endl
	     << "        -C         use half-Cauchy priors (default: false)" << endl
	     << "        -P #       width of half-Cauchy priors (default: 0.5)" << endl
	     << "        -q         run quietly (default: false)" << endl
	     << "        -h         print help" << endl
	     << endl;
	exit(errCode);
}


std::vector<double>
getEigenvector(const TMatrixT<double>& eigenvectors,
               const unsigned int& col)
{
	std::vector<double> eigenvector(eigenvectors.GetNrows());
	for(int r = 0; r < eigenvectors.GetNrows(); ++r) {
		eigenvector[r] = eigenvectors[r][col];
	}
	return eigenvector;
}


std::vector<double>
getParsForPoint(const std::vector<double>& eigenvector,
                const std::vector<double>& maximum,
                const double& eigenvalue,
                const int& p)
{
	std::vector<double> newPoint(maximum.size(), 0.);
	for(unsigned int r = 0; r < maximum.size(); ++r) {
		newPoint[r] = maximum[r] + eigenvector[r] * p / 50 * std::sqrt(eigenvalue);
	}
	return newPoint;
}


bool
checkNorm(const std::vector<double>& eigenvector)
{
	double norm = 0.;
	for(unsigned int r = 0; r < eigenvector.size(); ++r) {
		norm += eigenvector[r] * eigenvector[r];
	}
	if(std::fabs(norm - 1.) < 1e-14) {
		return true;
	} else {
		return false;
	}
}


void
normalizeVector(std::vector<double>& eigenvector)
{
	double norm = 0.;
	for (unsigned int r = 0; r < eigenvector.size(); ++r) {
		norm += eigenvector[r] * eigenvector[r];
	}
	for (unsigned int r = 0; r < eigenvector.size(); ++r) {
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
	string       ampDirName          = ".";                    // decay amplitude directory name
	bool         useRootAmps         = false;                  // if true .root amplitude files are read
	//bool         useRootAmps         = true;                   // if true .root amplitude files are read
	string       inFileName          = "";                     // input filename
	string       outFileName         = "eigenSlices.root";     // output filename
	bool         useNormalizedAmps   = false;                  // if true normalized amplitudes are used
	string       normIntFileName     = "";                     // file with normalization integrals
	string       accIntFileName      = "";                     // file with acceptance integrals
	unsigned int numbAccEvents       = 0;                      // number of events used for acceptance integrals
	bool         cauchyPriors        = false;
	double       cauchyWidth         = 0.5;
	bool         quiet               = false;
	extern char* optarg;
	// extern int optind;
	int c;
	while ((c = getopt(argc, argv, "d:Rf:o:Nn:a:A:CP:qh")) != -1)
		switch (c) {
		case 'd':
			ampDirName = optarg;
			break;
		case 'R':
			useRootAmps = true;
			break;
		case 'f':
			inFileName = optarg;
			break;
		case 'o':
			outFileName = optarg;
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
		case 'C':
			cauchyPriors = true;
			break;
		case 'P':
			cauchyWidth = atof(optarg);
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
	if (inFileName.length() <= 1) {
		printErr << "no input file name specified. Aborting..." << endl;
		usage(progName, 1);
	}
	// report parameters
	printInfo << "running " << progName << " with the following parameters:" << endl;
	cout << "    path to amplitude directory .................... '" << ampDirName               << "'" << endl
	     << "    use .root amplitude files ...................... "  << yesNo(useRootAmps)       << endl
	     << "    path to input file ............................. '" << inFileName               << "'" << endl
	     << "    path to output file ............................ '" << outFileName              << "'" << endl
	     << "    use normalization .............................. "  << yesNo(useNormalizedAmps) << endl
	     << "        path to file with normalization integral ... '" << normIntFileName          << "'" << endl
	     << "        path to file with acceptance integral ...... '" << accIntFileName           << "'" << endl
	     << "        number of acceptance norm. events .......... "  << numbAccEvents            << endl
	     << "    use half-Cauchy priors ......................... '" << yesNo(cauchyPriors)      << "'" << endl;
	if(cauchyPriors) {
		cout << "    width of cauchy priors.......................... "  << cauchyWidth << endl;
	}
	cout << "    quiet .......................................... '" << yesNo(quiet)             << "'" << endl;

	TFile* inFile = TFile::Open(inFileName.c_str(), "READ");
	if(not inFile || inFile->IsZombie()) {
		printErr << "could not open input file '" << inFileName << "'. Aborting..." << endl;
		return 1;
	}

	TTree* inTree = 0;
	inFile->GetObject(valTreeName.c_str(), inTree);
	if(not inTree) {
		printErr << "could not find result tree '" << valTreeName << "' in input file '" << inFileName << "'. Aborting..." << endl;
		return 1;
	}

	fitResult* result = 0;
	inTree->SetBranchAddress(valBranchName.c_str(), &result);

	if(inTree->GetEntries() != 1) {
		printErr << "result tree '" << valTreeName << "' has more than one entry, NOT IMPLEMENTED." << endl;
		return 1;
	}
	inTree->GetEntry(0);

	// create temporary file to store the wavelist
	char tempFileName[] = "XXXXXX";
	close(mkstemp(tempFileName));
	const string waveListFileName(tempFileName);
	ofstream waveListFile(waveListFileName.c_str());
	rpwa::partialWaveFitHelper::extractWaveList(*result, waveListFile);
	waveListFile.close();

	printInfo << "parameters extracted from input fit result" << endl;
	cout << "    mass bin centered at ........................... "  << result->massBinCenter() << " MeV/c^2" << endl
	     << "    path to temporary wave list file ............... '" << waveListFileName        << "'" << endl
	     << "    rank of spin density matrix .................... "  << result->rank()          << endl;

	// ---------------------------------------------------------------------------
	// setup likelihood function
	printInfo << "creating and setting up likelihood function" << endl;
	pwaLikelihood<complex<double> > L;
	if (quiet)
		L.setQuiet();
	L.useNormalizedAmps(useNormalizedAmps);
	if (cauchyPriors) {
		L.setPriorType(L.HALF_CAUCHY);
		L.setCauchyWidth(cauchyWidth);
	}
	L.init(result->rank(), result->massBinCenter(), waveListFileName, normIntFileName, accIntFileName,
	       ampDirName, numbAccEvents, useRootAmps);
	remove(waveListFileName.c_str());
	if (not quiet)
		cout << L << endl;
	const unsigned int nmbPar      = L.NDim();
	const unsigned int nmbParFixed = L.nmbParsFixed();

	unsigned int maxParNameLength = 0;  // maximum length of parameter names
	for(unsigned int i = 0; i < nmbPar; ++i) {
		const string& parName = L.parName(i);
		if (parName.length() > maxParNameLength)
			maxParNameLength = parName.length();
	}
	std::vector<double> pars(nmbPar);
	for(unsigned int i = 0; i < nmbPar; ++i) {
		const string& parName = L.parName(i);
		pars[i] = result->fitParameter(parName);
		printInfo << "setting parameter [" << setw(3) << i << "] " << setw(maxParNameLength) << parName << " = " << maxPrecisionAlign(pars[i]) << endl;
	}

	const double maxLogLikelihood = L.DoEval(pars.data());
	printInfo << "likelihood at maximum is " << maxPrecisionAlign(maxLogLikelihood) << "." << endl;

	const TMatrixT<double> covMatrixMinuit = result->fitParCovMatrix();
	TVectorT<double> eigenvaluesMinuit;
	TMatrixT<double> eigenvectorsMinuit;
	rpwa::partialWaveFitHelper::getEigenvectors(L, covMatrixMinuit, eigenvectorsMinuit, eigenvaluesMinuit);

	const TMatrixT<double> covMatrixAna = L.CovarianceMatrix(pars.data());
	TVectorT<double> eigenvaluesAna;
	TMatrixT<double> eigenvectorsAna;
	rpwa::partialWaveFitHelper::getEigenvectors(L, covMatrixAna, eigenvectorsAna, eigenvaluesAna);

	TFile* outFile = new TFile(outFileName.c_str(), "RECREATE");
	for(unsigned int par = 0; par < nmbPar-nmbParFixed; ++par) {
		std::vector<double> eigenvectorMinuit = getEigenvector(eigenvectorsMinuit, par);
		if (not checkNorm(eigenvectorMinuit)) {
			printWarn << "Minuit eigenvector " << par << " is not normalized. normalizing." << endl;
			normalizeVector(eigenvectorMinuit);
			// check again
			if (not checkNorm(eigenvectorMinuit)) {
				printErr << "could not normalize Minuit eigenvector! Aborting..." << endl;
				return 1;
			}
		}
		printInfo << "eigenvector          [" << setw(3) << par << "]: " << eigenvectorMinuit << endl;
		printInfo << "eigenvalue           [" << setw(3) << par << "]: " << eigenvaluesMinuit[par] << endl;

		std::vector<double> eigenvectorAna = getEigenvector(eigenvectorsAna, par);
		if (not checkNorm(eigenvectorAna)) {
			printWarn << "analytic eigenvector " << par << " is not normalized. normalizing." << endl;
			normalizeVector(eigenvectorAna);
			// check again
			if (not checkNorm(eigenvectorAna)) {
				printErr << "could not normalize analytic eigenvector! Aborting..." << endl;
				return 1;
			}
		}
		printInfo << "analytic eigenvector [" << setw(3) << par << "]: " << eigenvectorAna << endl;
		printInfo << "analytic eigenvalue  [" << setw(3) << par << "]: " << eigenvaluesAna[par] << endl;

		// The idea here is to scan the logLikelihood function from the
		// minimum in direction of the Eigenvectors.
		//
		// Around the minimum p the logLikelihood can be expanded as a
		// Taylor series:
		//
		// -logL(p+dp) = -logL(p) + dp^T * grad(-logL(x))|x=p + 1/2 * dp^T * hesse(-logL(x))|x=p * dp
		//
		// The second summand should be equal to zero at the minimum.
		// It is assumed that the Hessian matrix of -logL is the
		// inverse of the covariance matrix
		//
		// hesse(-logL(x))|x=p = C^-1
		//
		// Eigenvectors of the covariance matrix then are also
		// Eigenvectors of the Hessian matrix, albeit with inverse
		// Eigenvalues. If e is an Eigenvector of the covariance matrix
		// with Eigenvalue l
		//
		// C * e = l * e
		//
		// then
		//
		// hesse(-logL(x))|x=p * e = 1/l * e
		//
		// We now write dp = r * sqrt(l) * e where r in the following
		// is varied between -1 and +1. For a logLikelihood behaving
		// nicely we then should find the parabola
		//
		// -logL(p+dp) = -logL(p) + 1/2 * r^2
		//
		// This ideal parabola is drawn for a x between -sqrt(l) and
		// +sqrt(l) using x = r * sqrt(l).
		TGraph graphLikeli;
		for(int p = -50; p <= 50; ++p) {
			std::vector<double> newPars = getParsForPoint(eigenvectorMinuit, pars, eigenvaluesMinuit[par], p);
			const double likeli = L.DoEval(newPars.data()) - maxLogLikelihood;
			graphLikeli.SetPoint(p + 50, p / 50. * std::sqrt(eigenvaluesMinuit[par]), likeli);
		}

		const double lowerLimit = -1.2 * std::sqrt(eigenvaluesMinuit[par]);
		const double upperLimit =  1.2 * std::sqrt(eigenvaluesMinuit[par]);

		std::stringstream minuitParabolaFormula;
		// the magnitude of the eigenvector is 1. otherwise an
		// additional factor '|eigenvector|^2' would be required
		minuitParabolaFormula << "0.5 / " << eigenvaluesMinuit[par] << " * x*x";
		TF1 parabolaMinuit("minuitParabola", minuitParabolaFormula.str().c_str(), lowerLimit, upperLimit);

		std::stringstream anaParabolaFormula;
		// the magnitude of the eigenvector is 1. otherwise an
		// additional factor '|eigenvector|^2' would be required
		anaParabolaFormula << "0.5 / " << eigenvaluesAna[par] << " * x*x";
		TF1 parabolaAna("analyticParabola", anaParabolaFormula.str().c_str(), lowerLimit, upperLimit);

		std::stringstream canvasName;
		canvasName << "eigenvectorSlice" << par;
		TCanvas cnvs(canvasName.str().c_str());
		cnvs.cd();
		graphLikeli.Draw("A*");
		parabolaMinuit.Draw("Lsame");
		parabolaMinuit.SetLineColor(kBlue);
		parabolaAna.Draw("Lsame");
		parabolaAna.SetLineColor(kRed);
		cnvs.Write();
	}
	outFile->Close();
	printSucc << "slices successfully written to file '" << outFileName <<"'." << endl;

	printInfo << "setting Minuit line color to blue, analytical solution to red." << endl;

	return 0;
}
