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
// File and Version Information:
// $Id$
//
// Description:
//      fitting program for rootpwa
//      minimizes TPWALikelihood function
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

#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TComplex.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

#include "utilities.h"
#include "TFitBin.h"
#include "TFitResult.h"
#include "TPWALikelihood.h"


using namespace std;
using namespace ROOT::Math;


void
usage(const string& progName,
      const int     errCode = 0)
{
  cerr << "usage:" << endl
       << progName
       << " -l # -u # -w wavelist [-d amplitude directory -o outfile -S start value file -N -n normfile"
       << " [-a normfile] -r rank -M minimizer [-m algorithm] -q -h]" << endl
       << "    where:" << endl
       << "        -l #       lower edge of mass bin [MeV/c^2]" << endl
       << "        -u #       upper edge of mass bin [MeV/c^2]" << endl
       << "        -w file    path to wavelist file" << endl
       << "        -d dir     path to directory with decay amplitude files (default: '.')" << endl
       << "        -o file    path to output file (default: 'fitresult.root')" << endl
       << "        -S file    path to file with start values (default: none)" << endl
       << "        -N         use normalization of decay amplitudes (default: false)" << endl
       << "        -n file    path to normalization integral file (default: 'norm.int')" << endl
       << "        -a file    path to acceptance integral file (default: 'norm.int')" << endl
       << "        -r #       rank of spin density matrix (default: 1)" << endl
       << "        -M name    minimizer (default: Minuit2)" << endl
       << "        -m name    minimization algorithm (optional, default: Migrad)" << endl
       << "                   available minimizers: Minuit:      Migrad, Simplex, Minimize, Migrad_imp" << endl
       << "                                         Minuit2:     Migrad, Simplex, Combined, Scan Fumili" << endl
       << "                                         GSLMultiMin: ConjugateFR, ConjugatePR, BFGS, BFGS2, SteepestDescent" << endl
       << "                                         GSLMultiFit: -" << endl
       << "                                         GSLSimAn:    -" << endl
       << "                                         Linear:      Robust" << endl
       << "                                         Fumili:      -" << endl
       << "        -t #       minimizer tolerance (default: 1e-10)" << endl
       << "        -q         run quietly (default: false)" << endl
       << "        -h         print help" << endl
       << endl;
  exit(errCode);
}


int
main(int    argc,
     char** argv)
{
  // ---------------------------------------------------------------------------
  // internal parameters
  const string       valTreeName        = "pwa";
  const string       valBranchName      = "fitResult";
  const double       defaultStartValue  = 0.01;
  double             startValStep       = 0.0005;
  const unsigned int maxNmbOfIterations = 20000;
  const bool         runHesse           = false;
  const bool         runMinos           = false;

  // ---------------------------------------------------------------------------
  // parse command line options
  const string progName           = argv[0];
  double       massBinMin         = 0;                      // [MeV/c^2]
  double       massBinMax         = 0;                      // [MeV/c^2]
  string       waveListFileName   = "";                     // wavelist filename
  string       ampDirName         = ".";                    // decay amplitude directory name
  string       outFileName        = "fitresult.root";       // output filename
  string       startValFileName   = "";                     // file with start values
  bool         useStartVal        = false;                  // indicates whether there are valid start values
  bool         useNormalizedAmps  = false;                  // if true normalized amplitudes are used
  string       normIntFileName    = "";                     // file with normalization integrals
  string       accIntFileName     = "";                     // file with acceptance integrals
  unsigned int rank               = 1;                      // rank of fit
  string       minimizerType[2]   = {"Minuit2", "Migrad"};  // minimizer, minimization algorithm
  double       minimizerTolerance = 1e-10;                  // minimizer tolerance
  bool         quiet              = false;
  extern char* optarg;
  // extern int optind;
  int c;
  while ((c = getopt(argc, argv, "l:u:w:d:o:S:Nn:a:r:M:m:t:qh")) != -1)
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
    case 'd':
      ampDirName = optarg;
      break;
    case 'o':
      outFileName = optarg;
      break;
    case 'S':
      startValFileName = optarg;
      useStartVal      = 1;
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
    case 'r':
      rank = atoi(optarg);
      break;
    case 'M':
      minimizerType[0] = optarg;
      break;
    case 'm':
      minimizerType[1] = optarg;
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
    printWarn << "using default normalization integral file '" << normIntFileName << "'." << endl;
  }
  if (accIntFileName.length() <= 1) {
    accIntFileName = "norm.int";
    printWarn << "using default acceptance normalization integral file '" << accIntFileName << "'." << endl;
  }
  if (waveListFileName.length() <= 1) {
    printErr << "no wavelist file specified! aborting!" << endl;
    usage(progName, 1);
  }
  // report parameters
  printInfo << "running " << progName << " with the following parameters:" << endl;
  cout << "    mass bin [" <<  massBinMin << ", " <<  massBinMax << "] MeV/c^2" << endl
       << "    wave list file ......................... '" << waveListFileName << "'" << endl
       << "    output file ............................ '" << outFileName      << "'" << endl
       << "    file with start values ................. '" << startValFileName << "'" << endl
       << "    use normalization ...................... "  << useNormalizedAmps       << endl
       << "        file with normalization integral ... '" << normIntFileName  << "'" << endl
       << "        file with acceptance integral ...... '" << accIntFileName   << "'" << endl
       << "    rank of fit ............................ "  << rank                    << endl
       << "    minimizer .............................. "  << minimizerType[0] << ", " << minimizerType[1] << endl
       << "    quiet .................................. "  << quiet << endl;

  // ---------------------------------------------------------------------------
  // setup likelihood function
  printInfo << "creating and setting up likelihood function" << endl;
  TPWALikelihood<double> L;
  if (quiet)
    L.setQuiet();
  L.useNormalizedAmps(useNormalizedAmps);
  L.init(rank, waveListFileName, normIntFileName, accIntFileName, ampDirName);
  if (!quiet)
    cout << L << endl;
  const unsigned int nmbPar = L.NDim();
  
  // ---------------------------------------------------------------------------
  // setup minimizer
  printInfo << "creating and setting up minimizer " << minimizerType[0] << " using algorithm " << minimizerType[1] << endl;
  Minimizer* minimizer = Factory::CreateMinimizer(minimizerType[0], minimizerType[1]);
  if (!minimizer) { 
    printErr << "could not create minimizer! exiting!" << endl;
    throw;
  }
  minimizer->SetFunction(L);
  minimizer->SetPrintLevel((quiet) ? 0 : 3);

  // ---------------------------------------------------------------------------
  // read in TFitResult with start values
  printInfo << "reading start values from '" << startValFileName << "'." << endl;
  const double massBinCenter  = (massBinMin + massBinMax) / 2;
  TFitResult*  startFitResult = NULL;
  bool         hasStartVal    = false;
  TFile*       startValFile   = NULL;
  if (startValFileName.length() <= 2)
    printWarn << "start value file name '" << startValFileName << "' is invalid. using default start values." << endl;
  else {
    // open root file
    startValFile = TFile::Open(startValFileName.c_str(), "READ");
    if (!startValFile || startValFile->IsZombie())
      printWarn << "cannot open start value file '" << startValFileName << "'. using default start values." << endl;
    else {
      // get tree with start values
      TTree* tree;
      startValFile->GetObject(valTreeName.c_str(), tree);
      if (!tree)
	printWarn << "cannot find start value tree '"<< valTreeName << "' in file '" << startValFileName << "'." << endl;
      else {
	//startBin = new TFitBin();
	//tree->SetBranchAddress("fitbin", &startBin);
	startFitResult = new TFitResult();
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
	hasStartVal = true;
      }
    }
  }

  // ---------------------------------------------------------------------------
  // set start parameter values
  printInfo << "setting start values for " << nmbPar << " parameters." << endl
	    << "    parameter naming scheme is: V[rank index]_[IGJPCME][isobar spec]" << endl;
  unsigned int maxParNameLength = 0;       // maximum length of parameter names
  vector<bool> parIsFixed(nmbPar, false);  // memorizes state of variables; ROOT::Math::Minimizer has no corresponding accessor
  {
    bool success = true;
    for (unsigned int i = 0; i < nmbPar; ++i)
      if (L.parName(i).length() > maxParNameLength)
	maxParNameLength = L.parName(i).length();
    for (unsigned int i = 0; i < nmbPar; ++i) {
      double startVal      = defaultStartValue;
      const string parName = L.parName(i);
      if (hasStartVal) {
	// get parameter value from TFitResult
	assert(startFitResult);
	startVal = startFitResult->fitParameter(parName.c_str());
      }
      // check if parameter needs to be fixed because of threshold
      if ((L.parThreshold(i) == 0) || (L.parThreshold(i) < massBinCenter)) {
	if (!minimizer->SetVariable(i, parName, startVal, startValStep))
	  success = false;
	if (startVal == 0) {
	  cout << "    read start value 0 for parameter " << parName << ". using default start value." << endl;
	  startVal = defaultStartValue;
	}
	cout << "    setting parameter [" << setw(3) << i << "] "
	     << setw(maxParNameLength) << parName << " = " << maxPrecisionAlign(startVal) << endl;
      } else {
	if (!minimizer->SetFixedVariable(i, parName, 0.))  // fix this parameter to 0
	  success = false;
	cout << "    fixing parameter  [" << setw(3) << i << "] "
	     << setw(maxParNameLength) << parName << " = 0" << endl;
	parIsFixed[i] = true;
      }
      if (!success) {
	printErr << "something went wrong when setting log likelihood parameters! exiting." << endl;
	throw;
      }
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
  printInfo << "performing minimization." << endl;
  {
    minimizer->SetMaxIterations(maxNmbOfIterations);
    minimizer->SetTolerance    (minimizerTolerance);
    const bool success = minimizer->Minimize();
    if (success)
      printInfo << "minimization finished successfully." << endl;
    else
      printWarn << "minimization failed." << endl;
    if (runHesse) {
      printInfo << "calculating Hessian matrix." << endl;
      //success = minimizer->Hesse();  // comes only with ROOT 5.24+
      if (!success)
	printWarn << "calculation of Hessian matrix failed." << endl;
    }
    printInfo << "minimization stopped after " << minimizer->NCalls() << " function calls. minimizer status summary:" << endl
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
  }

  // ---------------------------------------------------------------------------
  // print results
  printInfo << "minimization result:" << endl;
  for (unsigned int i = 0; i< nmbPar; ++i) {
    cout << "    parameter [" << setw(3) << i << "] "
	 << setw(maxParNameLength) << L.parName(i) << " = ";
    if (parIsFixed[i])
      cout << minimizer->X()[i] << " (fixed)" << endl;
    else {
      cout << setw(12) << maxPrecisionAlign(minimizer->X()[i]) << " +- "
	   << setw(12) << maxPrecisionAlign(minimizer->Errors()[i]);
      if (runMinos && (i == 156)) {  // does not work for all parameters
	double minosErrLow = 0;
	double minosErrUp  = 0;
	const bool success = minimizer->GetMinosError(i, minosErrLow, minosErrUp);
	if (success)
	  cout << "    Minos: " << "[" << minosErrLow << ", +" << minosErrUp << "]" << endl;
      } else
	cout << endl;
    }
  }
  cout << "number of calls to likelihood function FdF() ... " << L.ncalls() << endl
       << "total time spent for likelihood calculation .... " << L.Ltime()  << " sec" << endl
       << "total time spent for normalization ............. " << L.Ntime()  << " sec" << endl;

  // ---------------------------------------------------------------------------
  // write out result
  printInfo << "writing result to '" << outFileName << "'." << endl;
  {
    // open output file and create tree for writing
    TFile* outFile = new TFile(outFileName.c_str(), "UPDATE");
    if (!outFile || outFile->IsZombie())
      printWarn << "cannot open output file '" << outFileName << "'. no results will be written." << endl;
    else {
      // check whether output tree already exists
      TTree*      tree;
      TFitBin*    result    = new TFitBin();
      TFitResult* fitResult = new TFitResult();
      outFile->GetObject(valTreeName.c_str(), tree);
      if (!tree) {
	printInfo << "file '" << outFileName << "' is empty. creating new tree '" << valTreeName << "' for PWA result." << endl;
	tree = new TTree(valTreeName.c_str(), valTreeName.c_str());
	tree->Branch("fitbin",              &result);  // depricated legacy branch
	tree->Branch(valBranchName.c_str(), &fitResult);
      } else {
	tree->SetBranchAddress("fitbin",              &result);
	tree->SetBranchAddress(valBranchName.c_str(), &fitResult);
      }

      { 
	// get data structures to construct TFitResult
	vector<complex<double> > prodAmps;                // production amplitudes
	vector<string>           prodAmpNames;            // names of production amplitudes used in fit
	vector<pair<int,int> >   fitParCovMatrixIndices;  // indices of fit parameters for real and imaginary part in covariance matrix matrix
	L.buildCAmps(minimizer->X(), prodAmps, fitParCovMatrixIndices, prodAmpNames, true);
	TMatrixT<double> fitParCovMatrix(nmbPar, nmbPar);  // covariance matrix of fit parameters
	for(unsigned int i = 0; i < nmbPar; ++i)
	  for(unsigned int j = 0; j < nmbPar; ++j)
	    fitParCovMatrix[i][j] = minimizer->CovMatrix(i,j);
	const unsigned int nmbWaves = L.nmbWaves() + 1;  // flat wave is not included in L.nmbWaves()
	TCMatrix normIntegral(nmbWaves, nmbWaves);  // normalization integral over full phase space without acceptance
	TCMatrix accIntegral (nmbWaves, nmbWaves);  // normalization integral over full phase space with acceptance
	L.getIntCMatrix(normIntegral, accIntegral);
	const int normNmbEvents = useNormalizedAmps ? 1 : L.nmbEvents();  // number of events to normalize to

	cout << "filling TFitResult:" << endl
	     << "    number of fit parameters ............... " << nmbPar                        << endl
	     << "    number of production amplitudes ........ " << prodAmps.size()               << endl
	     << "    number of production amplitude names ... " << prodAmpNames.size()           << endl
	     << "    number of wave names ................... " << nmbWaves                      << endl
	     << "    number of cov. matrix indices .......... " << fitParCovMatrixIndices.size() << endl
	     << "    dimension of covariance matrix ......... " << fitParCovMatrix.GetNrows() << " x " << fitParCovMatrix.GetNcols() << endl
	     << "    dimension of normalization matrix ...... " << normIntegral.nrows()       << " x " << normIntegral.ncols()       << endl
	     << "    dimension of acceptance matrix ......... " << accIntegral.nrows()        << " x " << accIntegral.ncols()        << endl;
	fitResult->fill(L.nmbEvents(),
			normNmbEvents,
			massBinCenter,
			minimizer->MinValue(),
			rank,
			prodAmps,
			prodAmpNames,
			fitParCovMatrix,
			fitParCovMatrixIndices,
			normIntegral);
      }

      if (1) { 
	// get data structures to construct TFitBin
	vector<TComplex>       prodAmplitudes;  // production amplitudes
	vector<pair<int,int> > indices;         // indices for error matrix access
	vector<TString>        waveNames;       // contains rank information 
	{
	  vector<complex<double> > V;
	  vector<string>           names;
	  L.buildCAmps(minimizer->X(), V, indices, names, true);
	  // convert to TComplex;
	  for (unsigned int i = 0; i < V.size(); ++i)
	    prodAmplitudes.push_back(TComplex(V[i].real(), V[i].imag()));
	  assert(indices.size() == prodAmplitudes.size());
	  for(unsigned int i = 0; i < names.size(); ++i)
	    waveNames.push_back(TString(names[i].c_str()));
	}
	vector<TString> waveTitles;  // without rank
	{
	  const vector<string>& titles = L.waveNames();
	  for(unsigned int i = 0; i < titles.size(); ++i)
	    waveTitles.push_back(TString(titles[i].c_str()));
	  waveTitles.push_back("flat");
	}
	// error matrix
	TMatrixD errMatrix(nmbPar, nmbPar);
	for(unsigned int i = 0; i < nmbPar; ++i)
	  for(unsigned int j = 0; j < nmbPar; ++j) {
	    errMatrix[i][j] = minimizer->CovMatrix(i,j);
	  }
	// normalixation integral and acceptance Matrix
	cout << " setting up integrals" << endl;
	const unsigned int n = waveTitles.size();
	TCMatrix integralMatrix(n, n);
	TCMatrix accMatrix     (n, n);
	L.getIntCMatrix(integralMatrix, accMatrix);
	//integralMatrix.Print();
	// representation of number of events depends on whether normalization was done
	const int nmbEvt = useNormalizedAmps ? 1 : L.nmbEvents();
	
	cout << "filling TFitBin:" << endl;
	cout << "    number of fit parameters ........... " << nmbPar                 << endl;
	cout << "    number of production amplitudes .... " << prodAmplitudes.size()  << endl;
	cout << "    number of indices .................. " << indices.size()         << endl;
	cout << "    number of wave names (with rank) ... " << waveNames.size()       << endl;
	cout << "    number of wave titles (w/o rank) ... " << waveTitles.size()      << endl;
	cout << "    dimension of error matrix .......... " << errMatrix.GetNrows()   << endl;
	cout << "    dimension of integral matrix ....... " << integralMatrix.nrows() << endl;
	cout << "    dimension of acceptance matrix ..... " << accMatrix.nrows()      << endl;
	result->fill(prodAmplitudes,
		     indices,
		     waveNames,
		     nmbEvt,
		     L.nmbEvents(), // raw number of data events
		     massBinCenter,
		     integralMatrix,
		     errMatrix,
		     minimizer->MinValue(),
		     rank);
	//result->PrintParameters();
      }

      // write result to file
      tree->Fill();
      tree->Write("", TObject::kOverwrite);
      outFile->Close();
    }
  }
  
  if (minimizer)
    delete minimizer;
  return 0;
}
