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
       << " -l # -u # -w wavelist [-o outfile -S start value file -N -n normfile [-a normfile]"
       << "-r rank -M minimizer [-m algorithm] -q -h]" << endl
       << "    where:" << endl
       << "        -l #       lower edge of mass bin [MeV/c^2]" << endl
       << "        -u #       upper edge of mass bin [MeV/c^2]" << endl
       << "        -w file    path to wavelist file" << endl
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
  double             defaultMinStep     = 0.0005;
  const unsigned int maxNmbOfIterations = 20000;
  const bool         runHesse           = false;
  const bool         runMinos           = false;

  // ---------------------------------------------------------------------------
  // parse command line options
  const string progName         = argv[0];
  double       binLowMass       = 0;                      // [MeV/c^2]
  double       binHighMass      = 0;                      // [MeV/c^2]
  string       waveListFileName = "";                     // wavelist filename
  string       outFileName      = "fitresult.root";       // output filename
  string       startValFileName = "";                     // file with start values
  bool         useStartVal      = false;                  // are there start values?
  bool         useNorm          = false;
  string       normIntFileName  = "";                     // file with normalization integrals
  string       accIntFileName   = "";                     // file with acceptance integrals
  unsigned int rank             = 1;                      // rank
  string       minimizerType[2] = {"Minuit2", "Migrad"};  // minimizer, minimization algorithm
  bool         quiet            = false;
  extern char* optarg;
  // extern int optind;
  int c;
  while ((c = getopt(argc, argv, "l:u:w:o:S:Nn:a:r:M:m:qh")) != -1)
    switch (c) {
    case 'l':
      binLowMass = atof(optarg);
      break;
    case 'u':
      binHighMass = atof(optarg);
      break;
    case 'w':
      waveListFileName = optarg;
      break;
    case 'o':
      outFileName = optarg;
      break;
    case 'S':
      startValFileName = optarg;
      useStartVal = 1;
      break;
    case 'N':
      useNorm = true;
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
    case 'q':
      quiet = true;
      break;
    case 'h':
      usage(progName);
      break;
    }
  if (normIntFileName.length() <= 1) {
    normIntFileName="norm.int";
    printWarn << "Using default normalization integral file '" << normIntFileName << "'." << endl;
  }
  if (accIntFileName.length() <= 1) {
    accIntFileName="norm.int";
    printWarn << "Using default acceptance normalization integral file '" << accIntFileName << "'." << endl;
  }
  if (waveListFileName.length() <= 1) {
    printErr << "No wavelist file specified! Aborting!" << endl;
    usage(progName, 1);
  }
  // report parameters
  printInfo << "Running " << progName << " with the following parameters:" << endl;
  cout << "    mass bin [" <<  binLowMass << ", " <<  binHighMass << "] MeV/c^2" << endl
       << "    wave list file ......................... '" << waveListFileName << "'" << endl
       << "    output file ............................ '" << outFileName      << "'" << endl
       << "    file with start values ................. '" << startValFileName << "'" << endl
       << "    use normalization ...................... "  << useNorm          << endl
       << "        file with normalization integral ... '" << normIntFileName  << "'" << endl
       << "        file with acceptance integral ...... '" << accIntFileName   << "'" << endl
       << "    rank ................................... "  << rank             << endl
       << "    minimizer .............................. "  << minimizerType[0] << ", " << minimizerType[1] << endl
       << "    quiet .................................. "  << quiet            << endl;

  // ---------------------------------------------------------------------------
  // setup likelihood function
  printInfo << "Creating and setting up likelihood function" << endl;
  TPWALikelihood<double> L;
  if (!quiet)
    L.SetQuiet();
  L.UseNormalizedAmps(useNorm);
  L.SetWavelist(waveListFileName);
  L.SetRank(rank);
  L.LoadIntegrals(normIntFileName, accIntFileName);
  L.LoadAmplitudes();
  const unsigned int nmbPar = L.NDim();
  
  // ---------------------------------------------------------------------------
  // setup minimizer
  printInfo << "Creating and setting up minimizer " << minimizerType[0] << " using algorithm " << minimizerType[1] << endl;
  Minimizer* minimizer = Factory::CreateMinimizer(minimizerType[0], minimizerType[1]);
  if (!minimizer) { 
    printErr << "Error creating minimizer. Exiting." << endl;
    throw;
  }
  minimizer->SetFunction(L);
  minimizer->SetPrintLevel((quiet) ? 0 : 3);

  // ---------------------------------------------------------------------------
  // read in object with start values
  printInfo << "Reading start values from '" << startValFileName << "'." << endl;
  double   binCenter    = 0.5 * (binLowMass + binHighMass);
  bool     hasStartVal  = false;
  TFile*   startValFile = NULL;
  TFitBin* startBin     = NULL;
  if (startValFileName.length() > 2) {
    // open root file
    startValFile = TFile::Open(startValFileName.c_str(), "READ");
    if (!startValFile || startValFile->IsZombie())
      printWarn << "Cannot open start value file '" << startValFileName << "'. Using default start values." << endl;
    else {
      // get tree with start values
      TTree* tree;
      startValFile->GetObject(valTreeName.c_str(), tree);
      if (!tree)
	printWarn << "Cannot find start value tree '"<< valTreeName << "' in file '" << startValFileName << "'." << endl;
      else {
	startBin = new TFitBin();
	tree->SetBranchAddress("fitbin", &startBin);
	//tree->SetBranchAddress(valBranchName.c_str(), &startBin);
	// find entry which is closest to mass bin center
	unsigned int iBest = 0;
	double mBest = 0;
	for (unsigned int i = 0; i < tree->GetEntriesFast(); ++i) {
	  tree->GetEntry(i);
	  if (fabs(binCenter - startBin->mass()) <= fabs(binCenter - mBest)) {
	    iBest = i;
	    mBest = startBin->mass();
	  }
	}  // end loop over TFitBins
	tree->GetEntry(iBest);
	hasStartVal = true;
      }
    }
  }  // endif startValFileName non-empty

  // ---------------------------------------------------------------------------
  // set start parameter values
  printInfo << "Setting start values for " << nmbPar << " parameters." << endl
	    << "    Parameter naming scheme is: V[rank index]_[IGJPCME][isobar spec]" << endl;
  unsigned int maxParNameLength = 0;
  for (unsigned int i = 0; i < nmbPar; ++i)
    if (L.parname(i).length() > maxParNameLength)
      maxParNameLength = L.parname(i).length();
  bool success = true;
  vector<bool> parIsFixed(nmbPar, false);  // memorizes state of variables; ROOT::Math::Minimizer has no corresponding accessor
  for (unsigned int i = 0; i < nmbPar; ++i) {
    double startVal = defaultStartValue;
    string parName = L.parname(i);
    if (hasStartVal){
      // look into start FitBin if we find a suitable parameter
      assert(startBin);
      startVal = startBin->getParameter(parName.c_str());
    }
    // check if parameter needs to be fixed because of threshold
    if ((L.parthreshold(i) == 0) || (L.parthreshold(i) < binHighMass)) {
      if (!minimizer->SetVariable(i, parName, startVal, defaultMinStep))
	success = false;
      if (startVal == 0) {
	cout << "    Read start value 0 for parameter " << parName << ". Setting default start value." << endl;
	startVal = defaultStartValue;
      }
      cout << "    Setting parameter [" << setw(3) << i << "] "
	   << setw(maxParNameLength) << parName << " = " << maxPrecisionAlign(startVal) << endl;
    } else {
      if (!minimizer->SetFixedVariable(i, parName, 0.))  // fix this parameter to 0
	success = false;
      cout << "    Fixing parameter  [" << setw(3) << i << "] "
	   << setw(maxParNameLength) << parName << " = 0" << endl;
      parIsFixed[i] = true;
    }
    if (!success) {
      printErr << "Something went wrong when setting minimizer parameters. Exiting." << endl;
      throw;
    }
  }
  // cleanup
  if(startValFile) {
    startValFile->Close();
    delete startValFile;
    startValFile = NULL;
  }

  // ---------------------------------------------------------------------------
  // find minimum of likelihood function
  printInfo << "Performing minimization." << endl;
  minimizer->SetMaxIterations(maxNmbOfIterations);
  success = minimizer->Minimize();
  if (!success)
    printWarn << "Minimization failed." << endl;
  if (runHesse) {
    printInfo << "Calculating Hessian matrix." << endl;
    //success = minimizer->Hesse();  // comes only with ROOT 5.24+
    if (!success)
      printWarn << "Calculation of Hessian matrix failed." << endl;
  }
  cout << "Minimization stopped after " << minimizer->NCalls() << " function calls." << endl
       << "    Total number of mimimzer parameters ................. " << minimizer->NDim()             << endl
       << "    Number of free minimzer parameters .................. " << minimizer->NFree()            << endl
       << "    Maximum allowed number of iterations ................ " << minimizer->MaxIterations()    << endl
       << "    Maximum allowed number of function calls ............ " << minimizer->MaxFunctionCalls() << endl
       << "    Minimizer status .................................... " << minimizer->Status()           << endl
       << "    Minimizer provides error and error matrix ........... " << minimizer->ProvidesError()    << endl
       << "    Minimizer has performed detailed error validation ... " << minimizer->IsValidError()     << endl
       << "    Estimated Distance to Minimum ....................... " << minimizer->Edm()              << endl
       << "    Statistical scale used for error calculation ........ " << minimizer->ErrorDef()         << endl
       << "    Minimizer strategy .................................. " << minimizer->Strategy()         << endl
       << "    Absolute tolerance .................................. " << minimizer->Tolerance()        << endl;

  // ---------------------------------------------------------------------------
  // print results
  printInfo << "Minimization result:" << endl;
  //double x[nmbPar];
  for (unsigned int i = 0; i< nmbPar; ++i) {
    //x[i] = minimizer->X()[i];  // step 1% from minimum for derivative
    cout << "    Parameter [" << setw(3) << i << "] "
	 << setw(maxParNameLength) << L.parname(i) << " = ";
    if (parIsFixed[i])
      cout << minimizer->X()[i] << " (fixed)" << endl;
    else {
      cout << setw(12) << maxPrecisionAlign(minimizer->X()[i]) << " +- "
	   << setw(12) << maxPrecisionAlign(minimizer->Errors()[i]);
      if (runMinos && (i == 156)) {  // does not work for all parameters
	double minosErrLow = 0;
	double minosErrUp  = 0;
	success = minimizer->GetMinosError(i, minosErrLow, minosErrUp);
	if (success)
	  cout << "    Minos: " << "[" << minosErrLow << ", +" << minosErrUp << "]" << endl;
      } else
	cout << endl;
    }
  }
  cout << "Number of calls to likelihood function FdF() ... " << L.ncalls() << endl
       << "Total time spent for likelihood calculation .... " << L.Ltime()  << endl
       << "Total time spent for normalization ............. " << L.Ntime()  << endl;

  // ---------------------------------------------------------------------------
  // write out result
  printInfo << "Writing result to '" << outFileName << "'." << endl;
  {
    // open output file and create tree for writing
    TFile* outFile = new TFile(outFileName.c_str(), "UPDATE");
    if (!outFile || outFile->IsZombie())
      printWarn << "Cannot open output file '" << outFileName << "'. No results will be written." << endl;
    else {
      // check whether output tree already exists
      TTree*      tree;
      TFitBin*    result    = new TFitBin();
      TFitResult* fitResult = new TFitResult();
      outFile->GetObject(valTreeName.c_str(), tree);
      if (!tree) {
	printInfo << "File '" << outFileName << "' is empty. Creating new tree '" << valTreeName << "' for PWA result." << endl;
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
	vector<string> waveNames = L.wavetitles();  // names of waves used in fit
	waveNames.push_back("flat");
	TMatrixT<double> fitParCovMatrix(nmbPar, nmbPar);  // covariance matrix of fit parameters
	for(unsigned int i = 0; i < nmbPar; ++i)
	  for(unsigned int j = 0; j < nmbPar; ++j)
	    fitParCovMatrix[i][j] = minimizer->CovMatrix(i,j);
	const unsigned int nmbWaves = waveNames.size();
	TCMatrix normIntegral(nmbWaves, nmbWaves);  // normalization integral over full phase space without acceptance
	TCMatrix accIntegral (nmbWaves, nmbWaves);  // normalization integral over full phase space with acceptance
	L.getIntCMatrix(normIntegral, accIntegral);
	const int normNmbEvents = useNorm ? 1 : L.nevents();  // number of events to normalize to

	cout << "Filling TFitResult:" << endl
	     << "    Number of fit parameters ............... " << nmbPar                        << endl
	     << "    Number of production amplitudes ........ " << prodAmps.size()               << endl
	     << "    Number of production amplitude names ... " << prodAmpNames.size()           << endl
	     << "    Number of wave names ................... " << waveNames.size()              << endl
	     << "    Number of cov. matrix indices .......... " << fitParCovMatrixIndices.size() << endl
	     << "    Dimension of covariance matrix ......... " << fitParCovMatrix.GetNrows()    << endl
	     << "    Dimension of normalization matrix ...... " << normIntegral.nrows()          << endl
	     << "    Dimension of acceptance matrix ......... " << accIntegral.nrows()           << endl;
	fitResult->fill(L.nevents(),
			normNmbEvents,
			binCenter,
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
	  vector<string> names;
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
	  const vector<string>& titles = L.wavetitles();
	  for(unsigned int i = 0; i < titles.size(); ++i)
	    waveTitles.push_back(TString(titles[i].c_str()));
	  waveTitles.push_back("flat");
	}
	// error matrix
	TMatrixD errMatrix(nmbPar, nmbPar);
	for(unsigned int i = 0; i < nmbPar; ++i)
	  for(unsigned int j = 0; j < nmbPar; ++j) {
	    errMatrix[i][j] = minimizer->CovMatrix(i,j);
	    //if(i==j)cout << "Sig"<< i << "=" << sqrt(errMatrix[i][i]) << endl;
	  }
	// normalixation integral and acceptance Matrix
	cout << " Setup integrals" << endl;
	unsigned int n = waveTitles.size();
	TCMatrix integralMatrix(n, n);
	TCMatrix accMatrix(n, n);
	L.getIntCMatrix(integralMatrix, accMatrix);
	//integralMatrix.Print();
	// representation of number of events depends on whether normalization was done
	int nmbEvt = useNorm ? 1 : L.nevents();
	
	cout << "Filling TFitBin:" << endl;
	cout << "    Number of fit parameters ........... " << nmbPar                 << endl;
	cout << "    Number of production amplitudes .... " << prodAmplitudes.size()  << endl;
	cout << "    Number of indices .................. " << indices.size()         << endl;
	cout << "    Number of wave names (with rank) ... " << waveNames.size()       << endl;
	cout << "    Number of wave titles (w/o rank) ... " << waveTitles.size()      << endl;
	cout << "    Dimension of error matrix .......... " << errMatrix.GetNrows()   << endl;
	cout << "    Dimension of integral matrix ....... " << integralMatrix.nrows() << endl;
	cout << "    Dimension of acceptance matrix ..... " << accMatrix.nrows()      << endl;
	result->fill(prodAmplitudes,
		     indices,
		     waveNames,
		     nmbEvt,
		     L.nevents(), // raw number of data events
		     binCenter,
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
