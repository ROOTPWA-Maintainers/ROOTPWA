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
//fitting program for rootpwa

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
#include "TPWALikelihoodC.h"


using namespace std;
using namespace ROOT::Math;


//int lineno = 1; // global variables needed for lex (not understood)
//char *progname;


void usage(const char* progName,
	   const int   errCode = 0)
{
  cerr << "usage:" << endl
       << progName << " -l # -u # [-i -h -q -N] [-n normfile -s start values -r rank] -w wavelist [-o outfile] " << endl
       << "    where:" << endl
       << "        -l #       lower edge of mass bin [MeV/c^2]" << endl
       << "        -u #       upper edge of mass bin [MeV/c^2]" << endl
       << "        -w file    path to wavelist file" << endl
       << "        -c file    path to constraint config file" << endl
       << "        -o file    path to output file (default: 'fitresult.root')" << endl
       << "        -a file    path to acceptance file" << endl
       << "        -n file    path to normalization file" << endl
       << "        -N         use normalization of decay amplitudes (default: false)" << endl
       << "        -S file    path to file with start values" <<endl
       << "        -r #       rank of spin density matrix (default: 1)" <<endl
       << "        -q         run quietly" << endl
       << "        -h         print help" << endl
       << endl;
  exit(errCode);
}


int main(int argc, char** argv)
{
  // ---------------------------------------------------------------------------
  // internal parameters
  const TString      valTreeName        = "pwa";
  const TString      valBranchName      = "fitbin";
  const double       defaultStartValue  = 0.01;
  double             defaultMinStep     = 0.0005;
  const unsigned int maxNmbOfIterations = 20000;
  const bool         runHesse           = false;
  const bool         runMinos           = false;
  const string       minimizerType[2]   = {"Minuit2", "Migrad"};
//   const string       minimizerType[2]   = {"Minuit2", "Simplex"};
//   const string       minimizerType[2]   = {"Minuit2", "Combined"};
//   const string       minimizerType[2]   = {"Minuit2", "Scan"};
//   const string       minimizerType[2]   = {"Minuit2", "Fumili"};
//   const string       minimizerType[2]   = {"GSLMultiMin", "ConjugateFR"};
//   const string       minimizerType[2]   = {"GSLMultiMin", "ConjugatePR"};
//   const string       minimizerType[2]   = {"GSLMultiMin", "BFGS"};
//   const string       minimizerType[2]   = {"GSLMultiMin", "BFGS2"};
//   const string       minimizerType[2]   = {"GSLMultiMin", "SteepestDescent"};
//   const string       minimizerType[2]   = {"GSLMultiFit", ""};
//   const string       minimizerType[2]   = {"GSLSimAn", ""};

  // ---------------------------------------------------------------------------
  // parse command line options
  char* progName = argv[0];
  double binLowMass = 0;
  double binHighMass = 0;
  TString waveListFileName = ""; // wavelist filename
  TString waveConstraintConfigFileName = ""; // constraint config file name
  TString outFileName = "fitresult.root"; // output filename
  TString startValFileName = ""; // file with start values
  TString normFileName; // file with normalization integrals
  TString accFileName;  // file with acceptance integrals
  bool quiet = false;
  bool useStartVal = false; // are there start values?
  bool useNorm = false;
  unsigned int rank = 1; // rank
  extern char* optarg;
  // extern int optind;
  int c;
  while ((c = getopt(argc, argv, "l:u:iw:o:c:a:S:n:Nr:qh")) != -1)
    switch (c) {
    case 'N':
      useNorm = true;
      break;
    case 'l':
      binLowMass = atof(optarg);
      break;
    case 'u':
      binHighMass = atof(optarg);
      break;
    case 'o':
      outFileName = optarg;
      break;
    case 'c':
       waveConstraintConfigFileName = optarg;
      break;
    case 'S':
      startValFileName = optarg;
      useStartVal = 1;
      break;
    case 'n':
      normFileName = optarg;
      break;
    case 'a':
      accFileName = optarg;
      break;
    case 'w':
      waveListFileName = optarg;
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
  if (normFileName.Length() <= 1) {
    normFileName="norm.int";
    printWarn << "Using default normalization integral file '" << normFileName << "'." << endl;
  }
  if (accFileName.Length() <= 1) {
    accFileName="norm.int";
    printWarn << "Using default acceptance normalization integral file '" << accFileName << "'." << endl;
  }
  if (waveListFileName.Length() <= 1) {
    printErr << "No wavelist file specified! Aborting!" << endl;
    usage(progName, 1);
  }

  // ---------------------------------------------------------------------------
  // setup likelihood function
  printInfo << "Creating and setting up likelihood function" << endl;
  TPWALikelihoodC L;
  if(quiet)L.SetQuiet();
  L.UseNormalizedAmps(useNorm);

  L.Init(waveListFileName,
	 rank,
	 normFileName,
	 accFileName,
	 20000,
	 waveConstraintConfigFileName);
  L.LoadAmplitudes();

  const unsigned int nmbPar = L.NDim();
  
  // 12 parameters + flat
  //double x[13]={0.52707,0.21068,-0.604365,0.17596,-0.216668,-0.0990815,-0.348459,0.208961,0,0,0,0,0};
  //double x[13]; for(int i=0;i<13;++i)x[i]=0.001;
  //string a[13]={"a","b","c","d","e","f","g","h","i","j","k","l","flat"};
  //cout << L.DoEval(x) << endl;
  //return 0;

  // ---------------------------------------------------------------------------
  // setup minimizer
  printInfo << "Creating and setting up minimizer " << minimizerType[0] << " " << minimizerType[1] << endl;
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
  if (startValFileName.Length() > 2) {
    // open root file
    startValFile = TFile::Open(startValFileName, "READ");
    if (!startValFile || startValFile->IsZombie())
      printWarn << "Cannot open start value file '" << startValFileName << "'. Using default start values." << endl;
    else {
      // get tree with start values
      TTree* tree;
      startValFile->GetObject(valTreeName, tree);
      if (!tree)
	printWarn << "Cannot find start value tree '"<< valTreeName << "' in file '" << startValFileName << "'." << endl;
      else {
	startBin = new TFitBin();
	tree->SetBranchAddress(valBranchName, &startBin);
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
  for (unsigned int i = 0; i < nmbPar; ++i){
    if (L.parname(i).length() > maxParNameLength)
     maxParNameLength = L.parname(i).length();
    cout<<L.parname(i)<<endl;
  }
  bool success = true;
  vector<bool> parIsFixed(nmbPar, false);  // memorizes state of variables; ROOT::Math::Minimizer has no corresponding accessor
  for (unsigned int i = 0; i < nmbPar; ++i) {
    cout << i << ": " ;
    cout.flush();
    double startVal = defaultStartValue;
    string parName = L.parname(i);
    cout<<parName<< endl;
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
  //L.SetMaxSampDL(1000);
  // Pre-Minimizaton:
  //minimizer->SetMaxIterations(10);
  //minimizer->Minimize();
  //L.SetMaxSampDL(20000);
  // Minimizaton:
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

  // check derivatives:
  // first numerical:
  /*
    double h=1E-8;
    double dxNum[nmbPar];
    double dxAna[nmbPar];
    for(unsigned int i=0; i<nmbPar;++i){
    x[i]+=h;
    double L2=L.DoEval(x);
    x[i]-=2*h;
    double L1=L.DoEval(x);
    dxNum[i]=0.5*(L2-L1)/h;
    x[i]+=h;
    }
    double F;
    L.FdF(x,F,dxAna);
    for(unsigned int i=0; i<nmbPar;++i){
    cout<< "dL/d"<<i<<"(num)="<<dxNum[i]<<endl;
    cout<< "dL/d"<<i<<"(ana)="<<dxAna[i]<<endl;
    }
  */

  // ---------------------------------------------------------------------------
  // write out result
  printInfo << "Writing result to '" << outFileName << "'." << endl;
  {
    // open output file and create tree for writing
    TFile* outFile = new TFile(outFileName, "UPDATE");
    if (!outFile || outFile->IsZombie())
      printWarn << "Cannot open output file '" << outFileName << "'. No results will be written." << endl;
    else {
      // check whether output tree already exists
      TTree* tree;
      TFitBin* result = new TFitBin();
      outFile->GetObject(valTreeName, tree);
      if (!tree) {
	printInfo << "File '" << outFileName << "' is empty. Creating new tree for PWA result." << endl;
	tree = new TTree(valTreeName, valTreeName);
	tree->Branch(valBranchName, &result);
      } else
	tree->SetBranchAddress(valBranchName, &result);

      // get data structures to construct result TFitBin
      vector<TComplex> prodAmplitudes;  // production amplitudes
      vector<pair<int,int> > indices;  // indices for error matrix access
      vector<TString> waveNames;  // contains rank information 
      vector<complex<double> > V;
      vector<string> names;
      
      // error matrix of fit parameters
      TMatrixD errMatrix(nmbPar, nmbPar);
      for(unsigned int i = 0; i < nmbPar; ++i){
	for(unsigned int j = 0; j < nmbPar; ++j) {
	  errMatrix[i][j] = minimizer->CovMatrix(i,j);
	  //if(i==j)cout << "Sig"<< i << "=" << sqrt(errMatrix[i][i]) << endl;
	}
      }
      TMatrixD covMatrix(2,2); // covariants of amplitudes
      // These two matrices can differ in the case when there are
      // constraints present such that the fit parameters do not 
      // coincide with amplitude Re/Im anymore.
      
      L.buildCAmps(minimizer->X(), V, indices, names, errMatrix,covMatrix,true);
      
      // convert to TComplex;
      for (unsigned int i = 0; i < V.size(); ++i)
	prodAmplitudes.push_back(TComplex(V[i].real(), V[i].imag()));
      assert(indices.size() == prodAmplitudes.size());
      for(unsigned int i = 0; i < names.size(); ++i)
	waveNames.push_back(TString(names[i].c_str()));
      
      vector<TString> waveTitles;  // without rank
      {
	const vector<string>& titles = L.wavetitles();
	for(unsigned int i = 0; i < titles.size(); ++i)
	  waveTitles.push_back(TString(titles[i].c_str()));
	waveTitles.push_back("flat");
      }
      
      //       // normalization integral
//       integral norm = L.normInt();
      // integral and acceptance Matrix
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
		   covMatrix,
		   minimizer->MinValue(),
		   rank);
      //result->PrintParameters();

      // write result to file
      TString binName = valBranchName;
      binName += binLowMass;
      binName += "_";
      binName += binHighMass;
      binName.ReplaceAll(" ","");
      tree->Fill();
      tree->Write("", TObject::kOverwrite);
      outFile->Close();
    }
  }
  
  if (minimizer)
    delete minimizer;
  return 0;
}
  

// dummy function needed since we link to but do not use minuit.o
int mnparm(int, string, double, double, double, double)
{
  printErr << "this is impossible" << endl;
  throw "aFit";
  return 0;
}
