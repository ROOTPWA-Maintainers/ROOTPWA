///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
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
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      generates tree with log likelihood values and its derivatives
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <map>
#include <list>

#include "TROOT.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"

#include "reportingUtils.hpp"
#include "TPWALikelihood.h"


using namespace std;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{
  cerr << "usage:" << endl
       << progName
       << " -w wavelist [-z # val -d amplitude directory -o outfile -s seed -N -n normfile"
       << " [-a normfile] -r rank -q -h]" << endl
       << "    where:" << endl
       << "        -w file    path to wavelist file" << endl
       << "        -z #       number of likelihood values (default: 100)" << endl
       << "        -d dir     path to directory with decay amplitude files (default: '.')" << endl
       << "        -o file    path to output file (default: 'testLikelihood.root')" << endl
       << "        -s #       seed for random start values (default: 1234567)" << endl
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
	// force loading of predefined std::map dictionary
	gROOT->ProcessLine("#include <map>");

	const string progName          = argv[0];
	unsigned int nmbLikelihoodVals = 100;                    // number of likelihood values to calculate
	int          startValSeed      = 1234567;
	string       waveListFileName  = "";                     // wavelist filename
	string       ampDirName        = ".";                    // decay amplitude directory name
  string       outFileName       = "testLikelihood.root";  // output filename
	bool         useNormalizedAmps = false;                  // if true normalized amplitudes are used
	string       normIntFileName   = "";                     // file with normalization integrals
	string       accIntFileName    = "";                     // file with acceptance integrals
	unsigned int numbAccEvents     = 0;                      // number of events used for acceptance integrals
	unsigned int rank              = 1;                      // rank of spin-density matrix
	bool         quiet             = false;
  extern char* optarg;
  // extern int optind;
  int c;
  while ((c = getopt(argc, argv, "z:w:d:o:s:Nn:a:A:r:qh")) != -1)
	  switch (c) {
	  case 'z':
		  nmbLikelihoodVals = atoi(optarg);
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
	  case 's':
		  startValSeed = atoi(optarg);
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
  cout << "    # of values to calculate ..................... "  << nmbLikelihoodVals       << endl
       << "    wave list file ............................... '" << waveListFileName << "'" << endl
       << "    output file .................................. '" << outFileName      << "'" << endl
       << "    seed for random parameters ................... "  << startValSeed            << endl
       << "    use normalization ............................ "  << useNormalizedAmps       << endl
       << "        file with normalization integral ......... '" << normIntFileName  << "'" << endl
       << "        file with acceptance integral ............ '" << accIntFileName   << "'" << endl
       << "        number of acceptance norm. events ........ "  << numbAccEvents           << endl
       << "    rank of spin-density matrix .................. "  << rank                    << endl
       << "    quiet ........................................ "  << ((quiet) ? "yes" : "no") << endl;

  // setup likelihood function
  printInfo << "creating and setting up likelihood function" << endl;
  TPWALikelihood<complex<double> > L;
	if (quiet)
		L.setQuiet();
	L.useNormalizedAmps(useNormalizedAmps);
#ifdef USE_CUDA
	//L.enableCuda(false);
	L.enableCuda(true);
#endif  
	L.init(rank, waveListFileName, normIntFileName, accIntFileName, ampDirName, numbAccEvents);
	if (!quiet)
		cout << L << endl;

  // create output file and tree
	TFile& outFile = *TFile::Open(outFileName.c_str(), "RECREATE");
  TTree  tree("testLikelihoodTree", "testLikelihoodTree");

  // define tree branches
	Double_t            logLikelihood;
	map<string, double> derivatives;
  tree.Branch("logLikelihood", &logLikelihood, "logLikelihood/D");
  tree.Branch("derivatives",   &derivatives);

  // setup random number generator and production amplitude map
  TRandom3     random(startValSeed);
  list<string> prodAmpNames;  // needed, to ensure that parameter values do not depend on ordering
	const unsigned int nmbPar = L.NDim();
  for (unsigned int i = 0; i < nmbPar; ++i)
	  prodAmpNames.push_back(L.parName(i));
  prodAmpNames.sort();

	// calculate some likelihood values and derivatives
  for (unsigned int i = 0; i < nmbLikelihoodVals; ++i) {
		// generate random production amplitudes independent of parameter order
	  double prodAmps[nmbPar];
	  bool   success = true;
	  for (list<string>::const_iterator j = prodAmpNames.begin(); j != prodAmpNames.end(); ++j) {
		  int parIndex = -1;
		  for (unsigned int k = 0; k < nmbPar; ++k)
			  if (L.parName(k) == *j) {
				  parIndex = k;
				  break;
			  }
		  if (parIndex < 0) {
			  cout << "cannot find parameter '" << *j << "'. skipping." << endl;
			  success = false;
			  continue;
		  }
		  prodAmps[parIndex] = random.Uniform(0, 10);
		  // cout << "    [" << *j << "] = [" << parIndex << "] = " << prodAmps[parIndex] << endl;
	  }
	  if (success) {
		  logLikelihood = L.DoEval(prodAmps);
		  // cout << "[" << i << "] = " << logLikelihood << endl;
		  double gradient[nmbPar];
		  L.Gradient(prodAmps, gradient);
		  derivatives.clear();
		  for (unsigned int j = 0; j < nmbPar; ++j) {
			  const string parName = L.parName(j);
			  derivatives[parName] = gradient[j];
		  }
		  tree.Fill();
	  }
	}

	tree.Write();
	outFile.Close();

	return 0;
}
