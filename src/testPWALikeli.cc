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
//      test likelihood calculation
//      program should be executed in amplitude directory
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include "TRandom3.h"

#include "TPWALikelihood.h"

#include "reportingUtils.hpp"


using namespace std;
using namespace rpwa;


int
main(int    argc,
     char** argv)
{
	using rpwa::cout;
	if (argc < 2) {
		printErr << "need wave list file name as argument. aborting." << endl;
		throw;
	}
	const bool         useNormalizedAmps = false;       // if true normalized amplitudes are used
	const unsigned int rank              = 2;           // rank of spin-density matrix
	const string       waveListFileName  = argv[1];     // wavelist filename
	const string       normIntFileName   = "norm.int";  // file with normalization integrals
	const string       accIntFileName    = "norm.int";  // file with acceptance integrals
	const string       ampDirName        = ".";         // decay amplitude directory name
	const unsigned int numbAccEvents     = 0;           // number of events used for acceptance integrals

	TPWALikelihood<complex<double> > L;
  L.useNormalizedAmps();
	L.useNormalizedAmps(useNormalizedAmps);
	L.init(rank, waveListFileName, normIntFileName, accIntFileName, ampDirName, numbAccEvents);
	printInfo << L;
	const unsigned int nmbPar = L.NDim();
	cout << "number of parameters = " << nmbPar << endl;

	// set random parameter values
	//double par[13]={0.52707,0.21068,-0.604365,0.17596,-0.216668,-0.0990815,-0.348459,0.208961,0,0,0,0,0};
  //double par[13]; for(int i=0;i<13;++i) par[i]=0.001;
  //string a[13]={"a","b","c","d","e","f","g","h","i","j","k","l","flat"};
	double par[nmbPar];
	gRandom->SetSeed(123456789);
  vector<unsigned int> parIndices  = L.orderedParIndices();
  for (unsigned int i = 0; i < parIndices.size(); ++i) {
	  const unsigned int parIndex = parIndices[i];
	  cout << "parameter[" << i << "] = " << L.parName(parIndex) << endl;
	  par[parIndex] = gRandom->Rndm();
  }

  {
	  // printInfo << "parameters:" << endl;
	  // for (unsigned int i = 0; i < nmbPar; ++i)
		//   cout << "par[" << i << "] = " << par[i] << endl;
	  printInfo << "production amplitudes:" << endl;
	  double                                          prodAmpFlat;
	  TPWALikelihood<complex<double> >::ampsArrayType prodAmps;
	  L.copyFromParArray(par, prodAmps, prodAmpFlat);
	  for (unsigned int iRank = 0; iRank < L.rank(); ++iRank)
		  for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			  for (unsigned int iWave = 0; iWave < L.nmbWaves(2 * iRefl - 1); ++iWave)
				  cout << "prodAmps[" << iRank << "][" << iRefl << "][" << iWave << "] = "
				       << maxPrecisionDouble(prodAmps[iRank][iRefl][iWave]) << endl;
	  printInfo << "parameters:" << endl;
	  double par2[nmbPar];
	  L.copyToParArray(prodAmps, prodAmpFlat, par2);
	  // for (unsigned int i = 0; i < nmbPar; ++i)
		//   cout << "par[" << i << "] = " << par[i] << endl;
	  for (unsigned int i = 0; i < nmbPar; ++i)
		  cout << "delta par[" << i << "] = " << maxPrecision(par[i] - par2[i]) << endl;
  }
	return 0;

  // calculate log likelihood
  const double logLikelihood = L.DoEval(par);
  printInfo << "log likelihood = " << maxPrecision(logLikelihood) << endl;

  // check derivatives
  // calculate numerical derivatives
  double numericalDerivates[nmbPar];
  {
	  const double delta = 1E-7;
	  for (unsigned int i = 0; i < nmbPar; ++i) {
		  const double oldParVal = par[i];
		  par[i] += delta;
		  const double logLikelihood2 = L.DoEval(par);
		  numericalDerivates[i] = (logLikelihood2 - logLikelihood) / delta;
		  par[i] = oldParVal;
	  }
  }
  // calculate analytical derivatives
  double analyticalDerivates[nmbPar];
  // double F;
  // L.FdF(par, F, analyticalDerivates);
  L.Gradient(par, analyticalDerivates);
  double maxDelta = 0;
  printInfo << "derivatives:" << endl;
  for (unsigned int i = 0; i < parIndices.size(); ++i) {
	  const unsigned int parIndex = parIndices[i];
	  const double       delta    = numericalDerivates[parIndex] - analyticalDerivates[parIndex];
	  if (abs(delta) > maxDelta)
		  maxDelta = abs(delta);
	  cout << "    dL / dPar[" << setw(nmbOfDigits(nmbPar - 1)) << parIndex << "]: numerical = "
	       << maxPrecisionAlign(numericalDerivates[parIndex])
	       << " vs. analytical = " << maxPrecisionAlign(analyticalDerivates[parIndex])
	       << "; (numer. - analyt.) = "
	       << maxPrecisionAlign(delta) << endl;
  }
  cout << "    maximum deviation = " << maxPrecision(maxDelta) << endl;
  
  return 0;
}
