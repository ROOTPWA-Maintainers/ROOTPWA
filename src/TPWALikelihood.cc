///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert, Boris Grube
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
//      Implementation of class TPWALikelihood
//      see TPWALikelihood.hh for details
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//      Boris Grube          Universe Cluster Munich
//
//
//-----------------------------------------------------------


#include <iomanip>
#include <fstream>
#include <cassert>

#include "TString.h"
#include "TSystem.h"
#include "TCMatrix.h"
#include "TStopwatch.h"

#include "utilities.h"


//#define USE_FDF


using namespace std;


template <typename T>
TPWALikelihood<T>::TPWALikelihood()
  : _nmbEvents(0),
    _rank(1),
    _nmbWaves(0),
    _nmbPars(0),
    _Ltime(0),
    _Ntime(0),
    _debug(true),
    _useNormalizedAmps(true)
{
  _nmbWavesRefl[0] = 0;
  _nmbWavesRefl[1] = 0;
  for (unsigned int i = 0; i < NMB_FUNCTIONCALLENUM; ++i)
    _nmbCalls[i] = 0;
  printInfo << "using DoEval()!!!" << endl;
}


template <typename T>
TPWALikelihood<T>::~TPWALikelihood()
{
  clearCache();
}


template <typename T>
void
TPWALikelihood<T>::FdF(const double* par,             // parameter array; reduced by rank conditions
		       double&       funcVal,         // function value
		       double*       gradient) const  // array of derivatives
{
  if (0) {  // run FdFX
    FdFX(par, funcVal, gradient);
    return;
  }

  ++(_nmbCalls[FDF]);

  // log consumed time
  TStopwatch timer;
  timer.Start(true);

  // copy arguments into parameter cache
  for (unsigned int i = 0; i < _nmbPars; ++i)
    _parCache[i] = par[i];

  // build complex production amplitudes from function parameters taking into account rank restrictions
  T            prodAmpFlat;
  complex<T>** prodAmps = NULL;
// !!! a possible optimization would be to create the array only once
// !!! and reuse it for each function call
  copyFromParArray(par, prodAmps, prodAmpFlat);
  const T prodAmpFlat2 = prodAmpFlat * prodAmpFlat;

  // array with likelihood derivative with resp to real and imaginary
  // parts of the production amplitudes although stored as complex
  // values, the dL are _not_ well defined complex numbers!
  T                   derivativeFlat = 0;
  vector2(complex<T>) derivatives(_rank, vector<complex<T> >(_nmbWaves, 0));

  // loop over events and calculate first term of log likelihood and derivatives with respect to parameters
  T logLikelihood = 0;
  for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
    T                   l = 0;                                        // likelihood for this event
    vector2(complex<T>) d(_rank, vector<complex<T> >(_nmbWaves, 0));  // likelihood derivative for this event
    for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks and reflectivities
      complex<T> ampPos = 0;  // positive reflectivity amplitude for this rank
      complex<T> ampNeg = 0;  // negative reflectivity amplitude for this rank
      for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {  // outer loop over waves
	// compute likelihood term
	const complex<T> amp  = prodAmps[iRank][iWave] * complex<T>(_decayAmps[iWave][iEvt]);
	const int        refl = _waveRefl[iWave];
	if (refl == -1)
	  ampNeg += amp;
	else
	  ampPos += amp;
	// compute derivatives
	for (unsigned int jWave = 0; jWave < _nmbWaves; ++jWave)  // inner loop over waves
	  // sum up amplitudes for current rank and only for waves with same reflectivity
	  if (refl == _waveRefl[jWave])
	    d[iRank][jWave] += amp;
      } // end outer loop over waves
      l += norm(ampPos);
      l += norm(ampNeg);
      assert(l >= 0);
      // loop again over waves for current rank
      for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave)
	d[iRank][iWave] *= conj(complex<T>(_decayAmps[iWave][iEvt]));
    }  // end loop over rank
    l             += prodAmpFlat2;
    logLikelihood -= TMath::Log(l);  // accumulate log likelihood
    // incorporate factor 2 / sigma
    const T factor = 2. / l;
    for (unsigned int iRank = 0; iRank < _rank; ++iRank)
      for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave)
 	derivatives[iRank][iWave] -= factor * d[iRank][iWave];
    derivativeFlat -= factor * prodAmpFlat;
  }  // end loop over events
 
  // log consumed time
  const double t1 = timer.RealTime();
  timer.Start(true);

  // compute normalization term of log likelihood and its derivatives with respect to parameters
  complex<T> normFactor  = 0;
  const T    nmbEvt      = (_useNormalizedAmps) ? 1 : (T)_nmbEvents;
  const T    twiceNmbEvt = 2 * nmbEvt;
  for (unsigned int iRank = 0; iRank < _rank; ++iRank)
    for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {    // outer loop over waves
      complex<T> normFactorDeriv = 0;
      for (unsigned int jWave = 0; jWave < _nmbWaves; ++jWave) {  // inner loop over waves
	if (_waveRefl[iWave] != _waveRefl[jWave])  // make sure that waves i and j have the same reflectivity
	  continue;
 	const complex<T> I = complex<T>(_accMatrix.element(iWave, jWave));
 	normFactor      += (prodAmps[iRank][iWave] * conj(prodAmps[iRank][jWave])) * I;
	normFactorDeriv += prodAmps[iRank][jWave] * conj(I);
      }
      derivatives[iRank][iWave] += normFactorDeriv * twiceNmbEvt;  // account for 2 * nmbEvents
    }
  // take care of flat wave
  normFactor.real() += prodAmpFlat2;
  derivativeFlat    += twiceNmbEvt * prodAmpFlat;

  // cleanup memory
  unsigned int prodAmpsDim[2];
  prodAmpsDim[0] = _rank;
  prodAmpsDim[1] = _nmbWaves;
  delete2DArray<complex<T> >(prodAmps, prodAmpsDim);

  // sort derivative results into output array and cache
  copyToParArray(derivatives, derivativeFlat, gradient);
  copyToParArray(derivatives, derivativeFlat, &(*(_derivCache.begin())));

  // log consumed time
  const double t2 = timer.RealTime();
  timer.Stop();
  _Ltime += t1;
  _Ntime += t2;
  
  if (_debug)
    printInfo << "log likelihood =  " << maxPrecisionAlign(logLikelihood) << ", "
	      << "normalization =  " << maxPrecisionAlign(normFactor.real()) << ", "
	      << "normalized likelihood = " << maxPrecisionAlign(logLikelihood + nmbEvt * normFactor.real()) << endl
	      << "    Time for likelihood = " << t1 << ", time for normalization = " << t2 << endl;

  // return likelihood value
  funcVal = logLikelihood + nmbEvt * normFactor.real();

  if (0) {  // compare to FdFX
    double fX;
    double dfX[_nmbPars];
    FdFX(par, fX, dfX);
    double maxDelta = 0;
    for (unsigned int i = 0; i < _nmbPars; ++i) {
      const double delta = (dfX[i] - gradient[i]) / gradient[i];
//       printInfo << "dfX[" << i << "] = " << maxPrecisionAlign(dfX[i])
// 		<< ", gradient[" << i << "] = " << maxPrecisionAlign(gradient[i])
// 		<< ", delta grad f[" << i << "] = " << maxPrecisionAlign(delta) << endl;
      if (fabs(delta) > fabs(maxDelta))
	maxDelta = delta;
    }
//     printInfo << "delta f        = " << maxPrecision(fX - logLikelihood) << endl;
//     printInfo << "delta f        = " << maxPrecision(fX - normFactor.real()) << endl;
    printInfo << "delta f = " << maxPrecisionAlign((fX - funcVal) / funcVal)
	       << ", !max delta grad f = " << maxPrecisionAlign(maxDelta) << endl;
  }
}


template <typename T>
void
TPWALikelihood<T>::FdFX(const double* par,             // parameter array; reduced by rank conditions
			double&       funcVal,         // function value
			double*       gradient) const  // array of derivatives
{
  ++(_nmbCalls[FDF]);
  
  // log consumed time
  TStopwatch timer;
  timer.Start(true);
  
  // copy arguments into parameter cache
  for (unsigned int i = 0; i < _nmbPars; ++i)
    _parCache[i] = par[i];

  // build complex production amplitudes from function parameters taking into account rank restrictions
  T                   prodAmpFlat;
  vector3(complex<T>) prodAmps;
  copyFromParArrayX(par, prodAmps, prodAmpFlat);
  const T prodAmpFlat2 = prodAmpFlat * prodAmpFlat;

  // array with likelihood derivative with resp to real and imaginary
  // parts of the production amplitudes although stored as complex
  // values, the dL are _not_ well defined complex numbers!
   T                   derivativeFlat = 0;
   vector3(complex<T>) derivatives(_rank, vector2(complex<T>)(2));
   for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
     derivatives[iRank][0].resize(_nmbWavesRefl[0], 0);
     derivatives[iRank][1].resize(_nmbWavesRefl[1], 0);
   }
  //!!! possible optimization: dL/dIm[Par_i] = 0 for all i

  // loop over events and calculate first term of log likelihood and derivatives with respect to parameters
  T logLikelihood = 0;
  for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
    T                   l = 0;                             // likelihood for this event
    vector3(complex<T>) d(_rank, vector2(complex<T>)(2));  // likelihood derivative for this event [rank][reflectivity][wave index]
    for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
      d[iRank][0].resize(_nmbWavesRefl[0], 0);
      d[iRank][1].resize(_nmbWavesRefl[1], 0);
    }
    for (unsigned int iRank = 0; iRank < _rank; ++iRank) {                    // incoherent sum over ranks
      complex<T> ampProdSum[2] = {0, 0};  // amplitude sum for negative/positive reflectivity for this rank
      for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {                      // incoherent sum over reflectivities
	for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
	  // compute likelihood term
	  const complex<T> amp = prodAmps[iRank][iRefl][iWave] * complex<T>(_decayAmpsX[iRefl][iWave][iEvt]);
	  ampProdSum[iRefl] += amp;
	  // compute derivatives
	  for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave)  // inner loop over waves
	    // sum up amplitudes for current rank and only for waves with same reflectivity
	    d[iRank][iRefl][jWave] += amp;
	}
      }
      l += norm(ampProdSum[1]);
      l += norm(ampProdSum[0]);
      assert(l >= 0);
      // loop again over waves for current rank
      for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
	for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
	  d[iRank][iRefl][iWave] *= conj(complex<T>(_decayAmpsX[iRefl][iWave][iEvt]));
    }  // end loop over rank
    l             += prodAmpFlat2;
    logLikelihood -= TMath::Log(l);  // accumulate log likelihood
    // incorporate factor 2 / sigma
    const T factor = 2. / l;
    for (unsigned int iRank = 0; iRank < _rank; ++iRank)
      for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
	for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
  	  derivatives[iRank][iRefl][iWave] -= factor * d[iRank][iRefl][iWave];
    derivativeFlat -= factor * prodAmpFlat;
  }  // end loop over events
 
  // log consumed time
  const double t1 = timer.RealTime();
  timer.Start(true);

  // compute normalization term of log likelihood and its derivatives with respect to parameters
  complex<T> normFactor  = 0;
  const T    nmbEvt      = (_useNormalizedAmps) ? 1 : (T)_nmbEvents;
  const T    twiceNmbEvt = 2 * nmbEvt;
  for (unsigned int iRank = 0; iRank < _rank; ++iRank)
    for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
      for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
	complex<T> normFactorDeriv = 0;
	for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave) { // inner loop over waves
	  const complex<T> I = complex<T>(_accMatrixX[iRefl][iWave][iRefl][jWave]);
 	  normFactor      += (prodAmps[iRank][iRefl][iWave] * conj(prodAmps[iRank][iRefl][jWave])) * I;
	  normFactorDeriv += prodAmps[iRank][iRefl][jWave] * conj(I);
	}
	derivatives[iRank][iRefl][iWave] += normFactorDeriv * twiceNmbEvt;  // account for 2 * nmbEvents
      }
  // take care of flat wave
  normFactor.real() += prodAmpFlat2;
  derivativeFlat    += twiceNmbEvt * prodAmpFlat;

  // sort derivative results into output array and cache
  copyToParArrayX(derivatives, derivativeFlat, gradient);
  copyToParArrayX(derivatives, derivativeFlat, &(*(_derivCache.begin())));

  // log consumed time
  const double t2 = timer.RealTime();
  timer.Stop();
  _Ltime += t1;
  _Ntime += t2;
  
  if (_debug)
    printInfo << "log likelihood =  " << maxPrecisionAlign(logLikelihood) << ", "
	      << "normalization =  " << maxPrecisionAlign(normFactor.real()) << ", "
	      << "normalized likelihood = " << maxPrecisionAlign(logLikelihood + nmbEvt * normFactor.real()) << endl
	      << "    Time for likelihood = " << t1 << ", time for normalization = " << t2 << endl;

  // return likelihood value
  funcVal = logLikelihood + nmbEvt * normFactor.real();
}


template <typename T>
double
TPWALikelihood<T>::DoEval(const double* par) const
{
  ++(_nmbCalls[DOEVAL]);

#ifdef USE_FDF

  // call FdF
  double logLikelihood;
  double gradient[_nmbPars];
  FdF(par, logLikelihood, gradient);
  return logLikelihood;

#else

  // log consumed time
  TStopwatch timer;
  timer.Start(true);

  // build complex production amplitudes from function parameters taking into account rank restrictions
  T            prodAmpFlat;
  complex<T>** prodAmps = NULL;
  copyFromParArray(par, prodAmps, prodAmpFlat);
  const T prodAmpFlat2 = prodAmpFlat * prodAmpFlat;

  // compute first term of log likelihood: sum over data
  T logLikelihood = 0;
  for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
    T l = 0;  // likelihood for this event
    for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks and reflectivities
      complex<T> ampPos = 0;  // positive reflectivity amplitude for this rank
      complex<T> ampNeg = 0;  // negative reflectivity amplitude for this rank
      for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {  // loop over waves
	const complex<T> amp  = prodAmps[iRank][iWave] * complex<T>(_decayAmps[iWave][iEvt]);
	if (_waveRefl[iWave] == -1)
	  ampNeg += amp;
	else
	  ampPos += amp;
      }
      l += norm(ampPos);
      l += norm(ampNeg);
      assert(l >= 0);
    }  // end loop over rank
    l             += prodAmpFlat2;
    logLikelihood -= TMath::Log(l);  // accumulate log likelihood
  }  // end loop over events
 
  // log consumed time
  const double t1 = timer.RealTime();
  timer.Start(true);

  // compute second term of log likelihood: normalization
  complex<T> normFactor = 0;
  const T    nmbEvt     = (_useNormalizedAmps) ? 1 : (T)_nmbEvents;
  for (unsigned int iRank = 0; iRank < _rank; ++iRank)
    for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave)      // outer loop
      for (unsigned int jWave = 0; jWave < _nmbWaves; ++jWave) {  // inner loop
	if (_waveRefl[iWave] != _waveRefl[jWave])  // make sure that waves i and j have the same reflectivity
	  continue;
 	const complex<T> I = complex<T>(_accMatrix.element(iWave, jWave));
	normFactor += (prodAmps[iRank][iWave] * conj(prodAmps[iRank][jWave])) * I;
      }
  // take care of flat wave
  normFactor.real() += prodAmpFlat2;

  // cleanup memory
  unsigned int prodAmpsDim[2];
  prodAmpsDim[0] = _rank;
  prodAmpsDim[1] = _nmbWaves;
  delete2DArray<complex<T> >(prodAmps, prodAmpsDim);

  // log consumed time
  const double t2 = timer.RealTime();
  timer.Stop();
  _Ltime += t1;
  _Ntime += t2;
  
  if (_debug)
    printInfo << "Log likelihood =  " << maxPrecisionAlign(logLikelihood) << ", "
	      << "normalization =  " << maxPrecisionAlign(normFactor.real()) << ", "
	      << "normalized log likelihood = " << maxPrecisionAlign(logLikelihood + nmbEvt * normFactor.real()) << endl
	      << "    Time for likelihood = " << t1 << ", time for normalization = " << t2 << endl;

  const double funcVal = logLikelihood + nmbEvt * normFactor.real();

  if (0) {  // compare to FdF
    double f;
    double df[_nmbPars];
    FdF(par, f, df);
    const double funcVal = logLikelihood + nmbEvt * normFactor.real();
    printInfo << "delta f = " << maxPrecision((f - funcVal) / funcVal) << endl;
  }

  // return log likelihood value
  return funcVal;

#endif  // USE_FDF
}


template <typename T>
double
TPWALikelihood<T>::DoDerivative(const double* par,
				unsigned int  derivativeIndex) const
{
  ++(_nmbCalls[DODERIVATIVE]);

  // check whether parameter is in cache
  bool samePar = true;
  for (unsigned int i = 0; i < _nmbPars; ++i)
    if (_parCache[i] != par[i]) {
      samePar = false;
      break;
    }
  if (samePar) {
    //cout << "using cached derivative! " << endl;
    return _derivCache[derivativeIndex];
  }
  // call FdF
  double logLikelihood;
  double gradient[_nmbPars];
  FdF(par, logLikelihood, gradient);
  return gradient[derivativeIndex];
}


// calculate derivatives with respect to parameters
template <typename T>
void
TPWALikelihood<T>::Gradient(const double* par,             // parameter array; reduced by rank conditions
			    double*       gradient) const  // array of derivatives
{
  ++(_nmbCalls[GRADIENT]);
  
#ifdef USE_FDF

  // check whether parameter is in cache
  bool samePar = true;
  for (unsigned int i = 0; i < _nmbPars; ++i)
    if (_parCache[i] != par[i]) {
      samePar = false;
      break;
    }
  if (samePar) {
    for (unsigned int i = 0; i < _nmbPars ; ++i)
      gradient[i] = _derivCache[i];
    return;
  }
  // call FdF
  double logLikelihood;
  FdF(par, logLikelihood, gradient);

#else

  // log consumed time
  TStopwatch timer;
  timer.Start(true);

  // copy arguments into parameter cache
  for (unsigned int i = 0; i < _nmbPars; ++i)
    _parCache[i] = par[i];

  // build complex production amplitudes from function parameters taking into account rank restrictions
  T            prodAmpFlat;
  complex<T>** prodAmps = NULL;
  copyFromParArray(par, prodAmps, prodAmpFlat);
  const T prodAmpFlat2 = prodAmpFlat * prodAmpFlat;

  // array with likelihood derivative with resp to real and imaginary
  // parts of the production amplitudes although stored as complex
  // values, the dL are _not_ well defined complex numbers!
  T                   derivativeFlat = 0;
  vector2(complex<T>) derivatives(_rank, vector<complex<T> >(_nmbWaves, 0));

  // compute derivative for first term of log likelihood
  for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
    T                   l = 0;                                        // likelihood for this event
    vector2(complex<T>) d(_rank, vector<complex<T> >(_nmbWaves, 0));  // likelihood derivative for this event
    for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks and reflectivities
      complex<T> ampPos = 0;  // positive reflectivity amplitude for this rank
      complex<T> ampNeg = 0;  // negative reflectivity amplitude for this rank
      for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {  // outer loop over waves
	const complex<T> amp  = prodAmps[iRank][iWave] * complex<T>(_decayAmps[iWave][iEvt]);
	const int        refl = _waveRefl[iWave];
	if (refl == -1)
	  ampNeg += amp;
	else
	  ampPos += amp;
	for (unsigned int jWave = 0; jWave < _nmbWaves; ++jWave)  // inner loop over waves
	  // sum up amplitudes for current rank and only for waves with same reflectivity
	  if (refl == _waveRefl[jWave])
	    d[iRank][jWave] += amp;
      } // end outer loop over waves
      l += norm(ampPos);
      l += norm(ampNeg);
      assert(l >= 0);
      // loop again over waves for current rank
      for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave)
	d[iRank][iWave] *= conj(complex<T>(_decayAmps[iWave][iEvt]));
    }  // end loop over rank
    l += prodAmpFlat2;
    // incorporate factor 2 / sigma
    const T factor = 2. / l;
    for (unsigned int iRank = 0; iRank < _rank; ++iRank)
      for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave)
	derivatives[iRank][iWave] -= factor * d[iRank][iWave];
    derivativeFlat -= factor * prodAmpFlat;
  }  // end loop over events
 
  // log consumed time
  const double t1 = timer.RealTime();
  timer.Start(true);

  // normalize derivatives
  const T nmbEvt      = (_useNormalizedAmps) ? 1 : (T)_nmbEvents;
  const T twiceNmbEvt = 2 * nmbEvt;
  for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
    for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {    // outer loop
      complex<T> normFactorDeriv = 0;
      for (unsigned int jWave = 0; jWave < _nmbWaves; ++jWave) {  // inner loop
	if (_waveRefl[iWave] != _waveRefl[jWave])  // make sure that waves i and j have the same reflectivity
	  continue;
 	const complex<T> I = complex<T>(_accMatrix.element(iWave, jWave));
	normFactorDeriv += prodAmps[iRank][jWave] * conj(I);
      }
      derivatives[iRank][iWave] += normFactorDeriv * twiceNmbEvt;  // account for 2 * nmbEvents
    }
  }
  // take care of flat wave
  derivativeFlat += twiceNmbEvt * prodAmpFlat;

  // cleanup memory
  unsigned int prodAmpsDim[2];
  prodAmpsDim[0] = _rank;
  prodAmpsDim[1] = _nmbWaves;
  delete2DArray<complex<T> >(prodAmps, prodAmpsDim);

  // return gradient values
  copyToParArray(derivatives, derivativeFlat, gradient);
  copyToParArray(derivatives, derivativeFlat, &(*(_derivCache.begin())));

  // log consumed time
  const double t2 = timer.RealTime();
  timer.Stop();
  _Ltime += t1;
  _Ntime += t2;

  if (0) {  // compare to FdF
    double f;
    double df[_nmbPars];
    FdF(par, f, df);
    double maxDelta = 0;
    for (unsigned int i = 0; i < _nmbPars; ++i) {
      const double delta = (df[i] - gradient[i]) / gradient[i];
      if (fabs(delta) > fabs(maxDelta))
	maxDelta = delta;
    }
    printInfo << "max delta df = " << maxPrecision(maxDelta) << endl;
  }

#endif  // USE_FDF
}


template <typename T>
void 
TPWALikelihood<T>::init(const unsigned int rank,
			const std::string& waveListFileName,
			const std::string& normIntFileName,
			const std::string& accIntFileName,
			const std::string& ampDirName)
{
  readWaveList(waveListFileName);
  buildParDataStruct(rank);
  readIntegrals(normIntFileName, accIntFileName);
  readDecayAmplitudes(ampDirName);
}


template <typename T>
void 
TPWALikelihood<T>::readWaveList(const string& waveListFileName)
{
  printInfo << "reading amplitude names and thresholds from wave list file '" << waveListFileName << "'." << endl;
  ifstream waveListFile(waveListFileName.c_str());
  if (!waveListFile) {
    printErr << "cannot open file '" << waveListFileName << "'. aborting." << endl;
    throw;
  }
  string line;
  int    lineNmb = 0;
  _waveNamesX.resize(2);
  while (getline(waveListFile, line)) {
    if (line[0] == '#')  // comments start with #
      continue;
    stringstream lineStream;
    lineStream.str(line);
    string waveName;
    if (lineStream >> waveName) {
      double threshold;
      // !!! it would be safer to make the threshold value in the wave list file mandatory
      if (!(lineStream >> threshold))
	threshold = 0;
      if (_debug)
	cout << "    reading line " << setw(3) << lineNmb + 1 << ": " << waveName << ",  threshold = " << setw(4) << threshold << " MeV/c^2" << endl;
      if (getReflectivity(waveName) > 0) {
	++_nmbWavesRefl[1];  // positive reflectivity
	_waveNamesX[1].push_back(waveName);
      } else {
	++_nmbWavesRefl[0];  // negative reflectivity
	_waveNamesX[0].push_back(waveName);
      }
      _waveNames.push_back(waveName);
      _waveThresholds.push_back(threshold);
    } else
      printWarn << "cannot parse line '" << line << "' in wave list file '" << waveListFileName << "'." << endl;
    ++lineNmb;
  }
  waveListFile.close();
  printInfo << "read " << lineNmb << " lines from wave list file '" << waveListFileName << "'." << endl;
  // calculate dimension
  _nmbWaves=_waveNames.size();
}


template <typename T>
void
TPWALikelihood<T>::buildParDataStruct(const unsigned int rank)
{
  if ((_waveNames.size() == 0) || (_waveThresholds.size() == 0)) {
    printErr << "no wave info. was readWaveList() executed successfully? aborting.";
    throw;
  }
  _rank = rank;
  // calculate dimension of function taking into account rank restrictions and flat wave
  _nmbPars = 0;
  for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
    const int nmbProdAmpsPos  = _nmbWavesRefl[1] - iRank;  // number non-zero production amplitudes in this rank with positive reflectivity
    int       nmbProdAmpsPosC = nmbProdAmpsPos - 1;        // number of complex-valued production amplitudes in this rank with positive reflectivity
    if (nmbProdAmpsPosC < 0)
      nmbProdAmpsPosC = 0;
    const int nmbProdAmpsNeg  = _nmbWavesRefl[0] - iRank;  // number non-zero production amplitudes in this rank with negative reflectivity
    int       nmbProdAmpsNegC = nmbProdAmpsNeg - 1;        // number of complex-valued production amplitudes in this rank with negative reflectivity
    if (nmbProdAmpsNegC < 0)
      nmbProdAmpsNegC = 0;
    _nmbPars += 2 * (nmbProdAmpsPosC + nmbProdAmpsNegC);  // 2 parameters for each complex production amplitude
    // 1 real production amplitude per rank and reflectivity
    if (nmbProdAmpsPos > 0)
      ++_nmbPars;  
    if (nmbProdAmpsNeg > 0)
      ++_nmbPars;
  }
  _nmbPars += 1;  // additonal flat wave
  printInfo << "dimension of likelihood function is " << _nmbPars << "." << endl;
  _parNames.resize     (_nmbPars,      "");
  _parThresholds.resize(_nmbPars,      0);
  _waveRefl.resize     (_nmbWaves + 1, 0);
  _parCache.resize     (_nmbPars,      0);
  _derivCache.resize   (_nmbPars,      0);
  _prodAmpToFuncParMapX.resize(_rank, vector2(protect__(pair<int, int>))(2));
  // build parameter names
  unsigned int parIndex = 0;
  for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
    _prodAmpToFuncParMapX[iRank][0].resize(_nmbWavesRefl[0], make_pair(0, 0));
    _prodAmpToFuncParMapX[iRank][1].resize(_nmbWavesRefl[1], make_pair(0, 0));
    unsigned int iWaveRefl[2] = {0, 0};  // indices for negative and positive reflectivity waves
    for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {
      _waveRefl[iWave] = getReflectivity(_waveNames[iWave]);
      const unsigned int refl = (_waveRefl[iWave] > 0) ? 1 : 0;
      if (iWaveRefl[refl] < iRank) {  // production amplitude is zero
	_prodAmpToFuncParMapX[iRank][refl][iWaveRefl[refl]] = make_pair(-1, -1);
// 	cout << _waveNames[iWave] << " [" << iRank << ", " << iWave << "] -> "
// 	     << " [" << iRank << ", " << refl << ", " << iWaveRefl[refl] << "] -> "
// 	     << _prodAmpToFuncParMapX[iRank][refl][iWaveRefl[refl]] << endl;
	++iWaveRefl[refl];
	continue;
      } else if (iWaveRefl[refl] == iRank) {  // production amplitude is real
	ostringstream parName;
	parName << "V" << iRank << "_" << _waveNames[iWave] << "_RE";
	_parNames     [parIndex] = parName.str();
	_parThresholds[parIndex] = _waveThresholds[iWave];
	_prodAmpToFuncParMapX[iRank][refl][iWaveRefl[refl]] = make_pair(parIndex, -1);
	++parIndex;
      } else {  // production amplitude is complex
	ostringstream parName;
	parName << "V" << iRank << "_" << _waveNames[iWave];
	_parNames     [parIndex] = parName.str() + "_RE";
	_parThresholds[parIndex] = _waveThresholds[iWave];
	++parIndex;
	_parNames     [parIndex] = parName.str() + "_IM";
	_parThresholds[parIndex] = _waveThresholds[iWave];
	_prodAmpToFuncParMapX[iRank][refl][iWaveRefl[refl]] = make_pair(parIndex - 1, parIndex);
	++parIndex;
      }
//       cout << _waveNames[iWave] << " [" << iRank << ", " << iWave << "] -> "
// 	   << " [" << iRank << ", " << refl << ", " << iWaveRefl[refl] << "] -> "
// 	   << _prodAmpToFuncParMapX[iRank][refl][iWaveRefl[refl]] << endl;
      ++iWaveRefl[refl];
    }  // end loop over waves
  }  // end loop over rank
  // flat wave
  _parNames     [parIndex]  = "V_flat";
  _parThresholds[parIndex]  = 0;
  _waveRefl     [_nmbWaves] = 1;

  if (0) {  // test copyFrom/ToParArrayX routines
    double* testPar = new double[_nmbPars];
    for (unsigned int i = 0; i < _nmbPars; ++i)
      testPar[i] = i + 1;
    vector3(complex<T>) testVal;
    T                   testFlatVal;
    copyFromParArrayX(testPar, testVal, testFlatVal);
    copyToParArrayX(testVal, testFlatVal, testPar);
    for (unsigned int i = 0; i < _nmbPars; ++i)
      cout << "!!!TEST " << setw(3) << i << ": " << testPar[i] << endl;
    delete[] testPar;
  }
}


// returns integral matrix reordered according to _waveNames array
template <typename T>
matrix<complex<double> >
TPWALikelihood<T>::reorderedIntegralMatrix(integral& integral) const
{
  // get original matrix and list of wave names
  const matrix<complex<double> > intMatrix    = integral.mat();
  const list<string>             intWaveNames = integral.files();
  // build index lookup-table
  vector<unsigned int> indexLookUp;  // lookup table: wave index -> index in normalization integral
  indexLookUp.resize(_nmbWaves, 0);
  for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {
    if (find(intWaveNames.begin(), intWaveNames.end(), _waveNames[iWave]) == intWaveNames.end()) {
      printErr << "wave " << _waveNames[iWave] << " is not in integral. aborting." << endl;
      throw;
    }
    indexLookUp[iWave] = integral.index(_waveNames[iWave]);  // not a const method!!! should be fixed
    if (_debug)
      cout << "    mapping wave [" << setw(3) << iWave << "] '" << _waveNames[iWave] << "' "
	   << "to index " << setw(3) << indexLookUp[iWave] << " in integral." << endl;
  }
  // create reordered matrix
  matrix<complex<double> > reorderedMatrix(intMatrix.nrows(), intMatrix.ncols());
  for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave)
    for (unsigned int jWave = 0; jWave < _nmbWaves; ++jWave)
      reorderedMatrix.element(iWave, jWave) = intMatrix.element(indexLookUp[iWave], indexLookUp[jWave]);
  return reorderedMatrix;
}


// returns integral matrix reordered according to _waveNames array
template <typename T>
void
TPWALikelihood<T>::reorderIntegralMatrixX(integral&            integral,
					  vector4(complex<T>)& reorderedMatrix) const
{
  // get original matrix and list of wave names
  const matrix<complex<double> > intMatrix    = integral.mat();
  const list<string>             intWaveNames = integral.files();
  // build index lookup-table
  vector2(unsigned int) indexLookUp;  // lookup table: [reflectivity][wave index] -> index in normalization integral
  indexLookUp.resize(2);
  for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {
    indexLookUp[iRefl].resize(_nmbWavesRefl[iRefl], 0);
    for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
      if (find(intWaveNames.begin(), intWaveNames.end(), _waveNamesX[iRefl][iWave]) == intWaveNames.end()) {
	printErr << "wave " << _waveNamesX[iRefl][iWave] << " is not in integral. Aborting." << endl;
	throw;
      }
      indexLookUp[iRefl][iWave] = integral.index(_waveNamesX[iRefl][iWave]);  // not a const method!!! should be fixed
      if (_debug)
	cout << "    mapping wave [" << setw(2) << (int)iRefl * 2 - 1 << ", " << setw(3) << iWave << "] '" << _waveNamesX[iRefl][iWave] << "' "
	     << "to index " << setw(3) << indexLookUp[iRefl][iWave] << " in integral." << endl;
    }
  }
  // create reordered matrix
  reorderedMatrix.clear();
  reorderedMatrix.resize(2);
  for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {
    reorderedMatrix[iRefl].resize(_nmbWavesRefl[iRefl]);
    for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
      reorderedMatrix[iRefl][iWave].resize(2);
      for (unsigned int jRefl = 0; jRefl < 2; ++jRefl) {
	reorderedMatrix[iRefl][iWave][jRefl].resize(_nmbWavesRefl[jRefl], 0);
	for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
	  reorderedMatrix[iRefl][iWave][jRefl][jWave] = intMatrix.element(indexLookUp[iRefl][iWave], indexLookUp[jRefl][jWave]);
// 	  cout << "[" << iRefl << "][" << iWave << "][" << jRefl << "][" << jWave << "] = "
// 	       << "[" << indexLookUp[iRefl][iWave] << "][" << indexLookUp[jRefl][jWave] << "] = "
// 	       << reorderedMatrix[iRefl][iWave][jRefl][jWave] << endl;
	}
      }
    }
  }
}


template <typename T>
void
TPWALikelihood<T>::readIntegrals(const string& normIntFileName,  // name of file with normalization integrals
				 const string& accIntFileName)   // name of file with acceptance integrals
{
  printInfo << "loading normalization integral from '" << normIntFileName << "'." << endl;
  ifstream intFile(normIntFileName.c_str());
  if (!intFile) {
    printErr << "cannot open file '" << normIntFileName << "'. aborting." << endl;
    throw;
  }
  integral normInt;
  // !!! integral.scan() performs no error checks!
  normInt.scan(intFile);
  intFile.close();
  _normMatrix  = reorderedIntegralMatrix(normInt);
  reorderIntegralMatrixX(normInt, _normMatrixX);

  printInfo << "loading acceptance integral from '" << accIntFileName << "'." << endl;
  intFile.open(accIntFileName.c_str());
  if (!intFile) {
    printErr << "cannot open file '" << accIntFileName << "'. exiting." << endl;
    throw;
  }
  integral accInt;
  // !!! integral.scan() performs no error checks!
  accInt.scan(intFile);
  intFile.close();
  //_accInt.events(100000); TODO: add possibility to rescale here!
  _accMatrix = reorderedIntegralMatrix(accInt);
  reorderIntegralMatrixX(accInt, _accMatrixX);
}


template <typename T>
void
TPWALikelihood<T>::readDecayAmplitudes(const string& ampDirName)
{
  // normalization integrals need to be loaded
  if (_normMatrix.nrows() == 0) {
    printErr << "normalization integrals have to be loaded before loading the amplitudes. aborting." << endl;
    throw;
  }
  printInfo << "loading amplitude data." << endl;
  
  // monitor memory usage
  ProcInfo_t infoBefore;
  if (_debug) {
    gSystem->GetProcInfo(&infoBefore);
    printInfo << "resident memory usage before loading amplitudes: " << infoBefore.fMemResident << "." << endl;
  }
  
  // clear cache
  clearCache();

  {
    // loop over amplitudes and read in data
    _decayAmpsX.resize(2);  // array dimension for reflectivity
    for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
      for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
	if (_debug)
	  cout << "    loading amplitude data for wave '" << _waveNamesX[iRefl][iWave] << "'." << endl;
	ifstream ampFile((ampDirName + "/" + _waveNamesX[iRefl][iWave]).c_str());
	if (!ampFile) {
	  printErr << "cannot open file '" << _waveNamesX[iRefl][iWave] << "'. aborting." << endl;
	  throw;
	}
	// get normalization
	const complex<double> normInt = _normMatrixX[iRefl][iWave][iRefl][iWave];
	_decayAmpsX[iRefl].push_back(vector<complex<double> >());
	_decayAmpsX[iRefl].back().reserve(_nmbEvents);  // number of events is known except for first wave
	complex<double> amp;
	while (ampFile.read((char*) &amp, sizeof(complex<double>))) {
	  if (_useNormalizedAmps)         // normalize data, if option is switched on
	    amp /= sqrt(normInt.real());  // rescale decay amplitude
	  _decayAmpsX[iRefl].back().push_back(amp);
	}
	if (_debug)
	  cout << "    read " << _decayAmpsX[iRefl].back().size() << " events from file '" << _waveNamesX[iRefl][iWave] << "'." << endl;
	_nmbEvents = _decayAmpsX[iRefl].back().size(); 
      }
    bool allWavesHaveSameNmbEvent = true;
    for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
      for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
	for (unsigned int jRefl = 0; jRefl <= iRefl; ++jRefl)
	  for (unsigned int jWave = 0; jWave <= ((jRefl == iRefl) ? iWave : _nmbWavesRefl[jRefl] - 1); ++jWave)
	    if (_decayAmpsX[iRefl][iWave].size() != _decayAmpsX[jRefl][jWave].size()) {
	      allWavesHaveSameNmbEvent = false;
	      break;
	    }
    if ((!allWavesHaveSameNmbEvent) || (_decayAmpsX[0][0].size() != _nmbEvents))
      printWarn << "amplitude files do not contain the same number of events for all waves." << endl;
    printInfo << "loaded " << _nmbEvents << " events into memory." << endl;
  }

  {
    // loop over amplitudes and read in data
    for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {
      if (_debug)
	cout << "    loading amplitude data for wave '" << _waveNames[iWave] << "'." << endl;
      ifstream ampFile((ampDirName + "/" + _waveNames[iWave]).c_str());
      if (!ampFile) {
	printErr << "cannot open file '" << _waveNames[iWave] << "'. aborting." << endl;
	throw;
      }
      // get normalization
      const complex<double> normInt = _normMatrix.element(iWave, iWave);
      _decayAmps.push_back(vector<complex<double> >());
      _decayAmps.back().reserve(_nmbEvents);  // number of events is known except for first wave
      complex<double> amp;
      while (ampFile.read((char*) &amp, sizeof(complex<double>))) {
	if (_useNormalizedAmps)         // normalize data, if option is switched on
	  amp /= sqrt(normInt.real());  // rescale decay amplitude
	_decayAmps.back().push_back(amp);
      }
      if (_debug)
	cout << "    read " << _decayAmps.back().size() << " events from file '" << _waveNames[iWave] << "'." << endl;
      _nmbEvents = _decayAmps.back().size(); 
    }
    bool allWavesHaveSameNmbEvent = true;
    for (unsigned int iWave = 1; iWave < _decayAmps.size(); ++iWave)
      if (_decayAmps[iWave - 1].size() != _decayAmps[iWave].size()) {
	allWavesHaveSameNmbEvent = false;
	break;
      }
    if ((!allWavesHaveSameNmbEvent) || (_decayAmps.back().size() != _nmbEvents))
      printWarn << "amplitude files do not contain the same number of events for all waves." << endl;
    printInfo << "loaded " << _nmbEvents << " events into memory." << endl;
  }

  // monitor memory usage
  if (_debug) {
    ProcInfo_t infoAfter;
    gSystem->GetProcInfo(&infoAfter);
    printInfo << "resident memory usage after loading amplitudes: " << infoAfter.fMemResident << ". "
	      << "memory used for amplitudes: " << infoAfter.fMemResident - infoBefore.fMemResident << "." << endl;
  }

  // rescale integrals, if necessary
  if (_useNormalizedAmps) {
    // matrices _normMatrix and _accMatrix are already normalized to number of Monte Carlo events
    printInfo << "rescaling integrals." << endl;
    // rescale normalization integral
    for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {
      const double norm_i = sqrt(_normMatrix.element(iWave, iWave).real());
      for (unsigned int jWave = 0; jWave < _nmbWaves; ++jWave) {
	if (iWave == jWave)
	  continue;  // set diagonal terms later so that norm_i,j stay unaffected
	const double norm_j = sqrt(_normMatrix.element(jWave, jWave).real());
	_normMatrix.element(iWave, jWave) /= norm_i * norm_j;
      }
    }
    // rescale acceptance integral
    for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {
      const double norm_i = sqrt(_normMatrix.element(iWave, iWave).real());
      for (unsigned int jWave = 0; jWave < _nmbWaves; ++jWave) {
	const double norm_j = sqrt(_normMatrix.element(jWave, jWave).real());
	_accMatrix.element(iWave, jWave) /= norm_i * norm_j;
// !!! magic do-nothing line; without it fit result changes
// this is probably due to the -O3 used in compilation which might break some things
// has no influence when using -O0
	cout << "";
      }
    }
    // set diagonal elements of normalization matrix
    for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave)
      _normMatrix.element(iWave, iWave) = 1;

    if (_debug)
      for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave)
	for (unsigned int jWave = 0; jWave < _nmbWaves; ++jWave)
	  cout << "    normalization matrix [" << setw(3) << iWave << ", " << setw(3) << jWave << "] = "
 	       << "("  << maxPrecisionAlign(_normMatrix.element(iWave, jWave).real())
 	       << ", " << maxPrecisionAlign(_normMatrix.element(iWave, jWave).imag()) << "), "
	       << "acceptance matrix [" << setw(3) << iWave << ", " << setw(3) << jWave << "] = "
 	       << "("  << maxPrecisionAlign(_accMatrix.element(iWave, jWave).real())
 	       << ", " << maxPrecisionAlign(_accMatrix.element(iWave, jWave).imag()) << ")"
	       << endl;

    // rescale normalization and acceptance integrals
    for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
      for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
	const double norm_i = sqrt(_normMatrixX[iRefl][iWave][iRefl][iWave].real());
	for (unsigned int jRefl = 0; jRefl < 2; ++jRefl)
	  for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
	    const double norm_j = sqrt(_normMatrixX[jRefl][jWave][jRefl][jWave].real());
	    if ((iRefl != jRefl) || (iWave != jWave))  // set diagonal terms later so that norm_i,j stay unaffected
	      _normMatrixX[iRefl][iWave][jRefl][jWave] /= norm_i * norm_j;
	    _accMatrixX[iRefl][iWave][jRefl][jWave] /= norm_i * norm_j;
	  }
      }
    // set diagonal elements of normalization matrix
    for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
      for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
	_normMatrixX[iRefl][iWave][iRefl][iWave] = 1;  // diagonal term

    if (_debug) {
      cout << endl << endl;
      for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
	for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
	  for (unsigned int jRefl = 0; jRefl < 2; ++jRefl)
	    for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
	      cout << "    normalization matrix [" << setw(2) << (int)iRefl * 2 - 1 << ", " << setw(3) << iWave << ", "
		   << setw(2) << (int)jRefl * 2 - 1 << ", " << setw(3) << jWave << "] = "
		   << "("  << maxPrecisionAlign(_normMatrixX[iRefl][iWave][jRefl][jWave].real())
		   << ", " << maxPrecisionAlign(_normMatrixX[iRefl][iWave][jRefl][jWave].imag()) << "), "
		   << "acceptance matrix [" << setw(3) << iWave << ", " << setw(3) << jWave << "] = "
		   << "("  << maxPrecisionAlign(_accMatrixX[iRefl][iWave][jRefl][jWave].real())
		   << ", " << maxPrecisionAlign(_accMatrixX[iRefl][iWave][jRefl][jWave].imag()) << ")"
		   << endl;
	    }
    }
  }  // end _useNormalizedAmps
}


template <typename T>
void
TPWALikelihood<T>::getIntCMatrix(TCMatrix& normMatrix,
				 TCMatrix& accMatrix)
{
  //normMatrix.ResizeTo(_nmbWaves,_nmbWaves);
  for (unsigned int i = 0; i < _nmbWaves; ++i)  // outer loop
    for (unsigned int j = 0; j < _nmbWaves; ++j) {  // inner loop
      normMatrix.set(i, j, _normMatrix.element(i, j));
      accMatrix.set (i, j, _accMatrix.element (i, j));
    }
  // add flat
  normMatrix.set(_nmbWaves, _nmbWaves, 1);
  accMatrix.set (_nmbWaves, _nmbWaves, 1);
}


// Complex valued Amplitudes and
// mapping of real and imaginary part of amplitudes
// in error matrix (needed by TFitBin)
// takes into account real-valued parameters
template <typename T>
void
TPWALikelihood<T>::buildCAmps(const double*             x,
			      vector<complex<double> >& V,
			      vector<pair<int,int> >&   indices,
			      vector<string>&           names,
			      const bool                withFlat)
{
  // build complex numbers from parameters
  // remember rank restrictions!
  V.clear();
  indices.clear();
  names.clear();
  //unsigned int namp=_rank*_nmbWaves;
  //if (withFlat)namp+=1;
  unsigned int parIndex = 0;
  for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
    unsigned int iWaveRefl[2] = {0, 0};  // indices for negative and positive reflectivity waves
    for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {
      int refl = (getReflectivity(_waveNames[iWave]) > 0 ? 1 : 0);
      double re, im;
      if (iWaveRefl[refl] == iRank) {  // real production amplitude
	indices.push_back(make_pair(parIndex, -1));
	re = x[parIndex++];
	im = 0;
	V.push_back(complex<double>(re, im));
	stringstream name;
	name << "V" << iRank << "_" << _waveNames[iWave];
	names.push_back(name.str());
      } else if (iWaveRefl[refl] > iRank) {  // complex production amplitude
	indices.push_back(make_pair(parIndex, parIndex + 1));
	re = x[parIndex++];
	im = x[parIndex++];
	V.push_back(complex<double>(re, im));
	stringstream name;
	name << "V" << iRank << "_" << _waveNames[iWave];
	names.push_back(name.str());
      }
      ++iWaveRefl[refl];
      //cout << "Wave" << iWave << "=" << V[iRank * _nmbWaves + iWave] << endl;
    }
  }
  if (withFlat){
    V.push_back(complex<double>(x[parIndex], 0));
    indices.push_back(make_pair(parIndex, -1));
    names.push_back("V_flat");
  }
}


template <typename T>
void
TPWALikelihood<T>::clearCache()
{
  for (unsigned int i = 0; i < _decayAmps.size(); ++i)
    _decayAmps[i].clear();
  _decayAmps.clear();
}


// depends on naming convention for waves!!!
// VR_IGJPCMEIso....
template <typename T>
int
TPWALikelihood<T>::getReflectivity(const TString& waveName) const
{
  int refl = 0;
  unsigned int reflIndex = 6;  // position of reflectivity in wave
  // check whether it is parameter or wave name
  if (waveName[0] == 'V')
    reflIndex = 9; 
  if (waveName[reflIndex] == '-')
    refl= -1;
  else if (waveName[reflIndex] == '+')
    refl= +1;
  else {
    printErr << "Cannot parse parameter/wave name '" << waveName << "'. Cannot not determine reflectivity. Aborting." << endl;
    throw;
  }
  if (_debug)
    printInfo << "Extracted reflectivity = " << refl << " from parameter name '" << waveName << "' (char position " << reflIndex << ")" << endl;
  return refl;
}


// copy values from array that corresponds to the function parameters
// to structure that corresponds to the complex production amplitudes
// taking into account rank restrictions
template <typename T>
void
TPWALikelihood<T>::copyFromParArray(const double*        inPar,             // input parameter array
				    vector2(complex<T>)& outVal,            // output values organized as 2D array of complex numbers with [rank][wave index]
				    T&                   outFlatVal) const  // output value corresponding to flat wave
{
  outVal.clear();
  outVal.resize(_rank, vector<complex<T> >(_nmbWaves, 0));
  unsigned int parIndex = 0;
  for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
    unsigned int iWaveRefl[2] = {0, 0};  // indices for negative and positive reflectivity waves
    for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {
      const unsigned int refl = (_waveRefl[iWave] > 0) ? 1 : 0;
      T re, im;
      if (iWaveRefl[refl] < iRank)          // production amplitude is zero
	re = im = 0;
      else if (iWaveRefl[refl] == iRank) {  // production amplitude is real
	re = inPar[parIndex++];
	im = 0;
      } else {                              // production amplitude is complex
	re = inPar[parIndex++];
	im = inPar[parIndex++];
      }
      outVal[iRank][iWave] = complex<T>(re, im);
      ++iWaveRefl[refl];
    }
  }
  outFlatVal = inPar[parIndex];
}


template <typename T>
void
TPWALikelihood<T>::copyFromParArray(const double* inPar,             // input parameter array
				    complex<T>**& outVal,            // output values organized as 2D array of complex numbers with [rank][wave index]
				    T&            outFlatVal) const  // output value corresponding to flat wave
{
  const unsigned int outValDim[2] = {_rank, _nmbWaves};
  allocate2DArray<complex<T> >(outVal, outValDim, 0);
  unsigned int parIndex = 0;
  for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
    unsigned int iWaveRefl[2] = {0, 0};  // indices for negative and positive reflectivity waves
    for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {
      const unsigned int refl = (_waveRefl[iWave] > 0) ? 1 : 0;
      T re, im;
      if (iWaveRefl[refl] < iRank)          // production amplitude is zero
	re = im = 0;
      else if (iWaveRefl[refl] == iRank) {  // production amplitude is real
	re = inPar[parIndex++];
	im = 0;
      } else {                              // production amplitude is complex
	re = inPar[parIndex++];
	im = inPar[parIndex++];
      }
      outVal[iRank][iWave] = complex<T>(re, im);
      ++iWaveRefl[refl];
    }
  }
  outFlatVal = inPar[parIndex];
}


template <typename T>
void
TPWALikelihood<T>::copyFromParArrayX(const double*        inPar,             // input parameter array
				     vector3(complex<T>)& outVal,            // output values organized as 3D array of complex numbers with [rank][reflectivity][wave index]
				     T&                   outFlatVal) const  // output value corresponding to flat wave
{
  // group parameters by rank and wave index only; keeps the order defined in wave list
  outVal.clear();
  outVal.resize(_rank, vector2(complex<T>)(2));
  unsigned int parIndex = 0;
  for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
    outVal[iRank][0].resize(_nmbWavesRefl[0], 0);
    outVal[iRank][1].resize(_nmbWavesRefl[1], 0);
    unsigned int iWaveRefl[2] = {0, 0};  // indices for negative and positive reflectivity waves
    for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {  // loop over all waves regardless of reflectivity
	const unsigned int refl = (_waveRefl[iWave] > 0) ? 1 : 0;
	T re, im;
	if (iWaveRefl[refl] < iRank)          // production amplitude is zero
	  re = im = 0;
	else if (iWaveRefl[refl] == iRank) {  // production amplitude is real
	  re = inPar[parIndex++];
	  im = 0;
	} else {                              // production amplitude is complex
	  re = inPar[parIndex++];
	  im = inPar[parIndex++];
	}
	outVal[iRank][refl][iWaveRefl[refl]] = complex<T>(re, im);
	//cout << "!!! [" << iRank << ", " << refl << ", " << iWaveRefl[refl] << "] = " << outVal[iRank][refl][iWaveRefl[refl]] << endl;
	++iWaveRefl[refl];
      }
    }
  outFlatVal = inPar[parIndex];

//   { // test
//     vector2(complex<T>) testVal;
//     copyFromParArray(inPar, testVal, outFlatVal);
//     for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
//       unsigned int waveIndex[2] = {0, 0};
//       for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {
// 	const unsigned int refl = (_waveRefl[iWave] > 0) ? 1 : 0;
// 	cout << "[" << iRank << ", " << setw(2) << iWave << "] = " << setw(3) << testVal[iRank][iWave] << " vs. "
// 	      << setw(3) << outVal[iRank][refl][waveIndex[refl]] << " at [" << iRank << ", " << refl << ", " << waveIndex[refl] << "]" << endl;
// 	++waveIndex[refl];
//       }
//     }
//   }

  // group parameters by rank, reflectivity and wave index
//   outVal.clear();
//   outVal.resize(_rank, vector2(complex<T>)(2));
//   unsigned int parIndex = 0;
//   for (unsigned int iRank = 0; iRank < _rank; ++iRank)
//     for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {
//       outVal[iRank][iRefl].resize(_nmbWavesRefl[iRefl], 0);
//       for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
// 	T re, im;
// 	if (iWave < iRank)          // production amplitude is zero
// 	  re = im = 0;
// 	else if (iWave == iRank) {  // production amplitude is real
// 	  re = inPar[parIndex++];
// 	  im = 0;
// 	} else {                    // production amplitude is complex
// 	  re = inPar[parIndex++];
// 	  im = inPar[parIndex++];
// 	}
// 	outVal[iRank][iRefl][iWave] = complex<T>(re, im);
//       }
//     }
//   outFlatVal = inPar[parIndex];
}


// copy values from structure that corresponds to complex
// production amplitudes to array that corresponds to function
// parameters taking into account rank restrictions
template <typename T>
void
TPWALikelihood<T>::copyToParArray(const vector2(complex<T>)& inVal,         // values corresponding to production amplitudes
				  const T                    inFlatVal,     // value corresponding to flat wave
				  double*                    outPar) const  // output parameter array
{
  unsigned int parIndex = 0;
  for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
    unsigned int iWaveRefl[2] = {0, 0};  // indices for negative and positive reflectivity waves
    for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {
      const unsigned int refl = (_waveRefl[iWave] > 0) ? 1 : 0;
      if (iWaveRefl[refl] == iRank)  // production amplitude is real
	outPar[parIndex++] = inVal[iRank][iWave].real();
      else if (iWaveRefl[refl] > iRank) {  // production amplitude is complex
	outPar[parIndex++] = inVal[iRank][iWave].real();
	outPar[parIndex++] = inVal[iRank][iWave].imag();
      }
      ++iWaveRefl[refl];
    }
  }
  outPar[parIndex] = inFlatVal;
}


template <typename T>
void
TPWALikelihood<T>::copyToParArrayX(const vector3(complex<T>)& inVal,         // values corresponding to production amplitudes
				   const T                    inFlatVal,     // value corresponding to flat wave
				   double*                    outPar) const  // output parameter array
{
  for (unsigned int iRank = 0; iRank < _rank; ++iRank)
    for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
      for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
	pair<int, int> parIndices = _prodAmpToFuncParMapX[iRank][iRefl][iWave];
// 	cout << " [" << iRank << ", " << iRefl << ", " << iWave << "] -> " << parIndices
// 	     << ": " << inVal[iRank][iRefl][iWave] << endl;
	if (parIndices.first >= 0) {   // real part
	  outPar[parIndices.first]  = inVal[iRank][iRefl][iWave].real();
// 	  cout << "!!! HERE: " << parIndices.first << ", " << inVal[iRank][iRefl][iWave].real() 
// 	       << " vs. " << outPar[parIndices.first] << endl;
	}
	if (parIndices.second >= 0)  // imaginary part
	  outPar[parIndices.second] = inVal[iRank][iRefl][iWave].imag();
      }
  outPar[_nmbPars - 1] = inFlatVal;
//   cout << "!!! outPar[0] = " << outPar[0] << " vs. " << inVal[0][1][0] << endl;
//   outPar[0] = inVal[0][1][0].real();
//   cout << "!!! outPar[0] = " << outPar[0] << " vs. " << inVal[0][1][0] << endl;
}


template <typename T>
ostream&
TPWALikelihood<T>::print(ostream& out) const
{
  out << "TPWALikelihood parameters:" << endl
      << "number of events ........................ " << _nmbEvents         << endl
      << "rank .................................... " << _rank              << endl
      << "number of waves ......................... " << _nmbWaves          << endl
      << "number of positive reflectivity waves ... " << _nmbWavesRefl[1]   << endl
      << "number of negative reflectivity waves ... " << _nmbWavesRefl[0]   << endl
      << "number of function parameters ........... " << _nmbPars           << endl
      << "print debug messages .................... " << _debug             << endl
      << "use normalized amplitudes ............... " << _useNormalizedAmps << endl
      << "list of waves: " << endl;
  for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave)
    out << "        " << setw(3) << iWave << " " << _waveNames[iWave] << "    "
	<< "reflectivity = " << _waveRefl[iWave] << "    "
	<< "threshold = " << _waveThresholds[iWave] << " MeV/c^2" << std::endl;
  out << "list of function parameters: " << endl;
  for (unsigned int iPar = 0; iPar < _nmbPars; ++iPar)
    out << "        " << setw(3) << iPar << " " << _parNames[iPar] << "    "
	<< "threshold = " << _parThresholds[iPar] << " MeV/c^2" << std::endl;
  return out;
}
