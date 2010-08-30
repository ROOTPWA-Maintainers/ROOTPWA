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

#include "boost/tuple/tuple.hpp"

#include "TString.h"
#include "TSystem.h"
#include "TCMatrix.h"
#include "TStopwatch.h"

#include "utilities.h"

#ifdef USE_CUDA
#include "../cuda/complex.cuh"
#include "../cuda/likelihoodInterface.cuh"
using namespace rpwa;
#endif


//#define USE_FDF


using namespace std;
using namespace boost;


template<typename T> bool TPWALikelihood<T>::_debug = true;


template<typename T>
TPWALikelihood<T>::TPWALikelihood()
	: _nmbEvents        (0),
	  _rank             (1),
	  _nmbWaves         (0),
	  _nmbPars          (0),
#ifdef USE_CUDA
	  _cudaEnabled      (false),
#endif
	  _useNormalizedAmps(true),
	  _numbAccEvents    (0),
	  _genCudaDiffHist  (false)
{
	_nmbWavesRefl[0] = 0;
	_nmbWavesRefl[1] = 0;
	resetFuncCallInfo();
#ifdef USE_FDF
	printInfo << "using FdF() to calculate likelihood" << endl;
#else
	printInfo << "using DoEval() to calculate likelihood" << endl;
#endif
}


template<typename T>
TPWALikelihood<T>::~TPWALikelihood()
{
	if (_cudaEnabled and _genCudaDiffHist) {
		_outFile->Write();
		_outFile->Close();
	}
	clear();
}


template<typename T>
void
TPWALikelihood<T>::FdF(const double* par,             // parameter array; reduced by rank conditions
                       double&       funcVal,         // function value
                       double*       gradient) const  // array of derivatives
{
	++(_funcCallInfo[FDF].nmbCalls);
  
	// timer for total time
	TStopwatch timerTot;
	timerTot.Start();
  
	// copy arguments into parameter cache
	for (unsigned int i = 0; i < _nmbPars; ++i)
		_parCache[i] = par[i];

	// build complex production amplitudes from function parameters taking into account rank restrictions
	T             prodAmpFlat;
	ampsArrayType prodAmps;
	copyFromParArray(par, prodAmps, prodAmpFlat);
	const T prodAmpFlat2 = prodAmpFlat * prodAmpFlat;

	// create array of likelihood derivative w.r.t. real and imaginary
	// parts of the production amplitudes
	// !NOTE! although stored as and constructed from complex values,
	// the dL themselves are _not_ well defined complex numbers!
	T                              derivativeFlat = 0;
	array<ampsArrayType::index, 3> derivShape     = {{ _rank, 2,
	                                                   max(_nmbWavesRefl[0], _nmbWavesRefl[1]) }};
	ampsArrayType                  derivatives(derivShape);

	// loop over events and calculate first term of log likelihood sum
	// as well as derivatives with respect to parameters
	TStopwatch timer;
	timer.Start();
	T logLikelihood = 0;
	for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
		T             likelihood = 0;          // likelihood for this event
		ampsArrayType d(derivShape);  // likelihood derivative for this event
		for (unsigned int iRank = 0; iRank < _rank; ++iRank) {                     // incoherent sum over ranks
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {                       // incoherent sum over reflectivities
				complex<T> ampProdSum = 0;                                             // amplitude sum for negative/positive reflectivity for this rank
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
					// compute likelihood term
					const complex<T> amp = prodAmps[iRank][iRefl][iWave]
						* complex<T>(_decayAmps[iEvt][iRefl][iWave]);
					ampProdSum += amp;
					// compute derivatives
					for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave)
						// sum up amplitudes for current rank and only for waves with same reflectivity
						d[iRank][iRefl][jWave] += amp;
				}
				likelihood += norm(ampProdSum);
			}
			assert(likelihood >= 0);
			// loop again over waves for current rank
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					d[iRank][iRefl][iWave] *= conj(complex<T>(_decayAmps[iEvt][iRefl][iWave]));
		}  // end loop over rank
		likelihood    += prodAmpFlat2;
		logLikelihood -= log(likelihood);  // accumulate log likelihood
		// incorporate factor 2 / sigma
		const T factor = 2. / likelihood;
		for (unsigned int iRank = 0; iRank < _rank; ++iRank)
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					derivatives[iRank][iRefl][iWave] -= factor * d[iRank][iRefl][iWave];
		derivativeFlat -= factor * prodAmpFlat;
	}  // end loop over events
	// log time needed for likelihood calculation
	timer.Stop();
	_funcCallInfo[FDF].funcTime += timer.RealTime();

	// compute normalization term of log likelihood and normalize derivatives w.r.t. parameters
	timer.Start();
	complex<T> normFactor  = 0;
	const T    nmbEvt      = (_useNormalizedAmps) ? 1 : (T)_nmbEvents;
	const T    twiceNmbEvt = 2 * nmbEvt;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				complex<T> normFactorDeriv = 0;
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave) {  // inner loop over waves with same reflectivity
					const complex<T> I = complex<T>(_accMatrix[iRefl][iWave][iRefl][jWave]);
					normFactor      += (prodAmps[iRank][iRefl][iWave]
					                    * conj(prodAmps[iRank][iRefl][jWave])) * I;
					normFactorDeriv += prodAmps[iRank][iRefl][jWave] * conj(I);
				}
				derivatives[iRank][iRefl][iWave] += normFactorDeriv * twiceNmbEvt;  // account for 2 * nmbEvents
			}
	// take care of flat wave
	normFactor.real() += prodAmpFlat2;
	derivativeFlat    += twiceNmbEvt * prodAmpFlat;
	// log time needed for normalization
	timer.Stop();
	_funcCallInfo[FDF].normTime += timer.RealTime();

	// sort derivative results into output array and cache
	copyToParArray(derivatives, derivativeFlat, gradient);
	copyToParArray(derivatives, derivativeFlat, toArray(_derivCache));

	// set function return value
	funcVal = logLikelihood + nmbEvt * normFactor.real();

	// log total consumed time
	timerTot.Stop();
	_funcCallInfo[FDF].totalTime += timerTot.RealTime();

	if (_debug)
		printInfo << "log likelihood =  " << maxPrecisionAlign(logLikelihood) << ", "
		          << "normalization =  " << maxPrecisionAlign(normFactor.real()) << ", "
		          << "normalized likelihood = "
		          << maxPrecisionAlign(logLikelihood + nmbEvt * normFactor.real()) << endl;
}


template<typename T>
double
TPWALikelihood<T>::DoEval(const double* par) const
{
	++(_funcCallInfo[DOEVAL].nmbCalls);

#ifdef USE_FDF

	// call FdF
	double logLikelihood;
	double gradient[_nmbPars];
	FdF(par, logLikelihood, gradient);
	return logLikelihood;

#else

	// timer for total time
	TStopwatch timerTot;
	timerTot.Start();
  
	// build complex production amplitudes from function parameters taking into account rank restrictions
	T             prodAmpFlat;
	ampsArrayType prodAmps;
	copyFromParArray(par, prodAmps, prodAmpFlat);
	const T prodAmpFlat2 = prodAmpFlat * prodAmpFlat;

	// loop over events and calculate real-data term of log likelihood
	TStopwatch timer;
	timer.Start();
	T logLikelihood = 0;
#ifdef USE_CUDA
	if (_cudaEnabled) {
		logLikelihood = cuda::likelihoodInterface<cuda::complex<T> >::logLikelihood
			(reinterpret_cast<cuda::complex<T>*>(prodAmps.data()),
			 prodAmps.num_elements(), prodAmpFlat, _rank);
	}
#endif
	if (not _cudaEnabled or _genCudaDiffHist)	{
		T logLikelihoodCpu = 0;
		for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
			T likelihood = 0;  // likelihood for this event
			for (unsigned int iRank = 0; iRank < _rank; ++iRank) {                   // incoherent sum over ranks
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {                     // incoherent sum over reflectivities
					complex<T> ampProdSum = 0;                                           // amplitude sum for negative/positive reflectivity for this rank
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)  // coherent sum over waves
						// compute likelihood term
						ampProdSum += prodAmps[iRank][iRefl][iWave] * complex<T>(_decayAmps[iEvt][iRefl][iWave]);
					likelihood += norm(ampProdSum);
				}
				assert(likelihood >= 0);
			}  // end loop over rank
			likelihood       += prodAmpFlat2;
			logLikelihoodCpu -= log(likelihood);  // accumulate log likelihood
		}  // end loop over events
		if (_cudaEnabled and _genCudaDiffHist) {  // fill difference histograms
			const T diffAbs = logLikelihoodCpu - logLikelihood;
			const T diffRel = 1 - logLikelihood / logLikelihoodCpu;
			_hLikelihoodDiffAbs->Fill(diffAbs);
			_hLikelihoodDiffRel->Fill(diffRel);
		}
		if (not _cudaEnabled)
			logLikelihood = logLikelihoodCpu;
	}
	// log time needed for likelihood calculation
	timer.Stop();
	_funcCallInfo[DOEVAL].funcTime += timer.RealTime();
	
	// compute normalization term of log likelihood
	timer.Start();
	complex<T> normFactor = 0;
	const T    nmbEvt     = (_useNormalizedAmps) ? 1 : (T)_nmbEvents;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave) {  // inner loop over waves with same reflectivity
					const complex<T> I = complex<T>(_accMatrix[iRefl][iWave][iRefl][jWave]);
					normFactor += (prodAmps[iRank][iRefl][iWave] * conj(prodAmps[iRank][iRefl][jWave])) * I;
				}
	// take care of flat wave
	normFactor.real() += prodAmpFlat2;
	// log time needed for normalization
	timer.Stop();
	_funcCallInfo[DOEVAL].normTime += timer.RealTime();
  
	// calculate and return log likelihood value
	const double funcVal = logLikelihood + nmbEvt * normFactor.real();

	// log total consumed time
	timerTot.Stop();
	_funcCallInfo[DOEVAL].totalTime += timerTot.RealTime();

	if (_debug)
		printInfo << "log likelihood =  " << maxPrecisionAlign(logLikelihood) << ", "
		          << "normalization =  " << maxPrecisionAlign(normFactor.real()) << ", "
		          << "normalized likelihood = "
		          << maxPrecisionAlign(logLikelihood + nmbEvt * normFactor.real()) << endl;
	
	return funcVal;

#endif  // USE_FDF
}


template<typename T>
double
TPWALikelihood<T>::DoDerivative(const double* par,
                                unsigned int  derivativeIndex) const
{
	++(_funcCallInfo[DODERIVATIVE].nmbCalls);

	// timer for total time
	TStopwatch timerTot;
	timerTot.Start();

	// check whether parameter is in cache
	bool samePar = true;
	for (unsigned int i = 0; i < _nmbPars; ++i)
		if (_parCache[i] != par[i]) {
			samePar = false;
			break;
		}
	if (samePar) {
		//cout << "using cached derivative! " << endl;
		timerTot.Stop();
		_funcCallInfo[DODERIVATIVE].totalTime += timerTot.RealTime();
		_funcCallInfo[DODERIVATIVE].funcTime   = _funcCallInfo[DODERIVATIVE].totalTime;
		return _derivCache[derivativeIndex];
	}
	// call FdF
	double logLikelihood;
	double gradient[_nmbPars];
	FdF(par, logLikelihood, gradient);
	return gradient[derivativeIndex];
}


// calculate derivatives with respect to parameters
template<typename T>
void
TPWALikelihood<T>::Gradient(const double* par,             // parameter array; reduced by rank conditions
                            double*       gradient) const  // array of derivatives
{
	++(_funcCallInfo[GRADIENT].nmbCalls);
  
	// timer for total time
	TStopwatch timerTot;
	timerTot.Start();

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
		timerTot.Stop();
		_funcCallInfo[GRADIENT].totalTime += timerTot.RealTime();
		_funcCallInfo[GRADIENT].funcTime   = _funcCallInfo[GRADIENT].totalTime;
		return;
	}
	// call FdF
	double logLikelihood;
	FdF(par, logLikelihood, gradient);

#else

	// copy arguments into parameter cache
	for (unsigned int i = 0; i < _nmbPars; ++i)
		_parCache[i] = par[i];

	// build complex production amplitudes from function parameters taking into account rank restrictions
	T             prodAmpFlat;
	ampsArrayType prodAmps;
	copyFromParArray(par, prodAmps, prodAmpFlat);

	// create array of likelihood derivative w.r.t. real and imaginary
	// parts of the production amplitudes
	// !NOTE! although stored as and constructed from complex values,
	// the dL themselves are _not_ well defined complex numbers!
	array<ampsArrayType::index, 3> derivShape = {{ _rank, 2, max(_nmbWavesRefl[0], _nmbWavesRefl[1]) }};
	ampsArrayType                  derivatives(derivShape);
	T                              derivativeFlat = 0;

	// compute derivative for first term of log likelihood
	TStopwatch timer;
	timer.Start();
#ifdef USE_CUDA
	if (_cudaEnabled) {
		cuda::likelihoodInterface<cuda::complex<T> >::logLikelihoodDeriv
			(reinterpret_cast<cuda::complex<T>*>(prodAmps.data()),
			 prodAmps.num_elements(), prodAmpFlat, _rank,
			 reinterpret_cast<cuda::complex<T>*>(derivatives.data()),
			 derivativeFlat);
	}
#endif
	if (not _cudaEnabled or _genCudaDiffHist)	{
		ampsArrayType derivativesCpu(derivShape);
		T             derivativeFlatCpu = 0;
		const T prodAmpFlat2 = prodAmpFlat * prodAmpFlat;
		for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
			T             likelihood = 0;          // likelihood for this event
			ampsArrayType derivative(derivShape);  // likelihood derivative for this event
			for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
					complex<T> ampProdSum = 0;  // amplitude sum for negative/positive reflectivity for this rank
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
						// compute likelihood term
						const complex<T> amp =   prodAmps[iRank][iRefl][iWave]
							                     * complex<T>(_decayAmps[iEvt][iRefl][iWave]);
						ampProdSum += amp;
					}
					likelihood += norm(ampProdSum);
					// set derivative term that is independent on derivative wave index
					for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave)
						// amplitude sums for current rank and for waves with same reflectivity
						derivative[iRank][iRefl][jWave] = ampProdSum;
				}
				assert(likelihood >= 0);
				// loop again over waves for current rank and multiply with complex conjugate
				// of decay amplitude of the wave with the derivative wave index
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
					for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave)
						derivative[iRank][iRefl][jWave] *= conj(complex<T>(_decayAmps[iEvt][iRefl][jWave]));
			}  // end loop over rank
			likelihood += prodAmpFlat2;
			// incorporate factor 2 / sigma
			const T factor = 2. / likelihood;
			for (unsigned int iRank = 0; iRank < _rank; ++iRank)
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
					for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave)
						derivativesCpu[iRank][iRefl][jWave] -= factor * derivative[iRank][iRefl][jWave];
			derivativeFlatCpu -= factor * prodAmpFlat;
		}  // end loop over events
		if (_cudaEnabled and _genCudaDiffHist) {  // fill difference histograms
			for (unsigned int iRank = 0; iRank < _rank; ++iRank)
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
						complex<T> diffAbs =   derivativesCpu[iRank][iRefl][iWave]
							                   - derivatives   [iRank][iRefl][iWave];
						complex<T> diffRel;
						diffRel.real() = 1 -   derivatives   [iRank][iRefl][iWave].real()
							                   / derivativesCpu[iRank][iRefl][iWave].real();
						diffRel.imag() = 1 -   derivatives   [iRank][iRefl][iWave].imag()
							                   / derivativesCpu[iRank][iRefl][iWave].imag();
						_hDerivDiffAbs[iRank][iRefl][iWave][0]->Fill(diffAbs.real());
						_hDerivDiffAbs[iRank][iRefl][iWave][1]->Fill(diffAbs.imag());
						_hDerivDiffRel[iRank][iRefl][iWave][0]->Fill(diffRel.real());
						_hDerivDiffRel[iRank][iRefl][iWave][1]->Fill(diffRel.imag());
						_hDerivDiffTotAbs[0]->Fill(diffAbs.real());
						_hDerivDiffTotAbs[1]->Fill(diffAbs.imag());
						_hDerivDiffTotRel[0]->Fill(diffRel.real());
						_hDerivDiffTotRel[1]->Fill(diffRel.imag());
					}
			const T diffFlatAbs = derivativeFlatCpu - derivativeFlat;
			const T diffFlatRel = 1 - derivativeFlat / derivativeFlatCpu;
			_hDerivDiffFlatAbs->Fill(diffFlatAbs);
			_hDerivDiffFlatRel->Fill(diffFlatRel);
		}
		if (not _cudaEnabled) {
			derivatives    = derivativesCpu;
			derivativeFlat = derivativeFlatCpu;
		}
	}
	// log time needed for gradient calculation
	timer.Stop();
	_funcCallInfo[GRADIENT].funcTime += timer.RealTime();

	// normalize derivatives w.r.t. parameters
	timer.Start();
	const T nmbEvt      = (_useNormalizedAmps) ? 1 : (T)_nmbEvents;
	const T twiceNmbEvt = 2 * nmbEvt;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				complex<T> normFactorDeriv = 0;
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave) {  // inner loop over waves with same reflectivity
					const complex<T> I = complex<T>(_accMatrix[iRefl][iWave][iRefl][jWave]);
					normFactorDeriv += prodAmps[iRank][iRefl][jWave] * conj(I);
				}
				derivatives[iRank][iRefl][iWave] += normFactorDeriv * twiceNmbEvt;  // account for 2 * nmbEvents
			}
	// take care of flat wave
	derivativeFlat += twiceNmbEvt * prodAmpFlat;
	// log time needed for normalization
	timer.Stop();
	_funcCallInfo[GRADIENT].normTime += timer.RealTime();

	// set return gradient values
	copyToParArray(derivatives, derivativeFlat, gradient);
	copyToParArray(derivatives, derivativeFlat, toArray(_derivCache));

	// log total consumed time
	timerTot.Stop();
	_funcCallInfo[GRADIENT].totalTime += timerTot.RealTime();

#endif  // USE_FDF
}


template<typename T>
void
#ifdef USE_CUDA
TPWALikelihood<T>::enableCuda(const bool enableCuda)
{
	_cudaEnabled = enableCuda;
}
#else
TPWALikelihood<T>::enableCuda(const bool) { }
#endif


template<typename T>
bool
TPWALikelihood<T>::cudaEnabled() const
{
#ifdef USE_CUDA
	return _cudaEnabled;
#else
	return false;
#endif
}


template<typename T>
void 
TPWALikelihood<T>::init(const unsigned int rank,
                        const std::string& waveListFileName,
                        const std::string& normIntFileName,
                        const std::string& accIntFileName,
                        const std::string& ampDirName,
                        const unsigned int numbAccEvents)
{
	_numbAccEvents=numbAccEvents;
	readWaveList(waveListFileName);
	buildParDataStruct(rank);
	readIntegrals(normIntFileName, accIntFileName);
	readDecayAmplitudes(ampDirName);
#ifdef USE_CUDA
	if (_cudaEnabled) {
		cuda::likelihoodInterface<cuda::complex<T> >::init
			(reinterpret_cast<cuda::complex<T>*>(_decayAmps.data()),
			 _decayAmps.num_elements(), _nmbEvents, _nmbWavesRefl, true);
		if (_cudaEnabled and _genCudaDiffHist) {
			if (not _outFile)
				_outFile = TFile::Open("cudaDiff.root", "RECREATE");
			_outFile->cd();
			_hLikelihoodDiffAbs = new TH1D
				("hLikelihoodDiffAbs", "hLikelihoodDiffAbs;L_{CPU} - L_{GPU};Count", 10000, -1e-7, 1e-7);
			_hLikelihoodDiffRel = new TH1D
				("hLikelihoodDiffRel", "hLikelihoodDiffRel;1 - L_{GPU} / L_{CPU};Count",
				 10000, -1e-11, 1e-11);
			_hDerivDiffAbs.resize(extents[_rank][2][max(_nmbWavesRefl[0], _nmbWavesRefl[1])][2]);
			_hDerivDiffRel.resize(extents[_rank][2][max(_nmbWavesRefl[0], _nmbWavesRefl[1])][2]);
			for (unsigned int iRank = 0; iRank < _rank; ++iRank)
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
						for (unsigned int reIm = 0; reIm < 2; ++reIm) {
							stringstream n;
							n << "_" << iRank << "_" << iRefl << "_" << iWave << "_"
							  << ((reIm == 0) ? "real" : "imag");
							_hDerivDiffAbs[iRank][iRefl][iWave][reIm] =
								new TH1D(("hDerivDiffAbs" + n.str()).c_str(),
								         ("hDerivDiffAbs" + n.str() + ";#nabla_{CPU} - #nabla_{GPU};Count").c_str(),
								         10000, -1e-11, 1e-11);
							_hDerivDiffRel[iRank][iRefl][iWave][reIm] =
								new TH1D(("hDerivDiffRel" + n.str()).c_str(),
								         ("hDerivDiffRel" + n.str() + ";1 - #nabla_{GPU} / #nabla_{CPU};Count").c_str(),
								         10000, -1e-11, 1e-11);
						}
			_hDerivDiffFlatAbs = new TH1D
				("hDerivDiffFlatAbs", "hDerivDiffFlatAbs;#nabla_{CPU} - #nabla_{GPU};Count",
				 10000, -1e-10, 1e-10);
			_hDerivDiffFlatRel = new TH1D
				("hDerivDiffFlatRel", "hDerivDiffFlatRel;1 - #nabla_{GPU} / #nabla_{CPU};Count",
				 10000, -1e-12, 1e-12);
			for (unsigned int reIm = 0; reIm < 2; ++reIm) {
				string n = "hDerivDiffTotAbs" + string((reIm == 0) ? "Real" : "Imag");
				_hDerivDiffTotAbs[reIm] = new TH1D
					(n.c_str(), (n + ";#nabla_{CPU} - #nabla_{GPU};Count").c_str(),
					 10000, -1e-10, 1e-10);
				n = "hDerivDiffTotRel" + string((reIm == 0) ? "Real" : "Imag");
				_hDerivDiffTotRel[reIm] = new TH1D
					(n.c_str(), (n + ";1 - #nabla_{GPU} / #nabla_{CPU};Count").c_str(),
					 10000, -1e-10, 1e-10);
			}
		}
	}
#endif
}


template<typename T>
void 
TPWALikelihood<T>::readWaveList(const string& waveListFileName)
{
	printInfo << "reading amplitude names and thresholds from wave list file "
	          << "'" << waveListFileName << "'." << endl;
	ifstream waveListFile(waveListFileName.c_str());
	if (not waveListFile) {
		printErr << "cannot open file '" << waveListFileName << "'. aborting." << endl;
		throw;
	}
	vector<string>       waveNames     [2];
	vector<double>       waveThresholds[2];
	vector<unsigned int> waveIndices   [2];
	unsigned int         countWave = 0;
	unsigned int         lineNmb   = 0;
	string               line;
	while (getline(waveListFile, line)) {
		if (line[0] == '#')  // comments start with #
			continue;
		stringstream lineStream;
		lineStream.str(line);
		string waveName;
		if (lineStream >> waveName) {
			double threshold;
			// !!! it would be safer to make the threshold value in the wave list file mandatory
			if (not (lineStream >> threshold))
				threshold = 0;
			if (_debug)
				printInfo << "reading line " << setw(3) << lineNmb + 1 << ": " << waveName<< ", "
				          << "threshold = " << setw(4) << threshold << " MeV/c^2" << endl;
			if (getReflectivity(waveName) > 0) {
				++_nmbWavesRefl[1];  // positive reflectivity
				waveNames     [1].push_back(waveName);
				waveThresholds[1].push_back(threshold);
				waveIndices   [1].push_back(countWave);
			} else {
				++_nmbWavesRefl[0];  // negative reflectivity
				waveNames     [0].push_back(waveName);
				waveThresholds[0].push_back(threshold);
				waveIndices   [0].push_back(countWave);
			}
			++countWave;
		} else
			printWarn << "cannot parse line '" << line << "' in wave list file "
			          << "'" << waveListFileName << "'" << endl;
		++lineNmb;
	}
	waveListFile.close();
	printInfo << "read " << lineNmb << " lines from wave list file "
	          << "'" << waveListFileName << "'" << endl;
	_nmbWaves = _nmbWavesRefl[0] + _nmbWavesRefl[1];
	_waveNames.resize     (extents[2][max(_nmbWavesRefl[0], _nmbWavesRefl[1])]);
	_waveThresholds.resize(extents[2][max(_nmbWavesRefl[0], _nmbWavesRefl[1])]);
	_waveToWaveList.resize(extents[2][max(_nmbWavesRefl[0], _nmbWavesRefl[1])]);
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
			_waveNames     [iRefl][iWave] = waveNames     [iRefl][iWave];
			_waveThresholds[iRefl][iWave] = waveThresholds[iRefl][iWave];
			_waveToWaveList[iRefl][iWave] = waveIndices   [iRefl][iWave];
		}
}


template<typename T>
void
TPWALikelihood<T>::buildParDataStruct(const unsigned int rank)
{
	if ((_nmbWavesRefl[0] + _nmbWavesRefl[1] == 0) or (_waveThresholds.size() == 0)) {
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
	_parNames.resize     (_nmbPars, "");
	_parThresholds.resize(_nmbPars, 0);
	_parCache.resize     (_nmbPars, 0);
	_derivCache.resize   (_nmbPars, 0);
	_prodAmpToFuncParMap.resize(extents[_rank][2][max(_nmbWavesRefl[0], _nmbWavesRefl[1])]);
	// build parameter names
	unsigned int parIndex = 0;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				ostringstream parName;
				if (iWave < iRank)  // production amplitude is zero
					_prodAmpToFuncParMap[iRank][iRefl][iWave] = make_tuple(-1, -1);
				else if (iWave == iRank) {  // production amplitude is real
					parName << "V" << iRank << "_" << _waveNames[iRefl][iWave] << "_RE";
					_parNames     [parIndex]                  = parName.str();
					_parThresholds[parIndex]                  = _waveThresholds[iRefl][iWave];
					_prodAmpToFuncParMap[iRank][iRefl][iWave] = make_tuple(parIndex, -1);
					++parIndex;
				} else {  // production amplitude is complex
					parName << "V" << iRank << "_" << _waveNames[iRefl][iWave];
					_parNames     [parIndex]                  = parName.str() + "_RE";
					_parThresholds[parIndex]                  = _waveThresholds[iRefl][iWave];
					_parNames     [parIndex + 1]              = parName.str() + "_IM";
					_parThresholds[parIndex + 1]              = _waveThresholds[iRefl][iWave];
					_prodAmpToFuncParMap[iRank][iRefl][iWave] = make_tuple(parIndex, parIndex + 1);
					parIndex += 2;
				}
			}
	// flat wave
	_parNames     [parIndex] = "V_flat";
	_parThresholds[parIndex] = 0;
}


// returns integral matrix reordered according to _waveNames array
template<typename T>
void
TPWALikelihood<T>::reorderIntegralMatrix(integral&            integral,
                                         normMatrixArrayType& reorderedMatrix) const
{
	// get original matrix and list of wave names
	const matrix<complex<double> > intMatrix    = integral.mat();
	const list<string>             intWaveNames = integral.files();
	// build index lookup-table [reflectivity][wave index] to index in normalization integral
	waveToIntMapType indexLookUp(extents[2][max(_nmbWavesRefl[0], _nmbWavesRefl[1])]);
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
			if (find(intWaveNames.begin(), intWaveNames.end(), _waveNames[iRefl][iWave])
			    == intWaveNames.end()) {
				printErr << "wave " << _waveNames[iRefl][iWave] << " is not in integral. aborting." << endl;
				throw;
			}
			indexLookUp[iRefl][iWave] = integral.index(_waveNames[iRefl][iWave]);
			if (_debug)
				printInfo << "    mapping wave [" << setw(2) << (int)iRefl * 2 - 1 << ", "
				          << setw(3) << iWave << "] '" << _waveNames[iRefl][iWave] << "' "
				          << "to index " << setw(3) << indexLookUp[iRefl][iWave] << " in integral." << endl;
		}
	// create reordered matrix
	reorderedMatrix.resize(extents[2][max(_nmbWavesRefl[0], _nmbWavesRefl[1])]
	                       [2][max(_nmbWavesRefl[0], _nmbWavesRefl[1])]);
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
			for (unsigned int jRefl = 0; jRefl < 2; ++jRefl)
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave)
					reorderedMatrix[iRefl][iWave][jRefl][jWave] = intMatrix.element(indexLookUp[iRefl][iWave],
					                                                                indexLookUp[jRefl][jWave]);
}


template<typename T>
void
TPWALikelihood<T>::readIntegrals(const string& normIntFileName,  // name of file with normalization integrals
                                 const string& accIntFileName)   // name of file with acceptance integrals
{
	printInfo << "loading normalization integral from '" << normIntFileName << "'" << endl;
	ifstream intFile(normIntFileName.c_str());
	if (not intFile) {
		printErr << "cannot open file '" << normIntFileName << "'. aborting." << endl;
		throw;
	}
	integral normInt;
	// !!! integral.scan() performs no error checks!
	normInt.scan(intFile);
	intFile.close();
	reorderIntegralMatrix(normInt, _normMatrix);

	printInfo << "loading acceptance integral from '" << accIntFileName << "'" << endl;
	intFile.open(accIntFileName.c_str());
	if (not intFile) {
		printErr << "cannot open file '" << accIntFileName << "'. exiting." << endl;
		throw;
	}
	integral accInt;
	// !!! integral.scan() performs no error checks!
	accInt.scan(intFile);
	intFile.close();
	if (_numbAccEvents != 0)
		accInt.events(_numbAccEvents); 
	reorderIntegralMatrix(accInt, _accMatrix);
}


template<typename T>
void
TPWALikelihood<T>::readDecayAmplitudes(const string& ampDirName)
{
	// check that normalization integrals are loaded
	if (_normMatrix.num_elements() == 0) {
		printErr << "normalization integrals have to be loaded before loading the amplitudes. aborting." << endl;
		throw;
	}
	printInfo << "loading amplitude data" << endl;

	clear();

	// loop over amplitudes and read in data
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
			if (_debug)
				printInfo << "loading amplitude data for wave '" << _waveNames[iRefl][iWave] << "'" << endl;
			ifstream ampFile((ampDirName + "/" + _waveNames[iRefl][iWave]).c_str());
			if (not ampFile) {
				printErr << "cannot open file '" << _waveNames[iRefl][iWave] << "'. aborting." << endl;
				throw;
			}
			// get normalization
			const complex<double>    normInt = _normMatrix[iRefl][iWave][iRefl][iWave];
			vector<complex<double> > amps;
			complex<double>          amp;
			amps.reserve(_nmbEvents);  // number of events is known except for first wave that is read in
			while (ampFile.read((char*) &amp, sizeof(complex<double>))) {
				if (_useNormalizedAmps)         // normalize data, if option is switched on
					amp /= sqrt(normInt.real());  // rescale decay amplitude
				amps.push_back(amp);
			}
			_nmbEvents = amps.size(); 
			// copy decay amplitudes into array that is indexed [event index][reflectivity][wave index]
			// this index scheme ensures a more linear memory access pattern in the likelihood function
			_decayAmps.resize(extents[_nmbEvents][2][max(_nmbWavesRefl[0], _nmbWavesRefl[1])]);
			for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt)
				_decayAmps[iEvt][iRefl][iWave] = amps[iEvt];
			if (_debug)
				printInfo << "read " << _nmbEvents << " events from file "
				          << "'" << _waveNames[iRefl][iWave] << "'" << endl;
		}
	printInfo << "loaded " << _nmbEvents << " events into memory" << endl;

	// rescale integrals, if necessary
	if (_useNormalizedAmps) {
		// matrices _normMatrix and _accMatrix are already normalized to number of Monte Carlo events
		printInfo << "rescaling integrals." << endl;
		// rescale normalization and acceptance integrals
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				const double norm_i = sqrt(_normMatrix[iRefl][iWave][iRefl][iWave].real());
				for (unsigned int jRefl = 0; jRefl < 2; ++jRefl)
					for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
						const double norm_j = sqrt(_normMatrix[jRefl][jWave][jRefl][jWave].real());
						if ((iRefl != jRefl) or (iWave != jWave))
							// set diagonal terms later so that norm_i,j stay unaffected
							_normMatrix[iRefl][iWave][jRefl][jWave] /= norm_i * norm_j;
						_accMatrix[iRefl][iWave][jRefl][jWave] /= norm_i * norm_j;
					}
			}
		// set diagonal elements of normalization matrix
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
				_normMatrix[iRefl][iWave][iRefl][iWave] = 1;  // diagonal term

		if (_debug) {
			printInfo << "normalized integral matrices" << endl;
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					for (unsigned int jRefl = 0; jRefl < 2; ++jRefl)
						for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
							cout << "    normalization matrix [" << setw(2) << (int)iRefl * 2 - 1 << ", "
							     << setw(3) << iWave << ", " << setw(2) << (int)jRefl * 2 - 1 << ", "
							     << setw(3) << jWave << "] = "
							     << "("  << maxPrecisionAlign(_normMatrix[iRefl][iWave][jRefl][jWave].real())
							     << ", " << maxPrecisionAlign(_normMatrix[iRefl][iWave][jRefl][jWave].imag())
							     << "), acceptance matrix [" << setw(3) << iWave << ", "
							     << setw(3) << jWave << "] = "
							     << "("  << maxPrecisionAlign(_accMatrix[iRefl][iWave][jRefl][jWave].real())
							     << ", " << maxPrecisionAlign(_accMatrix[iRefl][iWave][jRefl][jWave].imag()) << ")"
							     << endl;
						}
		}
	}  // _useNormalizedAmps
}


template<typename T>
void
TPWALikelihood<T>::getIntCMatrix(TCMatrix& normMatrix,
                                 TCMatrix& accMatrix) const
{
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
			for (unsigned int jRefl = 0; jRefl < 2; ++jRefl)
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
					const unsigned int iIndex = _waveToWaveList[iRefl][iWave];
					const unsigned int jIndex = _waveToWaveList[jRefl][jWave];
					normMatrix.set(iIndex, jIndex, _normMatrix[iRefl][iWave][jRefl][jWave]);
					accMatrix.set (iIndex, jIndex, _accMatrix [iRefl][iWave][jRefl][jWave]);
				}
	// add flat
	normMatrix.set(_nmbWaves, _nmbWaves, 1);
	accMatrix.set (_nmbWaves, _nmbWaves, 1);
}


// builds complex numbers from parameters
// maps real and imaginary part of amplitudes to error matrix
// for both rank restrictions are taken into account
template<typename T>
void
TPWALikelihood<T>::buildCAmps(const double*             inPar,
                              vector<complex<double> >& prodAmps,
                              vector<pair<int, int> >&  parIndices,
                              vector<string>&           prodAmpNames,
                              const bool                withFlat) const
{
	prodAmps.clear();
	parIndices.clear();
	prodAmpNames.clear();
	unsigned int parIndex = 0;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				double re, im;
				if (iWave < iRank)  // zero production amplitude
					continue;
				else if (iWave == iRank) {  // real production amplitude
					parIndices.push_back(make_pair(parIndex, -1));
					re = inPar[parIndex];
					im = 0;
					++parIndex;
				} else {  // complex production amplitude
					parIndices.push_back(make_pair(parIndex, parIndex + 1));
					re = inPar[parIndex];
					im = inPar[parIndex + 1];
					parIndex += 2;
				}
				prodAmps.push_back(complex<double>(re, im));
				stringstream prodAmpName;
				prodAmpName << "V" << iRank << "_" << _waveNames[iRefl][iWave];
				prodAmpNames.push_back(prodAmpName.str());
			}
	if (withFlat) {
		prodAmps.push_back(complex<double>(inPar[parIndex], 0));
		parIndices.push_back(make_pair(parIndex, -1));
		prodAmpNames.push_back("V_flat");
	}
}


template<typename T>
void
TPWALikelihood<T>::clear()
{
	_decayAmps.resize(extents[0][0][0]);
}


// depends on naming convention for waves!!!
// VR_IGJPCMEIso....
template<typename T>
int
TPWALikelihood<T>::getReflectivity(const TString& waveName)
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
		printErr << "cannot parse parameter/wave name '" << waveName << "'. "
		         << "cannot not determine reflectivity. aborting." << endl;
		throw;
	}
	if (_debug)
		printInfo << "extracted reflectivity = " << refl << " from parameter name "
		          << "'" << waveName << "' (char position " << reflIndex << ")" << endl;
	return refl;
}


// copy values from array that corresponds to the function parameters
// to structure that corresponds to the complex production amplitudes
// taking into account rank restrictions
template<typename T>
void
TPWALikelihood<T>::copyFromParArray(const double*  inPar,             // input parameter array
                                    ampsArrayType& outVal,            // array of complex output values [rank][reflectivity][wave index]
                                    T&             outFlatVal) const  // output value corresponding to flat wave
{
	// group parameters by rank and wave index only; keeps the order defined in wave list
	outVal.resize(extents[_rank][2][max(_nmbWavesRefl[0], _nmbWavesRefl[1])]);
	unsigned int parIndex = 0;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				T re, im;
				if (iWave < iRank)          // production amplitude is zero
					re = im = 0;
				else if (iWave == iRank) {  // production amplitude is real
					re = inPar[parIndex];
					im = 0;
					++parIndex;
				} else {                    // production amplitude is complex
					re = inPar[parIndex];
					im = inPar[parIndex + 1];
					parIndex += 2;
				}
				outVal[iRank][iRefl][iWave] = complex<T>(re, im);
			}
	outFlatVal = inPar[parIndex];
}


// copy values from structure that corresponds to complex
// production amplitudes to array that corresponds to function
// parameters taking into account rank restrictions
template<typename T>
void
TPWALikelihood<T>::copyToParArray(const ampsArrayType& inVal,         // values corresponding to production amplitudes
                                  const T              inFlatVal,     // value corresponding to flat wave
                                  double*              outPar) const  // output parameter array
{
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				tuple<int, int> parIndices = _prodAmpToFuncParMap[iRank][iRefl][iWave];
				if (get<0>(parIndices) >= 0)  // real part
					outPar[get<0>(parIndices)] = inVal[iRank][iRefl][iWave].real();
				if (get<1>(parIndices) >= 0)  // imaginary part
					outPar[get<1>(parIndices)] = inVal[iRank][iRefl][iWave].imag();
			}
	outPar[_nmbPars - 1] = inFlatVal;
}


template<typename T>
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
#ifdef USE_CUDA
	    << "use CUDA kernels ........................ " << _cudaEnabled       << endl
#endif	  
	    << "use normalized amplitudes ............... " << _useNormalizedAmps << endl
	    << "list of waves: " << endl;
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
			out << "        [" << setw(2) << sign(iRefl) << " " << setw(3) << iWave << "] "
			    << _waveNames[iRefl][iWave] << "    threshold = "
			    << _waveThresholds[iRefl][iWave] << " MeV/c^2" << endl;
	out << "list of function parameters: " << endl;
	for (unsigned int iPar = 0; iPar < _nmbPars; ++iPar)
		out << "        [" << setw(3) << iPar << "] " << _parNames[iPar] << "    "
		    << "threshold = " << _parThresholds[iPar] << " MeV/c^2" << endl;
	return out;
}


template<typename T>
ostream&
TPWALikelihood<T>::printFuncInfo(ostream& out) const
{
	const string funcNames[NMB_FUNCTIONCALLENUM] = {"FdF", "Gradient", "DoEval", "DoDerivative"};
	for (unsigned int i = 0; i < NMB_FUNCTIONCALLENUM; ++i)
		if (_funcCallInfo[i].nmbCalls > 0)
			out << "    " << _funcCallInfo[i].nmbCalls
			    << " calls to TPWALikelihood<T>::" << funcNames[i] << "()" << endl
			    << "    time spent in TPWALikelihood<T>::" << funcNames[i] << "(): " << endl
			    << "        total time ................. " << _funcCallInfo[i].totalTime << " sec" << endl
			    << "        return value calculation ... " << _funcCallInfo[i].funcTime  << " sec" << endl
			    << "        normalization .............. " << _funcCallInfo[i].normTime  << " sec" << endl;
	return out;
}


template<typename T>
vector<unsigned int>
TPWALikelihood<T>::orderedParIndices() const
{
	vector<unsigned int> orderedIndices;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int waveIndex = 0; waveIndex < _nmbWaves; ++waveIndex) {
			unsigned int iRefl, iWave;
			for (iRefl = 0; iRefl < 2; ++iRefl)
				for (iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					if (_waveToWaveList[iRefl][iWave] == waveIndex)
						goto found;
			printWarn << "indices are inconsistent. cannot find wave with index " << waveIndex
			          << " in wave list" << endl;
			continue;
      
		found:
      
			tuple<int, int> parIndices = _prodAmpToFuncParMap[iRank][iRefl][iWave];
			int parIndex = get<0>(parIndices);
			if (parIndex >= 0)
				orderedIndices.push_back(parIndex);
			parIndex = get<1>(parIndices);
			if (parIndex >= 0)
				orderedIndices.push_back(parIndex);
		}
	orderedIndices.push_back(_nmbPars - 1);  // flat wave
	if (orderedIndices.size() != _nmbPars)
		printWarn << "ordered list of parameter indices has inconsistent size "
		          << "(" << orderedIndices.size() << " vs. " << _nmbPars << ")" << endl;
	return orderedIndices;
}


template<typename T>
vector<string>
TPWALikelihood<T>::waveNames() const
{
	vector<string> names(_nmbWaves, "");
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
			const unsigned int index = _waveToWaveList[iRefl][iWave];
			names[index] = _waveNames[iRefl][iWave];
		}
	return names;
}


template<typename T>
void
TPWALikelihood<T>::resetFuncCallInfo() const
{
	for (unsigned int i = 0; i < NMB_FUNCTIONCALLENUM; ++i) {
		_funcCallInfo[i].nmbCalls   = 0;
		_funcCallInfo[i].funcTime   = 0;
		_funcCallInfo[i].normTime   = 0;
		_funcCallInfo[i].totalTime  = 0;
	}
}
