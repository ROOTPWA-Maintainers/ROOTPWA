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
//      see TPWALikelihood.h for details
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
#include <complex>
#include <cassert>

#include "TString.h"
#include "TSystem.h"
#include "TCMatrix.h"
#include "TStopwatch.h"
#include "TTree.h"

#include "reportingUtils.hpp"
#include "conversionUtils.hpp"
#include "fileUtils.hpp"
#ifdef USE_CUDA
#include "complex.cuh"
#include "likelihoodInterface.cuh"
#endif
#include "amplitudeTreeLeaf.h"
#include "TPWALikelihood.h"


// #define USE_FDF


using namespace std;
using namespace rpwa;
using namespace boost;
using namespace boost::accumulators;


template<typename complexT> bool TPWALikelihood<complexT>::_debug = true;


template<typename complexT>
TPWALikelihood<complexT>::TPWALikelihood()
	: _nmbEvents        (0),
	  _rank             (1),
	  _nmbWaves         (0),
	  _nmbWavesReflMax  (0),
	  _nmbPars          (0),
#ifdef USE_CUDA
	  _cudaEnabled      (false),
#endif
	  _useNormalizedAmps(true),
	  _numbAccEvents    (0)
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


template<typename complexT>
TPWALikelihood<complexT>::~TPWALikelihood()
{
	clear();
}


template<typename complexT>
void
TPWALikelihood<complexT>::FdF
(const double* par,             // parameter array; reduced by rank conditions
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
	value_type    prodAmpFlat;
	ampsArrayType prodAmps;
	copyFromParArray(par, prodAmps, prodAmpFlat);
	const value_type prodAmpFlat2 = prodAmpFlat * prodAmpFlat;

	// create array of likelihood derivative w.r.t. real and imaginary
	// parts of the production amplitudes
	// !NOTE! although stored as and constructed from complex values,
	// the dL themselves are _not_ well defined complex numbers!
	value_type                              derivativeFlat = 0;
	array<typename ampsArrayType::index, 3> derivShape     = {{ _rank, 2, _nmbWavesReflMax }};
	ampsArrayType                           derivatives(derivShape);

	// loop over events and calculate real-data term of log likelihood
	// as well as derivatives with respect to parameters
	TStopwatch timer;
	timer.Start();
	accumulator_set<value_type, stats<tag::sum(compensated)> > logLikelihoodAcc;
	accumulator_set<value_type, stats<tag::sum(compensated)> > derivativeFlatAcc;
	multi_array<accumulator_set<complexT, stats<tag::sum(compensated)> >, 3>
		derivativesAcc(derivShape);
	for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
		accumulator_set<value_type, stats<tag::sum(compensated)> > likelihoodAcc;
		ampsArrayType derivative(derivShape);  // likelihood derivative for this event
		for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
				accumulator_set<complexT, stats<tag::sum(compensated)> > ampProdAcc;
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)  // coherent sum over waves
					ampProdAcc(prodAmps[iRank][iRefl][iWave] * _decayAmps[iEvt][iRefl][iWave]);
				const complexT ampProdSum = sum(ampProdAcc);
				likelihoodAcc(norm(ampProdSum));
				// set derivative term that is independent on derivative wave index
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave)
					// amplitude sums for current rank and for waves with same reflectivity
					derivative[iRank][iRefl][jWave] = ampProdSum;
			}
			// loop again over waves for current rank and multiply with complex conjugate
			// of decay amplitude of the wave with the derivative wave index
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					derivative[iRank][iRefl][iWave] *= conj(_decayAmps[iEvt][iRefl][iWave]);
		}  // end loop over rank
		likelihoodAcc   (prodAmpFlat2            );
		logLikelihoodAcc(-log(sum(likelihoodAcc)));
		// incorporate factor 2 / sigma
		const value_type factor = 2. / sum(likelihoodAcc);
		for (unsigned int iRank = 0; iRank < _rank; ++iRank)
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					derivativesAcc[iRank][iRefl][iWave](-factor * derivative[iRank][iRefl][iWave]);
		derivativeFlatAcc(-factor * prodAmpFlat);
	}  // end loop over events
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
				derivatives[iRank][iRefl][iWave] = sum(derivativesAcc[iRank][iRefl][iWave]);
	derivativeFlat = sum(derivativeFlatAcc);
	// log time needed for likelihood calculation
	timer.Stop();
	_funcCallInfo[FDF].funcTime(timer.RealTime());

	// compute normalization term of log likelihood and normalize derivatives w.r.t. parameters
	timer.Start();
	accumulator_set<value_type, stats<tag::sum(compensated)> > normFactorAcc;
	const value_type nmbEvt      = (_useNormalizedAmps) ? 1 : _nmbEvents;
	const value_type twiceNmbEvt = 2 * nmbEvt;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				accumulator_set<complexT, stats<tag::sum(compensated)> > normFactorDerivAcc;
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave) {  // inner loop over waves with same reflectivity
					const complexT I = _accMatrix[iRefl][iWave][iRefl][jWave];
					normFactorAcc(real((prodAmps[iRank][iRefl][iWave] * conj(prodAmps[iRank][iRefl][jWave]))
					                   * I));
					normFactorDerivAcc(prodAmps[iRank][iRefl][jWave] * conj(I));
				}
				derivatives[iRank][iRefl][iWave] += sum(normFactorDerivAcc) * twiceNmbEvt;  // account for 2 * nmbEvents
			}
	// take care of flat wave
	normFactorAcc(prodAmpFlat2 * _totAcc);
	derivativeFlat += prodAmpFlat * twiceNmbEvt * _totAcc;
	// log time needed for normalization
	timer.Stop();
	_funcCallInfo[FDF].normTime(timer.RealTime());
	
	// sort derivative results into output array and cache
	copyToParArray(derivatives, derivativeFlat, gradient);
	copyToParArray(derivatives, derivativeFlat, toArray(_derivCache));
	
	// set function return value
	funcVal = sum(logLikelihoodAcc) + nmbEvt * sum(normFactorAcc);

	// log total consumed time
	timerTot.Stop();
	_funcCallInfo[FDF].totalTime(timerTot.RealTime());
	
	if (_debug)
		printInfo << "log likelihood =  "       << maxPrecisionAlign(sum(logLikelihoodAcc)) << ", "
		          << "normalization =  "        << maxPrecisionAlign(sum(normFactorAcc)   ) << ", "
		          << "normalized likelihood = " << maxPrecisionAlign(funcVal              ) << endl;
}


template<typename complexT>
double
TPWALikelihood<complexT>::DoEval(const double* par) const
{
	++(_funcCallInfo[DOEVAL].nmbCalls);
	
#ifdef USE_FDF

	// call FdF
	double logLikelihood;
	double gradient[_nmbPars];
	FdF(par, logLikelihood, gradient);
	return logLikelihood;

#else  // USE_FDF

	// timer for total time
	TStopwatch timerTot;
	timerTot.Start();
  
	// build complex production amplitudes from function parameters taking into account rank restrictions
	value_type    prodAmpFlat;
	ampsArrayType prodAmps;
	copyFromParArray(par, prodAmps, prodAmpFlat);
	const value_type prodAmpFlat2 = prodAmpFlat * prodAmpFlat;

	// loop over events and calculate real-data term of log likelihood
	TStopwatch timer;
	timer.Start();
	value_type logLikelihood = 0;
#ifdef USE_CUDA
	if (_cudaEnabled) {
		logLikelihood = cuda::likelihoodInterface<cuda::complex<value_type> >::logLikelihood
			(reinterpret_cast<cuda::complex<value_type>*>(prodAmps.data()),
			 prodAmps.num_elements(), prodAmpFlat, _rank);
	} else
#else  // USE_CUDA
	{
		accumulator_set<value_type, stats<tag::sum(compensated)> > logLikelihoodAcc;
		for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
			accumulator_set<value_type, stats<tag::sum(compensated)> > likelihoodAcc;
			for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
					accumulator_set<complexT, stats<tag::sum(compensated)> > ampProdAcc;
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
						ampProdAcc(prodAmps[iRank][iRefl][iWave] * _decayAmps[iEvt][iRefl][iWave]);
						// cout << "prodAmps[" << iRank << "][" << iRefl << "][" << iWave<< "] = "
						//      << maxPrecisionDouble(prodAmps[iRank][iRefl][iWave]) << "; "
						//      << "decayAmps[" << iEvt << "][" << iRefl << "][" << iWave << "] = "
						//      << maxPrecisionDouble(_decayAmps[iEvt][iRefl][iWave]) << endl;
					}
					likelihoodAcc(norm(sum(ampProdAcc)));
				}
			}
			likelihoodAcc(prodAmpFlat2);
			logLikelihoodAcc(-log(sum(likelihoodAcc)));
			// cout << endl;
		}
		logLikelihood = sum(logLikelihoodAcc);
	}
#endif  // USE_CUDA
	// log time needed for likelihood calculation
	timer.Stop();
	_funcCallInfo[DOEVAL].funcTime(timer.RealTime());
	
	// compute normalization term of log likelihood
	timer.Start();
	accumulator_set<value_type, stats<tag::sum(compensated)> > normFactorAcc;
	const value_type nmbEvt = (_useNormalizedAmps) ? 1 : _nmbEvents;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave)  // loop over waves with same reflectivity
					normFactorAcc(real((prodAmps[iRank][iRefl][iWave] * conj(prodAmps[iRank][iRefl][jWave]))
					                   * _accMatrix[iRefl][iWave][iRefl][jWave]));
	// take care of flat wave
	normFactorAcc(prodAmpFlat2 * _totAcc);
	// log time needed for normalization
	timer.Stop();
	_funcCallInfo[DOEVAL].normTime(timer.RealTime());
  
	// calculate and return log likelihood value
	const double funcVal = logLikelihood + nmbEvt * sum(normFactorAcc);

	// log total consumed time
	timerTot.Stop();
	_funcCallInfo[DOEVAL].totalTime(timerTot.RealTime());

	if (_debug)
		printInfo << "raw log likelihood =  "       << maxPrecisionAlign(logLikelihood     ) << ", "
		          << "normalization =  "            << maxPrecisionAlign(sum(normFactorAcc)) << ", "
		          << "normalized log likelihood = " << maxPrecisionAlign(funcVal           ) << endl;
	
	return funcVal;
	
#endif  // USE_FDF
}

	
template<typename complexT>
double
TPWALikelihood<complexT>::DoDerivative(const double* par,
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
		_funcCallInfo[DODERIVATIVE].totalTime(timerTot.RealTime());
		_funcCallInfo[DODERIVATIVE].funcTime (sum(_funcCallInfo[DODERIVATIVE].totalTime));
		return _derivCache[derivativeIndex];
	}
	// call FdF
	double logLikelihood;
	double gradient[_nmbPars];
	FdF(par, logLikelihood, gradient);
	return gradient[derivativeIndex];
}
 
 
// calculate derivatives with respect to parameters
template<typename complexT>
void
TPWALikelihood<complexT>::Gradient
(const double* par,             // parameter array; reduced by rank conditions
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
		_funcCallInfo[GRADIENT].totalTime(timerTot.RealTime());
		_funcCallInfo[GRADIENT].funcTime (sum(_funcCallInfo[GRADIENT].totalTime));
		return;
	}
	// call FdF
	double logLikelihood;
	FdF(par, logLikelihood, gradient);

#else  // USE_FDF

	// copy arguments into parameter cache
	for (unsigned int i = 0; i < _nmbPars; ++i)
		_parCache[i] = par[i];

	// build complex production amplitudes from function parameters taking into account rank restrictions
	value_type    prodAmpFlat;
	ampsArrayType prodAmps;
	copyFromParArray(par, prodAmps, prodAmpFlat);

	// create array of likelihood derivative w.r.t. real and imaginary
	// parts of the production amplitudes
	// !NOTE! although stored as and constructed from complex values,
	// the dL themselves are _not_ well defined complex numbers!
	value_type                              derivativeFlat = 0;
	array<typename ampsArrayType::index, 3> derivShape     = {{ _rank, 2, _nmbWavesReflMax }};
	ampsArrayType                           derivatives(derivShape);

	// loop over events and calculate derivatives with respect to parameters
	TStopwatch timer;
	timer.Start();
#ifdef USE_CUDA
	if (_cudaEnabled) {
		cuda::likelihoodInterface<cuda::complex<value_type> >::logLikelihoodDeriv
			(reinterpret_cast<cuda::complex<value_type>*>(prodAmps.data()),
			 prodAmps.num_elements(), prodAmpFlat, _rank,
			 reinterpret_cast<cuda::complex<value_type>*>(derivatives.data()),
			 derivativeFlat);
	} else
#else  // USE_CUDA
	{
		accumulator_set<value_type, stats<tag::sum(compensated)> > derivativeFlatAcc;
		multi_array<accumulator_set<complexT, stats<tag::sum(compensated)> >, 3>
			derivativesAcc(derivShape);
		const value_type prodAmpFlat2 = prodAmpFlat * prodAmpFlat;
		for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
			accumulator_set<value_type, stats<tag::sum(compensated)> > likelihoodAcc;
			ampsArrayType derivative(derivShape);  // likelihood derivatives for this event
			for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
					accumulator_set<complexT, stats<tag::sum(compensated)> > ampProdAcc;
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)  // coherent sum over waves
						ampProdAcc(prodAmps[iRank][iRefl][iWave] * _decayAmps[iEvt][iRefl][iWave]);
					const complexT ampProdSum = sum(ampProdAcc);
					likelihoodAcc(norm(ampProdSum));
					// set derivative term that is independent on derivative wave index
					for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave)
						// amplitude sums for current rank and for waves with same reflectivity
						derivative[iRank][iRefl][jWave] = ampProdSum;
				}
				// loop again over waves for current rank and multiply with complex conjugate
				// of decay amplitude of the wave with the derivative wave index
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
					for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave)
						derivative[iRank][iRefl][jWave] *= conj(_decayAmps[iEvt][iRefl][jWave]);
			}  // end loop over rank
			likelihoodAcc(prodAmpFlat2);
			// incorporate factor 2 / sigma
			const value_type factor = 2. / sum(likelihoodAcc);
			for (unsigned int iRank = 0; iRank < _rank; ++iRank)
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
					for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave)
						derivativesAcc[iRank][iRefl][jWave](-factor * derivative[iRank][iRefl][jWave]);
			derivativeFlatAcc(-factor * prodAmpFlat);
		}  // end loop over events
		for (unsigned int iRank = 0; iRank < _rank; ++iRank)
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					derivatives[iRank][iRefl][iWave] = sum(derivativesAcc[iRank][iRefl][iWave]);
		derivativeFlat = sum(derivativeFlatAcc);
	}
#endif  // USE_CUDA
	// log time needed for gradient calculation
	timer.Stop();
	_funcCallInfo[GRADIENT].funcTime(timer.RealTime());

	// normalize derivatives w.r.t. parameters
	timer.Start();
	const value_type nmbEvt      = (_useNormalizedAmps) ? 1 : _nmbEvents;
	const value_type twiceNmbEvt = 2 * nmbEvt;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				accumulator_set<complexT, stats<tag::sum(compensated)> > normFactorDerivAcc;
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave) {  // inner loop over waves with same reflectivity
					const complexT I = _accMatrix[iRefl][iWave][iRefl][jWave];
					normFactorDerivAcc(prodAmps[iRank][iRefl][jWave] * conj(I));
				}
				derivatives[iRank][iRefl][iWave] += sum(normFactorDerivAcc) * twiceNmbEvt;  // account for 2 * nmbEvents
			}
	// take care of flat wave
	derivativeFlat += prodAmpFlat * twiceNmbEvt * _totAcc;
	// log time needed for normalization
	timer.Stop();
	_funcCallInfo[GRADIENT].normTime(timer.RealTime());

	// set return gradient values
	copyToParArray(derivatives, derivativeFlat, gradient);
	copyToParArray(derivatives, derivativeFlat, toArray(_derivCache));

	// log total consumed time
	timerTot.Stop();
	_funcCallInfo[GRADIENT].totalTime(timerTot.RealTime());

#endif  // USE_FDF
}


template<typename complexT>
unsigned int
TPWALikelihood<complexT>::nmbWaves(const int reflectivity) const
{
	if (reflectivity == 0)
		return _nmbWaves;
	else if (reflectivity > 0)
		return _nmbWavesRefl[1];  // positive reflectivity
	else
		return _nmbWavesRefl[0];  // negative reflectivity
}


template<typename complexT>
void
#ifdef USE_CUDA
TPWALikelihood<complexT>::enableCuda(const bool enableCuda)
{
	_cudaEnabled = enableCuda;
}
#else
TPWALikelihood<complexT>::enableCuda(const bool) { }
#endif


template<typename complexT>
bool
TPWALikelihood<complexT>::cudaEnabled() const
{
#ifdef USE_CUDA
	return _cudaEnabled;
#else
	return false;
#endif
}


template<typename complexT>
void 
TPWALikelihood<complexT>::init(const unsigned int rank,
                               const std::string& waveListFileName,
                               const std::string& normIntFileName,
                               const std::string& accIntFileName,
                               const std::string& ampDirName,
                               const unsigned int numbAccEvents,
                               const bool         useRootAmps)
{
	_numbAccEvents = numbAccEvents;
	readWaveList(waveListFileName);
	buildParDataStruct(rank);
	readIntegrals(normIntFileName, accIntFileName);
	readDecayAmplitudes(ampDirName, useRootAmps);
#ifdef USE_CUDA
	if (_cudaEnabled)
		cuda::likelihoodInterface<cuda::complex<value_type> >::init
			(reinterpret_cast<cuda::complex<value_type>*>(_decayAmps.data()),
			 _decayAmps.num_elements(), _nmbEvents, _nmbWavesRefl, true);
#endif
}


template<typename complexT>
void 
TPWALikelihood<complexT>::readWaveList(const string& waveListFileName)
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
	_nmbWaves        = _nmbWavesRefl[0] + _nmbWavesRefl[1];
	_nmbWavesReflMax = max(_nmbWavesRefl[0], _nmbWavesRefl[1]);
	_waveNames.resize      (extents[2][_nmbWavesReflMax]);
	_waveThresholds.resize (extents[2][_nmbWavesReflMax]);
	_waveToWaveIndex.resize(extents[2][_nmbWavesReflMax]);
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
			_waveNames      [iRefl][iWave] = waveNames     [iRefl][iWave];
			_waveThresholds [iRefl][iWave] = waveThresholds[iRefl][iWave];
			_waveToWaveIndex[iRefl][iWave] = waveIndices   [iRefl][iWave];
		}
}


template<typename complexT>
void
TPWALikelihood<complexT>::buildParDataStruct(const unsigned int rank)
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
	_prodAmpToFuncParMap.resize(extents[_rank][2][_nmbWavesReflMax]);
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
template<typename complexT>
void
TPWALikelihood<complexT>::reorderIntegralMatrix(integral&            integral,
                                                normMatrixArrayType& reorderedMatrix) const
{
	// get original matrix and list of wave names
	const matrix<complex<double> > intMatrix    = integral.mat();
	list<string>                   intWaveNames = integral.files();
  // "int" saves filenames with path, this must be treated here
  for (list<string>::iterator it = intWaveNames.begin(); it != intWaveNames.end(); ++it) {
	  // find the slash, if not available -> take the first position of the string
	  const size_t slashPos = it->rfind('/');
	  if (slashPos != string::npos)
		  it->erase(0, slashPos + 1);
  }
	// build index lookup-table [reflectivity][wave index] to index in normalization integral
	waveToIntMapType indexLookUp(extents[2][_nmbWavesReflMax]);
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
			if (find(intWaveNames.begin(), intWaveNames.end(), _waveNames[iRefl][iWave])
			    == intWaveNames.end()) {
				printErr << "wave " << _waveNames[iRefl][iWave] << " is not in integral. aborting." << endl;
				throw;
			}
			indexLookUp[iRefl][iWave] = integral.index(_waveNames[iRefl][iWave]);
			if (_debug)
				printInfo << "    mapping wave [" << sign((int)iRefl * 2 - 1) << ", "
				          << setw(3) << iWave << "] '" << _waveNames[iRefl][iWave] << "' "
				          << "to index " << setw(3) << indexLookUp[iRefl][iWave] << " in integral." << endl;
		}
	// create reordered matrix
	reorderedMatrix.resize(extents[2][_nmbWavesReflMax][2][_nmbWavesReflMax]);
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
			for (unsigned int jRefl = 0; jRefl < 2; ++jRefl)
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
					const complex<double> val = intMatrix.element(indexLookUp[iRefl][iWave],
					                                              indexLookUp[jRefl][jWave]);
					reorderedMatrix[iRefl][iWave][jRefl][jWave] = complexT(val.real(), val.imag());
				}
}


// returns integral matrix reordered according to _waveNames array
template<typename complexT>
void
TPWALikelihood<complexT>::reorderIntegralMatrix(const normalizationIntegral& integral,
                                                normMatrixArrayType&         reorderedMatrix) const
{
	// create reordered matrix
	reorderedMatrix.resize(extents[2][_nmbWavesReflMax][2][_nmbWavesReflMax]);
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
			for (unsigned int jRefl = 0; jRefl < 2; ++jRefl)
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
					const complex<double> val = integral.element(_waveNames[iRefl][iWave],
					                                             _waveNames[jRefl][jWave]);
					reorderedMatrix[iRefl][iWave][jRefl][jWave] = complexT(val.real(), val.imag());
				}
}


template<typename complexT>
void
TPWALikelihood<complexT>::readIntegrals
(const string& normIntFileName,   // name of file with normalization integrals
 const string& accIntFileName,    // name of file with acceptance integrals
 const string& integralTKeyName)  // name of TKey which stores integral in .root file
{
	printInfo << "loading normalization integral from '" << normIntFileName << "'" << endl;
	const string normIntFileExt  = extensionFromPath(normIntFileName);
#if NORMALIZATIONINTEGRAL_ENABLED
	if (normIntFileExt == "root") {
		TFile* intFile  = TFile::Open(normIntFileName.c_str(), "READ");
		if (not intFile or intFile->IsZombie()) {
			printErr << "could open normalization integral file '" << normIntFileName << "'. "
			         << "aborting." << endl;
			throw;
		}
		normalizationIntegral* integral = 0;
		intFile->GetObject(integralTKeyName.c_str(), integral);
		if (not integral) {
			printErr << "cannot find integral object in TKey '" << integralTKeyName << "' in file "
			         << "'" << normIntFileName << "'. aborting." << endl;
			throw;
		}
		reorderIntegralMatrix(*integral, _normMatrix);
		intFile->Close();
	} else
#endif  // NORMALIZATIONINTEGRAL_ENABLED
		if (normIntFileExt == "int") {
		ifstream intFile(normIntFileName.c_str());
		if (not intFile) {
			printErr << "cannot open file '" << normIntFileName << "'. aborting." << endl;
			throw;
		}
		integral integral;
		// !!! integral.scan() performs no error checks!
		integral.scan(intFile);
		intFile.close();
		reorderIntegralMatrix(integral, _normMatrix);
	} else {
		printErr << "unknown file type '" << normIntFileName << "'. "
		         << "only .int and .root files are supported. aborting." << endl;
		throw;
	}

	printInfo << "loading acceptance integral from '" << accIntFileName << "'" << endl;
	const string accIntFileExt  = extensionFromPath(accIntFileName);
#if NORMALIZATIONINTEGRAL_ENABLED
	if (accIntFileExt == "root") {
		TFile* intFile  = TFile::Open(accIntFileName.c_str(), "READ");
		if (not intFile or intFile->IsZombie()) {
			printErr << "could open normalization integral file '" << accIntFileName << "'. "
			         << "aborting." << endl;
			throw;
		}
		normalizationIntegral* integral = 0;
		intFile->GetObject(integralTKeyName.c_str(), integral);
		if (not integral) {
			printErr << "cannot find integral object in TKey '" << integralTKeyName << "' in file "
			         << "'" << accIntFileName << "'. aborting." << endl;
			throw;
		}
		if (_numbAccEvents != 0) {
			_totAcc = ((double)integral->nmbEvents()) / (double)_numbAccEvents;
			printInfo << "total acceptance in this bin: " << _totAcc << endl;
			integral->setNmbEvents(_numbAccEvents);
		} else
			_totAcc = 1;
		reorderIntegralMatrix(*integral, _accMatrix);
		intFile->Close();
	} else
#endif  // NORMALIZATIONINTEGRAL_ENABLED
		if (accIntFileExt == "int") {
		ifstream intFile(accIntFileName.c_str());
		if (not intFile) {
			printErr << "cannot open file '" << accIntFileName << "'. aborting." << endl;
			throw;
		}
		integral integral;
		// !!! integral.scan() performs no error checks!
		integral.scan(intFile);
		intFile.close();
		if (_numbAccEvents != 0) {
			_totAcc = ((double)integral.nevents()) / (double)_numbAccEvents;
			printInfo << "total acceptance in this bin: " << _totAcc << endl;
			integral.events(_numbAccEvents);
		} else
			_totAcc = 1;
		reorderIntegralMatrix(integral, _accMatrix);
	} else {
		printErr << "unknown file type '" << accIntFileName << "'. "
		         << "only .int and .root files are supported. exiting." << endl;
		throw;
	}
}


template<typename complexT>
void
TPWALikelihood<complexT>::readDecayAmplitudes(const string& ampDirName,
                                              const bool    useRootAmps,
                                              const string& ampLeafName)
{
	// check that normalization integrals are loaded
	if (_normMatrix.num_elements() == 0) {
		printErr << "normalization integrals have to be loaded before loading the amplitudes. "
		         << "aborting." << endl;
		throw;
	}
	clear();

	printInfo << "loading amplitude data from " << ((useRootAmps) ? ".root" : ".amp")
	          << " files" << endl;
	// loop over amplitudes and read in data
	unsigned int nmbEvents = 0;
	bool         firstWave = true;
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
			// get normalization
			const complexT   normInt = _normMatrix[iRefl][iWave][iRefl][iWave];
			vector<complexT> amps;
			if (!firstWave)  // number of events is known except for first wave that is read in
				amps.reserve(nmbEvents);
			// read decay amplitudes
			string ampFilePath = ampDirName + "/" + _waveNames[iRefl][iWave];
#if AMPLITUDETREELEAF_ENABLED
			if (useRootAmps) {
				ampFilePath = changeFileExtension(ampFilePath, ".root");
				printInfo << "loading amplitude data from '" << ampFilePath << "'" << endl;
				// open amplitude file
				TFile* ampFile = TFile::Open(ampFilePath.c_str(), "READ");
				if (not ampFile or ampFile->IsZombie()) {
					printWarn << "cannot open amplitude file '" << ampFilePath << "'. skipping." << endl;
					continue;
				}
				// find amplitude tree
				TTree*       ampTree     = 0;
				const string ampTreeName = changeFileExtension(_waveNames[iRefl][iWave], ".amp");
				ampFile->GetObject(ampTreeName.c_str(), ampTree);
				if (not ampTree) {
					printWarn << "cannot find tree '" << ampTreeName << "' in file "
					          << "'" << ampFilePath << "'. skipping." << endl;
					continue;
				}
				// connect tree leaf
				amplitudeTreeLeaf* ampTreeLeaf = 0;
				ampTree->SetBranchAddress(ampLeafName.c_str(), &ampTreeLeaf);
				for (long int eventIndex = 0; eventIndex < ampTree->GetEntriesFast(); ++eventIndex) {
					ampTree->GetEntry(eventIndex);
					if (!ampTreeLeaf) {
						printWarn << "null pointer to amplitude leaf for event " << eventIndex << ". "
						          << "skipping." << endl;
						continue;
					}
					assert(ampTreeLeaf->nmbIncohSubAmps() == 1);
					complexT amp(ampTreeLeaf->incohSubAmp(0).real(), ampTreeLeaf->incohSubAmp(0).imag());
					if (_useNormalizedAmps)         // normalize data, if option is switched on
						amp /= sqrt(normInt.real());  // rescale decay amplitude
					amps.push_back(amp);
				}
			} else
#endif
			{
				printInfo << "loading amplitude data from '" << ampFilePath << "'" << endl;
				ifstream ampFile(ampFilePath.c_str());
				if (not ampFile) {
					printErr << "cannot open amplitude file '" << ampFilePath << "'. aborting." << endl;
					throw;
				}
				complexT amp;
				while (ampFile.read((char*)&amp, sizeof(complexT))) {
					if (_useNormalizedAmps)         // normalize data, if option is switched on
						amp /= sqrt(normInt.real());  // rescale decay amplitude
					amps.push_back(amp);
				}
			}
			if (firstWave)
				nmbEvents = _nmbEvents = amps.size();
			else {
				nmbEvents = amps.size();
				if (nmbEvents != _nmbEvents)
					printWarn << "size mismatch in amplitude files: this file contains " << nmbEvents
					          << " events, previous file had " << _nmbEvents << " events." << endl;
				nmbEvents = _nmbEvents;
			}
				
			// copy decay amplitudes into array that is indexed [event index][reflectivity][wave index]
			// this index scheme ensures a more linear memory access pattern in the likelihood function
			_decayAmps.resize(extents[_nmbEvents][2][_nmbWavesReflMax]);
			for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt)
				_decayAmps[iEvt][iRefl][iWave] = amps[iEvt];
			if (_debug)
				printInfo << "read " << _nmbEvents << " events from file "
				          << "'" << _waveNames[iRefl][iWave] << "'" << endl;
		}
	printInfo << "loaded decay amplitudes for " << _nmbEvents << " events into memory" << endl;

	// save phase space integrals
	_phaseSpaceIntegral.resize(extents[2][_nmbWavesReflMax]);
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
			_phaseSpaceIntegral[iRefl][iWave] = sqrt(_normMatrix[iRefl][iWave][iRefl][iWave].real());

	// rescale integrals, if necessary
	if (_useNormalizedAmps) {
		// matrices _normMatrix and _accMatrix are already normalized to number of Monte Carlo events
		printInfo << "rescaling integrals" << endl;
		// rescale normalization and acceptance integrals
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				const value_type norm_i = sqrt(_normMatrix[iRefl][iWave][iRefl][iWave].real());
				for (unsigned int jRefl = 0; jRefl < 2; ++jRefl)
					for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
						const value_type norm_j = sqrt(_normMatrix[jRefl][jWave][jRefl][jWave].real());
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
							cout << "    normalization matrix [" << sign((int)iRefl * 2 - 1) << ", "
							     << setw(3) << iWave << ", " << sign((int)jRefl * 2 - 1) << ", "
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


template<typename complexT>
void
TPWALikelihood<complexT>::getIntegralMatrices(TCMatrix&       normMatrix,
                                              TCMatrix&       accMatrix,
                                              vector<double>& phaseSpaceIntegral) const
{
  phaseSpaceIntegral.clear();
  phaseSpaceIntegral.resize(_nmbWaves + 1, 0);
  unsigned int iIndex = 0;
  for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
			phaseSpaceIntegral[iIndex] = _phaseSpaceIntegral[iRefl][iWave];
			unsigned int jIndex = 0;
			for (unsigned int jRefl = 0; jRefl < 2; ++jRefl) {
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
					const complexT     normVal = _normMatrix[iRefl][iWave][jRefl][jWave];
					const complexT     accVal  = _accMatrix [iRefl][iWave][jRefl][jWave];
					normMatrix.set(iIndex, jIndex, complex<double>(normVal.real(), normVal.imag()));
					accMatrix.set (iIndex, jIndex, complex<double>(accVal.real(),  accVal.imag() ));
					++jIndex;
				}
			}
			++iIndex;
		}
  }
	// add flat
	normMatrix.set(_nmbWaves, _nmbWaves, 1);
	accMatrix.set (_nmbWaves, _nmbWaves, 1);
  phaseSpaceIntegral[_nmbWaves] = 1;
}


// builds complex numbers from parameters
// maps real and imaginary part of amplitudes to error matrix
// for both rank restrictions are taken into account
template<typename complexT>
void
TPWALikelihood<complexT>::buildProdAmpArrays(const double*             inPar,
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


template<typename complexT>
void
TPWALikelihood<complexT>::clear()
{
	_decayAmps.resize(extents[0][0][0]);
}


// depends on naming convention for waves!!!
// VR_IGJPCMEIso....
template<typename complexT>
int
TPWALikelihood<complexT>::getReflectivity(const TString& waveName)
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
template<typename complexT>
void
TPWALikelihood<complexT>::copyFromParArray
(const double*  inPar,             // input parameter array
 ampsArrayType& outVal,            // array of complex output values [rank][reflectivity][wave index]
 value_type&    outFlatVal) const  // output value corresponding to flat wave
{
	// group parameters by rank and wave index only; keeps the order defined in wave list
	outVal.resize(extents[_rank][2][max(_nmbWavesRefl[0], _nmbWavesRefl[1])]);
	unsigned int parIndex = 0;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				value_type re, im;
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
				outVal[iRank][iRefl][iWave] = complexT(re, im);
			}
	outFlatVal = inPar[parIndex];
}


// copy values from structure that corresponds to complex
// production amplitudes to array that corresponds to function
// parameters taking into account rank restrictions
template<typename complexT>
void
TPWALikelihood<complexT>::copyToParArray
(const ampsArrayType& inVal,         // values corresponding to production amplitudes
 const value_type     inFlatVal,     // value corresponding to flat wave
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


template<typename complexT>
ostream&
TPWALikelihood<complexT>::print(ostream& out) const
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
			out << "        [" << setw(2) << sign((int)iRefl * 2 - 1) << " " << setw(3) << iWave << "] "
			    << _waveNames[iRefl][iWave] << "    threshold = "
			    << _waveThresholds[iRefl][iWave] << " MeV/c^2" << endl;
	out << "list of function parameters: " << endl;
	for (unsigned int iPar = 0; iPar < _nmbPars; ++iPar)
		out << "        [" << setw(3) << iPar << "] " << _parNames[iPar] << "    "
		    << "threshold = " << _parThresholds[iPar] << " MeV/c^2" << endl;
	return out;
}


template<typename complexT>
ostream&
TPWALikelihood<complexT>::printFuncInfo(ostream& out) const
{
	const string funcNames[NMB_FUNCTIONCALLENUM] = {"FdF", "Gradient", "DoEval", "DoDerivative"};
	for (unsigned int i = 0; i < NMB_FUNCTIONCALLENUM; ++i)
		if (_funcCallInfo[i].nmbCalls > 0)
			out << "    " << _funcCallInfo[i].nmbCalls
			    << " calls to TPWALikelihood<complexT>::" << funcNames[i] << "()" << endl
			    << "    time spent in TPWALikelihood<complexT>::" << funcNames[i] << "(): " << endl
			    << "        total time ................. " << sum(_funcCallInfo[i].totalTime) << " sec" << endl
			    << "        return value calculation ... " << sum(_funcCallInfo[i].funcTime ) << " sec" << endl
			    << "        normalization .............. " << sum(_funcCallInfo[i].normTime ) << " sec" << endl;
	return out;
}


template<typename complexT>
vector<unsigned int>
TPWALikelihood<complexT>::orderedParIndices() const
{
	vector<unsigned int> orderedIndices;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int waveIndex = 0; waveIndex < _nmbWaves; ++waveIndex) {
			unsigned int iRefl, iWave;
			for (iRefl = 0; iRefl < 2; ++iRefl)
				for (iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					if (_waveToWaveIndex[iRefl][iWave] == waveIndex)
						goto found;
			printWarn << "indices are inconsistent. cannot find wave with index " << waveIndex
			          << " in wave list" << endl;
			continue;
      
		found:
      
			tuple<int, int> parIndices = _prodAmpToFuncParMap[iRank][iRefl][iWave];
			int             parIndex   = get<0>(parIndices);
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


template<typename complexT>
vector<string>
TPWALikelihood<complexT>::waveNames() const
{
	vector<string> names(_nmbWaves, "");
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
			const unsigned int index = _waveToWaveIndex[iRefl][iWave];
			names[index] = _waveNames[iRefl][iWave];
		}
	return names;
}


template<typename complexT>
void
TPWALikelihood<complexT>::resetFuncCallInfo() const
{
	for (unsigned int i = 0; i < NMB_FUNCTIONCALLENUM; ++i) {
		_funcCallInfo[i].nmbCalls  = 0;
		// boost accumulators do not have a reset function
		_funcCallInfo[i].funcTime  = typename functionCallInfo::timeAccType();
		_funcCallInfo[i].normTime  = typename functionCallInfo::timeAccType();
		_funcCallInfo[i].totalTime = typename functionCallInfo::timeAccType();
	}
}


// explicit specializations
template class TPWALikelihood<complex<float > >;
template class TPWALikelihood<complex<double> >;
