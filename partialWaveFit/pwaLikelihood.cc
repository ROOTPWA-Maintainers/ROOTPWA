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
//
// Description:
//      Implementation of class pwaLikelihood
//      see pwaLikelihood.h for details
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//      Boris Grube          Universe Cluster Munich
//
//
//-----------------------------------------------------------


#include "pwaLikelihood.h"

#include <cassert>
#include <complex>
#include <fstream>
#include <iomanip>
#include <limits>

#include "TStopwatch.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include "amplitudeMetadata.h"
#include "amplitudeTreeLeaf.h"
#include "complexMatrix.h"
#include "conversionUtils.hpp"
#include "eventMetadata.h"
#include "fileUtils.hpp"
#include "reportingUtils.hpp"
#ifdef USE_CUDA
#include "complex.cuh"
#include "likelihoodInterface.cuh"
#endif


// #define USE_FDF


using namespace std;
using namespace rpwa;
using namespace boost;
using namespace boost::accumulators;
namespace bt = boost::tuples;


template<typename complexT> bool pwaLikelihood<complexT>::_debug = true;


namespace {

	double cauchyFunction(const double x, const double gamma)
	{
		if(x < 0.) {
			printWarn << "got negative argument." << endl;
		}
		return 1. / (1. + ((x*x) / (gamma*gamma)));
	}

	double cauchyFunctionDerivative(const double x, const double gamma)
	{
		if(x < 0.) {
			printWarn << "got negative argument." << endl;
		}
		return (-2. * x * gamma*gamma) / ((gamma*gamma + x*x) * (gamma*gamma + x*x));
	}

	double cauchyFunctionSecondDerivative(const double x, const double gamma)
	{
		if(x < 0.) {
			printWarn << "got negative argument." << endl;
		}
		const double x2 = x*x;
		const double gamma2 = gamma*gamma;
		return gamma2 / ((gamma2 + x2) * (gamma2 + x2)) * (-2. + 8. * x2 / (gamma2 + x2));
	}

}


template<typename complexT>
string
pwaLikelihood<complexT>::fitParameter::parName() const
{
	ostringstream parName;
	parName << "V";
	if (_waveName != "flat")
		parName << _rank;
	parName << "_" << _waveName;
	if (_realPart) {
		parName << "_RE";
	} else {
		parName << "_IM";
	}
	return parName.str();
}


template<typename complexT>
pwaLikelihood<complexT>::pwaLikelihood()
	: _nmbEvents        (0),
	  _rank             (1),
	  _nmbWaves         (0),
	  _nmbWavesReflMax  (0),
	  _nmbPars          (0),
	  _initialized      (false),
	  _normIntAdded     (false),
	  _accIntAdded      (false),
	  _initFinished     (false),
#ifdef USE_CUDA
	  _cudaEnabled      (false),
#endif
	  _useNormalizedAmps(true),
	  _priorType        (FLAT),
	  _cauchyWidth      (0.5),
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
pwaLikelihood<complexT>::~pwaLikelihood()
{
	clear();
}


template<typename complexT>
void
pwaLikelihood<complexT>::FdF
(const double* par,             // parameter array; reduced by rank conditions
 double&       funcVal,         // function value
 double*       gradient) const  // array of derivatives
{
	if (not _initFinished) {
		printErr << "pwaLikelihood::finishInit has not been called. Aborting..." << endl;
		throw;
	}
	++(_funcCallInfo[FDF].nmbCalls);

	// timer for total time
	TStopwatch timerTot;
	timerTot.Start();

	// copy arguments into parameter cache
	for (unsigned int i = 0; i < _nmbPars; ++i)
		_parCache[i] = par[i];

	// build complex production amplitudes from function parameters taking into account rank restrictions
	value_type        prodAmpFlat;
	prodAmpsArrayType prodAmps;
	copyFromParArray(par, prodAmps, prodAmpFlat);
	const value_type prodAmpFlat2 = prodAmpFlat * prodAmpFlat;

	// create array of likelihood derivatives w.r.t. real and imaginary
	// parts of the production amplitudes
	// !NOTE! although stored as and constructed from complex values,
	// the dL themselves are _not_ well defined complex numbers!
	value_type                                         derivativeFlat = 0;
	boost::array<typename prodAmpsArrayType::index, 3> derivShape     = {{ _rank, 2, _nmbWavesReflMax }};
	prodAmpsArrayType                                  derivatives(derivShape);

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
		prodAmpsArrayType derivative(derivShape);  // likelihood derivatives for this event
		for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
				accumulator_set<complexT, stats<tag::sum(compensated)> > ampProdAcc;
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
					ampProdAcc(prodAmps[iRank][iRefl][iWave] * _decayAmps[iRefl][iEvt][iWave]);
				}
				const complexT ampProdSum = sum(ampProdAcc);
				likelihoodAcc(norm(ampProdSum));
				// set derivative term that is independent on derivative wave index
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					// amplitude sums for current rank and for waves with same reflectivity
					derivative[iRank][iRefl][iWave] = ampProdSum;
			}
			// loop again over waves for current rank and multiply with complex conjugate
			// of decay amplitude of the wave with the derivative wave index
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					derivative[iRank][iRefl][iWave] *= conj(_decayAmps[iRefl][iEvt][iWave]);
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

	double priorValue = 0.;
	switch(_priorType)
	{
		case FLAT:
			break;
		case HALF_CAUCHY:
			for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
						const double r = abs(prodAmps[iRank][iRefl][iWave]);
						const double cauchyFunctionValue = cauchyFunction(r, _cauchyWidth);
						const double factor = (1./(r*cauchyFunctionValue)) * cauchyFunctionDerivative(r, _cauchyWidth);
						complexT derivative(factor*prodAmps[iRank][iRefl][iWave].real(), factor*prodAmps[iRank][iRefl][iWave].imag());
						priorValue -= log(cauchyFunctionValue);
						derivatives[iRank][iRefl][iWave] -= derivative;
					}
				}
			}
			break;
	}

	// sort derivative results into output array and cache
	copyToParArray(derivatives, derivativeFlat, gradient);
	copyToParArray(derivatives, derivativeFlat, _derivCache.data());

	// calculate log likelihood value
	funcVal = sum(logLikelihoodAcc) + nmbEvt * sum(normFactorAcc) + priorValue;

	// log total consumed time
	timerTot.Stop();
	_funcCallInfo[FDF].totalTime(timerTot.RealTime());

	if (_debug)
		printDebug << "raw log likelihood = "        << maxPrecisionAlign(sum(logLikelihoodAcc)) << ", "
		           << "normalization = "             << maxPrecisionAlign(sum(normFactorAcc)   ) << ", "
		           << "prior = "                     << maxPrecisionAlign(priorValue           ) << ", "
		           << "normalized log likelihood = " << maxPrecisionAlign(funcVal              ) << endl;
}


template<typename complexT>
double
pwaLikelihood<complexT>::DoEval(const double* par) const
{
	if (not _initFinished) {
		printErr << "pwaLikelihood::finishInit has not been called. Aborting..." << endl;
		throw;
	}
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
	value_type        prodAmpFlat;
	prodAmpsArrayType prodAmps;
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
#endif
	{
		accumulator_set<value_type, stats<tag::sum(compensated)> > logLikelihoodAcc;
		for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
			accumulator_set<value_type, stats<tag::sum(compensated)> > likelihoodAcc;
			for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
					accumulator_set<complexT, stats<tag::sum(compensated)> > ampProdAcc;
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
						ampProdAcc(prodAmps[iRank][iRefl][iWave] * _decayAmps[iRefl][iEvt][iWave]);
					}
					const complexT ampProdSum = sum(ampProdAcc);
					likelihoodAcc(norm(ampProdSum));
				}
			}  // end loop over rank
			likelihoodAcc   (prodAmpFlat2            );
			logLikelihoodAcc(-log(sum(likelihoodAcc)));
		}  // end loop over events
		logLikelihood = sum(logLikelihoodAcc);
	}
	// log time needed for likelihood calculation
	timer.Stop();
	_funcCallInfo[DOEVAL].funcTime(timer.RealTime());

	// compute normalization term of log likelihood
	timer.Start();
	accumulator_set<value_type, stats<tag::sum(compensated)> > normFactorAcc;
	const value_type nmbEvt = (_useNormalizedAmps) ? 1 : _nmbEvents;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave) {  // inner loop over waves with same reflectivity
					const complexT I = _accMatrix[iRefl][iWave][iRefl][jWave];
					normFactorAcc(real((prodAmps[iRank][iRefl][iWave] * conj(prodAmps[iRank][iRefl][jWave]))
					                   * I));
				}
			}
	// take care of flat wave
	normFactorAcc(prodAmpFlat2 * _totAcc);
	// log time needed for normalization
	timer.Stop();
	_funcCallInfo[DOEVAL].normTime(timer.RealTime());

	double priorValue = 0.;
	switch(_priorType)
	{
		case FLAT:
			break;
		case HALF_CAUCHY:
			for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
						const double r = abs(prodAmps[iRank][iRefl][iWave]);
						const double cauchyFunctionValue = cauchyFunction(r, _cauchyWidth);
						priorValue -= log(cauchyFunctionValue);
					}
				}
			}
			break;
	}

	// calculate log likelihood value
	const double funcVal = logLikelihood + nmbEvt * sum(normFactorAcc) + priorValue;

	// log total consumed time
	timerTot.Stop();
	_funcCallInfo[DOEVAL].totalTime(timerTot.RealTime());

	if (_debug)
		printDebug << "raw log likelihood = "        << maxPrecisionAlign(logLikelihood     ) << ", "
		           << "normalization = "             << maxPrecisionAlign(sum(normFactorAcc)) << ", "
		           << "prior = "                     << maxPrecisionAlign(priorValue        ) << ", "
		           << "normalized log likelihood = " << maxPrecisionAlign(funcVal           ) << endl;

	return funcVal;

#endif  // USE_FDF
}


template<typename complexT>
double
pwaLikelihood<complexT>::DoDerivative(const double*      par,
                                      const unsigned int derivativeIndex) const
{
	if (not _initFinished) {
		printErr << "pwaLikelihood::finishInit has not been called. Aborting..." << endl;
		throw;
	}
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
	timerTot.Stop();
	_funcCallInfo[DODERIVATIVE].totalTime(timerTot.RealTime());
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
template<typename complexT>
void
pwaLikelihood<complexT>::Gradient
(const double* par,             // parameter array; reduced by rank conditions
 double*       gradient) const  // array of derivatives
{
	if (not _initFinished) {
		printErr << "pwaLikelihood::finishInit has not been called. Aborting..." << endl;
		throw;
	}
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
	timerTot.Stop();
	_funcCallInfo[GRADIENT].totalTime(timerTot.RealTime());
	if (samePar) {
		for (unsigned int i = 0; i < _nmbPars ; ++i)
			gradient[i] = _derivCache[i];
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
	value_type        prodAmpFlat;
	prodAmpsArrayType prodAmps;
	copyFromParArray(par, prodAmps, prodAmpFlat);

	// create array of likelihood derivatives w.r.t. real and imaginary
	// parts of the production amplitudes
	// !NOTE! although stored as and constructed from complex values,
	// the dL themselves are _not_ well defined complex numbers!
	value_type                                         derivativeFlat = 0;
	boost::array<typename prodAmpsArrayType::index, 3> derivShape     = {{ _rank, 2, _nmbWavesReflMax }};
	prodAmpsArrayType                                  derivatives(derivShape);

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
#endif
	{
		accumulator_set<value_type, stats<tag::sum(compensated)> > derivativeFlatAcc;
		multi_array<accumulator_set<complexT, stats<tag::sum(compensated)> >, 3>
			derivativesAcc(derivShape);
		const value_type prodAmpFlat2 = prodAmpFlat * prodAmpFlat;
		for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
			accumulator_set<value_type, stats<tag::sum(compensated)> > likelihoodAcc;
			prodAmpsArrayType derivative(derivShape);  // likelihood derivatives for this event
			for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
					accumulator_set<complexT, stats<tag::sum(compensated)> > ampProdAcc;
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
						ampProdAcc(prodAmps[iRank][iRefl][iWave] * _decayAmps[iRefl][iEvt][iWave]);
					}
					const complexT ampProdSum = sum(ampProdAcc);
					likelihoodAcc(norm(ampProdSum));
					// set derivative term that is independent on derivative wave index
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
						// amplitude sums for current rank and for waves with same reflectivity
						derivative[iRank][iRefl][iWave] = ampProdSum;
				}
				// loop again over waves for current rank and multiply with complex conjugate
				// of decay amplitude of the wave with the derivative wave index
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
						derivative[iRank][iRefl][iWave] *= conj(_decayAmps[iRefl][iEvt][iWave]);
			}  // end loop over rank
			likelihoodAcc(prodAmpFlat2);
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
	}
	// log time needed for likelihood calculation
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

	switch(_priorType)
	{
		case FLAT:
			break;
		case HALF_CAUCHY:
			for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
						const double r = abs(prodAmps[iRank][iRefl][iWave]);
						const double cauchyFunctionValue = cauchyFunction(r, _cauchyWidth);
						const double factor = (1./(r*cauchyFunctionValue)) * cauchyFunctionDerivative(r, _cauchyWidth);
						complexT derivative(factor*prodAmps[iRank][iRefl][iWave].real(), factor*prodAmps[iRank][iRefl][iWave].imag());
						derivatives[iRank][iRefl][iWave] -= derivative;
					}
				}
			}
			break;
	}

	// sort derivative results into output array and cache
	copyToParArray(derivatives, derivativeFlat, gradient);
	copyToParArray(derivatives, derivativeFlat, _derivCache.data());

	// log total consumed time
	timerTot.Stop();
	_funcCallInfo[GRADIENT].totalTime(timerTot.RealTime());

#endif  // USE_FDF
}


// calculate Hessian with respect to parameters
template<typename complexT>
TMatrixT<double>
pwaLikelihood<complexT>::Hessian
(const double* par) const  // parameter array; reduced by rank conditions
{
	if (not _initFinished) {
		printErr << "pwaLikelihood::finishInit has not been called. Aborting..." << endl;
		throw;
	}
	++(_funcCallInfo[HESSIAN].nmbCalls);

	// timer for total time
	TStopwatch timerTot;
	timerTot.Start();

	// build complex production amplitudes from function parameters taking into account rank restrictions
	value_type        prodAmpFlat;
	prodAmpsArrayType prodAmps;
	copyFromParArray(par, prodAmps, prodAmpFlat);
	const value_type prodAmpFlat2 = prodAmpFlat * prodAmpFlat;

	// create array to store likelihood hessian w.r.t. real and imaginary
	// parts of the production amplitudes
	value_type                                         hessianFlat  = 0;       // term in Hessian matrix where we derive twice w.r.t. the flat wave
	boost::array<typename prodAmpsArrayType::index, 3> derivShape   = {{ _rank, 2, _nmbWavesReflMax }};
	prodAmpsArrayType                                  flatTerms(derivShape);  // array for terms where we first derive w.r.t to
	                                                                           // a non-flat term and then w.r.t. the flat wave
	boost::array<typename prodAmpsArrayType::index, 7> hessianShape = {{ _rank, 2, _nmbWavesReflMax, _rank, 2, _nmbWavesReflMax, 3 }};
	boost::multi_array<value_type, 7>                  hessian(hessianShape);  // array to store components for Hessian matrix

	// loop over events and calculate second derivatives with respect to
	// parameters for the raw likelihood part
	TStopwatch timer;
	timer.Start();
	accumulator_set<value_type, stats<tag::sum(compensated)> > hessianFlatAcc;
	multi_array<accumulator_set<complexT, stats<tag::sum(compensated)> >, 3>
		flatTermsAcc(derivShape);
	multi_array<accumulator_set<value_type, stats<tag::sum(compensated)> >, 7>
		hessianAcc(hessianShape);
	for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
		accumulator_set<value_type, stats<tag::sum(compensated)> > likelihoodAcc;
		prodAmpsArrayType derivative(derivShape);  // likelihood derivatives for this event
		for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
				accumulator_set<complexT, stats<tag::sum(compensated)> > ampProdAcc;
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
					ampProdAcc(prodAmps[iRank][iRefl][iWave] * _decayAmps[iRefl][iEvt][iWave]);
				}
				const complexT ampProdSum = sum(ampProdAcc);
				likelihoodAcc(norm(ampProdSum));
				// set derivative term that is independent on derivative wave index
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					// amplitude sums for current rank and for waves with same reflectivity
					derivative[iRank][iRefl][iWave] = ampProdSum;
			}
			// loop again over waves for current rank and multiply with complex conjugate
			// of decay amplitude of the wave with the derivative wave index
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					derivative[iRank][iRefl][iWave] *= conj(_decayAmps[iRefl][iEvt][iWave]);
		}  // end loop over rank
		likelihoodAcc(prodAmpFlat2);
		// incorporate factor 2 / sigma
		const value_type factor  = 2. / sum(likelihoodAcc);
		const value_type factor2 = factor*factor;
		for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
					for (unsigned int jRank = 0; jRank < _rank; ++jRank) {
						for (unsigned int jRefl = 0; jRefl < 2; ++jRefl) {
							for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
								// last array index 0 indicates derivative w.r.t. real part of first prodAmp and real part of the second prodAmp
								hessianAcc[iRank][iRefl][iWave][jRank][jRefl][jWave][0](factor2 * derivative[jRank][jRefl][jWave].real() * derivative[iRank][iRefl][iWave].real());
								// last array index 1 indicates derivative w.r.t. real part of first prodAmp and imaginary part of the second prodAmp
								hessianAcc[iRank][iRefl][iWave][jRank][jRefl][jWave][1](factor2 * derivative[jRank][jRefl][jWave].imag() * derivative[iRank][iRefl][iWave].real());
								// last array index 2 indicates derivative w.r.t. imaginary part of first prodAmp and imaginary part of the second prodAmp
								hessianAcc[iRank][iRefl][iWave][jRank][jRefl][jWave][2](factor2 * derivative[jRank][jRefl][jWave].imag() * derivative[iRank][iRefl][iWave].imag());
								if(iRank == jRank and iRefl == jRefl) {
									const complexT uPrime = conj(_decayAmps[jRefl][iEvt][jWave]) * _decayAmps[iRefl][iEvt][iWave];
									hessianAcc[iRank][iRefl][iWave][jRank][jRefl][jWave][0](-factor * uPrime.real());
									hessianAcc[iRank][iRefl][iWave][jRank][jRefl][jWave][1](-factor * uPrime.imag());
									hessianAcc[iRank][iRefl][iWave][jRank][jRefl][jWave][2](-factor * uPrime.real());
								}
							}
						}
					}
					flatTermsAcc[iRank][iRefl][iWave](factor2 * prodAmpFlat * derivative[iRank][iRefl][iWave]);  // calculate terms where we first derive w.r.t. real/imag part
					                                                                                             // of a prodAmp and then w.r.t. the flat wave
				}
			}
		}
		hessianFlatAcc(factor2 * prodAmpFlat2 - factor);
	}  // end loop over events
	for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				for (unsigned int jRank = 0; jRank < _rank; ++jRank) {
					for (unsigned int jRefl = 0; jRefl < 2; ++jRefl) {
						for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
							hessian[iRank][iRefl][iWave][jRank][jRefl][jWave][0] = sum(hessianAcc[iRank][iRefl][iWave][jRank][jRefl][jWave][0]);
							hessian[iRank][iRefl][iWave][jRank][jRefl][jWave][1] = sum(hessianAcc[iRank][iRefl][iWave][jRank][jRefl][jWave][1]);
							hessian[iRank][iRefl][iWave][jRank][jRefl][jWave][2] = sum(hessianAcc[iRank][iRefl][iWave][jRank][jRefl][jWave][2]);
						}
					}
				}
				flatTerms[iRank][iRefl][iWave] = sum(flatTermsAcc[iRank][iRefl][iWave]);
			}
		}
	}
	hessianFlat = sum(hessianFlatAcc);
	// log time needed for calculation of second derivatives of raw likelhood part
	timer.Stop();
	_funcCallInfo[HESSIAN].funcTime(timer.RealTime());

	// normalize second derivatives w.r.t. parameters
	timer.Start();
	const value_type nmbEvt      = (_useNormalizedAmps) ? 1 : _nmbEvents;
	const value_type twiceNmbEvt = 2 * nmbEvt;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave) {
					const complexT I = _accMatrix[iRefl][iWave][iRefl][jWave];
					hessian[iRank][iRefl][iWave][iRank][iRefl][jWave][0] += I.real() * twiceNmbEvt;
					hessian[iRank][iRefl][iWave][iRank][iRefl][jWave][1] += I.imag() * twiceNmbEvt;
					hessian[iRank][iRefl][iWave][iRank][iRefl][jWave][2] += I.real() * twiceNmbEvt;
				}
	hessianFlat += twiceNmbEvt * _totAcc;
	// log time needed for normalization
	timer.Stop();
	_funcCallInfo[HESSIAN].normTime(timer.RealTime());

	switch(_priorType)
	{
		case FLAT:
			break;
		case HALF_CAUCHY:
			for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
						const double r = abs(prodAmps[iRank][iRefl][iWave]);

						const double cauchyFunctionValue                 = cauchyFunction                (r, _cauchyWidth);
						const double cauchyFunctionDerivativeValue       = cauchyFunctionDerivative      (r, _cauchyWidth);
						const double cauchyFunctionSecondDerivativeValue = cauchyFunctionSecondDerivative(r, _cauchyWidth);

						const double a1 = cauchyFunctionSecondDerivativeValue / (cauchyFunctionValue * r*r);
						const double a2 = cauchyFunctionDerivativeValue*cauchyFunctionDerivativeValue / (cauchyFunctionValue*cauchyFunctionValue * r*r);
						const double a3 = cauchyFunctionDerivativeValue / (cauchyFunctionValue * r);
						const double a4 = cauchyFunctionDerivativeValue / (cauchyFunctionValue * r*r*r);
						const double a124 = a1 - a2 - a4;

						const double re = prodAmps[iRank][iRefl][iWave].real();
						const double im = prodAmps[iRank][iRefl][iWave].imag();

						hessian[iRank][iRefl][iWave][iRank][iRefl][iWave][0] -= (a124) * re*re + a3;
						hessian[iRank][iRefl][iWave][iRank][iRefl][iWave][1] -= (a124) * re*im;
						hessian[iRank][iRefl][iWave][iRank][iRefl][iWave][2] -= (a124) * im*im + a3;
					}
				}
			}
			break;
	}

	TMatrixT<double> hessianMatrix(_nmbPars, _nmbPars);
	for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				bt::tuple<int, int> parIndices1 = _prodAmpToFuncParMap[iRank][iRefl][iWave];
				const int r1 = get<0>(parIndices1);
				const int i1 = get<1>(parIndices1);

				if (r1 < 0)
					continue;
				if (_parameters[r1].fixed())
					continue;
				assert(i1 < 0 || not _parameters[i1].fixed());

				for (unsigned int jRank = 0; jRank < _rank; ++jRank) {
					for (unsigned int jRefl = 0; jRefl < 2; ++jRefl) {
						for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
							bt::tuple<int, int> parIndices2 = _prodAmpToFuncParMap[jRank][jRefl][jWave];
							const int r2 = get<0>(parIndices2);
							const int i2 = get<1>(parIndices2);

							if (r2 < 0)
								continue;
							if (_parameters[r2].fixed())
								continue;
							assert(i2 < 0 || not _parameters[i2].fixed());

							hessianMatrix[r1][r2] = hessian[iRank][iRefl][iWave][jRank][jRefl][jWave][0];
							if (i2 >= 0) // real/imaginary derivative
								hessianMatrix[r1][i2] = hessian[iRank][iRefl][iWave][jRank][jRefl][jWave][1];
							if (i1 >= 0) // imaginary/real derivative
								hessianMatrix[i1][r2] = hessian[jRank][jRefl][jWave][iRank][iRefl][iWave][1];
							if (i1 >= 0 && i2 >= 0) // imaginary/imaginary derivative
								hessianMatrix[i1][i2] = hessian[iRank][iRefl][iWave][jRank][jRefl][jWave][2];
						}
					}
				}

				// second derivatives w.r.t. flat wave
				hessianMatrix[r1][_nmbPars - 1] = flatTerms[iRank][iRefl][iWave].real();
				hessianMatrix[_nmbPars - 1][r1] = flatTerms[iRank][iRefl][iWave].real();
				if (i1 >= 0) {
					hessianMatrix[i1][_nmbPars - 1] = flatTerms[iRank][iRefl][iWave].imag();
					hessianMatrix[_nmbPars - 1][i1] = flatTerms[iRank][iRefl][iWave].imag();
				}
			}
		}
	}
	hessianMatrix[_nmbPars - 1][_nmbPars - 1] = hessianFlat;  // enter flat/flat term

	// log total consumed time
	timerTot.Stop();
	_funcCallInfo[HESSIAN].totalTime(timerTot.RealTime());

	return hessianMatrix;
}


/// calculates eigenvectors/-values of Hessian
template<typename complexT>
vector<pair<TVectorT<double>, double> >
pwaLikelihood<complexT>::HessianEigenVectors(const TMatrixT<double>& hessian) const
{
	// reduce Hessian matrix by removing the rows and columns corresponding
	// to fixed parameters
	TMatrixT<double> hessianRed(_nmbPars - _nmbParsFixed, _nmbPars - _nmbParsFixed);
	{
		unsigned int iSkip = 0;
		for (unsigned int i = 0; i < _nmbPars; ++i) {
			if (_parameters[i].fixed()) {
				iSkip++;
				continue;
			}

			unsigned int jSkip = 0;
			for (unsigned int j = 0; j < _nmbPars; ++j) {
				if (_parameters[j].fixed()) {
					jSkip++;
					continue;
				}

				hessianRed[i - iSkip][j - jSkip] = hessian[i][j];
			}
		}
	}

	TVectorT<double> eigenValuesRed;
	TMatrixT<double> eigenVectorsRed = hessianRed.EigenVectors(eigenValuesRed);

	vector<pair<TVectorT<double>, double> > eigen;
	{
		for (unsigned int i = 0; i < _nmbPars - _nmbParsFixed; ++i) {
			const TVectorT<double> eigenVectorRed(TMatrixTColumn_const<double>(eigenVectorsRed, i));
			TVectorT<double> eigenVector(_nmbPars);

			unsigned int jSkip = 0;
			for (unsigned int j = 0; j < _nmbPars; ++j) {
				if (_parameters[j].fixed()) {
					jSkip++;
					eigenVector[j] = 0.;
				} else {
					eigenVector[j] = eigenVectorRed[j - jSkip];
				}
			}

			eigen.push_back(make_pair(eigenVector, eigenValuesRed[i]));
		}
	}
	return eigen;
}


/// calculates covariance matrix of function at point defined by par
template<typename complexT>
TMatrixT<double>
pwaLikelihood<complexT>::CovarianceMatrix(const double* par) const
{
	const TMatrixT<double> hessian = Hessian(par);
	return CovarianceMatrix(hessian);
}


/// turns hessian into covariance matrix
template<typename complexT>
TMatrixT<double>
pwaLikelihood<complexT>::CovarianceMatrix(const TMatrixT<double>& hessian) const
{
	// reduce Hessian matrix by removing the rows and columns corresponding
	// to fixed parameters
	TMatrixT<double> hessianRed(_nmbPars - _nmbParsFixed, _nmbPars - _nmbParsFixed);
	{
		unsigned int iSkip = 0;
		for (unsigned int i = 0; i < _nmbPars; ++i) {
			if (_parameters[i].fixed()) {
				iSkip++;
				continue;
			}

			unsigned int jSkip = 0;
			for (unsigned int j = 0; j < _nmbPars; ++j) {
				if (_parameters[j].fixed()) {
					jSkip++;
					continue;
				}

				hessianRed[i - iSkip][j - jSkip] = hessian[i][j];
			}
		}
	}

	// invert Hessian to get covariance matrix
	TMatrixT<double> covarianceRed(TMatrixT<double>::kInverted, hessianRed);

	// blow up covariance matrix by adding the rows and columns corresponding
	// to fixed parameters
	TMatrixT<double> covariance(_nmbPars, _nmbPars);
	{
		unsigned int iSkip = 0;
		for (unsigned int i = 0; i < _nmbPars; ++i) {
			if (_parameters[i].fixed()) {
				iSkip++;
				continue;
			}

			unsigned int jSkip = 0;
			for (unsigned int j = 0; j < _nmbPars; ++j) {
				if (_parameters[j].fixed()) {
					jSkip++;
					continue;
				}

				covariance[i][j] = covarianceRed[i - iSkip][j - jSkip];
			}
		}
	}

	return covariance;
}


template<typename complexT>
vector<double>
pwaLikelihood<complexT>::CorrectParamSigns(const double* par) const
{
	// build complex production amplitudes from function parameters taking into account rank restrictions
	value_type        prodAmpFlat;
	prodAmpsArrayType prodAmps;
	copyFromParArray(par, prodAmps, prodAmpFlat);

	if (prodAmpFlat < 0) {
		prodAmpFlat *= -1;  // flip sign of flat wave if it is negative
	}

	for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
			if (prodAmps[iRank][iRefl][iRank].real() < 0 and prodAmps[iRank][iRefl][iRank].imag() == 0) {
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
					prodAmps[iRank][iRefl][iWave] *= -1;  // flip sign of coherent waves if anchorwave is negative
				}
			}
		}
	}

	vector<double> returnParams(_nmbPars);
	copyToParArray(prodAmps, prodAmpFlat, returnParams.data());
	return returnParams;
}


template<typename complexT>
unsigned int
pwaLikelihood<complexT>::nmbWaves(const int reflectivity) const
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
pwaLikelihood<complexT>::enableCuda(const bool enableCuda)
{
	_cudaEnabled = enableCuda;
}
#else
pwaLikelihood<complexT>::enableCuda(const bool) { }
#endif


template<typename complexT>
bool
pwaLikelihood<complexT>::cudaEnabled() const
{
#ifdef USE_CUDA
	return _cudaEnabled;
#else
	return false;
#endif
}


template<typename complexT>
bool
pwaLikelihood<complexT>::init(const vector<waveDescThresType>& waveDescThresType,
                              const unsigned int               rank,
                              const double                     massBinCenter)
{
	if (not readWaveList(waveDescThresType))
		return false;
	if (not buildParDataStruct(rank, massBinCenter))
		return false;

	_initialized = true;
	return true;
}


template<typename complexT>
bool
pwaLikelihood<complexT>::addNormIntegral(const ampIntegralMatrix& normMatrix)
{
	if (not _initialized) {
		printErr << "likelihood not initialized!" << endl
		         << "call pwaLikelihood::init() before pwaLikelihood::addNormIntegral()" << endl
		         << "Aborting..." << endl;
		return false;
	}

	_numbAccEvents = normMatrix.nmbEvents();

	reorderIntegralMatrix(normMatrix, _normMatrix);

	_normIntAdded = true;
	return true;
}


template<typename complexT>
bool
pwaLikelihood<complexT>::addAccIntegral(ampIntegralMatrix& accMatrix,
                                        const unsigned int accEventsOverride)
{
	if (not _normIntAdded) {
		printErr << "normalization integral not added before acceptance integral!" << endl
		         << "call pwaLikelihood::addNormIntegral() before pwaLikelihood::addAccIntegral()" << endl
		         << "Aborting..." << endl;
		return false;
	}

	_numbAccEvents = accEventsOverride==0 ? _numbAccEvents : accEventsOverride;
	_totAcc = (double) accMatrix.nmbEvents() / _numbAccEvents;
	accMatrix.setNmbEvents(_numbAccEvents);
	printInfo << "total acceptance in this bin: " << _totAcc << endl;

	reorderIntegralMatrix(accMatrix, _accMatrix);

	_accIntAdded = true;
	return true;
}


template<typename complexT>
bool
pwaLikelihood<complexT>::addAmplitude(const vector<const amplitudeMetadata*>& ampMetas)
{
	if (not _accIntAdded) {
		printErr << "no acceptance integral found. "
		         << "Call pwaLikelihood::addAccIntegral() before pwaLikelihood::addAmplitude(). "
		         << "Aborting..." << endl;
		return false;
	}
	if (ampMetas.size() == 0) {
		printErr << "no amplitudeMetadata given. Aborting..." << endl;
		return false;
	}
	for (size_t iAmpMeta = 0; iAmpMeta < ampMetas.size(); ++iAmpMeta) {
		if (not ampMetas[iAmpMeta]) {
			printErr << "amplitudeMetadata metas[" << iAmpMeta << "] not valid. Aborting..." << endl;
			return false;
		}
	}
	const string& waveName = ampMetas[0]->objectBaseName();
	for (size_t iAmpMeta = 0; iAmpMeta < ampMetas.size(); ++iAmpMeta) {
		if (ampMetas[iAmpMeta]->objectBaseName() != waveName) {
			printErr << "wave names in different amplitude metadatas do not match. Aborting..." << endl;
			return false;
		}
	}
	if (_waveParams.count(waveName) == 0) {
		printErr << "requested to add decay amplitudes for wave '" << waveName << "' which is not in wavelist. Aborting..." << endl;
		return false;
	}

	// counting all events
	const bool onTheFlyBinning = not _eventFileProperties.empty();
	size_t totalEvents = 0;
	for (size_t iAmpMeta = 0; iAmpMeta < ampMetas.size(); ++iAmpMeta) {
		if (onTheFlyBinning) {
			const vector<eventMetadata>& evtMetas = ampMetas[iAmpMeta]->eventMetadata();
			for (size_t iEvtMeta = 0; iEvtMeta < evtMetas.size(); ++iEvtMeta) {
				totalEvents += _eventFileProperties[evtMetas[iEvtMeta].contentHash()].second.size();
			}
		} else {
			totalEvents += ampMetas[iAmpMeta]->amplitudeTree()->GetEntriesFast();
		}
	}
	if (totalEvents == 0) {
		printErr << "no events to load. Aborting..." << endl;
		return false;
	}

	{
		// check ordering of event files
		const bool eventFileHashOrderWasEmpty = _eventFileHashOrder.empty();
		size_t runningEventFileIndex = 0;
		for(size_t iAmpMeta = 0; iAmpMeta < ampMetas.size(); ++iAmpMeta) {
			const amplitudeMetadata* ampMeta = ampMetas[iAmpMeta];
			for(size_t iEvtMeta = 0; iEvtMeta < ampMeta->eventMetadata().size(); ++iEvtMeta) {
				const string& eventFileHash = ampMeta->eventMetadata()[iEvtMeta].contentHash();
				if(eventFileHashOrderWasEmpty) {
					_eventFileHashOrder.push_back(eventFileHash);
				} else {
					if (runningEventFileIndex >= _eventFileHashOrder.size() or _eventFileHashOrder[runningEventFileIndex] != eventFileHash) {
						printErr << "order of event files differs between corresponding amplitude files. Aborting..." << endl;
						return false;
					}
					++runningEventFileIndex;
				}
				if(onTheFlyBinning) {
					if(_eventFileProperties.count(eventFileHash) == 0) {
						printErr << "event file hash from amplitude metadata not in "
						         << "on-the-fly binning information. Aborting..." << endl;
						return false;
					}
				}
			}
		}
	}

	vector<complexT> amps(totalEvents);
	size_t eventCount = 0; // Running count for event number over all single files
	for (size_t iAmpMeta = 0; iAmpMeta < ampMetas.size(); ++iAmpMeta) {
		const amplitudeMetadata* ampMeta = ampMetas[iAmpMeta];
		// connect tree leaf
		amplitudeTreeLeaf* ampTreeLeaf = 0;
		ampMeta->amplitudeTree()->SetBranchAddress(amplitudeMetadata::amplitudeLeafName.c_str(), &ampTreeLeaf);
		if (not ampTreeLeaf) {
			printWarn << "null pointer to amplitude leaf. Aborting..." << endl;
			return false;
		}

		if (onTheFlyBinning) {
			size_t skipEvents  = 0;
			for(size_t iEvtMeta = 0; iEvtMeta < ampMeta->eventMetadata().size(); ++iEvtMeta) {
				const string& eventFileHash = ampMeta->eventMetadata()[iEvtMeta].contentHash();
				const vector<size_t>& entriesInBin = _eventFileProperties[eventFileHash].second;
				for(size_t iEvent = 0; iEvent < entriesInBin.size(); ++iEvent, ++eventCount) {
					ampMeta->amplitudeTree()->GetEntry(skipEvents + entriesInBin[iEvent]);
					assert(ampTreeLeaf->nmbIncohSubAmps() == 1);
					complexT amp(ampTreeLeaf->incohSubAmp(0).real(), ampTreeLeaf->incohSubAmp(0).imag());
					amps[eventCount] = amp;
				}
				skipEvents += _eventFileProperties[eventFileHash].first;
			}
		} else {
			for(long iEvent = 0; iEvent < ampMeta->amplitudeTree()->GetEntriesFast(); ++iEvent, ++eventCount) {
				ampMeta->amplitudeTree()->GetEntry(iEvent);
				assert(ampTreeLeaf->nmbIncohSubAmps() == 1);
				complexT amp(ampTreeLeaf->incohSubAmp(0).real(), ampTreeLeaf->incohSubAmp(0).imag());
				amps[eventCount] = amp;
			}
		}
	}

	if (_nmbEvents == 0) {
		// first amplitude file read
		_nmbEvents = totalEvents;
		_decayAmps[0].resize(extents[_nmbEvents][_nmbWavesRefl[0]]);
		_decayAmps[1].resize(extents[_nmbEvents][_nmbWavesRefl[1]]);
	}
	if (totalEvents != _nmbEvents) {
		printWarn << "size mismatch in amplitude files: this file contains " << totalEvents
		          << " events, previous file had " << _nmbEvents << " events." << endl;
		return false;
	}

	const unsigned int refl = _waveParams[waveName].first;
	const unsigned int waveIndex = _waveParams[waveName].second;

	// get normalization
	const complexT normInt = _normMatrix[refl][waveIndex][refl][waveIndex];

	// copy decay amplitudes into array that is indexed [event index][reflectivity][wave index]
	// this index scheme ensures a more linear memory access pattern in the likelihood function
	for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
		if (_useNormalizedAmps) {  // normalize data, if option is switched on
			if (normInt == (value_type)0. && amps[iEvt] != (value_type)0.) {
				printErr << "normalization integral for wave '" << waveName << "' is zero, but the amplitude is not. Aborting...";
				return false;
			}
			if (normInt != (value_type)0.)
				amps[iEvt] /= sqrt(normInt.real());  // rescale decay amplitude
		}
		_decayAmps[refl][iEvt][waveIndex] = amps[iEvt];
	}

	_waveAmpAdded[refl][waveIndex] = true; // note that this amplitude has been added to the likelihood

	printInfo << "read decay amplitudes of " << totalEvents << " events for wave '" << waveName << "' into memory" << endl;
	return true;
}


template<typename complexT>
bool
pwaLikelihood<complexT>::finishInit()
{
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
			if (not _waveAmpAdded[iRefl][iWave]) {
				printErr << "not all expected amplitudes were added to the likelihood ('" << _waveNames[iRefl][iWave] << "'). Aborting..." << endl;
				return false;
			}
		}
	}

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
				value_type norm_i = sqrt(_normMatrix[iRefl][iWave][iRefl][iWave].real());
				if(norm_i==0)norm_i=1;
				for (unsigned int jRefl = 0; jRefl < 2; ++jRefl)
					for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
						value_type norm_j = sqrt(_normMatrix[jRefl][jWave][jRefl][jWave].real());
						// protect against empty amplitudes (which will not contribute anyway)
						if(norm_j==0)norm_j=1;
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
			printDebug << "normalized integral matrices" << endl;
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

#ifdef USE_CUDA
	if (_cudaEnabled)
		cuda::likelihoodInterface<cuda::complex<value_type> >::init
			(reinterpret_cast<cuda::complex<value_type>*>(_decayAmps.data()),
			 _decayAmps.num_elements(), _nmbEvents, _nmbWavesRefl, true);
#endif

	_initFinished = true;
	return true;
}


template<typename complexT>
bool
pwaLikelihood<complexT>::setOnTheFlyBinning(const map<string, pair<double, double> >& binningMap,
                                            const vector<const eventMetadata*>&       evtMetas)
{
	for (size_t metaIndex = 0; metaIndex < evtMetas.size(); ++metaIndex) {
		const eventMetadata* evtMeta = evtMetas[metaIndex];
		if (not evtMeta) {
			printErr << "eventMetadata not set. Aborting..." << endl;
			return false;
		}
		TTree* evtTree = evtMeta->eventTree();
		if(not evtTree) {
			printErr << "event tree invalid. Aborting..." << endl;
			return false;
		}
		const string& evtHash = evtMeta->contentHash();
		map<string, double> binningVariables;
		typedef map<string, pair<double, double> >::const_iterator it_type;
		for(it_type iterator = binningMap.begin(); iterator != binningMap.end(); ++iterator) {
			const string& additionalVar = iterator->first;
			binningVariables[additionalVar] = 0.;
			evtTree->SetBranchAddress(additionalVar.c_str(), &binningVariables[additionalVar]);
		}
		vector<size_t> eventIndices;
		for (long eventIndex = 0; eventIndex < evtTree->GetEntriesFast(); ++eventIndex) {
			evtTree->GetEntry(eventIndex);
			bool useEvent = true;
			for(it_type iterator = binningMap.begin(); iterator != binningMap.end(); ++iterator) {
				const string& additionalVar = iterator->first;
				const double& lowerBound = iterator->second.first;
				const double& upperBound = iterator->second.second;
				const double& binVarValue = binningVariables[additionalVar];
				if (binVarValue < lowerBound or binVarValue >= upperBound) {
					useEvent = false;
					break;
				}
			}
			if (useEvent) {
				eventIndices.push_back(eventIndex);
			}
		}
		printInfo << "found " << eventIndices.size() << " of " << evtMeta->eventTree()->GetEntriesFast() << " in the event file to be in given bin." << endl;
		_eventFileProperties[evtHash] = pair<size_t, vector<size_t> >(evtMeta->eventTree()->GetEntriesFast(), eventIndices);
	}
	return true;
}


template<typename complexT>
bool
pwaLikelihood<complexT>::readWaveList(const vector<waveDescThresType>& waveDescThres)
{
	vector<string> waveNames[2];
	vector<double> waveThres[2];
	for (unsigned int i = 0; i < waveDescThres.size(); ++i) {
		const waveDescription& waveDesc = boost::get<1>(waveDescThres[i]);
		isobarDecayTopologyPtr decayTopo;
		waveDesc.constructDecayTopology(decayTopo);
		isobarDecayVertexPtr   decayVertex = decayTopo->XIsobarDecayVertex();
		particlePtr            X           = decayVertex->parent();

		const int    refl      = X->reflectivity();
		const string waveName  = boost::get<0>(waveDescThres[i]);
		const double threshold = boost::get<2>(waveDescThres[i]);

		assert(refl == -1 || refl == +1);
		if (refl > 0) {
			++_nmbWavesRefl[1];  // positive reflectivity
			waveNames      [1].push_back(waveName);
			waveThres      [1].push_back(threshold);
		} else {
			++_nmbWavesRefl[0];  // negative reflectivity
			waveNames      [0].push_back(waveName);
			waveThres      [0].push_back(threshold);
		}
	}
	_nmbWaves        = _nmbWavesRefl[0] + _nmbWavesRefl[1];
	_nmbWavesReflMax = max(_nmbWavesRefl[0], _nmbWavesRefl[1]);
	_waveNames.resize     (extents[2][_nmbWavesReflMax]);
	_waveThresholds.resize(extents[2][_nmbWavesReflMax]);
	_waveAmpAdded.resize  (extents[2][_nmbWavesReflMax]);
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
			_waveNames     [iRefl][iWave] = waveNames[iRefl][iWave];
			_waveThresholds[iRefl][iWave] = waveThres[iRefl][iWave];
			_waveAmpAdded  [iRefl][iWave] = false;
			_waveParams.insert(pair<string, pair<unsigned int, unsigned int> >(waveNames[iRefl][iWave], pair<unsigned int, unsigned int>(iRefl, iWave)));
		}
	return true;
}


template<typename complexT>
bool
pwaLikelihood<complexT>::buildParDataStruct(const unsigned int rank,
                                            const double       massBinCenter)
{
	if ((_nmbWavesRefl[0] + _nmbWavesRefl[1] == 0) or (_waveThresholds.size() == 0)) {
		printErr << "no wave info. was readWaveList() executed successfully? Aborting..." << endl;
		return false;
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
	_nmbParsFixed = 0;
	_parameters.resize(_nmbPars);
	_parCache.resize  (_nmbPars, 0);
	_derivCache.resize(_nmbPars, 0);
	_prodAmpToFuncParMap.resize(extents[_rank][2][_nmbWavesReflMax]);
	// build parameter names
	unsigned int parIndex = 0;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				bool fixed = false;
				if (_waveThresholds[iRefl][iWave] != 0 && _waveThresholds[iRefl][iWave] >= massBinCenter) {
					fixed = true;
				}
				if (iWave < iRank)  // production amplitude is zero
					_prodAmpToFuncParMap[iRank][iRefl][iWave] = bt::make_tuple(-1, -1);
				else if (iWave == iRank) {  // production amplitude is real
					_parameters[parIndex] = fitParameter(_waveNames[iRefl][iWave], iRank, _waveThresholds[iRefl][iWave],
					                                     fixed, true);
					if (fixed) {
						++_nmbParsFixed;
					}
					_prodAmpToFuncParMap[iRank][iRefl][iWave] = bt::make_tuple(parIndex, -1);
					++parIndex;
				} else {  // production amplitude is complex
					_parameters[parIndex  ] = fitParameter(_waveNames[iRefl][iWave], iRank, _waveThresholds[iRefl][iWave],
					                                       fixed, true);
					_parameters[parIndex+1] = fitParameter(_waveNames[iRefl][iWave], iRank, _waveThresholds[iRefl][iWave],
					                                       fixed, false);
					if (fixed) {
						_nmbParsFixed += 2;
					}
					_prodAmpToFuncParMap[iRank][iRefl][iWave] = bt::make_tuple(parIndex, parIndex + 1);
					parIndex += 2;
				}
			}
	// flat wave
	_parameters[parIndex] = fitParameter("flat", 0, 0., false, true);
	return true;
}


// returns integral matrix reordered according to _waveNames array
template<typename complexT>
void
pwaLikelihood<complexT>::reorderIntegralMatrix(const ampIntegralMatrix& integral,
                                               normMatrixArrayType&     reorderedMatrix) const
{
	// create reordered matrix
	reorderedMatrix.resize(extents[2][_nmbWavesReflMax][2][_nmbWavesReflMax]);
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
			for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave) {
				const complex<double> val = integral.element(_waveNames[iRefl][iWave],
				                                             _waveNames[iRefl][jWave]);
				reorderedMatrix[iRefl][iWave][iRefl][jWave] = complexT(val.real(), val.imag());
			}
}


template<typename complexT>
void
pwaLikelihood<complexT>::getIntegralMatrices(complexMatrix&  normMatrix,
                                             complexMatrix&  accMatrix,
                                             vector<double>& phaseSpaceIntegral,
                                             const bool      withFlat) const
{
	const unsigned int size = _nmbWaves + (withFlat ? 1 : 0);
	normMatrix.resizeTo(size, size);
	accMatrix.resizeTo(size, size);
	phaseSpaceIntegral.resize(size, 0);
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
	if (withFlat) {
		// set unused entries to 0
		for (unsigned int i = 0; i < normMatrix.nCols(); ++i) {
			normMatrix.set(_nmbWaves, i, 0);
			normMatrix.set(i, _nmbWaves, 0);
			accMatrix.set (_nmbWaves, i, 0);
			accMatrix.set (i, _nmbWaves, 0);
		}
		// add flat
		normMatrix.set(_nmbWaves, _nmbWaves, 1.);
		accMatrix.set (_nmbWaves, _nmbWaves, _totAcc);
		phaseSpaceIntegral[_nmbWaves] = 1.;
	}
}


// builds complex numbers from parameters
// maps real and imaginary part of amplitudes to error matrix
// for both rank restrictions are taken into account
template<typename complexT>
void
pwaLikelihood<complexT>::buildProdAmpArrays(const double*             inPar,
                                            vector<complex<double> >& prodAmps,
                                            vector<pair<int, int> >&  parIndices,
                                            vector<string>&           prodAmpNames,
                                            const bool                withFlat) const
{
	prodAmps.clear();
	parIndices.clear();
	prodAmpNames.clear();
	unsigned int parIndex = 0;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {
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
		}
	}
	if (withFlat) {
		prodAmps.push_back(complex<double>(inPar[parIndex], 0));
		parIndices.push_back(make_pair(parIndex, -1));
		prodAmpNames.push_back("V_flat");
	}
}


template<typename complexT>
void
pwaLikelihood<complexT>::clear()
{
	_decayAmps[0].resize(extents[0][0]);
	_decayAmps[1].resize(extents[0][0]);
}


// copy values from array that corresponds to the function parameters
// to structure that corresponds to the complex production amplitudes
// taking into account rank restrictions
template<typename complexT>
void
pwaLikelihood<complexT>::copyFromParArray
(const double*      inPar,             // input parameter array
 prodAmpsArrayType& outVal,            // array of complex output values [rank][reflectivity][wave index]
 value_type&        outFlatVal) const  // output value corresponding to flat wave
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
pwaLikelihood<complexT>::copyToParArray
(const prodAmpsArrayType& inVal,         // values corresponding to production amplitudes
 const value_type         inFlatVal,     // value corresponding to flat wave
 double*                  outPar) const  // output parameter array
{
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				bt::tuple<int, int> parIndices = _prodAmpToFuncParMap[iRank][iRefl][iWave];
				if (bt::get<0>(parIndices) >= 0)  // real part
					outPar[bt::get<0>(parIndices)] = inVal[iRank][iRefl][iWave].real();
				if (bt::get<1>(parIndices) >= 0)  // imaginary part
					outPar[bt::get<1>(parIndices)] = inVal[iRank][iRefl][iWave].imag();
			}
	outPar[_nmbPars - 1] = inFlatVal;
}


template<typename complexT>
ostream&
pwaLikelihood<complexT>::print(ostream& out) const
{
	out << "pwaLikelihood parameters:" << endl
	    << "number of events ........................ " << _nmbEvents         << endl
	    << "rank .................................... " << _rank              << endl
	    << "number of waves ......................... " << _nmbWaves          << endl
	    << "number of positive reflectivity waves ... " << _nmbWavesRefl[1]   << endl
	    << "number of negative reflectivity waves ... " << _nmbWavesRefl[0]   << endl
	    << "number of function parameters ........... " << _nmbPars           << endl
	    << "number of fixed function parameters ..... " << _nmbParsFixed      << endl
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
			    << _waveThresholds[iRefl][iWave] << " GeV/c^2" << endl;
	out << "list of function parameters: " << endl;
	for (unsigned int iPar = 0; iPar < _nmbPars; ++iPar)
		out << "        [" << setw(3) << iPar << "] " << _parameters[iPar].parName() << "    "
		    << "threshold = " << _parameters[iPar].threshold() << " GeV/c^2" << endl;
	return out;
}


template<typename complexT>
ostream&
pwaLikelihood<complexT>::printFuncInfo(ostream& out) const
{
	const string funcNames[NMB_FUNCTIONCALLENUM] = {"FdF", "Gradient", "DoEval", "DoDerivative", "Hessian"};
	for (unsigned int i = 0; i < NMB_FUNCTIONCALLENUM; ++i)
		if (_funcCallInfo[i].nmbCalls > 0)
			out << "    " << _funcCallInfo[i].nmbCalls
			    << " calls to pwaLikelihood<complexT>::" << funcNames[i] << "()" << endl
			    << "    time spent in pwaLikelihood<complexT>::" << funcNames[i] << "(): " << endl
			    << "        total time ................. " << sum(_funcCallInfo[i].totalTime) << " sec" << endl
			    << "        return value calculation ... " << sum(_funcCallInfo[i].funcTime ) << " sec" << endl
			    << "        normalization .............. " << sum(_funcCallInfo[i].normTime ) << " sec" << endl;
	return out;
}


template<typename complexT>
void
pwaLikelihood<complexT>::resetFuncCallInfo() const
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
template class rpwa::pwaLikelihood<complex<float > >;
template class rpwa::pwaLikelihood<complex<double> >;
