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
//
// Description:
//      test Likelihood
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <complex>
#include <cassert>
#include <ctime>
#include <algorithm>

#define BOOST_DISABLE_ASSERTS
#include "boost/multi_array.hpp"

#include "TRandom3.h"

#include "sumAccumulators.hpp"
#include "reportingUtils.hpp"
#include "arrayUtils.hpp"

#include "complex.cuh"
#include "likelihoodInterface.cuh"


using namespace std;
using namespace boost;
using namespace boost::accumulators;
using namespace rpwa;


typedef multi_array<complex<double>, 3> ampsArrayType;  // host array type of production and decay amplitudes


template<typename T>
inline
T
cascadeSum(const std::vector<T>& data,
           const std::size_t     nmbOfSumsAtEachStage = 2)
{
	const std::vector<T>* sumsPrev       = &data;
	std::vector<T>*       sumsNext       = 0;
	bool                  firstIteration = true;
	do {
		std::size_t nmbSumsNext = sumsPrev->size() / nmbOfSumsAtEachStage;
		if (sumsPrev->size() % nmbOfSumsAtEachStage > 0)
			++nmbSumsNext;
		sumsNext = new std::vector<T>(nmbSumsNext, T());
		for (std::size_t i = 0; i < nmbSumsNext; ++i)
			for (std::size_t j = i * nmbOfSumsAtEachStage;
			     j < min((i + 1) * nmbOfSumsAtEachStage, sumsPrev->size()); ++j)
				(*sumsNext)[i] += (*sumsPrev)[j];
		// cleanup
		if (firstIteration)
			firstIteration = false;
		else
			delete sumsPrev;
		// prepare for next iteration
		sumsPrev = sumsNext;
	} while (sumsPrev->size() > 1);
	const T sum = (*sumsPrev)[0];
	delete sumsPrev;
	return sum;
}


template<typename complexT>
void
generateData(const unsigned int             nmbEvents,
             const unsigned int             rank,
             const unsigned int             nmbWavesRefl[2],
             multi_array<complexT, 3>&      decayAmps,
             multi_array<complexT, 3>&      prodAmps,
             typename complexT::value_type& prodAmpFlat,
             const unsigned int             seed = 123456789)
{
	TRandom3 random(seed);

	// set decay amplitudes
	decayAmps.resize(extents[nmbEvents][2][max(nmbWavesRefl[0], nmbWavesRefl[1])]);
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave)
			for (unsigned int iEvt = 0; iEvt < nmbEvents; ++iEvt) {
				//const complexT::value_type val = iEvt * 1000 + iRefl * 100 + iWave;
				//decayAmps[iEvt][iRefl][iWave] = complexT(val, val + 0.5);
				decayAmps[iEvt][iRefl][iWave] = complexT(random.Uniform(0, 1), random.Uniform(0, 1));
			}

	// set production amplitudes
	prodAmps.resize(extents[rank][2][max(nmbWavesRefl[0], nmbWavesRefl[1])]);
	//unsigned int parIndex = 0;
	for (unsigned int iRank = 0; iRank < rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {
				//const complexT::value_type val = iRank * 1000 + iRefl * 100 + iWave;
				//prodAmps[iRank][iRefl][iWave] = complexT(val, val + 0.5);
				prodAmps[iRank][iRefl][iWave] = complexT(random.Uniform(0, 1), random.Uniform(0, 1));
			}
	//prodAmpFlat = parIndex;
	prodAmpFlat = random.Uniform(0, 1);
}


template<typename complexT>
typename complexT::value_type
logLikelihoodMultiArray(const multi_array<complexT, 3>&      decayAmps,
                        const multi_array<complexT, 3>&      prodAmps,
                        const typename complexT::value_type& prodAmpFlat,
                        const unsigned int                   nmbEvents,
                        const unsigned int                   rank,
                        const unsigned int                   nmbWavesRefl[2])
{
	typedef typename complexT::value_type T;
	const T prodAmpFlat2 = prodAmpFlat * prodAmpFlat;
	// loop over events and calculate real-data term of log likelihood
	accumulator_set<T, stats<tag::sum(compensated)> > logLikelihoodAcc;
	for (unsigned int iEvt = 0; iEvt < nmbEvents; ++iEvt) {
		accumulator_set<T, stats<tag::sum(compensated)> > likelihoodAcc;
		for (unsigned int iRank = 0; iRank < rank; ++iRank) {  // incoherent sum over ranks
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
				accumulator_set<complexT, stats<tag::sum(compensated)> > ampProdAcc;
				for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave)  // coherent sum over waves
					ampProdAcc(prodAmps[iRank][iRefl][iWave] * decayAmps[iEvt][iRefl][iWave]);
				likelihoodAcc(norm(sum(ampProdAcc)));
			}
		}
		likelihoodAcc(prodAmpFlat2);
		logLikelihoodAcc(-log(sum(likelihoodAcc)));
	}

	return sum(logLikelihoodAcc);
}


template<typename complexT>
typename complexT::value_type
runLogLikelihoodMultiArray(const unsigned int nmbRepitions,
                           const unsigned int nmbEvents,
                           const unsigned int rank,
                           const unsigned int nmbWavesRefl[2],
                           double&            elapsedTime)
{
	multi_array<complexT, 3>      decayAmps;
	multi_array<complexT, 3>      prodAmps;
	typename complexT::value_type prodAmpFlat;
	generateData(nmbEvents, rank, nmbWavesRefl, decayAmps, prodAmps, prodAmpFlat);

	// call function
	clock_t                       start, finish;  // rough estimation of run time
	typename complexT::value_type logLikelihood = 0;
	start = clock();
	for (unsigned int i = 0; i < nmbRepitions; ++i)
		logLikelihood = logLikelihoodMultiArray(decayAmps, prodAmps, prodAmpFlat,
		                                        nmbEvents, rank, nmbWavesRefl);
	finish = clock();
	elapsedTime = ((double)(finish - start)) / CLOCKS_PER_SEC;  // [sec]

	return logLikelihood;
}


template<typename complexT>
void
logLikelihoodDerivMultiArray(const multi_array<complexT, 3>&      decayAmps,
                             const multi_array<complexT, 3>&      prodAmps,
                             const typename complexT::value_type& prodAmpFlat,
                             const unsigned int                   nmbEvents,
                             const unsigned int                   rank,
                             const unsigned int                   nmbWavesRefl[2],
                             multi_array<complexT, 3>&            derivatives,
                             typename complexT::value_type&       derivativeFlat)
{
	typedef typename complexT::value_type T;

	// create array of likelihood derivative w.r.t. real and imaginary
	// parts of the production amplitudes
	// !NOTE! although stored as and constructed from complex values,
	// the dL themselves are _not_ well defined complex numbers!
	derivatives.resize(extents[rank][2][max(nmbWavesRefl[0], nmbWavesRefl[1])]);

	// loop over events and calculate derivatives with respect to parameters
	accumulator_set<T, stats<tag::sum(compensated)> > derivativeFlatAcc;
	multi_array<accumulator_set<complexT, stats<tag::sum(compensated)> >, 3>
			derivativesAcc(extents[rank][2][max(nmbWavesRefl[0], nmbWavesRefl[1])]);
	const T prodAmpFlat2 = prodAmpFlat * prodAmpFlat;
	for (unsigned int iEvt = 0; iEvt < nmbEvents; ++iEvt) {
		accumulator_set<T, stats<tag::sum(compensated)> > likelihoodAcc;
		multi_array<complexT, 3> derivative(extents[rank][2][max(nmbWavesRefl[0], nmbWavesRefl[1])]);  // likelihood derivative for this event
		for (unsigned int iRank = 0; iRank < rank; ++iRank) {  // incoherent sum over ranks
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
				accumulator_set<complexT, stats<tag::sum(compensated)> > ampProdAcc;
				for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave)  // coherent sum over waves
					ampProdAcc(prodAmps[iRank][iRefl][iWave] * decayAmps[iEvt][iRefl][iWave]);
				const complexT ampProdSum = sum(ampProdAcc);
				likelihoodAcc(norm(ampProdSum));
				// set derivative term that is independent on derivative wave index
				for (unsigned int jWave = 0; jWave < nmbWavesRefl[iRefl]; ++jWave)
					// amplitude sums for current rank and for waves with same reflectivity
					derivative[iRank][iRefl][jWave] = ampProdSum;
			}
			// loop again over waves for current rank and multiply with complex conjugate
			// of decay amplitude of the wave with the derivative wave index
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int jWave = 0; jWave < nmbWavesRefl[iRefl]; ++jWave)
					derivative[iRank][iRefl][jWave] *= conj(complexT(decayAmps[iEvt][iRefl][jWave]));
		}  // end loop over rank
		likelihoodAcc(prodAmpFlat2);
		// incorporate factor 2 / sigma
		const T factor = 2. / sum(likelihoodAcc);
		for (unsigned int iRank = 0; iRank < rank; ++iRank)
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int jWave = 0; jWave < nmbWavesRefl[iRefl]; ++jWave)
					derivativesAcc[iRank][iRefl][jWave](-factor * derivative[iRank][iRefl][jWave]);
		derivativeFlatAcc(-factor * prodAmpFlat);
	}  // end loop over events
	for (unsigned int iRank = 0; iRank < rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave)
				derivatives[iRank][iRefl][iWave] = sum(derivativesAcc[iRank][iRefl][iWave]);
	derivativeFlat = sum(derivativeFlatAcc);
}


template<typename complexT>
void
runLogLikelihoodDerivMultiArray(const unsigned int             nmbRepitions,
                                const unsigned int             nmbEvents,
                                const unsigned int             rank,
                                const unsigned int             nmbWavesRefl[2],
                                multi_array<complexT, 3>&      derivatives,
                                typename complexT::value_type& derivativeFlat,
                                double&                        elapsedTime)
{
	multi_array<complexT, 3>      decayAmps;
	multi_array<complexT, 3>      prodAmps;
	typename complexT::value_type prodAmpFlat;
	generateData(nmbEvents, rank, nmbWavesRefl, decayAmps, prodAmps, prodAmpFlat);

	// call function
	clock_t start, finish;  // rough estimation of run time
	start = clock();
	for (unsigned int i = 0; i < nmbRepitions; ++i) {
		logLikelihoodDerivMultiArray(decayAmps, prodAmps, prodAmpFlat,
		                             nmbEvents, rank, nmbWavesRefl, derivatives, derivativeFlat);
	}
	finish = clock();
	elapsedTime = ((double)(finish - start)) / CLOCKS_PER_SEC;  // [sec]
	// for (unsigned int iRank = 0; iRank < rank; ++iRank)
	// 	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
	// 		for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {
	// 			complexT derivative = derivatives[iRank][iRefl][iWave];
	// 			cout  << "    derivative[" << iRank << "][" << iRefl << "][" << iWave << "] = "
	// 			      << "(" << maxPrecision(derivative.real()) << ", "
	// 			      << maxPrecision(derivative.real()) << ")" << endl;
	// 		}
}


template<typename complexT>
typename complexT::value_type
logLikelihoodPseudoArray(const complexT*                      decayAmps,
                         const complexT*                      prodAmps,
                         const typename complexT::value_type& prodAmpFlat,
                         const unsigned int                   nmbEvents,
                         const unsigned int                   rank,
                         const unsigned int                   nmbWavesRefl[2])
{
	typedef typename complexT::value_type T;
	const T            prodAmpFlat2   = prodAmpFlat * prodAmpFlat;
	const unsigned int prodAmpDim [3] = {rank,      2, max(nmbWavesRefl[0], nmbWavesRefl[1])};
	const unsigned int decayAmpDim[3] = {nmbEvents, 2, max(nmbWavesRefl[0], nmbWavesRefl[1])};
	// loop over events and calculate first term of log likelihood
	accumulator_set<T, stats<tag::sum(compensated)> > logLikelihoodAcc;
	for (unsigned int iEvt = 0; iEvt < nmbEvents; ++iEvt) {
		accumulator_set<T, stats<tag::sum(compensated)> > likelihoodAcc;
		for (unsigned int iRank = 0; iRank < rank; ++iRank) {  // incoherent sum over ranks
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
				accumulator_set<complexT, stats<tag::sum(compensated)> > ampProdAcc;
				for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
					// compute likelihood term
					const unsigned int prodAmpIndices [3] = {iRank, iRefl, iWave};
					const unsigned int decayAmpIndices[3] = {iEvt,  iRefl, iWave};
					ampProdAcc(  prodAmps [indicesToOffset<unsigned int>(prodAmpIndices,  prodAmpDim,  3)]
					           * decayAmps[indicesToOffset<unsigned int>(decayAmpIndices, decayAmpDim, 3)]);
				}
				likelihoodAcc(norm(sum(ampProdAcc)));
			}
		}
		likelihoodAcc(prodAmpFlat2);
		logLikelihoodAcc(-log(sum(likelihoodAcc)));
	}
  
	return sum(logLikelihoodAcc);
}


double
runLogLikelihoodPseudoArray(const unsigned int nmbRepitions,
                            const unsigned int nmbEvents,
                            const unsigned int rank,
                            const unsigned int nmbWavesRefl[2],
                            double&            elapsedTime)
{
	ampsArrayType decayAmps;
	ampsArrayType prodAmps;
	double        prodAmpFlat;
	generateData(nmbEvents, rank, nmbWavesRefl, decayAmps, prodAmps, prodAmpFlat);

	// call function
	clock_t start, finish;  // rough estimation of run time
	double  logLikelihood = 0;
	start = clock();
	for (unsigned int i = 0; i < nmbRepitions; ++i)
		if (0)
			logLikelihood = logLikelihoodPseudoArray<complex<double> >
				(decayAmps.data(), prodAmps.data(), prodAmpFlat, nmbEvents, rank, nmbWavesRefl);
		else
			logLikelihood = logLikelihoodPseudoArray<cuda::complex<double> >
				(reinterpret_cast<cuda::complex<double>*>(decayAmps.data()),
				 reinterpret_cast<cuda::complex<double>*>(prodAmps.data()),
				 prodAmpFlat, nmbEvents, rank, nmbWavesRefl);
	finish = clock();
	elapsedTime = ((double)(finish - start)) / CLOCKS_PER_SEC;  // [sec]

	return logLikelihood;
}


double
runLogLikelihoodCuda(const unsigned int nmbRepitions,
                     const unsigned int nmbEvents,
                     const unsigned int rank,
                     const unsigned int nmbWavesRefl[2],
                     double&            elapsedTime)
{
	ampsArrayType decayAmps;
	ampsArrayType prodAmps;
	double        prodAmpFlat;
	generateData(nmbEvents, rank, nmbWavesRefl, decayAmps, prodAmps, prodAmpFlat);

	// initialize CUDA environment
	cuda::likelihoodInterface<cuda::complex<double> >::init
		(reinterpret_cast<cuda::complex<double>*>(decayAmps.data()),
		 decayAmps.num_elements(), nmbEvents, nmbWavesRefl, true);
	//cuda::likelihoodInterface<cuda::complex<double> >::setDebug(true);

	// call function
	clock_t start, finish;  // rough estimation of run time
	double  logLikelihood = 0;
	start = clock();
	for (unsigned int i = 0; i < nmbRepitions; ++i)
		logLikelihood = cuda::likelihoodInterface<cuda::complex<double> >::logLikelihood
			(reinterpret_cast<cuda::complex<double>*>(prodAmps.data()),
			 prodAmps.num_elements(), prodAmpFlat, rank);
	finish = clock();
	elapsedTime = ((double)(finish - start)) / CLOCKS_PER_SEC;  // [sec]

	cuda::likelihoodInterface<cuda::complex<double> >::closeCudaDevice();
	return logLikelihood;
}


template<typename complexT>
void
runLogLikelihoodDerivCuda(const unsigned int             nmbRepitions,
                          const unsigned int             nmbEvents,
                          const unsigned int             rank,
                          const unsigned int             nmbWavesRefl[2],
                          multi_array<complexT, 3>&      derivatives,
                          typename complexT::value_type& derivativeFlat,
                          double&                        elapsedTime)
{
	multi_array<complexT, 3>      decayAmps;
	multi_array<complexT, 3>      prodAmps;
	typename complexT::value_type prodAmpFlat;
	generateData(nmbEvents, rank, nmbWavesRefl, decayAmps, prodAmps, prodAmpFlat);

	// initialize CUDA environment
	cuda::likelihoodInterface<cuda::complex<double> >::init
		(reinterpret_cast<cuda::complex<double>*>(decayAmps.data()),
		 decayAmps.num_elements(), nmbEvents, nmbWavesRefl, true);

	// call function
	clock_t start, finish;  // rough estimation of run time
	start = clock();
	derivatives.resize(extents[rank][2][max(nmbWavesRefl[0], nmbWavesRefl[1])]);
	for (unsigned int i = 0; i < nmbRepitions; ++i)
		cuda::likelihoodInterface<cuda::complex<double> >::logLikelihoodDeriv
			(reinterpret_cast<cuda::complex<double>*>(prodAmps.data()),
			 prodAmps.num_elements(), prodAmpFlat, rank,
			 reinterpret_cast<cuda::complex<double>*>(derivatives.data()), derivativeFlat);
	finish = clock();
	elapsedTime = ((double)(finish - start)) / CLOCKS_PER_SEC;  // [sec]

	cuda::likelihoodInterface<cuda::complex<double> >::closeCudaDevice();
}


template<typename T>
void
testLogLikelihoodDerivCuda(const unsigned int nmbRepitions,
                           const unsigned int _nmbEvents,
                           const unsigned int _rank,
                           const unsigned int _nmbWavesRefl[2])
{
	ampsArrayType _decayAmps;
	ampsArrayType prodAmps;
	double        prodAmpFlat;
	generateData(_nmbEvents, _rank, _nmbWavesRefl, _decayAmps, prodAmps, prodAmpFlat);
	const T prodAmpFlat2 = prodAmpFlat * prodAmpFlat;

	// initialize CUDA environment
	cuda::likelihoodInterface<cuda::complex<double> >& interface
		= cuda::likelihoodInterface<cuda::complex<double> >::instance();
	interface.setDebug(true);
	interface.init(reinterpret_cast<cuda::complex<double>*>(_decayAmps.data()),
	               _decayAmps.num_elements(), _nmbEvents, _nmbWavesRefl, true);

	for (unsigned int i = 0; i < nmbRepitions; ++i) {
		// create array of likelihood derivative w.r.t. real and imaginary
		// parts of the production amplitudes
		// !NOTE! although stored as and constructed from complex values,
		// the dL themselves are _not_ well defined complex numbers!
		array<ampsArrayType::index, 3> derivShape = {{ _rank, 2, max(_nmbWavesRefl[0], _nmbWavesRefl[1]) }};
		ampsArrayType                  derivatives(derivShape);
		T                              derivativeFlat = 0;

		// compute derivative for first term of log likelihood
#ifdef USE_CUDA
		// if (_useCuda) {
		ampsArrayType derivativesCuda(derivShape);
		T             derivativeFlatCuda = 0;
		cout << "#" << flush;
		cuda::likelihoodInterface<cuda::complex<T> >& interface
			= cuda::likelihoodInterface<cuda::complex<T> >::instance();
		interface.logLikelihoodDeriv(reinterpret_cast<cuda::complex<T>*>(prodAmps.data()),
		                             prodAmps.num_elements(), prodAmpFlat, _rank,
		                             reinterpret_cast<cuda::complex<T>*>(derivativesCuda.data()),
		                             derivativeFlatCuda);
		// logLikelihoodDerivMultiArray<complex<T> >(_decayAmps, prodAmps, prodAmpFlat, _nmbEvents, _rank,
		//                                           _nmbWavesRefl, derivativesCuda, derivativeFlatCuda);
		// } else
#endif
		{
			cout << "*" << flush;
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
							derivatives[iRank][iRefl][jWave] -= factor * derivative[iRank][iRefl][jWave];
				derivativeFlat -= factor * prodAmpFlat;
			}  // end loop over events
		}
	
		static complex<T> maxDiff = 0;
		for (unsigned int iRank = 0; iRank < _rank; ++iRank)
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
					complex<T> diff;
					diff.real() = 1 -   derivativesCuda[iRank][iRefl][iWave].real()
						/ derivatives[iRank][iRefl][iWave].real();
					diff.imag() = 1 -   derivativesCuda[iRank][iRefl][iWave].imag()
						/ derivatives[iRank][iRefl][iWave].imag();
					bool newMaxDiff = false;
					if (abs(diff.real()) > maxDiff.real()) {
						maxDiff.real() = abs(diff.real());
						newMaxDiff = true;
					}
					if (abs(diff.imag()) > maxDiff.imag()) {
						maxDiff.imag() = abs(diff.imag());
						newMaxDiff = true;
					}
					if (newMaxDiff)
						printInfo << "[" << iRank << "][" << iRefl << "][" << iWave << "]: " << maxDiff << "; "
						          << derivatives    [iRank][iRefl][iWave] << " vs. "
						          << derivativesCuda[iRank][iRefl][iWave] << endl;
				}
		static T maxDiffFlat = 0;
		const  T diffFlat    = abs(1 - derivativeFlatCuda / derivativeFlat);
		if (diffFlat > maxDiffFlat) {
			maxDiffFlat = diffFlat;
			printInfo << maxDiffFlat << endl;
		}

		for (unsigned int iRank = 0; iRank < _rank; ++iRank)
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					cout << "    [" << iRank << "][" << iRefl << "][" << iWave << "]: "
					     << maxPrecisionDouble(derivatives    [iRank][iRefl][iWave]) << " - "
					     << maxPrecisionDouble(derivativesCuda[iRank][iRefl][iWave]) << " = "
					     << maxPrecisionDouble(  derivatives[iRank][iRefl][iWave]
					                           - derivativesCuda[iRank][iRefl][iWave]) << endl;
	}

	interface.closeCudaDevice();
}


int
main(int    argc,
     char** argv)
{
	if (0) {
		const unsigned int nmbRepitions    = 2;
		const unsigned int nmbEvents       = 100000;
		const unsigned int nmbWavesRefl[2] = {2, 2};
		const unsigned int rank            = 1;
		
		printInfo << "running derivative test with parameters:"       << endl
		          << "    number of repitions ..................... " << nmbRepitions    << endl
		          << "    number of events ........................ " << nmbEvents       << endl
		          << "    rank of fit ............................. " << rank            << endl
		          << "    number of positive reflectivity waves ... " << nmbWavesRefl[1] << endl
		          << "    number of negative reflectivity waves ... " << nmbWavesRefl[0] << endl;
  
		testLogLikelihoodDerivCuda<double>(nmbRepitions, nmbEvents, rank, nmbWavesRefl);
	}


	if (0) {
		const unsigned int nmbRepitions    = 100;
		// setup parameters that roughly correspond to the pi- pi+ pi- PWA
		const unsigned int nmbEvents       = 10000;
		// 34 waves with positive, 7 waves with negative reflectivity, and flat wave: 42 in total
		const unsigned int nmbWavesRefl[2] = {7, 34};
		const unsigned int rank            = 2;

		printInfo << "running mock-up log likelihood calculation with parameters:"       << endl
		          << "    number of repitions ..................... " << nmbRepitions    << endl
		          << "    number of events ........................ " << nmbEvents       << endl
		          << "    rank of fit ............................. " << rank            << endl
		          << "    number of positive reflectivity waves ... " << nmbWavesRefl[1] << endl
		          << "    number of negative reflectivity waves ... " << nmbWavesRefl[0] << endl;
  
		double elapsedTime   [3];
		double logLikelihoods[3];

		//logLikelihoods[0]	= runLogLikelihoodMultiArray<cuda::complex<double> >
		logLikelihoods[0]	= runLogLikelihoodMultiArray<complex<double> >
			(nmbRepitions, nmbEvents, rank, nmbWavesRefl, elapsedTime[0]);
		logLikelihoods[1] = runLogLikelihoodPseudoArray
			(nmbRepitions, nmbEvents, rank, nmbWavesRefl, elapsedTime[1]);
		logLikelihoods[2] = runLogLikelihoodCuda
			(nmbRepitions, nmbEvents, rank, nmbWavesRefl, elapsedTime[2]);

		printInfo << "finished log likelihood calculation:" << endl
		          << "    elapsed time (multi_array) ...... " << elapsedTime[0] << " sec" << endl
		          << "    elapsed time (pseudoArray) ...... " << elapsedTime[1] << " sec" << endl
		          << "    elapsed time (CUDA) ............. " << elapsedTime[2] << " sec" << endl
		          << "    CUDA kernel time ................ "
		          << cuda::likelihoodInterface<cuda::complex<double> >::kernelTime() << " sec" << endl
		          << "    log(likelihood) (multi_array) ... " << maxPrecision(logLikelihoods[0]) << endl
		          << "    log(likelihood) (pseudoArray) ... " << maxPrecision(logLikelihoods[1]) << endl
		          << "    log(likelihood) (CUDA) .......... " << maxPrecision(logLikelihoods[2]) << endl
		          << "    delta[log(likelihood)] ........... "
		          << logLikelihoods[0] - logLikelihoods[2] << endl;
	}


	if (1) {
		const unsigned int nmbRepitions    = 100;
		// setup parameters that roughly correspond to the pi- pi+ pi- PWA
		const unsigned int nmbEvents       = 10000;
		// 34 waves with positive, 7 waves with negative reflectivity, and flat wave: 42 in total
		const unsigned int nmbWavesRefl[2] = {7, 34};
		const unsigned int rank            = 2;

		printInfo << "running mock-up log likelihood gradient calculation with parameters:" << endl
		          << "    number of repitions ..................... " << nmbRepitions    << endl
		          << "    number of events ........................ " << nmbEvents       << endl
		          << "    rank of fit ............................. " << rank            << endl
		          << "    number of positive reflectivity waves ... " << nmbWavesRefl[1] << endl
		          << "    number of negative reflectivity waves ... " << nmbWavesRefl[0] << endl;

		//typedef cuda::complex<double> complexT;
		typedef std::complex<double> complexT;
		double                   elapsedTime   [2];
		multi_array<complexT, 3> derivatives   [2];
		double                   derivativeFlat[2];

		runLogLikelihoodDerivMultiArray(nmbRepitions, nmbEvents, rank, nmbWavesRefl,
		                                derivatives[0], derivativeFlat[0], elapsedTime[0]);
		runLogLikelihoodDerivCuda(nmbRepitions, nmbEvents, rank, nmbWavesRefl,
		                          derivatives[1], derivativeFlat[1], elapsedTime[1]);

		complexT maxDiffAbs = 0;
		complexT maxDiffRel = 0;
		for (unsigned int iRank = 0; iRank < rank; ++iRank)
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {
					const complexT diffAbs =   derivatives[0][iRank][iRefl][iWave]
						                       - derivatives[1][iRank][iRefl][iWave];
					cout << "    [" << iRank << "][" << iRefl << "][" << iWave << "]: "
					     << maxPrecisionDouble(derivatives[0][iRank][iRefl][iWave]) << " - "
					     << maxPrecisionDouble(derivatives[1][iRank][iRefl][iWave]) << " = "
					     << maxPrecisionDouble(diffAbs) << endl;
					if (abs(diffAbs.real()) > maxDiffAbs.real())
						maxDiffAbs.real() = abs(diffAbs.real());
					if (abs(diffAbs.imag()) > maxDiffAbs.imag())
						maxDiffAbs.imag() = abs(diffAbs.imag());
					complexT diffRel;
					diffRel.real() = 1 -   derivatives[0][iRank][iRefl][iWave].real()
						/ derivatives[1][iRank][iRefl][iWave].real();
					diffRel.imag() = 1 -   derivatives[0][iRank][iRefl][iWave].imag()
						/ derivatives[1][iRank][iRefl][iWave].imag();
					if (abs(diffRel.real()) > maxDiffRel.real())
						maxDiffRel.real() = abs(diffRel.real());
					if (abs(diffRel.imag()) > maxDiffRel.imag())
						maxDiffRel.imag() = abs(diffRel.imag());
				}

		printInfo << "finished calculation of log likelihood derivatives:" << endl
		          << "    elapsed time (multi_array) ... " << elapsedTime[0]  << " sec" << endl
		          << "    elapsed time (CUDA) .......... " << elapsedTime[1]  << " sec" << endl
		          << "    max. absolute difference ..... " << maxDiffAbs                << endl
		          << "    max. relative difference ..... " << maxDiffRel                << endl
		          << "    flat deriv. (multi_array) .... " << maxPrecision(derivativeFlat[0]) << endl
		          << "    flat deriv. (CUDA) ........... " << maxPrecision(derivativeFlat[1]) << endl
		          << "    flat deriv. difference ....... "
		          << derivativeFlat[0] - derivativeFlat[1] << endl;
	}

	return 0;
}
