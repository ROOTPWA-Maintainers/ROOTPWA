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

#define BOOST_DISABLE_ASSERTS
#include "boost/multi_array.hpp"

#include "reportingUtils.hpp"
#include "nDimArrayUtils.hpp"
#include "nDimArrayUtils.hpp"

#include "complex.cuh"
#include "likelihoodInterface.cuh"


using namespace std;
using namespace boost;
using namespace rpwa;


typedef multi_array<complex<double>, 3> ampsArrayType;  // host array type of production and decay amplitudes


void
generateData(const unsigned int nmbEvents,
             const unsigned int rank,
             const unsigned int nmbWavesRefl[2],
             ampsArrayType&     decayAmps,
             ampsArrayType&     prodAmps,
             double&            prodAmpFlat)
{
	// set decay amplitudes
	decayAmps.resize(extents[nmbEvents][2][max(nmbWavesRefl[0], nmbWavesRefl[1])]);
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave)
			for (unsigned int iEvt = 0; iEvt < nmbEvents; ++iEvt) {
				const double val = iEvt * 1000 + iRefl * 100 + iWave;
				decayAmps[iEvt][iRefl][iWave] = complex<double>(val, val + 0.5);
			}

	// set production amplitudes
	prodAmps.resize(extents[rank][2][max(nmbWavesRefl[0], nmbWavesRefl[1])]);
	unsigned int parIndex = 0;
	for (unsigned int iRank = 0; iRank < rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {
				const double val = iRank * 1000 + iRefl * 100 + iWave;
				prodAmps[iRank][iRefl][iWave] = complex<double>(val, val + 0.5);
			}
	prodAmpFlat = parIndex;
}


template<typename complexT>
typename complexT::value_type
logLikelihoodMultiArray(const ampsArrayType&                 decayAmps,
                        const ampsArrayType&                 prodAmps,
                        const typename complexT::value_type& prodAmpFlat,
                        const unsigned int                   nmbEvents,
                        const unsigned int                   rank,
                        const unsigned int                   nmbWavesRefl[2])
{
	const typename complexT::value_type prodAmpFlat2 = prodAmpFlat * prodAmpFlat;
	// loop over events and calculate first term of log likelihood
	typename complexT::value_type logLikelihood = 0;
	for (unsigned int iEvt = 0; iEvt < nmbEvents; ++iEvt) {
		typename complexT::value_type likelihood = 0;  // likelihood for this event
		for (unsigned int iRank = 0; iRank < rank; ++iRank) {  // incoherent sum over ranks
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
				complexT ampProdSum = 0;  // amplitude sum for negative/positive reflectivity for this rank
				for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
					// compute likelihood term
					ampProdSum += prodAmps[iRank][iRefl][iWave] * decayAmps[iEvt][iRefl][iWave];
				}
				likelihood += norm(ampProdSum);
			}
			assert(likelihood >= 0);
		}  // end loop over rank
		likelihood    += prodAmpFlat2;
		logLikelihood -= log(likelihood);  // accumulate log likelihood
	}  // end loop over events

	return logLikelihood;
}


double
runLogLikelihoodMultiArray(const unsigned int nmbRepitions,
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
	double  logLikelihood;
	start = clock();
	for (unsigned int i = 0; i < nmbRepitions; ++i)
		logLikelihood = logLikelihoodMultiArray<complex<double> >
			(decayAmps, prodAmps, prodAmpFlat, nmbEvents, rank, nmbWavesRefl);
	finish = clock();
	elapsedTime = ((double)(finish - start)) / CLOCKS_PER_SEC;  // [sec]

	return logLikelihood;
}


template<typename complexT>
void
logLikelihoodDerivMultiArray(const ampsArrayType&                 decayAmps,
                             const ampsArrayType&                 prodAmps,
                             const typename complexT::value_type& prodAmpFlat,
                             const unsigned int                   nmbEvents,
                             const unsigned int                   rank,
                             const unsigned int                   nmbWavesRefl[2],
                             ampsArrayType&                       derivatives,
                             typename complexT::value_type&       derivativeFlat)
{
	// create array of likelihood derivative w.r.t. real and imaginary
	// parts of the production amplitudes
	// !NOTE! although stored as and constructed from complex values,
	// the dL themselves are _not_ well defined complex numbers!
	derivatives.resize(extents[rank][2][max(nmbWavesRefl[0], nmbWavesRefl[1])]);
	derivativeFlat = 0;

	// compute derivative for first term of log likelihood
	const typename complexT::value_type prodAmpFlat2 = prodAmpFlat * prodAmpFlat;
	for (unsigned int iEvt = 0; iEvt < nmbEvents; ++iEvt) {
		typename complexT::value_type likelihood = 0;  // likelihood for this event
		ampsArrayType                 derivative(extents[rank][2][max(nmbWavesRefl[0], nmbWavesRefl[1])]);  // likelihood derivative for this event
		for (unsigned int iRank = 0; iRank < rank; ++iRank) {  // incoherent sum over ranks
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
				complexT ampProdSum = 0;  // amplitude sum for negative/positive reflectivity for this rank
				for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
					// compute likelihood term
					const complexT amp =   prodAmps[iRank][iRefl][iWave]
						                   * complexT(decayAmps[iEvt][iRefl][iWave]);
					ampProdSum += amp;
				}
				likelihood += norm(ampProdSum);
				// set derivative term that is independent on derivative wave index
				for (unsigned int jWave = 0; jWave < nmbWavesRefl[iRefl]; ++jWave)
					// amplitude sums for current rank and for waves with same reflectivity
					derivative[iRank][iRefl][jWave] = ampProdSum;
			}
			assert(likelihood >= 0);
			// loop again over waves for current rank and multiply with complex conjugate
			// of decay amplitude of the wave with the derivative wave index
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int jWave = 0; jWave < nmbWavesRefl[iRefl]; ++jWave)
					derivative[iRank][iRefl][jWave] *= conj(complexT(decayAmps[iEvt][iRefl][jWave]));
		}  // end loop over rank
		likelihood += prodAmpFlat2;
		// incorporate factor 2 / sigma
		const typename complexT::value_type factor = 2. / likelihood;
		for (unsigned int iRank = 0; iRank < rank; ++iRank)
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int jWave = 0; jWave < nmbWavesRefl[iRefl]; ++jWave)
					derivatives[iRank][iRefl][jWave] -= factor * derivative[iRank][iRefl][jWave];
		derivativeFlat -= factor * prodAmpFlat;
	}  // end loop over events
}


void
runLogLikelihoodDerivMultiArray(const unsigned int nmbRepitions,
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
	clock_t       start, finish;  // rough estimation of run time
	ampsArrayType derivatives;
	double        derivativeFlat;
	start = clock();
	for (unsigned int i = 0; i < nmbRepitions; ++i)
		logLikelihoodDerivMultiArray<complex<double> >
			(decayAmps, prodAmps, prodAmpFlat, nmbEvents, rank, nmbWavesRefl, derivatives, derivativeFlat);
	finish = clock();
	elapsedTime = ((double)(finish - start)) / CLOCKS_PER_SEC;  // [sec]
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
	const typename complexT::value_type prodAmpFlat2   = prodAmpFlat * prodAmpFlat;
	const unsigned int                  prodAmpDim [3] = {rank,      2, max(nmbWavesRefl[0], nmbWavesRefl[1])};
	const unsigned int                  decayAmpDim[3] = {nmbEvents, 2, max(nmbWavesRefl[0], nmbWavesRefl[1])};
	// loop over events and calculate first term of log likelihood
	typename complexT::value_type logLikelihood = 0;
	for (unsigned int iEvt = 0; iEvt < nmbEvents; ++iEvt) {
		typename complexT::value_type likelihood = 0;  // likelihood for this event
		for (unsigned int iRank = 0; iRank < rank; ++iRank) {  // incoherent sum over ranks
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
				complexT ampProdSum = 0;  // amplitude sum for negative/positive reflectivity for this rank
				for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
					// compute likelihood term
					const unsigned int prodAmpIndices [3] = {iRank, iRefl, iWave};
					const unsigned int decayAmpIndices[3] = {iEvt,  iRefl, iWave};
					ampProdSum +=   prodAmps [indicesToOffset<unsigned int>(prodAmpIndices,  prodAmpDim,  3)]
						* decayAmps[indicesToOffset<unsigned int>(decayAmpIndices, decayAmpDim, 3)];
				}
				likelihood += norm(ampProdSum);
			}
			assert(likelihood >= 0);
		}  // end loop over rank
		likelihood    += prodAmpFlat2;
		logLikelihood -= log(likelihood);  // accumulate log likelihood
	}  // end loop over events
  
	return logLikelihood;
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
	double  logLikelihood;
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
	cuda::likelihoodInterface<cuda::complex<double> >& interface
		= cuda::likelihoodInterface<cuda::complex<double> >::instance();
	interface.setDebug(true);
	interface.init(reinterpret_cast<cuda::complex<double>*>(decayAmps.data()),
	               decayAmps.num_elements(), nmbEvents, nmbWavesRefl, true);

	// call function
	clock_t start, finish;  // rough estimation of run time
	double  logLikelihood;
	start = clock();
	for (unsigned int i = 0; i < nmbRepitions; ++i)
		logLikelihood = interface.logLikelihood
			(reinterpret_cast<cuda::complex<double>*>(prodAmps.data()),
			 prodAmps.num_elements(), prodAmpFlat, rank);
	finish = clock();
	elapsedTime = ((double)(finish - start)) / CLOCKS_PER_SEC;  // [sec]

	return logLikelihood;
}


int
main(int    argc,
     char** argv)
{
	const unsigned int nmbRepitions    = 10;
	// setup parameters that roughly correspond to the pi- pi+ pi- PWA
	const unsigned int nmbEvents       = 100000;
	// 34 waves with positive, 7 waves with negative reflectivity, and flat wave: 42 in total
	const unsigned int nmbWavesRefl[2] = {7, 34};
	const unsigned int rank            = 2;
  
	double elapsedTime   [3];
	double logLikelihoods[3];

	runLogLikelihoodDerivMultiArray(nmbRepitions, nmbEvents, rank, nmbWavesRefl, elapsedTime[0]);

	// logLikelihoods[0]
	// 	= runLogLikelihoodMultiArray (nmbRepitions, nmbEvents, rank, nmbWavesRefl, elapsedTime[0]);

	// logLikelihoods[1]
	// 	= runLogLikelihoodPseudoArray(nmbRepitions, nmbEvents, rank, nmbWavesRefl, elapsedTime[1]);

	// logLikelihoods[2]
	// 	= runLogLikelihoodCuda       (nmbRepitions, nmbEvents, rank, nmbWavesRefl, elapsedTime[2]);

	printInfo << "ran mock-up log likelihood calculation with parameters:"           << endl
	          << "    number of repitions ..................... " << nmbRepitions    << endl
	          << "    number of events ........................ " << nmbEvents       << endl
	          << "    rank of fit ............................. " << rank            << endl
	          << "    number of positive reflectivity waves ... " << nmbWavesRefl[1] << endl
	          << "    number of negative reflectivity waves ... " << nmbWavesRefl[0] << endl
	          << "    elapsed time (multiArray) ............... " << elapsedTime[0]  << " sec" << endl
	          << "    elapsed time (pseudo array) ............. " << elapsedTime[1]  << " sec" << endl
	          << "    elapsed time (CUDA) ..................... " << elapsedTime[2]  << " sec" << endl
	          << "    log(likelihood) (multiArray) ............ " << maxPrecision(logLikelihoods[0]) << endl
	          << "    log(likelihood) (pseudo array) .......... " << maxPrecision(logLikelihoods[1]) << endl
	          << "    log(likelihood) (CUDA) .................. " << maxPrecision(logLikelihoods[2]) << endl
	          << "    delta[log(likelihood)] .................. "
	          << maxPrecision(logLikelihoods[0] - logLikelihoods[2]) << endl;
  
	return 0;
}
