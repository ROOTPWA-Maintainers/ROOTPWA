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
#include "cudaLikelihoodInterface.h"


using namespace std;
using namespace boost;
using namespace rpwa;


typedef multi_array<std::complex<double>, 3> ampsArrayType;  // array for production and decay amplitudes


template<typename complexT, typename T>
T
logLikelihoodSumMultiArray(const ampsArrayType& decayAmps,
                           const ampsArrayType& prodAmps,
                           const T&             prodAmpFlat,
                           const unsigned int   nmbEvents,
                           const unsigned int   rank,
                           const unsigned int   nmbWavesRefl[2])
{
	const T prodAmpFlat2 = prodAmpFlat * prodAmpFlat;
	// loop over events and calculate first term of log likelihood
	T logLikelihood = 0;
	for (unsigned int iEvt = 0; iEvt < nmbEvents; ++iEvt) {
		T l = 0;  // likelihood for this event
		for (unsigned int iRank = 0; iRank < rank; ++iRank) {  // incoherent sum over ranks
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
				complexT ampProdSum = 0;  // amplitude sum for negative/positive reflectivity for this rank
				for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
					// compute likelihood term
					ampProdSum += prodAmps[iRank][iRefl][iWave] * decayAmps[iRefl][iWave][iEvt];
				}
				l += norm(ampProdSum);
			}
			assert(l >= 0);
		}  // end loop over rank
		l             += prodAmpFlat2;
		logLikelihood -= log(l);  // accumulate log likelihood
	}  // end loop over events

	return logLikelihood;
}


double
runLogLikelihoodSumMultiArray(const unsigned int nmbRepitions,
                              const unsigned int nmbEvents,
                              const unsigned int rank,
                              const unsigned int nmbWavesRefl[2],
                              double&            elapsedTime)
{
	// set decay amplitudes
	ampsArrayType decayAmps;
	decayAmps.resize(extents[2][max(nmbWavesRefl[0], nmbWavesRefl[1])][nmbEvents]);
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave)
			for (unsigned int iEvt = 0; iEvt < nmbEvents; ++iEvt) {
				const double val = iEvt * 1000 + iRefl * 100 + iWave;
				decayAmps[iRefl][iWave][iEvt] = std::complex<double>(val, val + 0.5);
			}

	// set production amplitudes
	ampsArrayType prodAmps;
	prodAmps.resize(extents[rank][2][max(nmbWavesRefl[0], nmbWavesRefl[1])]);
	unsigned int parIndex = 0;
	for (unsigned int iRank = 0; iRank < rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {
				const double val = iRank * 1000 + iRefl * 100 + iWave;
				prodAmps[iRank][iRefl][iWave] = std::complex<double>(val, val + 0.5);
			}
	const double prodAmpFlat = parIndex;

	// call function
	clock_t start, finish;  // rough estimation of run time
	double  logLikelihood;
	start = clock();
	for (unsigned int i = 0; i < nmbRepitions; ++i)
		logLikelihood = logLikelihoodSumMultiArray<std::complex<double>, double>
			(decayAmps, prodAmps, prodAmpFlat, nmbEvents, rank, nmbWavesRefl);
	finish = clock();
	elapsedTime = ((double)(finish - start)) / CLOCKS_PER_SEC;  // [sec]

	return logLikelihood;
}


template<typename complexT, typename T>
T
logLikelihoodSumPseudoArray(const complexT*    decayAmps,
                            const complexT*    prodAmps,
                            const T&           prodAmpFlat,
                            const unsigned int nmbEvents,
                            const unsigned int rank,
                            const unsigned int nmbWavesRefl[2])
{
	const T            prodAmpFlat2   = prodAmpFlat * prodAmpFlat;
	const unsigned int prodAmpDim [3] = {rank,      2, max(nmbWavesRefl[0], nmbWavesRefl[1])};
	const unsigned int decayAmpDim[3] = {2, max(nmbWavesRefl[0], nmbWavesRefl[1]), nmbEvents};
	// loop over events and calculate first term of log likelihood
	T logLikelihood = 0;
	for (unsigned int iEvt = 0; iEvt < nmbEvents; ++iEvt) {
		T l = 0;  // likelihood for this event
		for (unsigned int iRank = 0; iRank < rank; ++iRank) {  // incoherent sum over ranks
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
				complexT ampProdSum = 0;  // amplitude sum for negative/positive reflectivity for this rank
				for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
					// compute likelihood term
					const unsigned int prodAmpIndices [3] = {iRank, iRefl, iWave};
					const unsigned int decayAmpIndices[3] = {iRefl, iWave, iEvt};
					// ampProdSum +=   prodAmps [indicesToOffset<unsigned int>(prodAmpIndices,  prodAmpDim,  3)]
					// 	            * decayAmps[indicesToOffset<unsigned int>(decayAmpIndices, decayAmpDim, 3)];
					ampProdSum =   ampProdSum
						           +   prodAmps [indicesToOffset<unsigned int>(prodAmpIndices,  prodAmpDim,  3)]
					 	             * decayAmps[indicesToOffset<unsigned int>(decayAmpIndices, decayAmpDim, 3)];
					// cout << "    [" << iEvt << "]["  << iRank << "][" << iRefl << "][" << iWave << "] = "
					//      // << "(" << real(ampProdSum) << ", " << imag(ampProdSum) << ")"
					//      << "("
					//      << real(prodAmps [indicesToOffset<unsigned int>(prodAmpIndices,  prodAmpDim,  3)])
					//      << ", "
					//      << imag(prodAmps [indicesToOffset<unsigned int>(prodAmpIndices,  prodAmpDim,  3)])
					//      << "), ("
					//      << real(decayAmps[indicesToOffset<unsigned int>(decayAmpIndices, decayAmpDim, 3)])
					//      << ", "
					//      << imag(decayAmps[indicesToOffset<unsigned int>(decayAmpIndices, decayAmpDim, 3)])
					//      << ")" << endl;
				}
				l += norm(ampProdSum);
			}
			assert(l >= 0);
		}  // end loop over rank
		l             += prodAmpFlat2;
		logLikelihood -= log(l);  // accumulate log likelihood
	}  // end loop over events
  
	return logLikelihood;
}


double
runLogLikelihoodSumPseudoArray(const unsigned int nmbRepitions,
                               const unsigned int nmbEvents,
                               const unsigned int rank,
                               const unsigned int nmbWavesRefl[2],
                               double&            elapsedTime)
{
	// set decay amplitudes
	const unsigned int    decayAmpDim[3] = {2, max(nmbWavesRefl[0], nmbWavesRefl[1]), nmbEvents};
	std::complex<double>* decayAmps;
	const unsigned int    decayAmpsSize  = allocatePseudoNdimArray<std::complex<double>, unsigned int>
		(decayAmps, decayAmpDim, 3);
	printInfo << "size of decay amplitude array is " << decayAmpsSize / (1024. * 1024.) << " MiBytes; "
	          << 100 * nmbEvents * (nmbWavesRefl[0] + nmbWavesRefl[1]) * sizeof(std::complex<double>)
		/ (double)decayAmpsSize << " % used" << endl;
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave)
			for (unsigned int iEvt = 0; iEvt < nmbEvents; ++iEvt) {
				const unsigned int decayAmpIndices[3] = {iRefl, iWave, iEvt};
				const unsigned int decayAmpOffset     = indicesToOffset<unsigned int>(decayAmpIndices,
					                                                                    decayAmpDim, 3);
				const double       val                = iEvt * 1000 + iRefl * 100 + iWave;
				decayAmps[decayAmpOffset] = std::complex<double>(val, val + 0.5);
			}

	// set production amplitudes
	const unsigned int    prodAmpDim[3] = {rank, 2, max(nmbWavesRefl[0], nmbWavesRefl[1])};
	std::complex<double>* prodAmps;
	const unsigned int    prodAmpsSize  = allocatePseudoNdimArray<std::complex<double>, unsigned int>
		(prodAmps, prodAmpDim, 3);
	printInfo << "size of production amplitude array is " << prodAmpsSize << " bytes" << endl;
	unsigned int parIndex = 1;
	for (unsigned int iRank = 0; iRank < rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {
				const unsigned int prodAmpIndices[3] = {iRank, iRefl, iWave};
				const unsigned int prodAmpOffset     = indicesToOffset<unsigned int>(prodAmpIndices,
				                                                                     prodAmpDim, 3);
				const double       val               = iRank * 1000 + iRefl * 100 + iWave;
				prodAmps[prodAmpOffset] = std::complex<double>(val, val + 0.5);
			}
	const double prodAmpFlat = parIndex;

	// call function
	clock_t start, finish;  // rough estimation of run time
	double  logLikelihood;
	start = clock();
	for (unsigned int i = 0; i < nmbRepitions; ++i)
		logLikelihood = logLikelihoodSumPseudoArray<std::complex<double>, double>
			(decayAmps, prodAmps, prodAmpFlat, nmbEvents, rank, nmbWavesRefl);
	finish = clock();
	elapsedTime = ((double)(finish - start)) / CLOCKS_PER_SEC;  // [sec]

	delete[] decayAmps;
	delete[] prodAmps;

	return logLikelihood;
}


double
runLogLikelihoodSumCuda(const unsigned int nmbRepitions,
                        const unsigned int nmbEvents,
                        const unsigned int rank,
                        const unsigned int nmbWavesRefl[2],
                        double&            elapsedTime)
{
	// set decay amplitudes
	ampsArrayType decayAmps;
	decayAmps.resize(extents[2][max(nmbWavesRefl[0], nmbWavesRefl[1])][nmbEvents]);
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave)
			for (unsigned int iEvt = 0; iEvt < nmbEvents; ++iEvt) {
				const double val = iEvt * 1000 + iRefl * 100 + iWave;
				decayAmps[iRefl][iWave][iEvt] = std::complex<double>(val, val + 0.5);
			}

	// set production amplitudes
	ampsArrayType prodAmps;
	prodAmps.resize(extents[rank][2][max(nmbWavesRefl[0], nmbWavesRefl[1])]);
	unsigned int parIndex = 0;
	for (unsigned int iRank = 0; iRank < rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {
				const double val = iRank * 1000 + iRefl * 100 + iWave;
				prodAmps[iRank][iRefl][iWave] = std::complex<double>(val, val + 0.5);
			}
	const double prodAmpFlat = parIndex;

	if (1) {
		// test low-level access to multi_arrays
		const rpwa::complex<double>* prodAmpsPseudo = reinterpret_cast<rpwa::complex<double>*>(prodAmps.data());
		const unsigned int           prodAmpDim[3]  = {rank, 2, max(nmbWavesRefl[0], nmbWavesRefl[1])};
		for (unsigned int iRank = 0; iRank < rank; ++iRank)
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {
					const unsigned int prodAmpIndices[3] = {iRank, iRefl, iWave};
					if (
					    (   real(prodAmps[iRank][iRefl][iWave])
					     != real(prodAmpsPseudo[indicesToOffset<unsigned int>(prodAmpIndices, prodAmpDim, 3)]))
					    ||
					    (   imag(prodAmps[iRank][iRefl][iWave])
					     != imag(prodAmpsPseudo[indicesToOffset<unsigned int>(prodAmpIndices, prodAmpDim, 3)]))
					   )
						printWarn <<"boost[" << iRank << "][" << iRefl << "][" << iWave << "] = "
						          << prodAmps[iRank][iRefl][iWave]
						          << " vs. pseudo[" << iRank << "][" << iRefl << "][" << iWave << "] = ("
						          << real(prodAmpsPseudo[indicesToOffset<unsigned int>(prodAmpIndices, prodAmpDim, 3)])
						          << ", "
						          << imag(prodAmpsPseudo[indicesToOffset<unsigned int>(prodAmpIndices, prodAmpDim, 3)])
						          << ")" << endl;
				}
		const unsigned int size = prodAmps.num_elements() * sizeof(std::complex<double>);
		if (size != rank * 2 * max(nmbWavesRefl[0], nmbWavesRefl[1]) * sizeof(rpwa::complex<double>))
			printWarn << "boost size = " << size << " vs. pseudo size = "
			          << rank * 2 * max(nmbWavesRefl[0], nmbWavesRefl[1]) * sizeof(rpwa::complex<double>)
			          << " bytes" << endl;
	}

	// initialize CUDA environment
	rpwa::complex<double>* d_decayAmps;
	unsigned int           nmbBlocks;
	unsigned int           nmbThreadsPerBlock;
	initLogLikelihoodCuda(reinterpret_cast<rpwa::complex<double>*>(decayAmps.data()),
	                      decayAmps.num_elements() * sizeof(rpwa::complex<double>),
	                      d_decayAmps,
	                      nmbBlocks,
	                      nmbThreadsPerBlock);

	// call function
	clock_t start, finish;  // rough estimation of run time
	double  logLikelihood;
	start = clock();
	for (unsigned int i = 0; i < nmbRepitions; ++i)
		logLikelihood = sumLogLikelihoodCuda(reinterpret_cast<rpwa::complex<double>*>(prodAmps.data()),
		                                     prodAmps.num_elements() * sizeof(rpwa::complex<double>),
		                                     prodAmpFlat,
		                                     d_decayAmps,
		                                     nmbEvents,
		                                     rank,
		                                     nmbWavesRefl,
                                         nmbBlocks,
		                                     nmbThreadsPerBlock);
	finish = clock();
	elapsedTime = ((double)(finish - start)) / CLOCKS_PER_SEC;  // [sec]

	// cout << endl << "-----------------------------------------------------------------" << endl << endl;
	// const double foo = logLikelihoodSumPseudoArray<rpwa::complex<double>, double>
	// 		(reinterpret_cast<rpwa::complex<double>*>(decayAmps.data()),
	// 		 reinterpret_cast<rpwa::complex<double>*>(prodAmps.data()),
	// 		 prodAmpFlat, nmbEvents, rank, nmbWavesRefl);
	// printInfo << logLikelihood << " vs. " << foo << ": " << logLikelihood - foo << endl;

	//cutilSafeCall(cudaFree(d_decayAmps));
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
  
	double elapsedTime  [3];
	double logLikelihood[3];
	logLikelihood[0] = runLogLikelihoodSumMultiArray (nmbRepitions, nmbEvents, rank, nmbWavesRefl,
	                                                  elapsedTime[0]);

	logLikelihood[1] = runLogLikelihoodSumPseudoArray(nmbRepitions, nmbEvents, rank, nmbWavesRefl,
	                                                  elapsedTime[1]);

	logLikelihood[2] = runLogLikelihoodSumCuda       (nmbRepitions, nmbEvents, rank, nmbWavesRefl,
	                                                  elapsedTime[2]);

	printInfo << "ran mock-up log likelihood calculation with parameters:"            << endl
	          << "    number of repitions ..................... " << nmbRepitions     << endl
	          << "    number of events ........................ " << nmbEvents        << endl
	          << "    rank of fit ............................. " << rank             << endl
	          << "    number of positive reflectivity waves ... " << nmbWavesRefl[1]  << endl
	          << "    number of negative reflectivity waves ... " << nmbWavesRefl[0]  << endl
	          << "    elapsed time (multiArray) ............... " << elapsedTime[0]   << " sec" << endl
	          << "    elapsed time (pseudo array) ............. " << elapsedTime[1]   << " sec" << endl
	          << "    elapsed time (CUDA) ..................... " << elapsedTime[2]   << " sec" << endl
	          << "    log(likelihood) (multiArray) ............ " << logLikelihood[0] << endl
	          << "    log(likelihood) (pseudo array) .......... " << logLikelihood[1] << endl
	          << "    log(likelihood) (CUDA) .................. " << logLikelihood[2] << endl
	          << "    delta[log(likelihood)] .................. "
	          << logLikelihood[0] - logLikelihood[2] << endl;
  
	return 0;
}
