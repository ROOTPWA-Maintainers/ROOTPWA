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
//      CUDA kernels for summation of large arrays and for likelihood
//      calculation
//
//
// Author List:
//      Philipp Meyer          TUM            (original author)
//
//
//-------------------------------------------------------------------------



// cascadable kernel that computes sum of N values of an array
template<typename T>
__global__
void
sumKernel
(const T*           d_valArray,  // array with values to sum up
 const unsigned int nmbVals,     // total number of values all kernels have to process
 T*                 d_sumArray)  // output array of partial sums with one entry for each kernel
{
	const unsigned int threadId   = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int nmbThreads = gridDim.x * blockDim.x;
  
	// #pragma unroll 16
	T sum = T(0);
	for (unsigned int i = threadId; i < nmbVals; i += nmbThreads)
		sum += d_valArray[i];

	d_sumArray[threadId] = sum;
}


// kernel that calculates real-data term of log likelihood sum for a number of events
template<typename complexT>
__global__
void
sumLogLikelihoodKernel
(const complexT*                     d_prodAmps,           // 3-dim array of production amplitudes [iRank][iRefl][iWave]
 const typename complexT::value_type prodAmpFlat2,         // squared production amplitude of flat wave
 const complexT*                     d_decayAmps,          // 3-dim array of decay amplitudes [iRefl][iWave][iEvt]
 const unsigned int                  nmbEvents,            // total number of events all kernels have to process
 const unsigned int                  rank,                 // rank of spin-density matrix
 const unsigned int                  nmbWavesReflNeg,      // number waves with negative reflectivity
 const unsigned int                  nmbWavesReflPos,      // number waves with positive reflectivity
 const unsigned int                  nmbWavesMax,          // maximum extent of iWave index for production and decay amplitude arrays
 typename complexT::value_type*      d_logLikelihoodSums)  // output array of partial log likelihood sums with one entry for each kernel
{
	const unsigned int threadId   = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int nmbThreads = gridDim.x * blockDim.x;

	const unsigned int nmbWavesRefl[2] = {nmbWavesReflNeg, nmbWavesReflPos};
	// define extents of arrays
	const unsigned int prodAmpDim [3] = {rank, 2,           nmbWavesMax};
	const unsigned int decayAmpDim[3] = {2,    nmbWavesMax, nmbEvents};
	// loop over events and calculate real-data term of log likelihood
	typename complexT::value_type logLikelihoodSum = complexT::value_type(0);
	for (unsigned int iEvt = threadId; iEvt < nmbEvents; iEvt += nmbThreads) {
		typename complexT::value_type l = complexT::value_type(0);  // likelihood for this event
		for (unsigned int iRank = 0; iRank < rank; ++iRank) {  // incoherent sum over ranks
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
				complexT ampProdSum = complexT(0);  // amplitude sum for negative/positive reflectivity for this rank
				for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
					// compute likelihood term
					const unsigned int prodAmpIndices [3] = {iRank, iRefl, iWave};
					const unsigned int decayAmpIndices[3] = {iRefl, iWave, iEvt};
					ampProdSum +=   d_prodAmps [indicesToOffset<unsigned int>(prodAmpIndices,  prodAmpDim,  3)]
						            * d_decayAmps[indicesToOffset<unsigned int>(decayAmpIndices, decayAmpDim, 3)];
				}
				l += norm(ampProdSum);
			}
		}
		l                += prodAmpFlat2;
		logLikelihoodSum -= log(l);  // accumulate log likelihood
  }
	// write result
	d_logLikelihoodSums[threadId] = logLikelihoodSum;
}
