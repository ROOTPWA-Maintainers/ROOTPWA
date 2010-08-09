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
//      CUDA Kernels for likelihood and summation
//
//
// Author List:
//      Philipp Meyer          TUM            (original author)
//
//
//-------------------------------------------------------------------------


template<typename T>
__global__
void
SumCUDA (const T*           valArray,
         const unsigned int nmbVals,
         T*                 sumArray)
{
	const unsigned int threadId   = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int nmbThreads = gridDim.x * blockDim.x;
  
	// #pragma unroll 16
	T sum = T(0);
	for (unsigned int i = threadId; i < nmbVals; i += nmbThreads)
		sum += valArray[i];

	sumArray[threadId] = sum;
}


template<typename complexT, typename scalarT>
__global__
void
SumCUDAPseudoArray(const complexT*     decayAmps,
                   const complexT*     prodAmps,
                   const scalarT       prodAmpFlat,
                   const unsigned int  nmbEvents,
                   const unsigned int  rank,
                   // const unsigned int* nmbWavesRefl,
                   const unsigned int  nmbWavesRefl0,
                   const unsigned int  nmbWavesRefl1,
                   const unsigned int  nmbWavesReflMax,
                   scalarT*            logLikelihoodSums)
{
	const unsigned int threadId   = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int nmbThreads = gridDim.x * blockDim.x;

	const scalarT      prodAmpFlat2    = prodAmpFlat * prodAmpFlat;
	const unsigned int prodAmpDim [3]  = {rank, 2, nmbWavesReflMax};
	const unsigned int decayAmpDim[3]  = {2, nmbWavesReflMax, nmbEvents};
	const unsigned int nmbWavesRefl[2] = {nmbWavesRefl0, nmbWavesRefl1};
	// loop over events and calculate first term of log likelihood
	scalarT logLikelihood = scalarT(0);
	for (unsigned int iEvt = threadId; iEvt < nmbEvents; iEvt += nmbThreads) {
		scalarT l = scalarT(0);  // likelihood for this event
		for (unsigned int iRank = 0; iRank < rank; ++iRank) {  // incoherent sum over ranks
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
				complexT ampProdSum = complexT(0);  // amplitude sum for negative/positive reflectivity for this rank
				for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
					// compute likelihood term
					const unsigned int prodAmpIndices [3] = {iRank, iRefl, iWave};
					const unsigned int decayAmpIndices[3] = {iRefl, iWave, iEvt};
					ampProdSum =   ampProdSum
						           +   prodAmps [indicesToOffset<unsigned int>(prodAmpIndices,  prodAmpDim,  3)]
						             * decayAmps[indicesToOffset<unsigned int>(decayAmpIndices, decayAmpDim, 3)];
				}
				l += norm(ampProdSum);
			}
		}
		l             += prodAmpFlat2;
		logLikelihood -= log(l);  // accumulate log likelihood
  }

	logLikelihoodSums[threadId] = logLikelihood;
	//logLikelihoodSums[threadId] = threadId;
}


