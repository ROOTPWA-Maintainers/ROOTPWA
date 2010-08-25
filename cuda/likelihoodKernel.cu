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


namespace rpwa {

	namespace cuda {


		// generic cascadable kernel that computes sum of N values in an array
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


		// generic kernel that computes sum of N values in an array and
		// writes result to memory location defined by pointer and offset
		template<typename T>
		__global__
		void
		sumToMemCellKernel
		(const T*           d_valArray,    // array with values to sum up
		 const unsigned int nmbVals,       // total number of values all kernels have to process
		 T*                 d_sumArray,    // output array of partial sums
		 const unsigned int outputOffset)  // offset of memory location the result will be written to
		{
			const unsigned int threadId   = blockIdx.x * blockDim.x + threadIdx.x;
			const unsigned int nmbThreads = gridDim.x * blockDim.x;
  
			// #pragma unroll 16
			T sum = T(0);
			for (unsigned int i = threadId; i < nmbVals; i += nmbThreads)
				sum += d_valArray[i];

			d_sumArray[outputOffset] = sum;
		}


		// kernel that calculates real-data term of log likelihood sum for a number of events
		template<typename complexT>
		__global__
		void
		logLikelihoodKernel
		(const complexT*                     d_prodAmps,        // 3-dim. array of production amplitudes; [iRank][iRefl][iWave]
		 const typename complexT::value_type prodAmpFlat2,      // squared production amplitude of flat wave
		 const complexT*                     d_decayAmps,       // 3-dim. array of decay amplitudes; [iRefl][iWave][iEvt]
		 const unsigned int                  nmbEvents,         // total number of events all kernels have to process
		 const unsigned int                  rank,              // rank of spin-density matrix
		 const unsigned int                  nmbWavesReflNeg,   // number waves with negative reflectivity
		 const unsigned int                  nmbWavesReflPos,   // number waves with positive reflectivity
		 const unsigned int                  nmbWavesMax,       // maximum extent of iWave index for production and decay amplitude arrays
		 typename complexT::value_type*      d_logLikelihoods)  // output array of partial log likelihood sums with one entry for each kernel
		{
			const unsigned int threadId   = blockIdx.x * blockDim.x + threadIdx.x;
			const unsigned int nmbThreads = gridDim.x * blockDim.x;

			const unsigned int nmbWavesRefl[2] = {nmbWavesReflNeg, nmbWavesReflPos};
			// define extents of arrays
			const unsigned int prodAmpDim [3] = {rank, 2,           nmbWavesMax};
			const unsigned int decayAmpDim[3] = {2,    nmbWavesMax, nmbEvents};
			// loop over events and calculate real-data term of log likelihood
			typename complexT::value_type logLikelihood = 0;
			for (unsigned int iEvt = threadId; iEvt < nmbEvents; iEvt += nmbThreads) {
				typename complexT::value_type likelihood = 0;  // likelihood for this event
				for (unsigned int iRank = 0; iRank < rank; ++iRank) {  // incoherent sum over ranks
					for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
						complexT ampProdSum = complexT(0);  // amplitude sum for negative/positive reflectivity for this rank
						for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
							// compute likelihood term
							const unsigned int prodAmpIndices [3] = {iRank, iRefl, iWave};
							const unsigned int decayAmpIndices[3] = {iRefl, iWave, iEvt};
							ampProdSum +=
								  d_prodAmps [indicesToOffset<unsigned int>(prodAmpIndices,  prodAmpDim,  3)]
								* d_decayAmps[indicesToOffset<unsigned int>(decayAmpIndices, decayAmpDim, 3)];
						}
						likelihood += norm(ampProdSum);
					}
				}
				likelihood    += prodAmpFlat2;
				logLikelihood -= log(likelihood);  // accumulate log likelihood
			}
			// write result
			d_logLikelihoods[threadId] = logLikelihood;
		}


		// kernel that calculates the real-data term of the derivative of
		// the log likelihood that is independent of the indices of the
		// production amplitude w.r.t. which the derivative is taken
		template<typename complexT>
		__global__
		void
		logLikelihoodDerivFirstTermKernel
		(const complexT*                     d_prodAmps,       // 3-dim. array of production amplitudes; [iRank][iRefl][iWave]
		 const typename complexT::value_type prodAmpFlat2,     // squared production amplitude of flat wave
		 const complexT*                     d_decayAmps,      // 3-dim. array of decay amplitudes; [iRefl][iWave][iEvt]
		 const unsigned int                  nmbEvents,        // total number of events all kernels have to process
		 const unsigned int                  rank,             // rank of spin-density matrix
		 const unsigned int                  nmbWavesReflNeg,  // number waves with negative reflectivity
		 const unsigned int                  nmbWavesReflPos,  // number waves with positive reflectivity
		 const unsigned int                  nmbWavesMax,      // maximum extent of iWave index for production and decay amplitude arrays
		 complexT*                           d_derivTerms,     // 3-dim. output array of first derivative terms; [iRank][iRefl][iEvt]
		 typename complexT::value_type*      d_likelihoods)    // output array of likelihoods for each event
		{
			const unsigned int iEvt = blockIdx.x * blockDim.x + threadIdx.x;
			if (iEvt < nmbEvents) {
				const unsigned int nmbWavesRefl[2] = {nmbWavesReflNeg, nmbWavesReflPos};
				// define extents of arrays
				const unsigned int            prodAmpDim  [3] = {rank, 2,           nmbWavesMax};
				const unsigned int            decayAmpDim [3] = {2,    nmbWavesMax, nmbEvents};
				const unsigned int            derivTermDim[3] = {rank, 2,           nmbEvents};
				typename complexT::value_type likelihood      = 0;  // likelihood for this event
				for (unsigned int iRank = 0; iRank < rank; ++iRank) {  // incoherent sum over ranks
					for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
						complexT ampProdSum = complexT(0);  // amplitude sum for negative/positive reflectivity for this rank
						for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
							const unsigned int prodAmpIndices [3] = {iRank, iRefl, iWave};
							const unsigned int decayAmpIndices[3] = {iRefl, iWave, iEvt};
							ampProdSum +=
								  d_prodAmps [indicesToOffset<unsigned int>(prodAmpIndices,  prodAmpDim,  3)]
								* d_decayAmps[indicesToOffset<unsigned int>(decayAmpIndices, decayAmpDim, 3)];
						}
						likelihood += norm(ampProdSum);
						// write derivative term
						const unsigned int derivTermIndices[3] = {iRank, iRefl, iEvt};
						d_derivTerms[indicesToOffset<unsigned int>(derivTermIndices, derivTermDim, 3)]
							= ampProdSum;
					}
				}
				likelihood += prodAmpFlat2;
				// write likelihood
				d_likelihoods[iEvt] = likelihood;
			}
		}


		// kernel that operates on the output of logLikelihoodDerivFirstTermKernel
		// and calculates the real-data derivative sum of the log likelihood for a
		// particular production amplitude
		template<typename complexT>
		__global__
		void
		logLikelihoodDerivKernel
		(const complexT*                      d_decayAmps,       // 3-dim. array of decay amplitudes; [iRefl][iWave][iEvt]
		 const complexT*                      d_derivTerms,      // precalculated 3-dim. array of first derivative terms; [iRank][iRefl][iEvt]
		 const typename complexT::value_type* d_likelihoods,     // precalculated array of likelihoods for each event
		 const unsigned int                   nmbEvents,         // total number of events all kernels have to process
		 const unsigned int                   rank,              // rank of spin-density matrix
		 const unsigned int                   nmbWavesMax,       // maximum extent of iWave index for production and decay amplitude arrays
		 const unsigned int                   iRank,             // rank of production amplitude w.r.t. which derivative is taken
		 const unsigned int                   iRefl,             // reflectivity of production amplitude w.r.t. which derivative is taken
		 const unsigned int                   iWave,             // wave index of production amplitude w.r.t. which derivative is taken
		 complexT*                            d_derivativeSums)  // output array of partial derivative sums with one entry for each kernel
		{
			const unsigned int threadId   = blockIdx.x * blockDim.x + threadIdx.x;
			const unsigned int nmbThreads = gridDim.x * blockDim.x;

			// define extents of arrays
			const unsigned int decayAmpDim [3] = {2,    nmbWavesMax, nmbEvents};
			const unsigned int derivTermDim[3] = {rank, 2,           nmbEvents};
			// loop over events and calculate real-data term of derivative of log likelihood
			complexT derivativeSum = complexT(0);
			for (unsigned int iEvt = threadId; iEvt < nmbEvents; iEvt += nmbThreads) {
				// multiply derivative term 1 with with complex conjugate of
				// decay amplitude of the wave with the derivative wave index
				const unsigned int decayAmpIndices [3] = {iRefl, iWave, iEvt};
				const unsigned int derivTermIndices[3] = {iRank, iRefl, iEvt};
				const complexT     derivative          =
					  d_derivTerms[indicesToOffset<unsigned int>(derivTermIndices, derivTermDim, 3)]
					* conj(d_decayAmps[indicesToOffset<unsigned int>(decayAmpIndices, decayAmpDim, 3)]);
				// apply factor from derivative of log
				derivativeSum -= (2. / d_likelihoods[iEvt]) * derivative;
			}
			// write result
			d_derivativeSums[threadId] = derivativeSum;
		}


		template<typename complexT>
		__global__
		void
		logLikelihoodDerivKernelXXX
		(const complexT*                      d_decayAmps,       // 3-dim. array of decay amplitudes; [iRefl][iWave][iEvt]
		 const complexT*                      d_derivTerms,      // precalculated 3-dim. array of first derivative terms; [iRank][iRefl][iEvt]
		 const typename complexT::value_type* d_likelihoods,     // precalculated array of likelihoods for each event
		 const unsigned int                   nmbEvents,         // total number of events all kernels have to process
		 const unsigned int                   rank,              // rank of spin-density matrix
		 const unsigned int                   nmbWavesReflNeg,   // number waves with negative reflectivity
		 const unsigned int                   nmbWavesReflPos,   // number waves with positive reflectivity
		 const unsigned int                   nmbWavesMax,       // maximum extent of iWave index for production and decay amplitude arrays
		 complexT*                            d_derivativeSums)  // output array of partial derivative sums with one entry for each kernel
		{
			const unsigned int threadId   = blockIdx.x * blockDim.x + threadIdx.x;
			const unsigned int nmbThreads = gridDim.x * blockDim.x;

			const unsigned int nmbWavesRefl[2] = {nmbWavesReflNeg, nmbWavesReflPos};
			// define extents of arrays
			const unsigned int decayAmpDim [3] = {2,    nmbWavesMax, nmbEvents};
			const unsigned int derivTermDim[3] = {rank, 2,           nmbEvents};
			const unsigned int derivSumDim [4] = {rank, 2,           nmbWavesMax, nmbThreads};
			// loop over events and calculate real-data term of derivative of log likelihood
			for (unsigned int iRank = 0; iRank < rank; ++iRank)
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
					for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {
						complexT derivativeSum = complexT(0);
						for (unsigned int iEvt = threadId; iEvt < nmbEvents; iEvt += nmbThreads) {
							// multiply derivative term 1 with with complex conjugate of
							// decay amplitude of the wave with the derivative wave index
							const unsigned int decayAmpIndices [3] = {iRefl, iWave, iEvt};
							const unsigned int derivTermIndices[3] = {iRank, iRefl, iEvt};
							const complexT     derivative          =
								  d_derivTerms[indicesToOffset<unsigned int>(derivTermIndices, derivTermDim, 3)]
								* conj(d_decayAmps[indicesToOffset<unsigned int>(decayAmpIndices, decayAmpDim, 3)]);
							// apply factor from derivative of log
							derivativeSum -= (2. / d_likelihoods[iEvt]) * derivative;
						}
						// write result
						const unsigned int derivSumIndices[4] = {iRank, iRefl, iWave, threadId};
						d_derivativeSums[indicesToOffset<unsigned int>(derivSumIndices, derivSumDim, 4)]
							= derivativeSum;
					}
		}


		template<typename T>
		__global__
		void
		sumKernelXXX
		(const T*           d_derivatives,     // array with values to sum up
		 const unsigned int nmbDerivatives,    // total number of values all kernels have to process
		 const unsigned int rank,              // rank of spin-density matrix
		 const unsigned int nmbWavesReflNeg,   // number waves with negative reflectivity
		 const unsigned int nmbWavesReflPos,   // number waves with positive reflectivity
		 const unsigned int nmbWavesMax,       // maximum extent of iWave index for production and decay amplitude arrays
		 T*                 d_derivativeSums,  // output array of partial sums with one entry for each kernel
		 const unsigned int nmbSums)           // number of elements in output array
		{
			const unsigned int threadId = blockIdx.x * blockDim.x + threadIdx.x;
			if (threadId < nmbSums) {
				const unsigned int nmbWavesRefl[2] = {nmbWavesReflNeg, nmbWavesReflPos};
				// define extents of derivative array
				const unsigned int derivDim   [4] = {rank, 2, nmbWavesMax, nmbDerivatives};
				const unsigned int derivSumDim[4] = {rank, 2, nmbWavesMax, nmbSums};
				for (unsigned int iRank = 0; iRank < rank; ++iRank)
					for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
						for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {
							T derivativeSum = T(0);
							for (unsigned int i = threadId; i < nmbDerivatives; i += nmbSums) {
								const unsigned int derivIndices[4] = {iRank, iRefl, iWave, i};
								derivativeSum
									+= d_derivatives[indicesToOffset<unsigned int>(derivIndices, derivDim, 4)];
							}
							const unsigned int derivSumIndices[4] = {iRank, iRefl, iWave, threadId};
							d_derivativeSums[indicesToOffset<unsigned int>(derivSumIndices, derivSumDim, 4)]
								= derivativeSum;
						}
			}
		}


		// kernel that operates on the output of logLikelihoodDerivFirstTermKernel
		// and calculates the real-data derivative sum of the log likelihood for the flat wave
		template<typename complexT>
		__global__
		void
		logLikelihoodDerivFlatKernel
		(const typename complexT::value_type  prodAmpFlat,           // (real) production amplitude of flat wave
		 const typename complexT::value_type* d_likelihoods,         // precalculated array of likelihoods for each event
		 const unsigned int                   nmbEvents,             // total number of events all kernels have to process
		 typename complexT::value_type*       d_derivativeFlatSums)  // output array of partial derivative sums with one entry for each kernel
		{
			const unsigned int threadId   = blockIdx.x * blockDim.x + threadIdx.x;
			const unsigned int nmbThreads = gridDim.x * blockDim.x;

			// loop over events and calculate real-data term of derivative of log likelihood
			typename complexT::value_type derivativeFlatSum = 0;
			for (unsigned int iEvt = threadId; iEvt < nmbEvents; iEvt += nmbThreads)
				derivativeFlatSum -= (2. / d_likelihoods[iEvt]) * prodAmpFlat;
			// write result
			d_derivativeFlatSums[threadId] = derivativeFlatSum;
		}


	}  // namespace cuda

}  // namespace rpwa
