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
//      interface for CUDA likelihood functions
//
//
// Author List:
//      Philipp Meyer          TUM            (original author)
//
//
//-------------------------------------------------------------------------


// 'Standard' includes:
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include <cutil_math.h>
#include <ctime>
#include <cuda_runtime.h>
#include <iostream>
#include "nDimArrayUtils.hpp"
#include "reportingUtils.hpp"
#include <cutil_inline.h>

// my includes:
#include "cudaLikelihoodInterface.h" // includes my Complex datatype
using namespace rpwa;
using namespace std;


// kernel:
#include "cudaLikelihoodKernel.cu"


template<typename complexT>
void
PrepareCUDA(const complexT*    decayAmps,
            const unsigned int decay_mem_size,
            complexT*&         d_decayAmps,
            unsigned int&      nmbBlocks,
            unsigned int&      nmbThreadsPerBlock)
{
	cudaSetDevice(cutGetMaxGflopsDeviceId()); // selects the best gpu
	int dev;
	cudaGetDevice(&dev);
	struct cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, dev);
	nmbBlocks          = prop.multiProcessorCount;
	nmbThreadsPerBlock = prop.maxThreadsPerBlock;
	nmbThreadsPerBlock = 448;
  
	cutilSafeCall(cudaMalloc((void**)&d_decayAmps, decay_mem_size));
	cutilSafeCall(cudaMemcpy(d_decayAmps, decayAmps, decay_mem_size, cudaMemcpyHostToDevice));
}


template<typename complexT, typename T>
T
runCudaLogLikelihoodKernels
(const complexT*    prodAmps,
 const unsigned int prod_mem_size,
 const T            prodAmpFlat,
 const complexT*    d_decayAmps,
 const unsigned int nmbEvents,
 const unsigned int rank,
 const unsigned int nmbWavesRefl[2],
 const unsigned int nmbBlocks,
 const unsigned int nmbThreadsPerBlock)
{
	// copy production amplitudes to device
	complexT* d_prodAmps;
	{
		cutilSafeCall(cudaMalloc((void**)&d_prodAmps, prod_mem_size));
		cutilSafeCall(cudaMemcpy(d_prodAmps, prodAmps, prod_mem_size, cudaMemcpyHostToDevice));
	}

	// first summation stage
	T*                 d_logLikelihoodSums1;
	const unsigned int nmbElements1 = nmbThreadsPerBlock * nmbBlocks;
	{
		const	unsigned int size = sizeof(T) * nmbElements1;
		cutilSafeCall(cudaMalloc((void**)&d_logLikelihoodSums1, size));

		sumLogLikelihoodKernel<complexT, T><<<nmbBlocks, nmbThreadsPerBlock>>>
			(d_prodAmps, prodAmpFlat * prodAmpFlat, d_decayAmps, nmbEvents, rank,
			 nmbWavesRefl[0], nmbWavesRefl[1], max(nmbWavesRefl[0], nmbWavesRefl[1]),
			 d_logLikelihoodSums1);
	}

	// second summation stage
	T*                 d_logLikelihoodSums2;
	const unsigned int nmbElements2 = nmbThreadsPerBlock;
	{ 
		const	unsigned int size = sizeof(T) * nmbElements2;
		cutilSafeCall(cudaMalloc((void**)&d_logLikelihoodSums2, size));

		sumKernel<T><<<1, nmbThreadsPerBlock>>>(d_logLikelihoodSums1, nmbElements1, d_logLikelihoodSums2);
	}
	
	// third and last summation stage
	T*                 d_logLikelihoodSums3;
	const unsigned int nmbElements3 = 1;
	{
		const	unsigned int size = sizeof(T) * nmbElements3;
		cutilSafeCall(cudaMalloc((void**)&d_logLikelihoodSums3, size));
		sumKernel<T><<<1, 1>>>(d_logLikelihoodSums2, nmbElements2, d_logLikelihoodSums3);
	}

	// copy result to host
	T logLikelihoodSum;
	cutilSafeCall(cudaMemcpy(&logLikelihoodSum, d_logLikelihoodSums3,
	                         sizeof(T), cudaMemcpyDeviceToHost));

	cutilSafeCall(cudaFree(d_prodAmps          ));
	cutilSafeCall(cudaFree(d_logLikelihoodSums1));
	cutilSafeCall(cudaFree(d_logLikelihoodSums2));
	cutilSafeCall(cudaFree(d_logLikelihoodSums3));

	return logLikelihoodSum;
}


double
rpwa::sumLogLikelihoodCuda
(const rpwa::complex<double>* prodAmps,
 const unsigned int           prod_mem_size,
 const double                 prodAmpFlat,
 const rpwa::complex<double>* d_decayAmps,
 const unsigned int           nmbEvents,
 const unsigned int           rank,
 const unsigned int           nmbWavesRefl[2],
 const unsigned int           nmbBlocks,
 const unsigned int           nmbThreadsPerBlock)
{
	return runCudaLogLikelihoodKernels<rpwa::complex<double>, double>
		(prodAmps, prod_mem_size, prodAmpFlat, d_decayAmps, nmbEvents,
		 rank, nmbWavesRefl, nmbBlocks, nmbThreadsPerBlock);
}


void
rpwa::initLogLikelihoodCuda
(const rpwa::complex<double>* decayAmps,
 const unsigned int           decay_mem_size,
 rpwa::complex<double>*&      d_decayAmps,
 unsigned int&                nmbBlocks,
 unsigned int&                nmbThreadsPerBlock)
{
	PrepareCUDA<rpwa::complex<double> >(decayAmps, decay_mem_size, d_decayAmps,
	                                    nmbBlocks, nmbThreadsPerBlock);
}
