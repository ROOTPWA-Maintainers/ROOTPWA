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
//      Test CUDA program
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
#include "cuda_helper.h" // includes my Complex datatype
using namespace rpwa;
using namespace std;


// kernel:
#include "pwa_cuda_kernel.cu"


template<typename complexT, typename scalarT>
void
PrepareCUDA(complexT*          decayAmps,
            const unsigned int decay_mem_size,
            complexT**         d_decayAmps,
            unsigned int&      num_threads,
            unsigned int&      num_blocks)
{
	cudaSetDevice(cutGetMaxGflopsDeviceId()); // selects the best gpu
	int dev;
	cudaGetDevice(&dev);
	struct cudaDeviceProp prop;
	cudaGetDeviceProperties (&prop, dev);
	num_threads = prop.maxThreadsPerBlock;
	num_blocks  = prop.multiProcessorCount;
	// num_blocks  = 1;
	num_threads = 448;
  
	cudaMalloc((void**) d_decayAmps, decay_mem_size);
	cudaMemcpy(*d_decayAmps, decayAmps, decay_mem_size, cudaMemcpyHostToDevice);
}


template<typename complexT, typename scalarT>
scalarT
SumArrayCUDA(complexT*          prodAmps,
             const unsigned int prod_mem_size,
             const scalarT      prodAmpFlat,
             const unsigned int nmbEvents,
             const unsigned int rank,
             const unsigned int nmbWavesRefl[2],
             complexT*          d_decayAmps,
             unsigned int       num_threads,
             unsigned int       num_blocks)
{
	complexT* d_prodAmps;
	cudaMalloc((void**) &d_prodAmps, prod_mem_size);
	cudaMemcpy(d_prodAmps, prodAmps, prod_mem_size, cudaMemcpyHostToDevice);

	scalarT* d_preresult;
	int result_mem_size = sizeof(scalarT) * num_threads * num_blocks;
	cudaMalloc((void**) &d_preresult, result_mem_size);

	SumCUDAPseudoArray<complexT,scalarT><<<num_blocks, num_threads>>>
		(d_decayAmps, d_prodAmps, prodAmpFlat, nmbEvents, rank,
		 nmbWavesRefl[0], nmbWavesRefl[1], max(nmbWavesRefl[0], nmbWavesRefl[1]), d_preresult);

	scalarT* d_preresult2;
	cudaMalloc((void**) &d_preresult2, sizeof(scalarT) * num_threads);
	SumCUDA<scalarT><<<1,num_threads>>>(d_preresult, num_threads * num_blocks, d_preresult2);
	
	scalarT* d_result;
	cudaMalloc((void**) &d_result, sizeof(scalarT));
	SumCUDA<scalarT><<<1,1>>>(d_preresult2, num_threads, d_result);
  
	scalarT* result = (scalarT*) malloc(sizeof(scalarT));
	cudaMemcpy(result, d_result, sizeof(scalarT), cudaMemcpyDeviceToHost);

	cudaFree(d_result);
	//	cudaFree(d_decayAmps);
	cudaFree(d_prodAmps);
	cudaFree(d_preresult);
	cudaFree(d_preresult2);

	return *result;
}


double
rpwa::SumArrayCUDA2(complex<double>*   prodAmps,
                    const unsigned int prod_mem_size,
                    const double       prodAmpFlat,
                    const unsigned int nmbEvents,
                    const unsigned int rank,
                    const unsigned int nmbWavesRefl[2],
                    complex<double>*   d_decayAmps,
                    unsigned int       num_threads,
                    unsigned int       num_blocks)
{
	return SumArrayCUDA<complex<double>,double>
		(prodAmps, prod_mem_size, prodAmpFlat, nmbEvents, rank, nmbWavesRefl, d_decayAmps,
		 num_threads, num_blocks);
}


void
rpwa::PrepareCUDA2(complex<double>*   decayAmps,
                   const unsigned int decay_mem_size,
                   complex<double>**  d_decayAmps,
                   unsigned int&      num_threads,
                   unsigned int&      num_blocks)
{
	PrepareCUDA<complex<double>,double>(decayAmps, decay_mem_size, d_decayAmps,
	                                    num_threads, num_blocks);
}
