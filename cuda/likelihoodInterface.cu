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


#include <cuda.h>
#include <cutil_inline.h>

#include "nDimArrayUtils.hpp"
#include "reportingUtils.hpp"

#include "complex.cuh"
#include "likelihoodInterface.cuh"
#include "likelihoodKernel.cu"


using namespace std;
using namespace rpwa;
using namespace rpwa::cuda;


template<typename complexT> likelihoodInterface<complexT> likelihoodInterface<complexT>::_instance;

template<typename complexT> bool           likelihoodInterface<complexT>::_cudaInitialized    = false;
template<typename complexT> int            likelihoodInterface<complexT>::_nmbOfCudaDevices   = 0;
template<typename complexT> int            likelihoodInterface<complexT>::_cudaDeviceId       = -1;
template<typename complexT> cudaDeviceProp likelihoodInterface<complexT>::_cudaDeviceProp;
template<typename complexT> unsigned int   likelihoodInterface<complexT>::_nmbBlocks          = 0;
template<typename complexT> unsigned int   likelihoodInterface<complexT>::_nmbThreadsPerBlock = 0;
template<typename complexT> complexT*      likelihoodInterface<complexT>::_d_decayAmps        = 0;
template<typename complexT> unsigned int   likelihoodInterface<complexT>::_nmbEvents          = 0;
template<typename complexT> unsigned int   likelihoodInterface<complexT>::_nmbWavesRefl[2]    = {0, 0};
template<typename complexT> bool           likelihoodInterface<complexT>::_debug              = false;


template<typename complexT>
likelihoodInterface<complexT>::~likelihoodInterface()
{
	if (_d_decayAmps) {
		cutilSafeCall(cudaFree(_d_decayAmps));
		cutilSafeCall(cudaThreadExit());
	}
}


template<typename complexT>
unsigned int
likelihoodInterface<complexT>::totalDeviceMem()
{
	if (not _cudaInitialized) {
		printWarn << "cannot estimate total device memory. CUDA device not initialized." << endl;
		return 0;
	}
	return _cudaDeviceProp.totalGlobalMem;
}


template<typename complexT>
unsigned int
likelihoodInterface<complexT>::availableDeviceMem()
{
	// unfortunately the runtime API does not provide a call to get the
	// available memory on the device so we have to do it via the driver
	// API; sigh
	if (not _cudaInitialized) {
		printWarn << "cannot estimate available device memory. CUDA device not initialized." << endl;
		return 0;
	}
	if (cuInit(0) != CUDA_SUCCESS) {
		printWarn << "cuInit() failed." << endl;
		return 0;
	}
  CUdevice device;
  if (cuDeviceGet(&device, _cudaDeviceId) != CUDA_SUCCESS) {
	  printWarn << "cuDeviceGet() failed." << endl;
		return 0;
	}
  CUcontext context;
  if (cuCtxCreate(&context, 0, device) != CUDA_SUCCESS) {
	  printWarn << "cuCtxCreate() failed." << endl;
		return 0;
	}
	unsigned int free, total;
	if (cuMemGetInfo(&free, &total) != CUDA_SUCCESS) {
	  printWarn << "cuMemGetInfo() failed." << endl;
		return 0;
	}
  if (cuCtxDetach(context) != CUDA_SUCCESS)
	  printWarn << "cuCtxDetach() failed." << endl;
  return free;
}


template<typename complexT>
const struct cudaDeviceProp*
likelihoodInterface<complexT>::deviceProperties()
{
	if (not _cudaInitialized)
		return 0;
	return &_cudaDeviceProp;
}


template<typename complexT>
bool
likelihoodInterface<complexT>::init
(const complexT*    decayAmps,        // array of decay amplitudes; [iRefl][iWave][iEvt] or [iEvt][iRefl][iWave]
 const unsigned int nmbDecayAmps,     // number of elements in decay amplitude array
 const unsigned int nmbEvents,        // total number of events
 const unsigned int nmbWavesRefl[2],  // number of waves for each reflectivity
 const bool         reshuffleArray)   // if set devcay amplitude array is reshuffled from [iEvt][iRefl][iWave] to [iRefl][iWave][iEvt]
{
	if (not initCudaDevice()) {
		printWarn << "problems initializing CUDA device" << endl;
		return false;
	}
	if (_debug)
		printInfo << _instance;
	if (not loadDecayAmps(decayAmps, nmbDecayAmps, nmbEvents, nmbWavesRefl, reshuffleArray)) {
		printWarn << "problems loading decay amplitudes into CUDA device" << endl;
		return false;
	}
	return true;
}


template<typename complexT>
bool
likelihoodInterface<complexT>::initCudaDevice()
{
	_cudaInitialized = false;

  // get number of CUDA devices in system
  cutilSafeCall(cudaGetDeviceCount(&_nmbOfCudaDevices));
  if (_nmbOfCudaDevices == 0) {
    printWarn << "there are no CUDA devices in the system" << endl;
    return false;
  }
  printInfo << "found " << _nmbOfCudaDevices << " CUDA device(s)" << endl;

  // use most powerful GPU in system
  _cudaDeviceId = cutGetMaxGflopsDeviceId();
  cutilSafeCall(cudaGetDeviceProperties(&_cudaDeviceProp, _cudaDeviceId));
  printInfo << "using CUDA device[" << _cudaDeviceId << "]: '" << _cudaDeviceProp.name << "'" << endl;
  // fields for both major & minor fields are 9999, if device is not present
  if ((_cudaDeviceProp.major == 9999) and (_cudaDeviceProp.minor == 9999)) {
	  printWarn << "there is no CUDA device with ID " << _cudaDeviceId << endl;
	  return false;
  }
  cutilSafeCall(cudaSetDevice(_cudaDeviceId));

  // setup thread grid paramters
	_nmbBlocks          = _cudaDeviceProp.multiProcessorCount;
	//nmbThreadsPerBlock = _cudaDeviceProp.maxThreadsPerBlock;
	_nmbThreadsPerBlock = 448;
	printInfo << "using " << _nmbBlocks << " x " << _nmbThreadsPerBlock << " = "
	          << _nmbBlocks * _nmbThreadsPerBlock << " CUDA threads for likelihood calculation" << endl;

	_cudaInitialized = true;
	return true;
}


template<typename complexT>
bool
likelihoodInterface<complexT>::loadDecayAmps
(const complexT*    decayAmps,        // array of decay amplitudes; [iRefl][iWave][iEvt] or [iEvt][iRefl][iWave]
 const unsigned int nmbDecayAmps,     // number of elements in decay amplitude array
 const unsigned int nmbEvents,        // total number of events
 const unsigned int nmbWavesRefl[2],  // number of waves for each reflectivity
 const bool         reshuffleArray)   // if set devcay amplitude array is reshuffled from [iEvt][iRefl][iWave] to [iRefl][iWave][iEvt]
{
	if (not decayAmps) {
		printErr << "null pointer to decay amplitudes. aborting." << endl;
		throw;
	}
	if (not _cudaInitialized) {
		printWarn << "cannot load decay amplitudes. CUDA device is not initialized." << endl;
		return false;
	}

	_nmbEvents       = nmbEvents;
	_nmbWavesRefl[0] = nmbWavesRefl[0];
	_nmbWavesRefl[1] = nmbWavesRefl[1];

	complexT* h_decayAmps = 0;
	if (reshuffleArray) {
		// change memory layout of decay amplitudes from [iEvt][iRefl][iWave] to [iRefl][iWave][iEvt]
		h_decayAmps  = new complexT[nmbDecayAmps];
		const unsigned int decayAmpDimOld[3] = {_nmbEvents, 2, max(_nmbWavesRefl[0], _nmbWavesRefl[1])};
		const unsigned int decayAmpDimNew[3] = {2, max(_nmbWavesRefl[0], _nmbWavesRefl[1]), _nmbEvents};
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
				for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
					const unsigned int decayAmpIndicesOld[3] = {iEvt, iRefl, iWave};
					const unsigned int offsetOld             = indicesToOffset<unsigned int>
						                                           (decayAmpIndicesOld, decayAmpDimOld, 3);
					const unsigned int decayAmpIndicesNew[3] = {iRefl, iWave, iEvt};
					const unsigned int offsetNew             = indicesToOffset<unsigned int>
						                                           (decayAmpIndicesNew, decayAmpDimNew, 3);
					h_decayAmps[offsetNew].real() = real(decayAmps[offsetOld]);
					h_decayAmps[offsetNew].imag() = imag(decayAmps[offsetOld]);
				}
	}

	// copy decay amps to device memory
	const unsigned int size = nmbDecayAmps * sizeof(complexT);
	cutilSafeCall(cudaMalloc((void**)&_d_decayAmps, size));
	cutilSafeCall(cudaMemcpy(_d_decayAmps, (reshuffleArray) ? h_decayAmps : decayAmps,
	                         size, cudaMemcpyHostToDevice));
	printInfo << availableDeviceMem() / (1024. * 1024.) << " MiBytes left on CUDA device after loading "
	          << "decay amplitudes" << endl;
	return true;
}


template<typename complexT>
likelihoodInterface<complexT>::value_type
likelihoodInterface<complexT>::logLikelihood
(const complexT*    prodAmps,     // array of production amplitudes; [iRank][iRefl][iWave]
 const unsigned int nmbProdAmps,  // number of elements in production amplitude array
 const value_type   prodAmpFlat,  // (real) amplitude of flat wave
 const unsigned int rank)         // rank of spin-density matrix
{
	if (not prodAmps) {
		printErr << "null pointer to production amplitudes. aborting." << endl;
		throw;
	}

	// copy production amplitudes to device
	complexT* d_prodAmps;
	{
		const unsigned int size = nmbProdAmps * sizeof(complexT);
		cutilSafeCall(cudaMalloc((void**)&d_prodAmps, size));
		cutilSafeCall(cudaMemcpy(d_prodAmps, prodAmps, size, cudaMemcpyHostToDevice));
	}

	// first summation stage
	value_type*        d_logLikelihoods0;
	const unsigned int nmbElements0 = _nmbThreadsPerBlock * _nmbBlocks;
	cutilSafeCall(cudaMalloc((void**)&d_logLikelihoods0, sizeof(value_type) * nmbElements0));
	logLikelihoodKernel<complexT><<<_nmbBlocks, _nmbThreadsPerBlock>>>
		(d_prodAmps, prodAmpFlat * prodAmpFlat, _d_decayAmps, _nmbEvents, rank,
		 _nmbWavesRefl[0], _nmbWavesRefl[1], max(_nmbWavesRefl[0], _nmbWavesRefl[1]),
		 d_logLikelihoods0);
	//cutilSafeCall(cudaThreadSynchronize());
	// second summation stage
	value_type*        d_logLikelihoods1;
	const unsigned int nmbElements1 = _nmbThreadsPerBlock;
	cutilSafeCall(cudaMalloc((void**)&d_logLikelihoods1, sizeof(value_type) * nmbElements1));
	sumKernel<value_type><<<1, _nmbThreadsPerBlock>>>(d_logLikelihoods0, nmbElements0,
	                                                  d_logLikelihoods1);
	//cutilSafeCall(cudaThreadSynchronize());
	// third and last summation stage
	value_type* d_logLikelihoods2;
	cutilSafeCall(cudaMalloc((void**)&d_logLikelihoods2, sizeof(value_type)));
	sumKernel<value_type><<<1, 1>>>(d_logLikelihoods1, nmbElements1, d_logLikelihoods2);
	//cutilSafeCall(cudaThreadSynchronize());
	// copy result to host
	value_type logLikelihood;
	cutilSafeCall(cudaMemcpy(&logLikelihood, d_logLikelihoods2,
	                         sizeof(value_type), cudaMemcpyDeviceToHost));
	// cleanup
	cutilSafeCall(cudaFree(d_prodAmps       ));
	cutilSafeCall(cudaFree(d_logLikelihoods0));
	cutilSafeCall(cudaFree(d_logLikelihoods1));
	cutilSafeCall(cudaFree(d_logLikelihoods2));

	return logLikelihood;
}


//!!! this function does not work properly when -O3 is used in compilation
//!!!
//!!! the symptom is that from the _second_ invocation on, the function
//!!! calculates derivatives[0][0][0] wrongly. the reason seems to be
//!!! that the CUDA kernel is called with the iRank, iRefl, iWave from
//!!! the last kernel call of the previous function call;
//!!! testLikelihoodMockup.cc has a corresponding test case
//!!!
//!!! system specs: Ubuntu 10.04.1 LTS x86_64, CUDA 3.1, gcc version 4.3.4 (Ubuntu 4.3.4-10ubuntu1)
template<typename complexT>
likelihoodInterface<complexT>::value_type
likelihoodInterface<complexT>::logLikelihoodDeriv
(const complexT*    prodAmps,        // array of production amplitudes; [iRank][iRefl][iWave]
 const unsigned int nmbProdAmps,     // number of elements in production amplitude array
 const value_type   prodAmpFlat,     // (real) amplitude of flat wave
 const unsigned int rank,            // rank of spin-density matrix
 complexT*          derivatives,     // array of log likelihood derivatives; [iRank][iRefl][iWave]
 value_type&        derivativeFlat)  // log likelihood derivative of flat wave
{
	if (not prodAmps) {
		printErr << "null pointer to production amplitudes. aborting." << endl;
		throw;
	}

	// copy production amplitudes to device
	complexT* d_prodAmps;
	{
		const unsigned int size = nmbProdAmps * sizeof(complexT);
		cutilSafeCall(cudaMalloc((void**)&d_prodAmps, size));
		cutilSafeCall(cudaMemcpy(d_prodAmps, prodAmps, size, cudaMemcpyHostToDevice));
	}

	// first stage: precalculate derivative term and likelihoods for each event
	complexT*   d_derivTerms;
	value_type* d_likelihoods;
	{
		cutilSafeCall(cudaMalloc((void**)&d_derivTerms,  sizeof(complexT) * rank * 2 * _nmbEvents));
		cutilSafeCall(cudaMalloc((void**)&d_likelihoods, sizeof(value_type) * _nmbEvents));
		//!!! this logic does not handle well smaller number of events; in
		//!!! this case number of threads per block should be chosen much smaller
		const unsigned int nmbBlocks = _nmbEvents / _nmbThreadsPerBlock + 1;
		//printInfo << "nmbEvents = " << _nmbEvents << ", nmbBlocks = " << nmbBlocks << endl;
		logLikelihoodDerivFirstTermKernel<complexT><<<nmbBlocks, _nmbThreadsPerBlock>>>
			(d_prodAmps, prodAmpFlat * prodAmpFlat, _d_decayAmps, _nmbEvents, rank,
			 _nmbWavesRefl[0], _nmbWavesRefl[1], max(_nmbWavesRefl[0], _nmbWavesRefl[1]),
			 d_derivTerms, d_likelihoods);
		//cutilSafeCall(cudaThreadSynchronize());
	}

	// second stage: calculate derivative  sums for all production amplitudes
	complexT*          d_derivatives;
	const unsigned int derivativeDim[3] = {rank, 2, max(_nmbWavesRefl[0], _nmbWavesRefl[1])};
	const unsigned int nmbDerivElements = derivativeDim[0] * derivativeDim[1] * derivativeDim[2];
	cutilSafeCall(cudaMalloc((void**)&d_derivatives, sizeof(complexT) * nmbDerivElements));
	// cout << "!!!HERE1 d_derivatives = " << d_derivatives << endl;
	for (unsigned int iRank = 0; iRank < rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				// first summation stage
				complexT*          d_derivativeSums0;
				const unsigned int nmbElements0 = _nmbThreadsPerBlock * _nmbBlocks;
				cutilSafeCall(cudaMalloc((void**)&d_derivativeSums0, sizeof(complexT) * nmbElements0));
				// cout << "!!!HERE1 [" << iRank << "][" << iRefl << "][" << iWave << "]" << endl;
				logLikelihoodDerivKernel<complexT><<<_nmbBlocks, _nmbThreadsPerBlock>>>
					(_d_decayAmps, d_derivTerms, d_likelihoods, _nmbEvents, rank,
					 max(_nmbWavesRefl[0], _nmbWavesRefl[1]), iRank, iRefl, iWave,
					 d_derivativeSums0);
				//cutilSafeCall(cudaThreadSynchronize());
				// second summation stage
				complexT*          d_derivativeSums1;
				const unsigned int nmbElements1 = _nmbThreadsPerBlock;
				cutilSafeCall(cudaMalloc((void**)&d_derivativeSums1, sizeof(complexT) * nmbElements1));
				sumKernel<complexT><<<1, _nmbThreadsPerBlock>>>(d_derivativeSums0, nmbElements0,
				                                                d_derivativeSums1);
				//cutilSafeCall(cudaThreadSynchronize());
				cutilSafeCall(cudaFree(d_derivativeSums0));
				// third and last summation stage
				const unsigned int derivativeIndices[3] = {iRank, iRefl, iWave};
				unsigned int offset                     = indicesToOffset<unsigned int>(derivativeIndices,
					                                                                      derivativeDim, 3);
				sumToMemCellKernel<complexT><<<1, 1>>>(d_derivativeSums1, nmbElements1,
				                                       d_derivatives, offset);
				//cutilSafeCall(cudaThreadSynchronize());
				// cutilSafeCall(cudaFree(d_derivativeSums1));

				// cutilSafeCall(cudaMemcpy(derivatives, d_derivatives,
				//                          sizeof(complexT) * nmbDerivElements, cudaMemcpyDeviceToHost));
				// cout << "!!!HERE2 [" << iRank << "][" << iRefl << "][" << iWave << "]; [0] = "
				//      << "(" << maxPrecision(derivatives[0].real()) << ", "
				//      <<        maxPrecision(derivatives[0].imag()) << ")" << endl;
			}
	// copy result to host and cleanup
	cutilSafeCall(cudaMemcpy(derivatives, d_derivatives,
	                         sizeof(complexT) * nmbDerivElements, cudaMemcpyDeviceToHost));
	cutilSafeCall(cudaFree(d_derivatives));
	// cout << "!!!HERE2 derivatives[0] = "
	//      << "(" << maxPrecision(derivatives[0].real()) << ", "
	//      <<        maxPrecision(derivatives[0].imag()) << ")" << endl;
	
	// flat wave requires special treatment	
	{
		// first summation stage
		value_type*        d_derivativeSumsFlat0;
		const unsigned int nmbElements0 = _nmbThreadsPerBlock * _nmbBlocks;
		cutilSafeCall(cudaMalloc((void**)&d_derivativeSumsFlat0, sizeof(value_type) * nmbElements0));
		logLikelihoodDerivFlatKernel<complexT><<<_nmbBlocks, _nmbThreadsPerBlock>>>
			(prodAmpFlat, d_likelihoods, _nmbEvents, d_derivativeSumsFlat0);
		//cutilSafeCall(cudaThreadSynchronize());
		// second summation stage
		value_type*        d_derivativeSumsFlat1;
		const unsigned int nmbElements1 = _nmbThreadsPerBlock;
		cutilSafeCall(cudaMalloc((void**)&d_derivativeSumsFlat1, sizeof(value_type) * nmbElements1));
		sumKernel<value_type><<<1, _nmbThreadsPerBlock>>>(d_derivativeSumsFlat0, nmbElements0,
		                                                  d_derivativeSumsFlat1);
		//cutilSafeCall(cudaThreadSynchronize());
		// third and last summation stage
		value_type* d_derivativeSumsFlat2;
		cutilSafeCall(cudaMalloc((void**)&d_derivativeSumsFlat2, sizeof(value_type)));
		sumKernel<value_type><<<1, 1>>>(d_derivativeSumsFlat1, nmbElements1, d_derivativeSumsFlat2);
		//cutilSafeCall(cudaThreadSynchronize());
		// copy result to host
		cutilSafeCall(cudaMemcpy(&derivativeFlat, d_derivativeSumsFlat2,
		                         sizeof(value_type), cudaMemcpyDeviceToHost));
		// cleanup
		cutilSafeCall(cudaFree(d_derivativeSumsFlat0));
		cutilSafeCall(cudaFree(d_derivativeSumsFlat1));
		cutilSafeCall(cudaFree(d_derivativeSumsFlat2));
	}

	// cleanup
	cutilSafeCall(cudaFree(d_derivTerms ));
	cutilSafeCall(cudaFree(d_likelihoods));
	cutilSafeCall(cudaFree(d_prodAmps   ));

	return 0;
}


template<typename complexT>
ostream&
likelihoodInterface<complexT>::print(ostream& out)
{
  const unsigned int nGpuArchCoresPerSM[] = {1, 8, 32};  // from SDK/shared/inc/shrUtils.h

  if (not _cudaInitialized) {
	  printWarn << "CUDA device is not initialized." << endl;
	  return out;
  }
  
  // fields for both major & minor fields are 9999, if no CUDA capable devices are present
  if ((_cudaDeviceProp.major == 9999) and (_cudaDeviceProp.minor == 9999)) {
	  printWarn << "there is no CUDA device with ID " << _cudaDeviceId << endl;
	  return out;
  }
  out << "CUDA device[" << _cudaDeviceId << "]: '" << _cudaDeviceProp.name << "' properties:" << endl;
    
  // print info
  int driverVersion = 0;
  cutilSafeCall(cudaDriverGetVersion(&driverVersion));
  int runtimeVersion = 0;     
  cutilSafeCall(cudaRuntimeGetVersion(&runtimeVersion));
  out << "    driver version: .................................. " << driverVersion / 1000 << "." << driverVersion % 100 << endl
      << "    runtime version: ................................. " << runtimeVersion / 1000 << "." << runtimeVersion % 100 << endl
      << "    capability major revision number: ................ " << _cudaDeviceProp.major << endl
      << "    capability minor revision number: ................ " << _cudaDeviceProp.minor << endl
      << "    GPU clock frequency: ............................. " << _cudaDeviceProp.clockRate * 1e-6f << " GHz" << endl
      << "    number of multiprocessors: ....................... " << _cudaDeviceProp.multiProcessorCount << endl
      << "    number of cores: ................................. " << nGpuArchCoresPerSM[_cudaDeviceProp.major] * _cudaDeviceProp.multiProcessorCount << endl
      << "    warp size: ....................................... " << _cudaDeviceProp.warpSize << endl
      << "    maximum number of threads per block: ............. " << _cudaDeviceProp.maxThreadsPerBlock << endl
      << "    maximum block dimensions: ........................ " << _cudaDeviceProp.maxThreadsDim[0] << " x " << _cudaDeviceProp.maxThreadsDim[1]
                                                                   << " x " << _cudaDeviceProp.maxThreadsDim[2] << endl
      << "    maximum grid dimension ........................... " << _cudaDeviceProp.maxGridSize[0] << " x " << _cudaDeviceProp.maxGridSize[1]
                                                                   << " x " << _cudaDeviceProp.maxGridSize[2] << endl
      << "    total amount of global memory: ................... " << _cudaDeviceProp.totalGlobalMem / (1024. * 1024.) << " MiBytes" << endl
      << "    amount of available global memory: ............... " << availableDeviceMem() / (1024. * 1024.) << " MiBytes" << endl
      << "    total amount of constant memory: ................. " << _cudaDeviceProp.totalConstMem << " bytes" << endl 
      << "    total amount of shared memory per block: ......... " << _cudaDeviceProp.sharedMemPerBlock << " bytes" << endl
      << "    total number of registers available per block: ... " << _cudaDeviceProp.regsPerBlock << endl
      << "    maximum memory pitch: ............................ " << _cudaDeviceProp.memPitch << " bytes" << endl
      << "    texture alignment: ............................... " << _cudaDeviceProp.textureAlignment << " bytes" << endl
      << "    concurrent copy and execution: ................... " << ((_cudaDeviceProp.deviceOverlap)            ? "yes" : "no") << endl
      << "    run time limit on kernels: ....................... " << ((_cudaDeviceProp.kernelExecTimeoutEnabled) ? "yes" : "no") << endl
      << "    integrated: ...................................... " << ((_cudaDeviceProp.integrated)               ? "yes" : "no") << endl
      << "    support for host page-locked memory mapping: ..... " << ((_cudaDeviceProp.canMapHostMemory)         ? "yes" : "no") << endl
      << "    compute mode: .................................... ";
  if (_cudaDeviceProp.computeMode == cudaComputeModeDefault)
	  out << "default (multiple host threads can use this device simultaneously)";
  else if (_cudaDeviceProp.computeMode == cudaComputeModeExclusive)
	  out << "exclusive (only one host thread at a time can use this device)";
  else if (_cudaDeviceProp.computeMode == cudaComputeModeProhibited)
	  out << "prohibited (no host thread can use this device)";
  else
	  out << "unknown";
  out << endl;
  return out;
}


// explicit specializations
template class likelihoodInterface<cuda::complex<float > >;
template class likelihoodInterface<cuda::complex<double> >;
