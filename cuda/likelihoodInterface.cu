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
//      interface for CUDA likelihood functions
//
//
// Author List:
//      Philipp Meyer          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <helper_cuda.h>

#include "arrayUtils.hpp"
#include "reportingUtils.hpp"

#include "conversionUtils.hpp"
#include "complex.cuh"
#include "likelihoodInterface.cuh"
#include "likelihoodKernel.cu"


using namespace std;
using namespace boost;
using namespace rpwa;
using namespace rpwa::cuda;


template<typename complexT, typename T>
void
sumKernelCaller<complexT, T>::call(const T*           d_sumsPrev,
                                   const unsigned int nmbSumsPrev,
                                   T*                 d_sumsNext,
                                   const unsigned int nmbSumsNext,
                                   const bool         debug)
{
	unsigned int nmbBlocks, nmbThreadsPerBlock;
	if (debug)
		printInfo << "running sumKernel<T>" << endl;
	likelihoodInterface<complexT>::estimateOptimumKernelGrid(sumKernel<T>, nmbBlocks,
	                                                         nmbThreadsPerBlock, nmbSumsNext);
	likelihoodInterface<complexT>::startKernelTimer();
	sumKernel<T><<<nmbBlocks, nmbThreadsPerBlock>>>
		(d_sumsPrev, nmbSumsPrev, d_sumsNext, nmbSumsNext);
	likelihoodInterface<complexT>::stopKernelTimer();
}


template<typename complexT, typename T>
void
sumDerivativesKernelCaller<complexT, T>::call(const T*           d_sumsPrev,
                                              const unsigned int nmbSumsPrev,
                                              T*                 d_sumsNext,
                                              const unsigned int nmbSumsNext,
                                              const bool         debug)
{
	unsigned int nmbBlocks, nmbThreadsPerBlock;
	if (debug)
		printInfo << "running sumDerivativesKernel<T>" << endl;
	likelihoodInterface<complexT>::estimateOptimumKernelGrid(sumDerivativesKernel<T>, nmbBlocks,
	                                                         nmbThreadsPerBlock, nmbSumsNext);
	likelihoodInterface<complexT>::startKernelTimer();
	sumDerivativesKernel<T><<<nmbBlocks, nmbThreadsPerBlock>>>
		(d_sumsPrev, nmbSumsPrev, likelihoodInterface<complexT>::_rank,
		 likelihoodInterface<complexT>::_nmbWavesRefl[0],
		 likelihoodInterface<complexT>::_nmbWavesRefl[1],
		 likelihoodInterface<complexT>::_nmbWavesMax, d_sumsNext, nmbSumsNext);
	likelihoodInterface<complexT>::stopKernelTimer();
}


template<typename complexT> likelihoodInterface<complexT> likelihoodInterface<complexT>::_instance;

template<typename complexT> bool           likelihoodInterface<complexT>::_cudaInitialized  = false;
template<typename complexT> int            likelihoodInterface<complexT>::_nmbOfCudaDevices = 0;
template<typename complexT> int            likelihoodInterface<complexT>::_cudaDeviceId     = -1;
template<typename complexT> cudaDeviceProp likelihoodInterface<complexT>::_cudaDeviceProp;
template<typename complexT> complexT*      likelihoodInterface<complexT>::_d_decayAmps      = 0;
template<typename complexT> unsigned int   likelihoodInterface<complexT>::_nmbEvents        = 0;
template<typename complexT> unsigned int   likelihoodInterface<complexT>::_nmbWavesRefl[2]  = {0, 0};
template<typename complexT> unsigned int   likelihoodInterface<complexT>::_nmbWavesMax      = 0;
template<typename complexT> unsigned int   likelihoodInterface<complexT>::_rank             = 0;
template<typename complexT> timer          likelihoodInterface<complexT>::_timer;
template<typename complexT> double         likelihoodInterface<complexT>::_kernelTime       = 0;
template<typename complexT> bool           likelihoodInterface<complexT>::_debug            = false;


template<typename complexT>
likelihoodInterface<complexT>::~likelihoodInterface()
{
	closeCudaDevice();
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
  checkCudaErrors(cudaGetDeviceCount(&_nmbOfCudaDevices));
  if (_nmbOfCudaDevices == 0) {
    printWarn << "there are no CUDA devices in the system" << endl;
    return false;
  }
  printInfo << "found " << _nmbOfCudaDevices << " CUDA device(s)" << endl;

  // use most powerful GPU in system
  _cudaDeviceId = gpuGetMaxGflopsDeviceId();
  checkCudaErrors(cudaGetDeviceProperties(&_cudaDeviceProp, _cudaDeviceId));
  // fields for both major & minor fields are 9999, if device is not present
  if ((_cudaDeviceProp.major == 9999) and (_cudaDeviceProp.minor == 9999)) {
          printWarn << "there is no CUDA device with ID " << _cudaDeviceId << endl;
          return false;
  }
  checkCudaErrors(cudaSetDevice(_cudaDeviceId));
        _cudaInitialized = true;
  printInfo << "using CUDA device[" << _cudaDeviceId << "]: '" << _cudaDeviceProp.name << "' "
            << availableDeviceMem() << " bytes available memory" << endl;
	return true;
}


template<typename complexT>
void
likelihoodInterface<complexT>::closeCudaDevice()
{
  _nmbOfCudaDevices = 0;
  _nmbEvents        = 0;
  _nmbWavesRefl[0]  = 0;
  _nmbWavesRefl[1]  = 0;
  _nmbWavesMax      = 0;
  _rank             = 0;
        if (_d_decayAmps) {
                checkCudaErrors(cudaFree(_d_decayAmps));
                _d_decayAmps = 0;
        }
        checkCudaErrors(cudaThreadExit());
        if (_cudaInitialized)
                printInfo << "closing CUDA device[" << _cudaDeviceId << "]: '" << _cudaDeviceProp.name << "' "
                          << availableDeviceMem() << " bytes available memory" << endl;
	_cudaInitialized = false;
  _cudaDeviceId    = -1;
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
		printErr << "null pointer to decay amplitudes. Aborting..." << endl;
		throw;
	}
	if (not _cudaInitialized) {
		printWarn << "cannot load decay amplitudes. CUDA device is not initialized." << endl;
		return false;
	}

	_nmbEvents       = nmbEvents;
	_nmbWavesRefl[0] = nmbWavesRefl[0];
	_nmbWavesRefl[1] = nmbWavesRefl[1];
	_nmbWavesMax     = max(_nmbWavesRefl[0], _nmbWavesRefl[1]);

	complexT* h_decayAmps = 0;
	if (reshuffleArray) {
		// change memory layout of decay amplitudes from [iEvt][iRefl][iWave] to [iRefl][iWave][iEvt]
		h_decayAmps  = new complexT[nmbDecayAmps];
		const unsigned int decayAmpDimOld[3] = {_nmbEvents, 2, _nmbWavesMax};
		const unsigned int decayAmpDimNew[3] = {2, _nmbWavesMax, _nmbEvents};
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
        checkCudaErrors(cudaMalloc((void**)&_d_decayAmps, size));
        checkCudaErrors(cudaMemcpy(_d_decayAmps, (reshuffleArray) ? h_decayAmps : decayAmps,
                                 size, cudaMemcpyHostToDevice));
        printInfo << availableDeviceMem() / (1024. * 1024.) << " MiBytes left on CUDA device after loading "
                  << "decay amplitudes" << endl;
	return true;
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
	if (not _cudaInitialized) {
		printWarn << "cannot estimate available device memory. CUDA device not initialized." << endl;
                return 0;
        }
        size_t free, total;
        checkCudaErrors(cudaMemGetInfo(&free, &total));
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
template<typename kernelT>
ostream&
likelihoodInterface<complexT>::printKernelAttributes(ostream& out,
                                                     kernelT* kernel)
{
        if (not kernel)
                out << "null pointer to kernel" << endl;
        struct cudaFuncAttributes       kernelAttr;
        checkCudaErrors(cudaFuncGetAttributes(&kernelAttr, kernel));
        out << "kernel attributes:" << endl
            << "    compiled for architecture ........ " << kernelAttr.binaryVersion / 10. << endl
            << "    PTX virtual architecture ......... " << kernelAttr.binaryVersion / 10. << endl
	    << "    amount of constant memory ........ " << kernelAttr.constSizeBytes  << " bytes" << endl
	    << "    amount of local memory ........... " << kernelAttr.localSizeBytes  << " bytes" << endl
	    << "    amount of shared memory .......... " << kernelAttr.sharedSizeBytes << " bytes" << endl
	    << "    number of registers used ......... " << kernelAttr.numRegs            << endl
	    << "    maximum # of threads per block ... " << kernelAttr.maxThreadsPerBlock << endl;
	return out;
}


template<typename complexT>
template<typename kernelT>
void
likelihoodInterface<complexT>::estimateOptimumKernelGrid(kernelT*           kernel,
                                                         unsigned int&      nmbBlocks,
                                                         unsigned int&      nmbThreadsPerBlock,
                                                         const unsigned int minNmbThreads)
{
	nmbBlocks = nmbThreadsPerBlock = 0;
        if (not kernel)
                printWarn << "null pointer to kernel" << endl;
        struct cudaFuncAttributes       kernelAttr;
        checkCudaErrors(cudaFuncGetAttributes(&kernelAttr, kernel));
        unsigned int       maxNmbThreadsPerBlock = kernelAttr.maxThreadsPerBlock;
        const unsigned int warpSize              = _cudaDeviceProp.warpSize;
        // make sure maxNmbThreadsPerBlock is multiple of warp size
	if (maxNmbThreadsPerBlock > warpSize)
		maxNmbThreadsPerBlock = (maxNmbThreadsPerBlock / warpSize) * warpSize;
	const unsigned int maxNmbBlocks = _cudaDeviceProp.multiProcessorCount;
	if (minNmbThreads == 0) {
		// special case: use maximum occupancy grid
		nmbBlocks          = maxNmbBlocks;
		nmbThreadsPerBlock = maxNmbThreadsPerBlock;
	} else if (minNmbThreads / maxNmbThreadsPerBlock >= maxNmbBlocks) {
		// create needed number of blocks with max. number of threads per block
		nmbBlocks = minNmbThreads / maxNmbThreadsPerBlock;
		if (minNmbThreads % maxNmbThreadsPerBlock > 0)
			++nmbBlocks;
		nmbThreadsPerBlock = maxNmbThreadsPerBlock;
	} else {
		// there are not enough threads for maximum occupancy grid
		// run with lower number of threads per block and try to maximize number of blocks
		nmbThreadsPerBlock = minNmbThreads / maxNmbBlocks;
		if (minNmbThreads % maxNmbBlocks > 0)
			++nmbThreadsPerBlock;
		// make sure nmbThreadsPerBlock is multiple of warp size
		if (nmbThreadsPerBlock > warpSize)
			nmbThreadsPerBlock = (nmbThreadsPerBlock / warpSize) * warpSize;
		else
			nmbThreadsPerBlock = warpSize;
		nmbBlocks = minNmbThreads / nmbThreadsPerBlock;
		if (minNmbThreads % nmbThreadsPerBlock > 0)
			++nmbBlocks;
	}
	if (_debug) {
		printInfo;
		if (minNmbThreads != 0) {
			cout << "optimum grid for " << minNmbThreads << " requested threads is "
			     << "[" << nmbBlocks << " blocks x " << nmbThreadsPerBlock << " threads] = "
			     << nmbBlocks * nmbThreadsPerBlock << " threads. "
			     << "maximum occupancy grid for this kernel would be ";
		} else
			cout << "setting maximum occupancy grid for this kernel ";
		cout << "[" << maxNmbBlocks << " blocks x " << maxNmbThreadsPerBlock << " threads] = "
		     << maxNmbBlocks * maxNmbThreadsPerBlock << " threads" << endl;
	}
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
		printErr << "null pointer to production amplitudes. Aborting..." << endl;
		throw;
	}
	
	const unsigned int initialAvailMem = availableDeviceMem();

	// copy production amplitudes to device
        complexT* d_prodAmps = 0;
        {
                const unsigned int size = nmbProdAmps * sizeof(complexT);
                checkCudaErrors(cudaMalloc((void**)&d_prodAmps, size));
                checkCudaErrors(cudaMemcpy(d_prodAmps, prodAmps, size, cudaMemcpyHostToDevice));
        }

        // first summation stage
	value_type*  d_logLikelihoodSumsPrev = 0;
	unsigned int nmbSumsPrev             = 0;
	{
		// run kernel with maximum occupancy grid
		unsigned int nmbBlocks, nmbThreadsPerBlock;
		if (_debug)
                        printInfo << "running logLikelihoodKernel<complexT>" << endl;
                estimateOptimumKernelGrid(logLikelihoodKernel<complexT>, nmbBlocks, nmbThreadsPerBlock);
                nmbSumsPrev = nmbBlocks * nmbThreadsPerBlock;
                checkCudaErrors(cudaMalloc((void**)&d_logLikelihoodSumsPrev, sizeof(value_type) * nmbSumsPrev));
                startKernelTimer();
                logLikelihoodKernel<complexT><<<nmbBlocks, nmbThreadsPerBlock>>>
                        (d_prodAmps, prodAmpFlat * prodAmpFlat, _d_decayAmps, _nmbEvents, rank,
                         _nmbWavesRefl[0], _nmbWavesRefl[1], _nmbWavesMax, d_logLikelihoodSumsPrev);
                stopKernelTimer();
                checkCudaErrors(cudaFree(d_prodAmps));
        }
        // cascaded summation of log likelihoods
        cascadedKernelSum<sumKernelCaller<complexT, value_type> >
                (6, d_logLikelihoodSumsPrev, nmbSumsPrev);
        // copy result to host and cleanup
        value_type logLikelihood;
        checkCudaErrors(cudaMemcpy(&logLikelihood, d_logLikelihoodSumsPrev,
                                 sizeof(value_type), cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaFree(d_logLikelihoodSumsPrev));

        // check for memory leaks
        if (availableDeviceMem() != initialAvailMem)
		printWarn << "potential CUDA device memory leak: memory before function call = "
		          << initialAvailMem << " bytes, after = " << availableDeviceMem() << " bytes" << endl;
		
	return logLikelihood;
}


template<typename complexT>
void
likelihoodInterface<complexT>::logLikelihoodDeriv
(const complexT*    prodAmps,        // array of production amplitudes; [iRank][iRefl][iWave]
 const unsigned int nmbProdAmps,     // number of elements in production amplitude array
 const value_type   prodAmpFlat,     // (real) amplitude of flat wave
 const unsigned int rank,            // rank of spin-density matrix
 complexT*          derivatives,     // array of log likelihood derivatives; [iRank][iRefl][iWave]
 value_type&        derivativeFlat)  // log likelihood derivative of flat wave
{
	if (not prodAmps) {
		printErr << "null pointer to production amplitudes. Aborting..." << endl;
		throw;
	}
	if (not derivatives) {
		printErr << "null pointer to derivatives array. Aborting..." << endl;
		throw;
	}

	const unsigned int initialAvailMem = availableDeviceMem();

	// copy production amplitudes to device
        complexT* d_prodAmps = 0;
        {
                const unsigned int size = nmbProdAmps * sizeof(complexT);
                checkCudaErrors(cudaMalloc((void**)&d_prodAmps, size));
                checkCudaErrors(cudaMemcpy(d_prodAmps, prodAmps, size, cudaMemcpyHostToDevice));
        }

        // first stage: precalculate derivative term and likelihoods for each event
        complexT*   d_derivTerms  = 0;
        value_type* d_likelihoods = 0;
        {
                checkCudaErrors(cudaMalloc((void**)&d_derivTerms,  sizeof(complexT) * rank * 2 * _nmbEvents));
                checkCudaErrors(cudaMalloc((void**)&d_likelihoods, sizeof(value_type) * _nmbEvents));
                unsigned int nmbBlocks, nmbThreadsPerBlock;
                if (_debug)
                        printInfo << "running logLikelihoodDerivFirstTermKernel<complexT>" << endl;
		estimateOptimumKernelGrid(logLikelihoodDerivFirstTermKernel<complexT>, nmbBlocks,
		                          nmbThreadsPerBlock, _nmbEvents);
		startKernelTimer();
		logLikelihoodDerivFirstTermKernel<complexT><<<nmbBlocks, nmbThreadsPerBlock>>>
                        (d_prodAmps, prodAmpFlat * prodAmpFlat, _d_decayAmps, _nmbEvents, rank,
                         _nmbWavesRefl[0], _nmbWavesRefl[1], _nmbWavesMax, d_derivTerms, d_likelihoods);
                stopKernelTimer();
                checkCudaErrors(cudaFree(d_prodAmps));
        }

        // second stage: sum derivatives
	{
		// first summation stage
		complexT*          d_derivativeSumsPrev = 0;
		unsigned int       nmbSumsPrev          = 0;
		const unsigned int nmbDerivElements     = rank * 2 * _nmbWavesMax;
		{
			// run kernel with maximum occupancy grid
			unsigned int nmbBlocks, nmbThreadsPerBlock;
			if (_debug)
                                printInfo << "running logLikelihoodDerivKernel<complexT>" << endl;
                        estimateOptimumKernelGrid(logLikelihoodDerivKernel<complexT>, nmbBlocks, nmbThreadsPerBlock);
                        nmbSumsPrev = nmbBlocks * nmbThreadsPerBlock;
                        checkCudaErrors(cudaMalloc((void**)&d_derivativeSumsPrev,
                                                 sizeof(complexT) * nmbDerivElements * nmbSumsPrev));
                        startKernelTimer();
                        logLikelihoodDerivKernel<complexT><<<nmbBlocks, nmbThreadsPerBlock>>>
				(_d_decayAmps, d_derivTerms, d_likelihoods, _nmbEvents, rank, _nmbWavesRefl[0],
				 _nmbWavesRefl[1], _nmbWavesMax, d_derivativeSumsPrev);
			stopKernelTimer();
		}
		// cascaded summation of derivatives
		_rank = rank;  // needed for kernel caller
                cascadedKernelSum<sumDerivativesKernelCaller<complexT, complexT> >
                        (6, d_derivativeSumsPrev, nmbSumsPrev, nmbDerivElements);
                // copy result to host and cleanup
                checkCudaErrors(cudaMemcpy(derivatives, d_derivativeSumsPrev,
                                         sizeof(complexT) * nmbDerivElements, cudaMemcpyDeviceToHost));
                checkCudaErrors(cudaFree(d_derivativeSumsPrev));
        }
        
        // flat wave requires special treatment 
	{
		// first summation stage
		value_type*  d_derivativeSumsPrev = 0;
		// value_type*  d_derivativeSumsNext = 0;
		unsigned int nmbSumsPrev          = 0;
		{
			// run kernel with maximum occupancy grid
			unsigned int nmbBlocks, nmbThreadsPerBlock;
			if (_debug)
				printInfo << "running logLikelihoodDerivFlatKernel<complexT>" << endl;
                        estimateOptimumKernelGrid(logLikelihoodDerivFlatKernel<complexT>, nmbBlocks,
                                                  nmbThreadsPerBlock);
                        nmbSumsPrev = nmbBlocks * nmbThreadsPerBlock;
                        checkCudaErrors(cudaMalloc((void**)&d_derivativeSumsPrev, sizeof(value_type) * nmbSumsPrev));
                        startKernelTimer();
                        logLikelihoodDerivFlatKernel<complexT><<<nmbBlocks, nmbThreadsPerBlock>>>
                                (prodAmpFlat, d_likelihoods, _nmbEvents, d_derivativeSumsPrev);
			stopKernelTimer();
		}
		// cascaded summation of flat derivatives
                cascadedKernelSum<sumKernelCaller<complexT, value_type> >
                        (6, d_derivativeSumsPrev, nmbSumsPrev);
                // copy result to host and cleanup
                checkCudaErrors(cudaMemcpy(&derivativeFlat, d_derivativeSumsPrev,
                                         sizeof(value_type), cudaMemcpyDeviceToHost));
                checkCudaErrors(cudaFree(d_derivativeSumsPrev));
        }

        // cleanup
        checkCudaErrors(cudaFree(d_derivTerms ));
        checkCudaErrors(cudaFree(d_likelihoods));

        // check for memory leaks
        if (availableDeviceMem() != initialAvailMem)
		printWarn << "potential CUDA device memory leak: memory before function call = "
		          << initialAvailMem << " bytes, after = " << availableDeviceMem() << " bytes" << endl;
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
  checkCudaErrors(cudaDriverGetVersion(&driverVersion));
  int runtimeVersion = 0;     
  checkCudaErrors(cudaRuntimeGetVersion(&runtimeVersion));
  out << "    driver version: ........................................ "
      << driverVersion / 1000 << "." << driverVersion % 100 << endl
      << "    runtime version: ....................................... "
      << runtimeVersion / 1000 << "." << runtimeVersion % 100 << endl
      << "    capability major revision number: ...................... "
      << _cudaDeviceProp.major << endl
      << "    capability minor revision number: ...................... "
      << _cudaDeviceProp.minor << endl
      << "    GPU clock frequency: ................................... "
      << _cudaDeviceProp.clockRate * 1e-6f << " GHz" << endl
      << "    number of multiprocessors: ............................. "
      << _cudaDeviceProp.multiProcessorCount << endl
      << "    number of cores: ....................................... "
      << nGpuArchCoresPerSM[_cudaDeviceProp.major] * _cudaDeviceProp.multiProcessorCount << endl
      << "    warp size: ............................................. "
      << _cudaDeviceProp.warpSize << endl
      << "    maximum number of threads per block: ................... "
      << _cudaDeviceProp.maxThreadsPerBlock << endl
      << "    maximum block dimensions: .............................. "
      << _cudaDeviceProp.maxThreadsDim[0] << " x " << _cudaDeviceProp.maxThreadsDim[1]
      << " x " << _cudaDeviceProp.maxThreadsDim[2] << endl
      << "    maximum grid dimension ................................. "
      << _cudaDeviceProp.maxGridSize[0] << " x " << _cudaDeviceProp.maxGridSize[1]
      << " x " << _cudaDeviceProp.maxGridSize[2] << endl
      << "    total amount of global memory: ......................... "
      << _cudaDeviceProp.totalGlobalMem / (1024. * 1024.) << " MiBytes" << endl
      << "    amount of available global memory: ..................... "
      << availableDeviceMem() / (1024. * 1024.) << " MiBytes" << endl
      << "    total amount of constant memory: ....................... "
      << _cudaDeviceProp.totalConstMem << " bytes" << endl 
      << "    total amount of shared memory per block: ............... "
      << _cudaDeviceProp.sharedMemPerBlock << " bytes" << endl
      << "    total number of registers available per block: ......... "
      << _cudaDeviceProp.regsPerBlock << endl
      << "    maximum 1-dim texture size: ............................ "
      << _cudaDeviceProp.maxTexture1D << " pixels" << endl
      << "    maximum 2-dim texture size: ............................ "
      << _cudaDeviceProp.maxTexture2D[0] << " x " << _cudaDeviceProp.maxTexture2D[1]
      << " pixels" << endl
      << "    maximum 3-dim texture size: ............................ "
      << _cudaDeviceProp.maxTexture3D[0] << " x " << _cudaDeviceProp.maxTexture3D[1]
      << " x " << _cudaDeviceProp.maxTexture3D[2] << " pixels" << endl
      << "    alignment required for textures: ....................... "
      << _cudaDeviceProp.textureAlignment << " bytes" << endl
      << "    alignment required for surfaces: ....................... "
      << _cudaDeviceProp.surfaceAlignment << " bytes" << endl
      << "    maximum memory pitch: .................................. "
      << _cudaDeviceProp.memPitch << " bytes" << endl
      << "    support for concurrent execution of multiple kernels ... "
      << yesNo(_cudaDeviceProp.concurrentKernels) << endl
      << "    support for ECC memory protection: ..................... "
      << yesNo(_cudaDeviceProp.ECCEnabled) << endl
      << "    concurrent copy and execution: ......................... "
      << yesNo(_cudaDeviceProp.deviceOverlap) << endl
      << "    run time limit on kernels: ............................. "
      << yesNo(_cudaDeviceProp.kernelExecTimeoutEnabled) << endl
      << "    device is an integrated component ...................... "
      << yesNo(_cudaDeviceProp.integrated) << endl
      << "    support for host page-locked memory mapping: ........... "
      << yesNo(_cudaDeviceProp.canMapHostMemory) << endl
      << "    compute mode: .......................................... ";
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


template<typename complexT>
template<class kernelCaller>
void
likelihoodInterface<complexT>::cascadedKernelSum
(const unsigned int                  nmbOfSumsAtEachStage,
 typename kernelCaller::value_type*& d_sumsPrev,
 unsigned int                        nmbSumsPrev,
 const unsigned int                  sumElementSize)
{ 
	typename kernelCaller::value_type* d_sumsNext = 0;
	do {
		// allocate device array for new sums
                unsigned int nmbSumsNext = nmbSumsPrev / nmbOfSumsAtEachStage;
                if (nmbSumsNext <= nmbOfSumsAtEachStage)
                        nmbSumsNext = 1;
                checkCudaErrors
                        (cudaMalloc((void**)&d_sumsNext,
                                    sizeof(typename kernelCaller::value_type) * sumElementSize * nmbSumsNext));
                // run kernel
                kernelCaller::call(d_sumsPrev, nmbSumsPrev, d_sumsNext, nmbSumsNext, _debug);
                // cleanup
                checkCudaErrors(cudaFree(d_sumsPrev));
                // prepare for next iteration
                d_sumsPrev  = d_sumsNext;
                nmbSumsPrev = nmbSumsNext;
	} while (nmbSumsPrev > nmbOfSumsAtEachStage);
}


template<typename complexT>
void
likelihoodInterface<complexT>::stopKernelTimer()
{
        checkCudaErrors(cudaThreadSynchronize());
        _kernelTime += _timer.elapsed();
}


// explicit specializations
template class likelihoodInterface<cuda::complex<float > >;
template class likelihoodInterface<cuda::complex<double> >;
