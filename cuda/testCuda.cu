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
//      test some CUDA stuff
//      http://blog.icare3d.org/2010/06/cuda-dynamic-template-parameters-22.html
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>

#include <cutil_inline.h>

#include "complex.hpp"
#include "textureReader.hpp"


using namespace std;


bool
printCudaDeviceInfo(const int deviceId)
{
  const unsigned int nGpuArchCoresPerSM[] = {1, 8, 32};  // from SDK/shared/inc/shrUtils.h

  cudaDeviceProp deviceProp;
  cutilSafeCall(cudaGetDeviceProperties(&deviceProp, deviceId));
  if (deviceId == 0) {
    // fields for both major & minor fields are 9999, if no CUDA capable devices are present
    if ((deviceProp.major == 9999) and (deviceProp.minor == 9999)) {
      cout << "there is no CUDA device" << endl;
      return false;
    }
  }
  cout << "CUDA device[" << deviceId << "]: '" << deviceProp.name << "'" << endl;
    
  // print info
  int driverVersion = 0;
  cutilSafeCall(cudaDriverGetVersion(&driverVersion));
  int runtimeVersion = 0;     
  cutilSafeCall(cudaRuntimeGetVersion(&runtimeVersion));
  cout << "    driver version: .................................. " << driverVersion / 1000 << "." << driverVersion % 100 << endl
       << "    runtime version: ................................. " << runtimeVersion / 1000 << "." << runtimeVersion % 100 << endl
       << "    capability major revision number: ................ " << deviceProp.major << endl
       << "    capability minor revision number: ................ " << deviceProp.minor << endl
       << "    GPU clock frequency: ............................. " << deviceProp.clockRate * 1e-6f << " GHz" << endl
       << "    number of multiprocessors: ....................... " << deviceProp.multiProcessorCount << endl
       << "    number of cores: ................................. " << nGpuArchCoresPerSM[deviceProp.major] * deviceProp.multiProcessorCount << endl
       << "    warp size: ....................................... " << deviceProp.warpSize << endl
       << "    maximum number of threads per block: ............. " << deviceProp.maxThreadsPerBlock << endl
       << "    maximum block dimensions: ........................ " << deviceProp.maxThreadsDim[0] << " x " << deviceProp.maxThreadsDim[1]
                                                                    << " x " << deviceProp.maxThreadsDim[2] << endl
       << "    maximum grid dimension ........................... " << deviceProp.maxGridSize[0] << " x " << deviceProp.maxGridSize[1]
                                                                    << " x " << deviceProp.maxGridSize[2] << endl
       << "    total amount of global memory: ................... " << deviceProp.totalGlobalMem / (1024. * 1024. * 1024.) << " GiBytes" << endl
       << "    total amount of constant memory: ................. " << deviceProp.totalConstMem << " bytes" << endl 
       << "    total amount of shared memory per block: ......... " << deviceProp.sharedMemPerBlock << " bytes" << endl
       << "    total number of registers available per block: ... " << deviceProp.regsPerBlock << endl
       << "    maximum memory pitch: ............................ " << deviceProp.memPitch << " bytes" << endl
       << "    texture alignment: ............................... " << deviceProp.textureAlignment << " bytes" << endl
       << "    concurrent copy and execution: ................... " << ((deviceProp.deviceOverlap)            ? "yes" : "no") << endl
       << "    run time limit on kernels: ....................... " << ((deviceProp.kernelExecTimeoutEnabled) ? "yes" : "no") << endl
       << "    integrated: ...................................... " << ((deviceProp.integrated)               ? "yes" : "no") << endl
       << "    support for host page-locked memory mapping: ..... " << ((deviceProp.canMapHostMemory)         ? "yes" : "no") << endl
       << "    compute mode: .................................... " << ((deviceProp.computeMode == cudaComputeModeDefault) ?
									"default (multiple host threads can use this device simultaneously)" :
									(deviceProp.computeMode == cudaComputeModeExclusive) ?
									"exclusive (only one host thread at a time can use this device)" :
									(deviceProp.computeMode == cudaComputeModeProhibited) ?
									"prohibited (no host thread can use this device)" :
									"unknown") << endl;
  return true;
}


template<typename T, typename textureReaderT>
__global__ void
sumTextureMemKernel(T*                 outData,  // pointer to device output data in global memory
		    const unsigned int nmbElementsPerThread)
{
  const unsigned int threadId   = blockDim.x * blockIdx.x + threadIdx.x;
  const unsigned int nmbThreads = gridDim.x * blockDim.x;
  T                  sum        = 0;
#pragma unroll 1
  for (int i = 0; i < nmbElementsPerThread; ++i)
    sum += textureReaderT::fetch((i * nmbThreads ) + threadId);
  outData[threadId] = sum;
}


template<typename T>
__global__ void
sumGlobalMemKernel(const T*           inData,   // pointer to device input data in global memory
		   T*                 outData,  // pointer to device output data in global memory
		   const unsigned int nmbElementsPerThread)
{
  const unsigned int threadId   = blockDim.x * blockIdx.x + threadIdx.x;
  const unsigned int nmbThreads = gridDim.x * blockDim.x;
  T                  sum        = 0;
#pragma unroll 1
  for (int i = 0; i < nmbElementsPerThread; ++i)
    sum += inData[(i * nmbThreads ) + threadId];
  outData[threadId] = sum;
}


template<typename T>
bool
verifySumKernel(const T*           inData,
		const T*           outData,  // output of GPU kernel
		const unsigned int nmbBlocks,
		const unsigned int nmbThreadsPerBlock,
		const unsigned int nmbElementsPerThread)
{
  if (not inData or not outData) {
    cout << "null pointer for data" << endl;
    return false;
  }
  const unsigned int nmbThreads = nmbBlocks * nmbThreadsPerBlock;
  T                  data[nmbThreads];
  bool               success = true;
  for (unsigned int threadId = 0; threadId < nmbThreads; ++threadId) {
    T sum = 0;
    for (int i = 0; i < nmbElementsPerThread; ++i)
      sum += inData[(i * nmbThreads ) + threadId];
    data[threadId] = sum;
    if (data[threadId] != outData[threadId]) {
      cout << "(CPU[" << threadId << "] = " << data[threadId]    << ") != "
	   << "(GPU[" << threadId << "] = " << outData[threadId] << "); "
	   << "delta = " << data[threadId] - outData[threadId] << endl;
      success = false;
    }
  }
  return success;
}


template<typename T, typename textureReaderT>
struct sumTextureMemKernelCaller {
  
  typedef T                                     value_type;
  typedef typename textureReaderT::texture_type texture_type;

  static void call(const unsigned int nmbBlocks,
		   const unsigned int nmbThreadsPerBlock,
		   const T*,
		   T*                 deviceOutData,
		   const unsigned int nmbElementsPerThread)
  {
    sumTextureMemKernel<T, textureReaderT><<< nmbBlocks, nmbThreadsPerBlock >>>
      (deviceOutData, nmbElementsPerThread);
  }

  static bool verify(const T*           hostInData,
		     const T*           hostOutData,  // output of GPU kernel
		     const unsigned int nmbBlocks,
		     const unsigned int nmbThreadsPerBlock,
		     const unsigned int nmbElementsPerThread)
  {
    return verifySumKernel(hostInData, hostOutData, nmbBlocks,
			   nmbThreadsPerBlock, nmbElementsPerThread);
  }

};


template<typename T>
struct sumGlobalMemKernelCaller {
  
  typedef T value_type;

  static void call(const unsigned int nmbBlocks,
		   const unsigned int nmbThreadsPerBlock,
		   const T*           deviceInData,
		   T*                 deviceOutData,
		   const unsigned int nmbElementsPerThread)
  {
    sumGlobalMemKernel<T><<< nmbBlocks, nmbThreadsPerBlock >>>
      (deviceInData, deviceOutData, nmbElementsPerThread);
  }

  static bool verify(const T*           hostInData,
		     const T*           hostOutData,  // output of GPU kernel
		     const unsigned int nmbBlocks,
		     const unsigned int nmbThreadsPerBlock,
		     const unsigned int nmbElementsPerThread)
  {
    return verifySumKernel(hostInData, hostOutData, nmbBlocks,
			   nmbThreadsPerBlock, nmbElementsPerThread);
  }

};


template<typename T, typename textureReaderT, typename kernelCallerT>
void runKernel(const unsigned int nmbBlocks,
	       const unsigned int nmbThreadsPerBlock,
	       const unsigned int nmbIterations)
{
  // create maximum sized texture (512 MB)
  unsigned int       nmbElements          = ((unsigned int)1 << 29) / sizeof(T);
  const unsigned int nmbThreads           = nmbBlocks * nmbThreadsPerBlock;
  const unsigned int nmbElementsPerThread = nmbElements / nmbThreads;
  nmbElements = nmbElementsPerThread * nmbThreads;
  const unsigned int dataSizeIn  = nmbElements * sizeof(T);
  const unsigned int dataSizeOut = nmbThreads  * sizeof(T);

  // create and initalize host arrays
  cout << "allocating " << dataSizeIn / (1024. * 1024.) << " MiBytes in global memory "
       << "(" << nmbElementsPerThread << " data elements per thread)" << endl;
  T* hostInData  = (T*) malloc(dataSizeIn );
  T* hostOutData = (T*) malloc(dataSizeOut);
  srand (time(NULL));
  for (unsigned int i = 0; i < nmbElements; ++i)
    hostInData[i] = (T)rand();

  // create device arrays and copy host data to device
  T* deviceInData;
  T* deviceOutData;
  cutilSafeCall(cudaMalloc((void**) &deviceInData,  dataSizeIn ));
  cutilSafeCall(cudaMalloc((void**) &deviceOutData, dataSizeOut));
  cutilSafeCall(cudaMemcpy(deviceInData, hostInData, dataSizeIn, cudaMemcpyHostToDevice));

  // bind texture
  textureReaderT::bindTexture(deviceInData, dataSizeIn);

  // dry-run kernel first to avoid any setup and caching effects
  kernelCallerT::call(nmbBlocks, nmbThreadsPerBlock, deviceInData,
		      deviceOutData, nmbElementsPerThread);

  // setup and start timer
  cudaEvent_t start, end;
  cutilSafeCall(cudaEventCreate(&start));
  cutilSafeCall(cudaEventCreate(&end  ));
  cutilSafeCall(cudaEventRecord(start, 0));

  // run kernel
  for (unsigned int iteration = 0; iteration < nmbIterations; ++iteration)
    kernelCallerT::call(nmbBlocks, nmbThreadsPerBlock, deviceInData,
			deviceOutData, nmbElementsPerThread);

  // stop timer
  cutilSafeCall(cudaEventRecord(end, 0));
  cutilSafeCall(cudaEventSynchronize(end));

  // calculate and report bandwidth
  float elapsedTime;
  cutilSafeCall(cudaEventElapsedTime(&elapsedTime, start, end));
  elapsedTime /= nmbIterations * 1000;  // [sec] per iteration
  const float bandwidth = dataSizeIn / elapsedTime;
  cout << "processed " << dataSizeIn / (1024. * 1024.) << " MiBytes in "
       << elapsedTime * 1000 << " msec; "
       << "total throughput: " << bandwidth / (1024 * 1024 * 1024) << " GiByte/sec" << endl;
  cutilSafeCall(cudaEventDestroy(start));
  cutilSafeCall(cudaEventDestroy(end  ));

  // copy kernel output data to host
  cutilSafeCall(cudaMemcpy(hostOutData, deviceOutData, dataSizeOut, cudaMemcpyDeviceToHost));

  // test data
  if (kernelCallerT::verify(hostInData, hostOutData, nmbBlocks,
			    nmbThreadsPerBlock, nmbElementsPerThread))
    cout << "verification successful" << endl;
  else
    cout << "verification failed" << endl;

  // unbind texture
  textureReaderT::unbindTexture();

  // cleanup memory
  free(hostInData );
  free(hostOutData);
  cutilSafeCall(cudaFree(deviceInData ));
  cutilSafeCall(cudaFree(deviceOutData));
  cudaThreadExit();
}


int main(int    argc,
	 char** argv) 
{
  // get number of CUDA devices in system
  int deviceCount = 0;
  cutilSafeCall(cudaGetDeviceCount(&deviceCount));
  if (deviceCount == 0) {
    cout << "there is no CUDA device" << endl;
    return 0;
  }
  cout << "found " << deviceCount << " CUDA device(s)" << endl;

  // print info for all CUDA devices in system
  for (int deviceId = 0; deviceId < deviceCount; ++deviceId)
    printCudaDeviceInfo(deviceId);
  
  // use most powerful GPU in system
  const int deviceId = cutGetMaxGflopsDeviceId();
  cudaDeviceProp deviceProp;
  cutilSafeCall(cudaGetDeviceProperties(&deviceProp, deviceId));
  cout << "using CUDA device[" << deviceId << "]: '" << deviceProp.name << "'" << endl;
  cutilSafeCall(cudaSetDevice(deviceId));

  // create maximum number of threads for all blocks
  const unsigned int nmbBlocks          = deviceProp.multiProcessorCount;
  const unsigned int nmbThreadsPerBlock = deviceProp.maxThreadsPerBlock;
  const unsigned int nmbIterations      = 100;
  cout << "using grid (" << nmbBlocks << " blocks) x "
       << "(" << nmbThreadsPerBlock << " threads per block)" << endl
       << "running " << nmbIterations << " kernel iterations" << endl;
  
  // run kernels
  runKernel<complex<float2, float>, floatComplexTextureReader,
     sumGlobalMemKernelCaller<complex<float2, float> > >
    (nmbBlocks, nmbThreadsPerBlock, nmbIterations);

  runKernel<complex<double2, double>, doubleComplexTextureReader,
    sumGlobalMemKernelCaller<complex<double2, double> > >
  (nmbBlocks, nmbThreadsPerBlock, nmbIterations);

  runKernel<complex<float2, float>, floatComplexTextureReader,
    sumTextureMemKernelCaller<complex<float2, float>, floatComplexTextureReader> >
  (nmbBlocks, nmbThreadsPerBlock, nmbIterations);

  runKernel<complex<double2, double>, doubleComplexTextureReader,
    sumTextureMemKernelCaller<complex<double2, double>, doubleComplexTextureReader> >
  (nmbBlocks, nmbThreadsPerBlock, nmbIterations);
}
