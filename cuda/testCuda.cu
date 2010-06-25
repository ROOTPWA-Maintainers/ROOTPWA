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
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <iostream>

#include <cuda_runtime.h>
#include <cutil_inline.h>


using namespace std;


bool printCudaDeviceInfo(const int deviceId)
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


// see http://forums.nvidia.com/lofiversion/index.php?t73978.html
// and http://forums.nvidia.com/index.php?showtopic=108787
template<typename storageT, typename scalarT>
struct __align__(16) testComplex {
  
  storageT _complex;

  // accessors for real and imaginary components
  inline __host__ __device__ scalarT& real() { return _complex.x; };
  inline __host__ __device__ scalarT& imag() { return _complex.y; };

  inline __host__ __device__ const scalarT& real() const { return _complex.x; };
  inline __host__ __device__ const scalarT& imag() const { return _complex.y; };

  testComplex(const scalarT& re = scalarT(),
  	      const scalarT& im = scalarT())
  {
    _complex.x = re;
    _complex.y = im;
  }
  testComplex(const testComplex& b)
  {
    _complex.x = b._complex.x;
    _complex.y = b._complex.y;
  }
  
  // assignment operators
  inline __host__ __device__ testComplex<storageT, scalarT>& operator =(const testComplex& a)
  {
    _complex.x = a._complex.x;
    _complex.y = a._complex.y;
    return *this;
  };
  inline __host__ __device__ testComplex<storageT, scalarT>& operator =(const scalarT (&a)[2])
  {
    _complex.x = a[0];
    _complex.y = a[1];
    return *this;
  };
  inline __host__ __device__ testComplex<storageT, scalarT>& operator =(const scalarT& a)
  {
    _complex.x = a;
    _complex.y = 0;
    return *this;
  };

  // basic arithmetic operations with complex numbers
  inline __host__ __device__ testComplex<storageT, scalarT>& operator +=(const testComplex<storageT, scalarT>& b)
  {
    this->_complex.x += b._complex.x;
    this->_complex.y += b._complex.y;
    return *this;
  }
  inline __host__ __device__ testComplex<storageT, scalarT> operator +(const testComplex<storageT, scalarT>& b) const
  {
    testComplex<storageT, scalarT> result = {{_complex.x + b._complex.x, _complex.y  + b._complex.y}};
    return result;
  }
  inline __host__ __device__ testComplex<storageT, scalarT> operator -(const testComplex<storageT, scalarT>& b) const
  {
    testComplex<storageT, scalarT> result = {{_complex.x - b._complex.x, _complex.y  - b._complex.y}};
    return result;
  }
  inline __host__ __device__ testComplex<storageT, scalarT> operator *(const testComplex<storageT, scalarT>& b) const
  {
    testComplex<storageT, scalarT> result = {{_complex.x * b._complex.x - _complex.y * b._complex.y,
					      _complex.y * b._complex.x + _complex.x * b._complex.y}};
    return result;
  }
  inline __host__ __device__ testComplex<storageT, scalarT> operator /(const testComplex<storageT, scalarT>& b) const
  {
    scalarT denom = (b._complex.x * b._complex.x + b._complex.y * b._complex.y );
    testComplex<storageT, scalarT> result = {{(_complex.x * b._complex.x + _complex.y * b._complex.y ) / denom,
					      (_complex.y * b._complex.x - _complex.x * b._complex.y ) / denom}};
    return result;
  }


  // basic arithmetic operations with scalars
  inline __host__ __device__ testComplex<storageT, scalarT> operator +(const scalarT& b) const
  {
    testComplex<storageT, scalarT> result = {{_complex.x + b, _complex.y}};
    return result;
  }
  inline __host__ __device__ testComplex<storageT, scalarT> operator -(const scalarT& b) const
  {
    testComplex<storageT, scalarT> result = {{_complex.x - b, _complex.y}};
    return result;
  }
  inline __host__ __device__ testComplex<storageT, scalarT> operator *(const scalarT& b) const
  {
    testComplex<storageT, scalarT> result = {{_complex.x * b, _complex.y * b}};
    return result;
  }
  inline __host__ __device__ testComplex<storageT, scalarT> operator /(const scalarT& b) const
  {
    testComplex<storageT, scalarT> result = {{_complex.x / b, _complex.y / b}};
    return result;
  }

  // negate a complex number
  inline __host__ __device__ testComplex<storageT, scalarT> operator -() const
  {
    testComplex<storageT, scalarT> result = {{-_complex.x, -_complex.y}};
    return result;
  }

  // complex conjugate
  inline __host__ __device__ testComplex<storageT, scalarT> conj() const
  {
    testComplex<storageT, scalarT> result = {{_complex.x, -_complex.y}};
    return result;
  }

  // complex norm = Re^2 + Im^2
  inline __host__ __device__ scalarT norm() const
  {
    return _complex.x * _complex.x + _complex.y * _complex.y;
  }

  // complex absolute value
  inline __host__ __device__ scalarT abs() const
  {
    return sqrt(norm());
  }

  // complex phase angle
  inline __host__ __device__ testComplex<storageT, scalarT> arg() const
  {
    return atan2(_complex.y, _complex.x);
  }

  // a possible alternative to a testComplex constructor
  static __host__ __device__ testComplex<storageT, scalarT> maketestComplex(const scalarT& a, const scalarT& b)
  {
    testComplex<storageT, scalarT> res;
    res.real() = a;
    res.imag() = b;
    return res;
  }

  // constants: 0, 1, i
  static __host__ __device__ const testComplex<storageT, scalarT> zero()
  {
    return maketestComplex((scalarT)0, (scalarT)0);
  }
  static __host__ __device__ const testComplex<storageT, scalarT> one()
  {
    return maketestComplex((scalarT)1, (scalarT)0);
  }
  static __host__ __device__ const testComplex<storageT, scalarT> i()
  {
    return maketestComplex((scalarT)0, (scalarT)1);
  }
};

template<typename storageT, typename scalarT>
__host__ __device__ testComplex<storageT, scalarT>
operator +(const scalarT&                        a,
	   const testComplex<storageT, scalarT>& b)
{
  testComplex<storageT, scalarT> result = {{a + b.value.x, b.value.y}};
  return result;
}
template<typename storageT, typename scalarT>
__host__ __device__ testComplex<storageT, scalarT>
operator -(const scalarT&                        a,
	   const testComplex<storageT, scalarT>& b)
{
  testComplex<storageT, scalarT> result = {{a - b.value.x, -b.value.y}};
  return result;
}
template<typename storageT, typename scalarT>
__host__ __device__ testComplex<storageT, scalarT>
operator *(const scalarT&                        a,
	   const testComplex<storageT, scalarT>& b)
{
  testComplex<storageT, scalarT> result = {{a * b._complex.x, a * b._complex.y}};
  return result;
}
template<typename storageT, typename scalarT>
__host__ __device__ testComplex<storageT, scalarT>
operator /(const scalarT&                        a,
	   const testComplex<storageT, scalarT>& b)
{
  scalarT denom = b._complex.x * b._complex.x + b._complex.y * b._complex.y;
  testComplex<storageT, scalarT> result = {{a * b._complex.x / denom, -a * b._complex.y / denom}};
  return result;
}


template<typename T>
struct __align__(16) complexStorage {
  T x;
  T y;
};


template<typename T>
__global__ void
testGlobMemReadBandwidthKernel(const T*           inData,       // pointer to input data in global memory
			       T*                 outData,      // pointer to output data in global memory
			       const unsigned int nmbElements)  // number of data elements for each thread
{
  const unsigned int threadId   = blockDim.x * blockIdx.x + threadIdx.x;
  const unsigned int nmbThreads = gridDim.x * blockDim.x;
  T sum = 0;
  // #pragma unroll 1
  for (unsigned int i = 0; i < nmbElements; ++i)
    sum += inData[(i * nmbThreads ) + threadId];
  outData[threadId] = sum;
}


template<typename T>
void
testGlobMemBandwidth(const unsigned int nmbBlocks            = 30,
		     const unsigned int nmbThreadsPerBlock   = 512,
		     const unsigned int nmbElementsPerThread = 10000)
{
  cout << "testing global memory bandwidth" << endl;

  // allocate host memory
  const unsigned int nmbElements = nmbBlocks * nmbThreadsPerBlock * nmbElementsPerThread;
  const unsigned int memSizeIn   = nmbElements * sizeof(T);
  const unsigned int memSizeOut  = nmbBlocks * nmbThreadsPerBlock * sizeof(T);
  cout << "allocating " << memSizeIn / (1024. * 1024.) << " MiBytes of global memory for " << nmbElements << " data elements:" << endl
       << "(" << nmbBlocks << " blocks) x (" << nmbThreadsPerBlock << " threads per block) x "
       << "(" << nmbElementsPerThread <<  " number data elements per thread)"<< endl;
  T* hostInData  = (T*) malloc(memSizeIn);
  T* hostOutData = (T*) malloc(memSizeOut);

  // initalize host memory
  for (unsigned int i = 0; i < nmbElements; ++i)
    hostInData[i] = (T)i;

  // allocate device memory and copy host data to device
  T* deviceInData;
  cutilSafeCall(cudaMalloc((void**) &deviceInData, memSizeIn));
  cutilSafeCall(cudaMemcpy(deviceInData, hostInData, memSizeIn, cudaMemcpyHostToDevice));
  T* deviceOutData;
  cutilSafeCall(cudaMalloc((void**) &deviceOutData, memSizeOut));

  // initialize timer
  unsigned int timer = 0;
  cutilCheckError(cutCreateTimer(&timer));
  cutilCheckError(cutStartTimer (timer));
  cudaEvent_t start, stop;
  cutilSafeCall(cudaEventCreate(&start));
  cutilSafeCall(cudaEventCreate(&stop));
  cutilSafeCall(cudaEventRecord(start, 0));

  // execute test kernel
  testGlobMemReadBandwidthKernel<T><<< nmbBlocks, nmbThreadsPerBlock >>>(deviceInData, deviceOutData, nmbElementsPerThread);
  cutilCheckMsg("there were errors executing testGlobMemBandwidthKernel");  // check if kernel execution generated an error
  
  // calculate memory bandwidth
  cutilSafeCall(cudaEventRecord(stop, 0));
  cutilSafeCall(cudaEventSynchronize(stop));
  cutilSafeThreadSync();
  cutilCheckError(cutStopTimer(timer));
  const float elapsedTimeCutTimer = cutGetTimerValue(timer) / 1000.;  // [sec]
  float elapsedTimeEvent;  // [sec]
  cutilSafeCall(cudaEventElapsedTime(&elapsedTimeEvent, start, stop));
  elapsedTimeEvent /= 1000.;
  cout << "CUT timer time: " << elapsedTimeCutTimer * 1000 << " msec" << endl
       << "global memory read bandwidth: "<< memSizeIn / (1024 * 1024 * 1024 * elapsedTimeCutTimer) << " GiByte/sec" << endl
       << "CUDA event time: " << elapsedTimeEvent * 1000 << " msec" << endl
       << "global memory read bandwidth: "<< memSizeIn / (1024 * 1024 * 1024 * elapsedTimeEvent) << " GiByte/sec" << endl;
  cutilSafeCall(cudaEventDestroy(start));
  cutilSafeCall(cudaEventDestroy(stop));
  cutilCheckError(cutDeleteTimer(timer));

  // copy kernel output data to host
  cutilSafeCall(cudaMemcpy(hostOutData, deviceOutData, memSizeOut, cudaMemcpyDeviceToHost));

  // cleanup memory
  free(hostInData);
  free(hostOutData);
  cutilSafeCall(cudaFree(deviceInData));
  cutilSafeCall(cudaFree(deviceOutData));
  cudaThreadExit();

  testComplex<double2, double> a = 1.;
  testComplex<double2, double> b(2., 3.);
  //testComplex<double2, double> c = {(double)3., (double)4.};
  // = (2. * a + b * a - a) / 2.;
}


template<typename T>
__global__ void
testGlobMemComplexReadBandwidthKernel(const T*           inDataRe,     // pointer to real parts in global memory
				      const T*           inDataIm,     // pointer to imaginary parts in global memory
				      T*                 outData,      // pointer to output data in global memory
				      const unsigned int nmbElements)  // number of data elements for each thread
{
  const unsigned int threadId   = blockDim.x * blockIdx.x + threadIdx.x;
  const unsigned int nmbThreads = gridDim.x * blockDim.x;
  T sumRe = 0;
  T sumIm = 0;
  // #pragma unroll 1
  for (unsigned int i = 0; i < nmbElements; ++i) {
    sumRe += inDataRe[(i * nmbThreads ) + threadId];
    sumIm += inDataIm[(i * nmbThreads ) + threadId];
  }
  outData[threadId] = sumRe + sumIm;
}


template<typename T>
void
testGlobMemComplexBandwidth(const unsigned int nmbBlocks            = 30,
			    const unsigned int nmbThreadsPerBlock   = 512,
			    const unsigned int nmbElementsPerThread = 10000)
{
  cout << "testing global memory bandwidth for complex SoA" << endl;

  // allocate host memory
  const unsigned int nmbElements = nmbBlocks * nmbThreadsPerBlock * nmbElementsPerThread;
  const unsigned int memSizeIn   = nmbElements * sizeof(T);
  const unsigned int memSizeOut  = nmbBlocks * nmbThreadsPerBlock * sizeof(T);
  cout << "allocating " << 2 * memSizeIn / (1024. * 1024.) << " MiBytes of global memory for " << nmbElements << " data elements:" << endl
       << "(" << nmbBlocks << " blocks) x (" << nmbThreadsPerBlock << " threads per block) x "
       << "(" << nmbElementsPerThread <<  " number data elements per thread)"<< endl;
  T* hostInDataRe = (T*) malloc(memSizeIn);
  T* hostInDataIm = (T*) malloc(memSizeIn);
  T* hostOutData  = (T*) malloc(memSizeOut);

  // initalize host memory
  for (unsigned int i = 0; i < nmbElements; ++i) {
    hostInDataRe[i] = (T)i;
    hostInDataIm[i] = (T)i;
  }

  // allocate device memory and copy host data to device
  T* deviceInDataRe;
  cutilSafeCall(cudaMalloc((void**) &deviceInDataRe, memSizeIn));
  cutilSafeCall(cudaMemcpy(deviceInDataRe, hostInDataRe, memSizeIn, cudaMemcpyHostToDevice));
  T* deviceInDataIm;
  cutilSafeCall(cudaMalloc((void**) &deviceInDataIm, memSizeIn));
  cutilSafeCall(cudaMemcpy(deviceInDataIm, hostInDataIm, memSizeIn, cudaMemcpyHostToDevice));
  T* deviceOutData;
  cutilSafeCall(cudaMalloc((void**) &deviceOutData, memSizeOut));

  // initialize timer
  unsigned int timer = 0;
  cutilCheckError(cutCreateTimer(&timer));
  cutilCheckError(cutStartTimer (timer));
  cudaEvent_t start, stop;
  cutilSafeCall(cudaEventCreate(&start));
  cutilSafeCall(cudaEventCreate(&stop));
  cutilSafeCall(cudaEventRecord(start, 0));

  // execute test kernel
  testGlobMemComplexReadBandwidthKernel<T><<< nmbBlocks, nmbThreadsPerBlock >>>(deviceInDataRe, deviceInDataIm, deviceOutData, nmbElementsPerThread);
  cutilCheckMsg("there were errors executing testGlobMemBandwidthKernel");  // check if kernel execution generated an error
  
  // calculate memory bandwidth
  cutilSafeCall(cudaEventRecord(stop, 0));
  cutilSafeCall(cudaEventSynchronize(stop));
  cutilSafeThreadSync();
  cutilCheckError(cutStopTimer(timer));
  const float elapsedTimeCutTimer = cutGetTimerValue(timer) / 1000.;  // [sec]
  float elapsedTimeEvent;  // [sec]
  cutilSafeCall(cudaEventElapsedTime(&elapsedTimeEvent, start, stop));
  elapsedTimeEvent /= 1000.;
  cout << "CUT timer time: " << elapsedTimeCutTimer * 1000 << " msec" << endl
       << "global memory read bandwidth: "<< 2 * memSizeIn / (1024 * 1024 * 1024 * elapsedTimeCutTimer) << " GiByte/sec" << endl
       << "CUDA event time: " << elapsedTimeEvent * 1000 << " msec" << endl
       << "global memory read bandwidth: "<< 2 * memSizeIn / (1024 * 1024 * 1024 * elapsedTimeEvent) << " GiByte/sec" << endl;
  cutilSafeCall(cudaEventDestroy(start));
  cutilSafeCall(cudaEventDestroy(stop));
  cutilCheckError(cutDeleteTimer(timer));

  // copy kernel output data to host
  cutilSafeCall(cudaMemcpy(hostOutData, deviceOutData, memSizeOut, cudaMemcpyDeviceToHost));

  // cleanup memory
  free(hostInDataRe);
  free(hostInDataIm);
  free(hostOutData);
  cutilSafeCall(cudaFree(deviceInDataRe));
  cutilSafeCall(cudaFree(deviceInDataIm));
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
  
  // take most powerful GPU
  const int deviceId = cutGetMaxGflopsDeviceId();
  cudaDeviceProp deviceProp;
  cutilSafeCall(cudaGetDeviceProperties(&deviceProp, deviceId));
  cout << "using CUDA device[" << deviceId << "]: '" << deviceProp.name << "'" << endl;
  cutilSafeCall(cudaSetDevice(deviceId));

  const unsigned int nmbBlocks            = deviceProp.multiProcessorCount;
  const unsigned int nmbThreadsPerBlock   = deviceProp.maxThreadsPerBlock;
  const unsigned int nmbElementsPerThread = 5000;
  // cout << "testing global memory read bandwidth for int:" << endl;
  // testGlobMemBandwidth<int>(nmbBlocks, nmbThreadsPerBlock, nmbElementsPerThread);
  // cout << "testing global memory read bandwidth for float:" << endl;
  // testGlobMemBandwidth<float>(nmbBlocks, nmbThreadsPerBlock, nmbElementsPerThread);
  cout << "testing global memory read bandwidth for double:" << endl;
  testGlobMemBandwidth<double>(nmbBlocks, nmbThreadsPerBlock, nmbElementsPerThread);
  cout << endl;
  cout << "testing global memory read bandwidth for complex<float2, float>:" << endl;
  testGlobMemBandwidth<testComplex<float2, float> >(nmbBlocks, nmbThreadsPerBlock, nmbElementsPerThread);
  cout << endl;
  cout << "testing global memory read bandwidth for complex<double2, double>:" << endl;
  testGlobMemBandwidth<testComplex<double2, double> >(nmbBlocks, nmbThreadsPerBlock, nmbElementsPerThread);
  cout << endl;
  cout << "testing global memory read bandwidth for complex<struct<double>, double>:" << endl;
  testGlobMemBandwidth<testComplex<complexStorage<double>, double> >(nmbBlocks, nmbThreadsPerBlock, nmbElementsPerThread);
  cout << endl;
  cout << "testing global memory read bandwidth for float re[], im[]:" << endl;
  testGlobMemComplexBandwidth<float>(nmbBlocks, nmbThreadsPerBlock, nmbElementsPerThread);
  cout << endl;
  cout << "testing global memory read bandwidth for double re[], im[]:" << endl;
  testGlobMemComplexBandwidth<double>(nmbBlocks, nmbThreadsPerBlock, nmbElementsPerThread);
}
