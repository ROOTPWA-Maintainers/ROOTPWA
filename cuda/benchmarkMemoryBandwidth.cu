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
//      benchmarks coalesced memory access for global and texture memory
//      using various data types
// 
//      based on MisterAnderson42's post in the nvidia forum
//      http://forums.nvidia.com/index.php?s=&showtopic=52806&view=findpost&p=292058
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


#define SHARED_MEM_BLOCK_SIZE 128


texture<float,  1, cudaReadModeElementType> textureFloat;
texture<float2, 1, cudaReadModeElementType> textureFloat2;
texture<float4, 1, cudaReadModeElementType> textureFloat4;


// global memory kernels
template<typename T>
__global__ void
copyGlobalMemKernel(T* inData,
		    T* outData,
		    T)
{
  const unsigned int threadId = threadIdx.x + blockIdx.x * blockDim.x;
  outData[threadId] = inData[threadId];
}

template<typename T>
__global__ void
writeOnlyGlobalMemKernel(T* inData,
			 T* outData,
			 T  c)
{
  const unsigned int threadId = threadIdx.x + blockIdx.x * blockDim.x;
  outData[threadId] = c;
}

template<typename T>
__global__ void
readOnlyGlobalMemKernel(T* inData,
			T* outData,
			T)
{
  const unsigned int threadId = threadIdx.x + blockIdx.x * blockDim.x;
  __shared__ T shared[SHARED_MEM_BLOCK_SIZE];
  shared[threadIdx.x] = inData[threadId];
  *((float *)(&shared[(threadIdx.x + 1) & (SHARED_MEM_BLOCK_SIZE - 1)])) += 1.0;
}


// texture memory kernels
__global__ void
copyTextureMemFloatKernel(float* inData,
			  float* outData,
			  float)
{
  const unsigned int threadId = threadIdx.x + blockIdx.x * blockDim.x;
  outData[threadId] = tex1Dfetch(textureFloat, threadId);
}

__global__ void
copyTextureMemFloat2Kernel(float2* inData,
			   float2* outData,
			   float2)
{
  const unsigned int threadId = threadIdx.x + blockIdx.x * blockDim.x;
  outData[threadId] = tex1Dfetch(textureFloat2, threadId);
}

__global__ void
copyTextureMemFloat4Kernel(float4* inData,
			   float4* outData,
			   float4)
{
  const unsigned int threadId = threadIdx.x + blockIdx.x * blockDim.x;
  outData[threadId] = tex1Dfetch(textureFloat4, threadId);
}

__global__ void
readOnlyTextureMemFloatKernel(float* inData,
			      float* outData,
			      float)
{
  const unsigned int threadId = threadIdx.x + blockIdx.x * blockDim.x;
  __shared__ float shared[SHARED_MEM_BLOCK_SIZE];
  shared[threadIdx.x] = tex1Dfetch(textureFloat, threadId);
  *((float *)(&shared[(threadIdx.x + 1) & (SHARED_MEM_BLOCK_SIZE - 1)])) += 1.0;
}
 
__global__ void
readOnlyTextureMemFloat2Kernel(float2* inData,
			       float2* outData,
			       float2)
{
  const unsigned int threadId = threadIdx.x + blockIdx.x * blockDim.x;
  __shared__ float2 shared[SHARED_MEM_BLOCK_SIZE];
  shared[threadIdx.x] = tex1Dfetch(textureFloat2, threadId);
  *((float *)(&shared[(threadIdx.x + 1) & (SHARED_MEM_BLOCK_SIZE - 1)])) += 1.0;
}
 
__global__ void
readOnlyTextureMemFloat4Kernel(float4* inData, 
			       float4* outData,
			       float4)
{
  const unsigned int threadId = threadIdx.x + blockIdx.x * blockDim.x;
  __shared__ float4 shared[SHARED_MEM_BLOCK_SIZE];
  shared[threadIdx.x] = tex1Dfetch(textureFloat4, threadId);
  *((float *)(&shared[(threadIdx.x + 1) & (SHARED_MEM_BLOCK_SIZE - 1)])) += 1.0;
}


#define BENCHMARK(kernel, elementType, value, kernelName, nmbTransfersPerElement) \
  {									\
    kernel<<< grid, threads >>>((elementType *)deviceInData, (elementType *)deviceOutData, value); \
									\
    cudaEvent_t start, end;						\
    cutilSafeCall(cudaEventCreate(&start));				\
    cutilSafeCall(cudaEventCreate(&end));				\
    cutilSafeCall(cudaEventRecord(start, 0));				\
									\
    for (int i = 0; i < nmbIterations; ++i)				\
      kernel<<< grid, threads >>>((elementType *)deviceInData, (elementType *)deviceOutData, value); \
									\
    cutilSafeCall(cudaEventRecord(end, 0));				\
    cutilSafeCall(cudaEventSynchronize(end));				\
    float runTime;							\
    cutilSafeCall(cudaEventElapsedTime(&runTime, start, end));		\
    runTime /= float(nmbIterations) * 1000;				\
    const float bandwidth = nmbTransfersPerElement * nmbElements * sizeof(elementType) / runTime; \
    cout << kernelName << " bandwidth: " << bandwidth / (1024 * 1024 * 1024) << " GiByte/sec" << endl; \
    cutilSafeCall(cudaEventDestroy(start));				\
    cutilSafeCall(cudaEventDestroy(end));				\
  }


int main()
{

  const unsigned int nmbElements   = 1 << 22;
  const unsigned int nmbThreads    = SHARED_MEM_BLOCK_SIZE;
  const unsigned int nmbIterations = 10000;

  float4 *deviceInData, *deviceOutData;
  cutilSafeCall(cudaMalloc((void**)&deviceInData,  sizeof(float4) * nmbElements));
  cutilSafeCall(cudaMalloc((void**)&deviceOutData, sizeof(float4) * nmbElements));
  cutilSafeCall(cudaBindTexture(0, textureFloat,  deviceInData, sizeof(float ) * nmbElements));
  cutilSafeCall(cudaBindTexture(0, textureFloat2, deviceInData, sizeof(float2) * nmbElements));
  cutilSafeCall(cudaBindTexture(0, textureFloat4, deviceInData, sizeof(float4) * nmbElements));
	
  dim3 threads(nmbThreads, 1, 1);
  dim3 grid(nmbElements / nmbThreads, 1, 1);

  BENCHMARK(copyGlobalMemKernel<float  >, float,   0.0f,                                "copyGlobalMemKernel<float>",   2);
  BENCHMARK(copyGlobalMemKernel<float2 >, float2,  make_float2(0.0f, 0.0f),             "copyGlobalMemKernel<float2>",  2);
  BENCHMARK(copyGlobalMemKernel<float4 >, float4,  make_float4(0.0f, 0.0f, 0.0f, 0.0f), "copyGlobalMemKernel<float4>",  2);
  BENCHMARK(copyGlobalMemKernel<double >, double,  0.0,                                 "copyGlobalMemKernel<double>",  2);
  BENCHMARK(copyGlobalMemKernel<double2>, double2, make_double2(0.0, 0.0),              "copyGlobalMemKernel<double2>", 2);
	
  cout << endl;
  BENCHMARK(writeOnlyGlobalMemKernel<float>,   float,   0.0f,                                "writeOnlyGlobalMemKernel<float>",   1);
  BENCHMARK(writeOnlyGlobalMemKernel<float2>,  float2,  make_float2(0.0f, 0.0f),             "writeOnlyGlobalMemKernel<float2>",  1);
  BENCHMARK(writeOnlyGlobalMemKernel<float4>,  float4,  make_float4(0.0f, 0.0f, 0.0f, 0.0f), "writeOnlyGlobalMemKernel<float4>",  1);
  BENCHMARK(writeOnlyGlobalMemKernel<double>,  double,  0.0,                                 "writeOnlyGlobalMemKernel<double>",  1);
  BENCHMARK(writeOnlyGlobalMemKernel<double2>, double2, make_double2(0.0, 0.0),              "writeOnlyGlobalMemKernel<double2>", 1);
	
  cout << endl;
  BENCHMARK(readOnlyGlobalMemKernel<float>,   float,   0.0f,                                "readOnlyGlobalMemKernel<float>",   1);
  BENCHMARK(readOnlyGlobalMemKernel<float2>,  float2,  make_float2(0.0f, 0.0f),             "readOnlyGlobalMemKernel<float2>",  1);
  BENCHMARK(readOnlyGlobalMemKernel<float4>,  float4,  make_float4(0.0f, 0.0f, 0.0f, 0.0f), "readOnlyGlobalMemKernel<float4>",  1);
  BENCHMARK(readOnlyGlobalMemKernel<double>,  double,  0.0,                                 "readOnlyGlobalMemKernel<double>",  1);
  BENCHMARK(readOnlyGlobalMemKernel<double2>, double2, make_double2(0.0, 0.0),              "readOnlyGlobalMemKernel<double2>", 1);
	
  cout << endl;
  BENCHMARK(copyTextureMemFloatKernel,  float, 0.0f,                                 "copyTextureMemKernel<float>",  2);
  BENCHMARK(copyTextureMemFloat2Kernel, float2, make_float2(0.0f, 0.0f),             "copyTextureMemKernel<float2>", 2);
  BENCHMARK(copyTextureMemFloat4Kernel, float4, make_float4(0.0f, 0.0f, 0.0f, 0.0f), "copyTextureMemKernel<float4>", 2);
	
  cout << endl;
  BENCHMARK(readOnlyTextureMemFloatKernel,  float,  0.0f,                                "readOnlyTextureMemKernel<float>",  1);
  BENCHMARK(readOnlyTextureMemFloat2Kernel, float2, make_float2(0.0f, 0.0f),             "readOnlyTextureMemKernel<float2>", 1);
  BENCHMARK(readOnlyTextureMemFloat4Kernel, float4, make_float4(0.0f, 0.0f, 0.0f, 0.0f), "readOnlyTextureMemKernel<float4>", 1);

}
