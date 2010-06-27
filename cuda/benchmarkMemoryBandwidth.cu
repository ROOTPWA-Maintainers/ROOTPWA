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
#include <cutil_math.h>  // operators for float2 and float4


#define SHARED_MEM_BLOCK_SIZE 128  // a power of two larger than half-warp size


using namespace std;


// cutil_math.h does not define anything for double2
inline __host__ __device__
void
operator += (double2&       a,
	     const double2& b)
{
    a.x += b.x;
    a.y += b.y;
}


// global memory kernels
template<typename T>
__global__ void
copyGlobalMemKernel(const T*           inData,       // pointer to input data in global memory
		    T*                 outData,      // pointer to output data in global memory
		    const unsigned int nmbElements)  // number of data elements for each thread
{
  const unsigned int threadId   = blockDim.x * blockIdx.x + threadIdx.x;
  const unsigned int nmbThreads = gridDim.x * blockDim.x;
  // #pragma unroll 1
  for (unsigned int i = 0; i < nmbElements; ++i) {
    const unsigned int index = (i * nmbThreads ) + threadId;
    outData[index] = inData[index];
  }
}

template<typename T>
__global__ void
writeOnlyGlobalMemKernel(const T*,
			 T*                 outData,      // pointer to output data in global memory
			 const unsigned int nmbElements)  // number of data elements for each thread
{
  const unsigned int threadId   = blockDim.x * blockIdx.x + threadIdx.x;
  const unsigned int nmbThreads = gridDim.x * blockDim.x;
  T                  val;  // some dummy value
  // #pragma unroll 1
  for (unsigned int i = 0; i < nmbElements; ++i)
    outData[(i * nmbThreads ) + threadId] = val;
}

template<typename T>
__global__ void
readOnlyGlobalMemKernel(const T*           inData,       // pointer to input data in global memory
			T*,
			const unsigned int nmbElements)  // number of data elements for each thread
{
  __shared__ T       sharedMem[SHARED_MEM_BLOCK_SIZE];  // dummy shared memory to store values
  const unsigned int threadId   = blockDim.x * blockIdx.x + threadIdx.x;
  const unsigned int nmbThreads = gridDim.x * blockDim.x;
  T                  val;
  // #pragma unroll 1
  for (unsigned int i = 0; i < nmbElements; ++i)
    val += inData[(i * nmbThreads ) + threadId];  // dummy operation to prevent removal of operations by compiler
  const unsigned int index = threadIdx.x & (SHARED_MEM_BLOCK_SIZE - 1);
  sharedMem[index] = val;                  // dummy operation to prevent removal of operations by compiler
  *((float *)(&sharedMem[index])) += 1.0;  // dummy operation to prevent removal of operations by compiler
}


// texture memory kernels

// since texture references are implicit static file scope variables, they cannot be passed as kernel parameters
// see http://forums.nvidia.com/lofiversion/index.php?t70630.html
// use helper structures to work around this
texture<float,  1, cudaReadModeElementType> floatTexture;
struct floatTextureReader {
  __device__ float operator ()(const int index) const { return tex1Dfetch(floatTexture, index); }
};
texture<float2, 1, cudaReadModeElementType> float2Texture;
struct float2TextureReader {
  __device__ float2 operator ()(const int index) const { return tex1Dfetch(float2Texture, index); }
};
texture<float4, 1, cudaReadModeElementType> float4Texture;
struct float4TextureReader {
  __device__ float4 operator ()(const int index) const { return tex1Dfetch(float4Texture, index); }
};

template<typename T, typename TTextureReader>  // T is limited to 1-, 2-, or 4-component signed or unsigned 8-, 16-, or 32-bit integers, or 32-bit floats
__global__ void
copyTextureMemKernel(T*                 outData,      // pointer to output data in global memory
		     const unsigned int nmbElements)  // number of data elements for each thread
{
  const unsigned int threadId   = blockDim.x * blockIdx.x + threadIdx.x;
  const unsigned int nmbThreads = gridDim.x * blockDim.x;
  TTextureReader     texFetch;
  // #pragma unroll 1
  for (unsigned int i = 0; i < nmbElements; ++i) {
    const unsigned int index = (i * nmbThreads ) + threadId;
    outData[index] = texFetch(index);
  }
}

template<typename T, typename TTextureReader>  // T is limited to 1-, 2-, or 4-component signed or unsigned 8-, 16-, or 32-bit integers, or 32-bit floats
__global__ void
readOnlyTextureMemKernel(T*,
			 const unsigned int nmbElements)  // number of data elements for each thread
{
  __shared__ T       sharedMem[SHARED_MEM_BLOCK_SIZE];  // dummy shared memory to store values
  const unsigned int threadId   = blockDim.x * blockIdx.x + threadIdx.x;
  const unsigned int nmbThreads = gridDim.x * blockDim.x;
  TTextureReader     texFetch;
  T                  val;
  // #pragma unroll 1
  for (unsigned int i = 0; i < nmbElements; ++i)
    val += texFetch((i * nmbThreads ) + threadId);  // dummy operation to prevent removal of operations by compiler
  const unsigned int index = threadIdx.x & (SHARED_MEM_BLOCK_SIZE - 1);
  sharedMem[index] = val;                  // dummy operation to prevent removal of operations by compiler
  *((float *)(&sharedMem[index])) += 1.0;  // dummy operation to prevent removal of operations by compiler
}
 
// __global__ void
// readOnlyTextureMemFloat4Kernel(float4* inData, 
// 			       float4* outData,
// 			       float4)
// {
//   const unsigned int threadId = threadIdx.x + blockIdx.x * blockDim.x;
//   __shared__ float4 shared[SHARED_MEM_BLOCK_SIZE];
//   shared[threadIdx.x] = tex1Dfetch(textureFloat4, threadId);
//   *((float *)(&shared[(threadIdx.x + 1) & (SHARED_MEM_BLOCK_SIZE - 1)])) += 1.0;
// }


#define BENCHMARK_GLOB_MEM(kernel, elementType, nmbTransfersPerElement)	\
  {									\
    /* dry-run kernel first to avoid any setup and caching effects*/	\
    kernel<elementType><<< nmbBlocks, nmbThreadsPerBlock >>>((elementType*)deviceInData, (elementType*)deviceOutData, nmbElementsPerThread); \
    /* setup and start timer */						\
    cudaEvent_t start, end;						\
    cutilSafeCall(cudaEventCreate(&start));				\
    cutilSafeCall(cudaEventCreate(&end));				\
    cutilSafeCall(cudaEventRecord(start, 0));				\
    /* run kernel */							\
    for (unsigned int iteration = 0; iteration < nmbIterations; ++iteration) \
      kernel<elementType><<< nmbBlocks, nmbThreadsPerBlock >>>((elementType*)deviceInData, (elementType*)deviceOutData, nmbElementsPerThread); \
    /* stop timer */							\
    cutilSafeCall(cudaEventRecord(end, 0));				\
    cutilSafeCall(cudaEventSynchronize(end));				\
    /* calculate and report bandwidth */				\
    float runTime;							\
    cutilSafeCall(cudaEventElapsedTime(&runTime, start, end));		\
    runTime /= nmbIterations * 1000;  /* [sec] per iteration */		\
    const float bandwidth = nmbTransfersPerElement * nmbElements * sizeof(elementType) / runTime; \
    cout << #kernel << "<" << #elementType << "> bandwidth: " << bandwidth / (1024 * 1024 * 1024) << " GiByte/sec" << endl; \
    cutilSafeCall(cudaEventDestroy(start));				\
    cutilSafeCall(cudaEventDestroy(end));				\
  }


#define BENCHMARK_TEX_MEM(kernel, elementType, nmbTransfersPerElement)	\
  {									\
    /* bind texture to input data */								\
    cutilSafeCall(cudaBindTexture(0, elementType ## Texture,  deviceInData, sizeof(elementType) * nmbElements)); \
    /* dry-run kernel first to avoid any setup and caching effects*/	\
    kernel<elementType, elementType ## TextureReader><<< nmbBlocks, nmbThreadsPerBlock >>>((elementType*)deviceOutData, nmbElementsPerThread); \
    /* setup and start timer */						\
    cudaEvent_t start, end;						\
    cutilSafeCall(cudaEventCreate(&start));				\
    cutilSafeCall(cudaEventCreate(&end));				\
    cutilSafeCall(cudaEventRecord(start, 0));				\
    /* run kernel */							\
    for (unsigned int iteration = 0; iteration < nmbIterations; ++iteration) \
      kernel<elementType, elementType ## TextureReader><<< nmbBlocks, nmbThreadsPerBlock >>>((elementType*)deviceOutData, nmbElementsPerThread); \
    /* stop timer */							\
    cutilSafeCall(cudaEventRecord(end, 0));				\
    cutilSafeCall(cudaEventSynchronize(end));				\
    /* unbind texture */						\
    cutilSafeCall(cudaUnbindTexture(elementType ## Texture));		\
    /* calculate and report bandwidth */				\
    float runTime;							\
    cutilSafeCall(cudaEventElapsedTime(&runTime, start, end));		\
    runTime /= nmbIterations * 1000;  /* [sec] per iteration */		\
    const float bandwidth = nmbTransfersPerElement * nmbElements * sizeof(elementType) / runTime; \
    cout << #kernel << "<" << #elementType << "> bandwidth: " << bandwidth / (1024 * 1024 * 1024) << " GiByte/sec" << endl; \
    cutilSafeCall(cudaEventDestroy(start));				\
    cutilSafeCall(cudaEventDestroy(end));				\
  }


int main()
{

  // use most powerful GPU in sytem
  const int deviceId = cutGetMaxGflopsDeviceId();
  cudaDeviceProp deviceProp;
  cutilSafeCall(cudaGetDeviceProperties(&deviceProp, deviceId));
  cout << "using CUDA device[" << deviceId << "]: '" << deviceProp.name << "'" << endl;
  cutilSafeCall(cudaSetDevice(deviceId));

  // create maximum number of threads for all blocks
  const unsigned int nmbBlocks            = deviceProp.multiProcessorCount;
  const unsigned int nmbThreadsPerBlock   = deviceProp.maxThreadsPerBlock;
  const unsigned int nmbIterations        = 100;
  cout << "using grid (" << nmbBlocks << " blocks) x (" << nmbThreadsPerBlock << " threads per block)" << endl
       << "running " << nmbIterations << " kernel iterations" << endl;

  // create largest device data arrays possible with 16 bytes words
  const unsigned int nmbElementsPerThread = (deviceProp.totalGlobalMem - deviceProp.totalGlobalMem / 10)  // leave 10 % margin
                                            / (nmbBlocks * nmbThreadsPerBlock * 2 * sizeof(float4));
  const unsigned int nmbElements          = nmbBlocks * nmbThreadsPerBlock * nmbElementsPerThread;
  float4*            deviceInData;
  float4*            deviceOutData;
  cout << "allocating 2 * " << sizeof(float4) * nmbElements / (1024. * 1024.) << " MiBytes in global memory "
       << "(" << nmbElementsPerThread << " data elements per thread)" << endl;
  cutilSafeCall(cudaMalloc((void**) &deviceInData,  sizeof(float4) * nmbElements));
  cutilSafeCall(cudaMalloc((void**) &deviceOutData, sizeof(float4) * nmbElements));
  
  if (0) {
    cout << "--------------------------------------------------------------------------------" << endl
	 << "running global memory copy benchmarks" << endl;
    BENCHMARK_GLOB_MEM(copyGlobalMemKernel, float,   2);
    BENCHMARK_GLOB_MEM(copyGlobalMemKernel, float2,  2);
    BENCHMARK_GLOB_MEM(copyGlobalMemKernel, float4,  2);
    BENCHMARK_GLOB_MEM(copyGlobalMemKernel, double,  2);
    BENCHMARK_GLOB_MEM(copyGlobalMemKernel, double2, 2);
  }

  if (0) {
    cout << "--------------------------------------------------------------------------------" << endl
	 << "running global memory write-only benchmarks" << endl;
    BENCHMARK_GLOB_MEM(writeOnlyGlobalMemKernel, float,   1);
    BENCHMARK_GLOB_MEM(writeOnlyGlobalMemKernel, float2,  1);
    BENCHMARK_GLOB_MEM(writeOnlyGlobalMemKernel, float4,  1);
    BENCHMARK_GLOB_MEM(writeOnlyGlobalMemKernel, double,  1);
    BENCHMARK_GLOB_MEM(writeOnlyGlobalMemKernel, double2, 1);
  }

  if (1) {
    cout << "--------------------------------------------------------------------------------" << endl
	 << "running global memory read-only benchmarks" << endl;
    BENCHMARK_GLOB_MEM(readOnlyGlobalMemKernel, float,   1);
    BENCHMARK_GLOB_MEM(readOnlyGlobalMemKernel, float2,  1);
    BENCHMARK_GLOB_MEM(readOnlyGlobalMemKernel, float4,  1);
    BENCHMARK_GLOB_MEM(readOnlyGlobalMemKernel, double,  1);
    BENCHMARK_GLOB_MEM(readOnlyGlobalMemKernel, double2, 1);
  }

  if (1) {
    cout << "--------------------------------------------------------------------------------" << endl
	 << "running texture memory copy benchmarks" << endl;
    BENCHMARK_TEX_MEM(copyTextureMemKernel, float,  2);
    BENCHMARK_TEX_MEM(copyTextureMemKernel, float2, 2);
    BENCHMARK_TEX_MEM(copyTextureMemKernel, float4, 2);
  }
	
  if (1) {
    cout << "--------------------------------------------------------------------------------" << endl
	 << "running texture memory read-only benchmarks" << endl;
    BENCHMARK_TEX_MEM(readOnlyTextureMemKernel, float,  1);
    BENCHMARK_TEX_MEM(readOnlyTextureMemKernel, float2, 1);
    BENCHMARK_TEX_MEM(readOnlyTextureMemKernel, float4, 1);
  }

}
