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
//      benchmarks coalesced memory access for global and texture memory
//      using various data types of up to 16 bytes
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
#include <algorithm>

#include <cuda_runtime.h>
#include <cutil_inline.h>
#include <cutil_math.h>  // operators for float2 and float4

#include "reportingUtils.hpp"
#include "complexTest.cuh"
#include "textureReader.cuh"


#define SHARED_MEM_BLOCK_SIZE 128  // a power of two larger than half-warp size


using namespace std;
using namespace rpwa;


// increasing ALIGN from 8 to 16 doubles throughput for double scalars
// in copy and write operations, but leaves read-only bandwidth unchanged
// effective bandwitdth for floats drops to half
// explicit specializations with correct alignment
template struct ALIGN( 8) cuda::complexStruct<float >;
template struct ALIGN(16) cuda::complexStruct<double>;


typedef cuda::complexTest<float2,                      float > float2Complex;
typedef cuda::complexTest<double2,                     double> double2Complex;
typedef cuda::complexTest<cuda::complexStruct<float >, float > floatStructComplex;
typedef cuda::complexTest<cuda::complexStruct<double>, double> doubleStructComplex;


// cutil_math.h does not define any operators for double2
inline
HOST_DEVICE
void
operator += (double2&       a,
             const double2& b)
{
	a.x += b.x;
	a.y += b.y;
}


// global memory kernels
template<typename T>
__global__
void
copyGlobalMemKernel(const T*           inData,       // pointer to device input data in global memory
                    T*                 outData,      // pointer to device output data in global memory
                    const unsigned int nmbElements)  // number of data elements for each thread
{
	const unsigned int threadId   = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned int nmbThreads = gridDim.x * blockDim.x;
#pragma unroll 1
	for (unsigned int i = 0; i < nmbElements; ++i) {
		const unsigned int index = (i * nmbThreads ) + threadId;
		outData[index] = inData[index];
	}
}

template<typename T>
__global__
void
copyGlobalMem2Kernel(const T*           inData1,      // pointer to device input data in global memory
                     const T*           inData2,      // pointer to device input data in global memory
                     T*                 outData1,     // pointer to device output data in global memory
                     T*                 outData2,     // pointer to device output data in global memory
                     const unsigned int nmbElements)  // number of data elements for each thread
{
	const unsigned int threadId   = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned int nmbThreads = gridDim.x * blockDim.x;
#pragma unroll 1
	for (unsigned int i = 0; i < nmbElements; ++i) {
		const unsigned int index = (i * nmbThreads ) + threadId;
		outData1[index] = inData1[index];
		outData2[index] = inData2[index];
	}
}

template<typename T>
__global__
void
writeOnlyGlobalMemKernel(const T*,
                         T*                 outData,      // pointer to device output data in global memory
                         const unsigned int nmbElements)  // number of data elements for each thread
{
	const unsigned int threadId   = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned int nmbThreads = gridDim.x * blockDim.x;
	T                  val;  // some dummy value
#pragma unroll 1
	for (unsigned int i = 0; i < nmbElements; ++i)
		outData[(i * nmbThreads ) + threadId] = val;
}

template<typename T>
__global__
void
writeOnlyGlobalMem2Kernel(const T*,
                          const T*,
                          T*                 outData1,     // pointer to device output data in global memory
                          T*                 outData2,     // pointer to device output data in global memory
                          const unsigned int nmbElements)  // number of data elements for each thread
{
	const unsigned int threadId   = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned int nmbThreads = gridDim.x * blockDim.x;
	T                  val;  // some dummy value
#pragma unroll 1
	for (unsigned int i = 0; i < nmbElements; ++i) {
		outData1[(i * nmbThreads ) + threadId] = val;
		outData2[(i * nmbThreads ) + threadId] = val;
	}
}

template<typename T>
__global__
void
readOnlyGlobalMemKernel(const T*           inData,       // pointer to device input data in global memory
                        T*,
                        const unsigned int nmbElements)  // number of data elements for each thread
{
	__shared__ double  sharedMem[SHARED_MEM_BLOCK_SIZE][2];  // dummy shared memory to store 16 byte values
	const unsigned int threadId   = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned int nmbThreads = gridDim.x * blockDim.x;
	T                  val;
#pragma unroll 1
	for (unsigned int i = 0; i < nmbElements; ++i)
		val += inData[(i * nmbThreads ) + threadId];
	const unsigned int index = threadIdx.x & (SHARED_MEM_BLOCK_SIZE - 1);
	*((T*)(&sharedMem[index])) = val;  // dummy operation to prevent removal of operations by compiler
	sharedMem[index][0] += 1.0;        // dummy operation to prevent removal of operations by compiler
	sharedMem[index][1] += 1.0;        // dummy operation to prevent removal of operations by compiler
}

template<typename T>
__global__
void
readOnlyGlobalMem2Kernel(const T*           inData1,      // pointer to device input data in global memory
                         const T*           inData2,      // pointer to device input data in global memory
                         T*,
                         T*,
                         const unsigned int nmbElements)  // number of data elements for each thread
{
	__shared__ double  sharedMem[SHARED_MEM_BLOCK_SIZE][2];  // dummy shared memory to store 16 byte values
	const unsigned int threadId   = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned int nmbThreads = gridDim.x * blockDim.x;
	T                  val;
#pragma unroll 1
	for (unsigned int i = 0; i < nmbElements; ++i) {
		val += inData1[(i * nmbThreads ) + threadId];
		val += inData2[(i * nmbThreads ) + threadId];
	}
	const unsigned int index = threadIdx.x & (SHARED_MEM_BLOCK_SIZE - 1);
	*((T*)(&sharedMem[index])) = val;
	sharedMem[index][0] += 1.0;        // dummy operation to prevent removal of operations by compiler
	sharedMem[index][1] += 1.0;        // dummy operation to prevent removal of operations by compiler
}


template<typename T>
__global__
void
readOnlySameLocGlobalMemKernel(const T*           inData,       // pointer to device input data in global memory
                               T*,
                               const unsigned int nmbElements)  // number of data elements for each thread
{
	__shared__ double sharedMem[SHARED_MEM_BLOCK_SIZE][2];  // dummy shared memory to store 16 byte values
	T                 val;
#pragma unroll 1
	for (unsigned int i = 0; i < nmbElements; ++i)
		val += inData[i];
	const unsigned int index = threadIdx.x & (SHARED_MEM_BLOCK_SIZE - 1);
	*((T*)(&sharedMem[index])) = val;  // dummy operation to prevent removal of operations by compiler
	sharedMem[index][0] += 1.0;        // dummy operation to prevent removal of operations by compiler
	sharedMem[index][1] += 1.0;        // dummy operation to prevent removal of operations by compiler
}

// texture memory kernels
template<typename T, typename textureReaderT>  // T is limited to 1-, 2-, or 4-component signed or unsigned 8-, 16-, or 32-bit integers, or 32-bit floats
__global__
void
copyTextureMemKernel(T*                 outData,      // pointer to device output data in global memory
                     const unsigned int nmbElements)  // number of data elements for each thread
{
	const unsigned int threadId   = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned int nmbThreads = gridDim.x * blockDim.x;
#pragma unroll 1
	for (unsigned int i = 0; i < nmbElements; ++i) {
		const unsigned int index = (i * nmbThreads ) + threadId;
		outData[index] = textureReaderT::fetch(index);
	}
}

template<typename T, typename textureReaderT>  // T is limited to 1-, 2-, or 4-component signed or unsigned 8-, 16-, or 32-bit integers, or 32-bit floats
__global__
void
readOnlyTextureMemKernel(T*,
                         const unsigned int nmbElements)  // number of data elements for each thread
{
	__shared__ T       sharedMem[SHARED_MEM_BLOCK_SIZE];  // dummy shared memory to store values
	const unsigned int threadId   = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned int nmbThreads = gridDim.x * blockDim.x;
	T                  val;
#pragma unroll 1
	for (unsigned int i = 0; i < nmbElements; ++i)
		val += textureReaderT::fetch((i * nmbThreads ) + threadId);  // dummy operation to prevent removal of operations by compiler
	const unsigned int index = threadIdx.x & (SHARED_MEM_BLOCK_SIZE - 1);
	sharedMem[index] = val;                  // dummy operation to prevent removal of operations by compiler
	*((float *)(&sharedMem[index])) += 1.0;  // dummy operation to prevent removal of operations by compiler
}


// helper macros for different kernel types
#define GLOBAL_MEM_KERNEL(kernel, elementType)	  \
	{ \
		kernel<elementType><<< nmbBlocks, nmbThreadsPerBlock >>>((elementType*)deviceInData[0], \
		                                                         (elementType*)deviceOutData[0], \
		                                                         nmbElementsPerThread); \
	}
#define GLOBAL_MEM_2_KERNEL(kernel, elementType)	  \
	{ \
		kernel<elementType><<< nmbBlocks, nmbThreadsPerBlock >>>((elementType*)deviceInData [0], \
		                                                         (elementType*)deviceInData [1], \
		                                                         (elementType*)deviceOutData[0], \
		                                                         (elementType*)deviceOutData[1], \
		                                                         nmbElementsPerThread); \
	}
#define TEXTURE_MEM_KERNEL(kernel, elementType)	  \
	{ \
		kernel<elementType, cuda::elementType ## TextureReader><<< nmbBlocks, nmbThreadsPerBlock >>> \
			((elementType*)deviceOutData[0], nmbElementsPerThread); \
	}

#define BENCHMARK(memoryType, kernel, elementType, nmbTransfersPerElement) \
	{ \
		/* dry-run kernel first to avoid any setup and caching effects*/ \
		memoryType ## _KERNEL(kernel, elementType); \
		/* setup and start timer */ \
		cudaEvent_t start, end; \
		cutilSafeCall(cudaEventCreate(&start)); \
		cutilSafeCall(cudaEventCreate(&end)); \
		cutilSafeCall(cudaEventRecord(start, 0)); \
		/* run kernel */ \
		for (unsigned int iteration = 0; iteration < nmbIterations; ++iteration) \
			memoryType ## _KERNEL(kernel, elementType); \
		/* stop timer */ \
		cutilSafeCall(cudaEventRecord(end, 0)); \
		cutilSafeCall(cudaEventSynchronize(end)); \
		/* calculate and report bandwidth */ \
		float runTime; \
		cutilSafeCall(cudaEventElapsedTime(&runTime, start, end)); \
		runTime /= nmbIterations * 1000;  /* [sec] per iteration */ \
		const float dataSize  = nmbTransfersPerElement * nmbElements * sizeof(elementType); \
		const float bandwidth = dataSize / runTime; \
		printInfo << #kernel << "<" << #elementType << "> processed " \
		          << dataSize / (1024 * 1024) << " MiBytes in " \
		          << runTime * 1000 << " msec; total throughput: " \
		          << bandwidth / (1024 * 1024 * 1024) << " GiByte/sec" << endl; \
		cutilSafeCall(cudaEventDestroy(start)); \
		cutilSafeCall(cudaEventDestroy(end)); \
	}


int main()
{

	// use most powerful GPU in sytem
	const int deviceId = cutGetMaxGflopsDeviceId();
	cudaDeviceProp deviceProp;
	cutilSafeCall(cudaGetDeviceProperties(&deviceProp, deviceId));
	printInfo << "using CUDA device[" << deviceId << "]: '" << deviceProp.name << "'" << endl;
	cutilSafeCall(cudaSetDevice(deviceId));

	// create maximum number of threads for all blocks
	const unsigned int nmbBlocks          = deviceProp.multiProcessorCount;
	const unsigned int nmbThreadsPerBlock = deviceProp.maxThreadsPerBlock;
	const unsigned int nmbIterations      = 1000;
	printInfo << "using grid (" << nmbBlocks << " blocks) x "
	          << "(" << nmbThreadsPerBlock << " threads per block); "
	          << "running " << nmbIterations << " kernel iterations" << endl;

	// create largest data arrays possible with 16 bytes words
	// !!! this does not work for textures
	// const unsigned int nmbElementsPerThread = (deviceProp.totalGlobalMem - deviceProp.totalGlobalMem / 10)  // leave 10 % margin
	//                                           / (nmbBlocks * nmbThreadsPerBlock * 4 * sizeof(float4));
	// const unsigned int nmbElements          = nmbBlocks * nmbThreadsPerBlock * nmbElementsPerThread;
	// !!! somehow there seems to be a 512 MByte limit for textures
	unsigned int nmbElements =
	  min((unsigned long)1 << 29,
	      (unsigned long)(deviceProp.totalGlobalMem - deviceProp.totalGlobalMem / 8) / 4) / sizeof(float4);
	const unsigned int nmbElementsPerThread = nmbElements / (nmbBlocks * nmbThreadsPerBlock);
	nmbElements = nmbElementsPerThread * nmbBlocks * nmbThreadsPerBlock;
	float4* deviceInData [2];
	float4* deviceOutData[2];
	printInfo << "allocating 4 * " << sizeof(float4) * nmbElements / (1024. * 1024.)
	          << " MiBytes in global memory "
	          << "(" << nmbElementsPerThread << " data elements per thread)" << endl;
	cutilSafeCall(cudaMalloc((void**) &deviceInData [0], sizeof(float4) * nmbElements));
	cutilSafeCall(cudaMalloc((void**) &deviceInData [1], sizeof(float4) * nmbElements));
	cutilSafeCall(cudaMalloc((void**) &deviceOutData[0], sizeof(float4) * nmbElements));
	cutilSafeCall(cudaMalloc((void**) &deviceOutData[1], sizeof(float4) * nmbElements));
	// bind textures
	cutilSafeCall(cudaBindTexture(0, cuda::floatTexture,  deviceInData[0], sizeof(float ) * nmbElements));
	cutilSafeCall(cudaBindTexture(0, cuda::float2Texture, deviceInData[0], sizeof(float2) * nmbElements));
	cutilSafeCall(cudaBindTexture(0, cuda::float4Texture, deviceInData[0], sizeof(float4) * nmbElements));

	if (1) {
		printInfo << "running global memory copy benchmarks ------------------------------------" << endl;
		BENCHMARK(GLOBAL_MEM,   copyGlobalMemKernel,  float,               2);
		BENCHMARK(GLOBAL_MEM,   copyGlobalMemKernel,  float2,              2);
		BENCHMARK(GLOBAL_MEM,   copyGlobalMemKernel,  float4,              2);
		BENCHMARK(GLOBAL_MEM,   copyGlobalMemKernel,  double,              2);
		BENCHMARK(GLOBAL_MEM,   copyGlobalMemKernel,  double2,             2);
		BENCHMARK(GLOBAL_MEM,   copyGlobalMemKernel,  float2Complex,       2);
		BENCHMARK(GLOBAL_MEM,   copyGlobalMemKernel,  double2Complex,      2);
		BENCHMARK(GLOBAL_MEM,   copyGlobalMemKernel,  floatStructComplex,  2);
		BENCHMARK(GLOBAL_MEM,   copyGlobalMemKernel,  doubleStructComplex, 2);
		BENCHMARK(GLOBAL_MEM_2, copyGlobalMem2Kernel, float,               4);
		BENCHMARK(GLOBAL_MEM_2, copyGlobalMem2Kernel, double,              4);
	}

	if (1) {
		printInfo << "running global memory write-only benchmarks ------------------------------" << endl;
		BENCHMARK(GLOBAL_MEM,   writeOnlyGlobalMemKernel,  float,               1);
		BENCHMARK(GLOBAL_MEM,   writeOnlyGlobalMemKernel,  float2,              1);
		BENCHMARK(GLOBAL_MEM,   writeOnlyGlobalMemKernel,  float4,              1);
		BENCHMARK(GLOBAL_MEM,   writeOnlyGlobalMemKernel,  double,              1);
		BENCHMARK(GLOBAL_MEM,   writeOnlyGlobalMemKernel,  double2,             1);
		BENCHMARK(GLOBAL_MEM,   writeOnlyGlobalMemKernel,  float2Complex,       1);
		BENCHMARK(GLOBAL_MEM,   writeOnlyGlobalMemKernel,  double2Complex,      1);
		BENCHMARK(GLOBAL_MEM,   writeOnlyGlobalMemKernel,  floatStructComplex,  1);
		BENCHMARK(GLOBAL_MEM,   writeOnlyGlobalMemKernel,  doubleStructComplex, 1);
		BENCHMARK(GLOBAL_MEM_2, writeOnlyGlobalMem2Kernel, float,               2);
		BENCHMARK(GLOBAL_MEM_2, writeOnlyGlobalMem2Kernel, double,              2);
	}

	if (1) {
		printInfo << "running global memory read-only benchmarks -------------------------------" << endl;
		BENCHMARK(GLOBAL_MEM,   readOnlyGlobalMemKernel,  float,               1);
		BENCHMARK(GLOBAL_MEM,   readOnlyGlobalMemKernel,  float2,              1);
		BENCHMARK(GLOBAL_MEM,   readOnlyGlobalMemKernel,  float4,              1);
		BENCHMARK(GLOBAL_MEM,   readOnlyGlobalMemKernel,  double,              1);
		BENCHMARK(GLOBAL_MEM,   readOnlyGlobalMemKernel,  double2,             1);
		BENCHMARK(GLOBAL_MEM,   readOnlyGlobalMemKernel,  float2Complex,       1);
		BENCHMARK(GLOBAL_MEM,   readOnlyGlobalMemKernel,  double2Complex,      1);
		BENCHMARK(GLOBAL_MEM,   readOnlyGlobalMemKernel,  floatStructComplex,  1);
		BENCHMARK(GLOBAL_MEM,   readOnlyGlobalMemKernel,  doubleStructComplex, 1);
		BENCHMARK(GLOBAL_MEM_2, readOnlyGlobalMem2Kernel, float,               2);
		BENCHMARK(GLOBAL_MEM_2, readOnlyGlobalMem2Kernel, double,              2);
	}

	if (1) {
		printInfo << "running global memory read-only benchmarks (same memory locations) -------" << endl;
		BENCHMARK(GLOBAL_MEM, readOnlySameLocGlobalMemKernel, float,               1);
		BENCHMARK(GLOBAL_MEM, readOnlySameLocGlobalMemKernel, float2,              1);
		BENCHMARK(GLOBAL_MEM, readOnlySameLocGlobalMemKernel, float4,              1);
		BENCHMARK(GLOBAL_MEM, readOnlySameLocGlobalMemKernel, double,              1);
		BENCHMARK(GLOBAL_MEM, readOnlySameLocGlobalMemKernel, double2,             1);
		BENCHMARK(GLOBAL_MEM, readOnlySameLocGlobalMemKernel, float2Complex,       1);
		BENCHMARK(GLOBAL_MEM, readOnlySameLocGlobalMemKernel, double2Complex,      1);
		BENCHMARK(GLOBAL_MEM, readOnlySameLocGlobalMemKernel, floatStructComplex,  1);
		BENCHMARK(GLOBAL_MEM, readOnlySameLocGlobalMemKernel, doubleStructComplex, 1);
	}

	if (1) {
		printInfo << "running texture memory copy benchmarks -----------------------------------" << endl;
		BENCHMARK(TEXTURE_MEM, copyTextureMemKernel, float,  2);
		BENCHMARK(TEXTURE_MEM, copyTextureMemKernel, float2, 2);
		BENCHMARK(TEXTURE_MEM, copyTextureMemKernel, float4, 2);
	}

	if (1) {
		printInfo << "running texture memory read-only benchmarks ------------------------------" << endl;
		BENCHMARK(TEXTURE_MEM, readOnlyTextureMemKernel, float,  1);
		BENCHMARK(TEXTURE_MEM, readOnlyTextureMemKernel, float2, 1);
		BENCHMARK(TEXTURE_MEM, readOnlyTextureMemKernel, float4, 1);
	}

}
