#include <math.h>
#include <stdio.h>
#include <cuda.h>
#include <cutil_inline.h>
#include <time.h>
#include <iostream>
#include <unistd.h>
#include <typeinfo>
#include "complex.cuh"
#include <unistd.h>
#include <TFile.h>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>

#define imin(a,b) (a<b?a:b)

// Allocates an array with random double entries.
template<class T>
T RandomInit(T* data, int n)
{
    srand( time(NULL) );
    for (int i = 0; i < n; i++) {
      data[i] = (rand()%10000+1)/12. ;
    }
    return 0;
}

const int N = 50000000;
const int threadsPerBlock = 512;
const int blocksPerGrid = imin( 32, (N + threadsPerBlock-1) / threadsPerBlock );

__global__ void tsum( double* A, double* sum, int N)
{
  __shared__ double cache[threadsPerBlock];
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  int cacheIndex = threadIdx.x;
  
  double temp = 0;
  
  while (idx < N) {
    temp += A[idx];
    idx += blockDim.x * gridDim.x;
  }
  
  cache[cacheIndex] = temp;
  
  __syncthreads();
  
  int i = blockDim.x/2;
  while (i != 0) {
    if (cacheIndex < i) {
      cache[cacheIndex] += cache[cacheIndex + i];
    }
    __syncthreads();
    i /= 2;
  }
  if (cacheIndex == 0) {
    sum[blockIdx.x] = cache[0];
  }
}

int main(int argc, char *argv[]) {
  int N = atoi( argv[1] );
  float timer[10];
  cudaEvent_t start, stop;
  cutilSafeCall( cudaEventCreate(&start) );
  cutilSafeCall( cudaEventCreate(&stop)  );
  
  double *A, *partial_result, C, CPU_result;
  double *dev_A, *dev_partial_result;
  
  cutilSafeCall( cudaHostAlloc((void**)&A, N*sizeof(double), cudaHostAllocDefault) );
  cutilSafeCall( cudaHostAlloc((void**)&partial_result, N*sizeof(double), cudaHostAllocDefault) );
  
  cutilSafeCall( cudaMalloc((void**)&dev_A,N*sizeof(double)) );
  cutilSafeCall( cudaMalloc((void**)&dev_partial_result,blocksPerGrid*sizeof(double)) );
  
  timespec init_beg, init_end; 
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &init_beg);    
  RandomInit(A,N);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &init_end); 
  printf("%.3f ms required for initializing random numbers\n", ( ( (init_end.tv_sec - init_beg.tv_sec) ) + (init_end.tv_nsec - init_beg.tv_nsec) / 1e9 ) * 1e3 );     

  cutilSafeCall( cudaEventRecord( start, 0) );  
  cutilSafeCall( cudaMemcpy(dev_A,A,N*sizeof(double),cudaMemcpyHostToDevice) );
  cutilSafeCall( cudaEventRecord( stop, 0) );
  cutilSafeCall( cudaEventSynchronize( stop ) );
  cutilSafeCall( cudaEventElapsedTime( &timer[1], start, stop));
  printf("%.3f ms required for uploading\n",timer[1]);
  printf("%f GB/s Uploading bandwidth\n", (N*sizeof(double)/(timer[1]/1e3))/(1024*1024*1024) ); 
  
  cutilSafeCall( cudaEventRecord( start, 0) );  
  tsum<<<blocksPerGrid,threadsPerBlock>>>(dev_A,dev_partial_result,N);
  cutilSafeCall( cudaEventRecord( stop, 0) );
  cutilSafeCall( cudaEventSynchronize( stop ) );
  cutilSafeCall( cudaEventElapsedTime( &timer[0], start, stop));
  printf("--------------------------------%.3f ms required for calculation on device\n",timer[0]);
  printf("%f GB/s Calculation bandwidth\n", ((N+blocksPerGrid)*sizeof(double)/(timer[0]/1e3))/(1024*1024*1024) );   
  
  cutilSafeCall( cudaEventRecord( start, 0) );    
  cutilSafeCall( cudaMemcpy(partial_result,dev_partial_result,blocksPerGrid*sizeof(double),cudaMemcpyDeviceToHost) );
  cutilSafeCall( cudaEventRecord( stop, 0) );
  cutilSafeCall( cudaEventSynchronize( stop ) );
  cutilSafeCall( cudaEventElapsedTime( &timer[2], start, stop));
  printf("%.3f ms required for downloading\n",timer[2]);
  printf("%f GB/s Downloading bandwidth\n", (blocksPerGrid*sizeof(double)/(timer[2]/1e3))/(1024*1024*1024) ); 
  
  printf("%.3f ms required for whole process\n",timer[0]+timer[1]+timer[2]);
  
  C = 0;
  for (int i=0; i<blocksPerGrid; i++) {
    C += partial_result[i];
  }
  
  timespec cpu_calc_beg, cpu_calc_end; 
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_calc_beg);     
  CPU_result = 0;
  for (int i = 0; i<N; i++) {
    CPU_result += A[i];
  }
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_calc_end); 
  printf("--------------------------------%.3f ms required for calculation on host\n", ( ( (cpu_calc_end.tv_sec - cpu_calc_beg.tv_sec) ) + (cpu_calc_end.tv_nsec - cpu_calc_beg.tv_nsec) / 1e9 ) * 1e3 );     

  printf("GPU value = %.16f =?= \nCPU value = %.16f\n",C,CPU_result);
  
  cudaFree(dev_A);
  cudaFree(dev_partial_result);
  
  cutilSafeCall( cudaFreeHost(A) );
  cutilSafeCall( cudaFreeHost(partial_result) );
  
  cutilSafeCall( cudaEventDestroy(start) );
  cutilSafeCall( cudaEventDestroy(stop) );
}
