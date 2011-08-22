// square.cu for testing cuda applications


#include <math.h>
#include <stdio.h>
#include <cuda.h>
#include <cutil_inline.h>
#include <time.h>
#include <iostream>
#include <unistd.h>
#include <TFile.h>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <TROOT.h>
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include <TGraph.h>
#include <TLine.h>
#include <TH2.h>
#include <TH3F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TList.h>
#include <TTree.h>
#include <TRandom3.h>
#include "templtest_kernel.cu"
#include "factorial.hpp"
#include "clebschGordanCoeffnew.hpp"
#include "dFunctionnew.hpp"
#include "complex.cuh"
// using namespace rpwa;

#define N 1
#define coa 1
#define PI 3.1415926535897932384626433832795028841971693993751058209749
// cardinality of the cgc-calculation depending on a: N = 32a^6 - 72a^5 + 60a^4 - 22a^3 + 3a^2 is the total number of possible cgc's depending on a...
// explanation later on

// Functions
    void Cleanup(void);
    void RandomInit(double*, int);
    void cgcInitj(int*, int);
    void cgcInitm(int*, int);
    void cgcInittheta(double*, int);
    void cgcInitJ(int*, int);
    void cgcInitM(int*, int);
    void waiter(int);

// Kernel that executes on the CUDA device

// __global__ void square(const int *A, const int *B, const int *C, const int *D, const int *E, const int *F, double *result)
__global__ void square(const int *A, const int *B, const int *C, const double *D, cuda::complex<double> *result)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int facto = cufac(0);
  double cgctest1 			= cucgc(4,0,4,-4,4,-4);
  double cgctest2 			= cucgcwiki(8,0,0,0,4,0);
  double dfunctest 			= cudfunc(1,1,1,3.1415);
  double dfunc 				= cudfunc(1,1,1,PI/2);
  cuda::complex<double> Dfunctest1 	= cuDfunc(1,1,1,PI,PI/2,0);
  cuda::complex<double> Dfunctest2	= cuDfuncReflConj(4,2,4,1,-1,-0.564732,1.29395,0);
  cuda::complex<double> SHtest1		= cuSpherHarm(0,0,0,0);
//  cuda::complex<double> ampltest1	= cuHelAmplitude(PI/3,PI/5,PI/10,PI/4,1.5,1.27,4,2,4,4,4,0,4,0);
  double bumtest1			= cuBreakUpMomentum(1.5,1.27,0.13957);
  double bftest				= BFactor(2,bumtest1);
  int test 				= cuabs(5);
  bool test1 				= cuisOdd(3);
  int test2 				= cupmo(765672);
  cuda::complex<double> bwigdat		(0.183925,0.007331); 
  
//   if (idx<N) { result[idx] = cucgc(A[idx],B[idx],C[idx],D[idx],E[idx],F[idx]) ; }
//    if (idx<N) { result[idx] = cudfunc(A[idx],B[idx],C[idx],D[idx]) ; }
   if (idx<N) { result[idx] = 5.0*( cucgc(4,0,0,0,4,0)*cucgc(4,0,4,0,4,0)*cuDfuncReflConj(4,2,0,-1,1,-0.564732,1.29395,0)*cuDfuncConj(4,0,0,-1.6509,1.11757,0) + 
				    cucgc(4,2,0,0,4,2)*cucgc(4,0,4,2,4,2)*cuDfuncReflConj(4,2,2,-1,1,-0.564732,1.29395,0)*cuDfuncConj(4,2,0,-1.6509,1.11757,0) +
				    cucgc(4,4,0,0,4,4)*cucgc(4,0,4,4,4,4)*cuDfuncReflConj(4,2,4,-1,1,-0.564732,1.29395,0)*cuDfuncConj(4,4,0,-1.6509,1.11757,0) +
				    cucgc(4,-2,0,0,4,-2)*cucgc(4,0,4,-2,4,-2)*cuDfuncReflConj(4,2,-2,-1,1,-0.564732,1.29395,0)*cuDfuncConj(4,-2,0,-1.6509,1.11757,0) +
				    cucgc(4,-4,0,0,4,-4)*cucgc(4,0,4,-4,4,-4)*cuDfuncReflConj(4,2,-4,-1,1,-0.564732,1.29395,0)*cuDfuncConj(4,-4,0,-1.6509,1.11757,0)) * 
				    sqrt(BFactor(4,0.50514288800993423)) * 
				    (cucgc(0,0,0,0,0,0)*cucgc(4,0,0,0,4,0)*sqrt(BFactor(4,0.25846533761461093))*bwigdat); }
}


// main routine that executes on the host

int main(void)
{

  /// Device query
  cudaDeviceProp prop;
  int count;
  cutilSafeCall( cudaGetDeviceCount( &count ) );
  for (int i=0; i<count; i++) {
    cutilSafeCall( cudaGetDeviceProperties( &prop, i) );
    printf("Working on a %s\n", prop.name);
    printf("Compute capability: %d.%d\n", prop.major, prop.minor );
    printf("Clock rate: %d MHz\n", prop.clockRate/1024 );
    printf("Multiprocessor count: %d\n", prop.multiProcessorCount );   
    printf("Total global memory: %d MByte\n", prop.totalGlobalMem/(1024*1024) );
    printf("Threads in warp: %d\n", prop.warpSize );
    printf("Shared memory available per block: %d KB", prop.sharedMemPerBlock/1024 );
    
    
     printf("Major revision number:         %d\n",  prop.major);
    printf("Minor revision number:         %d\n",  prop.minor);
    printf("Name:                          %s\n",  prop.name);
    printf("Total global memory:           %u\n",  prop.totalGlobalMem);
    printf("Total shared memory per block: %u\n",  prop.sharedMemPerBlock);
    printf("Total registers per block:     %d\n",  prop.regsPerBlock);
    printf("Warp size:                     %d\n",  prop.warpSize);
    printf("Maximum memory pitch:          %u\n",  prop.memPitch);
    printf("Maximum threads per block:     %d\n",  prop.maxThreadsPerBlock);
    for (int i = 0; i < 3; ++i)
    printf("Maximum dimension %d of block:  %d\n", i, prop.maxThreadsDim[i]);
    for (int i = 0; i < 3; ++i)
    printf("Maximum dimension %d of grid:   %d\n", i, prop.maxGridSize[i]);
    printf("Clock rate:                    %d\n",  prop.clockRate);
    printf("Total constant memory:         %u\n",  prop.totalConstMem);
    printf("Texture alignment:             %u\n",  prop.textureAlignment);
    printf("Concurrent copy and execution: %s\n",  (prop.deviceOverlap ? "Yes" : "No"));
    printf("Number of multiprocessors:     %d\n",  prop.multiProcessorCount);

  }
  
  
  // Variables
    int* h_A;
    int* h_B;
    int* h_C;
    double* h_D;
/*
    int* h_E;
    int* h_F;
*/
    cuda::complex<double>* h_result;
    int* d_A;
    int* d_B;
    int* d_C;
    double* d_D;
/*
    int* d_E;
    int* d_F;
*/
    cuda::complex<double>* d_result;
    double* CPUpow;
    
    // Allocating memory for array on host to avoid segmentation fault
    CPUpow = (double *)malloc(N*sizeof(double));
      if(NULL == CPUpow) {
      printf("Error at malloc...?\n");
    }
    
  // Zeiten
    float t_c;		// check
    float t_htd;	// host to device
    float t_dth;	// device to host
    float t_cod;	// calculation on device
  
    printf("\n\n NOW RANDOMTEST: Exponentiating %i Random-Numbers\n\n", N);
    
     //Timestamp check with cuda
      cudaEvent_t start, stop;
      cutilSafeCall( cudaEventCreate(&start) );
      cutilSafeCall( cudaEventCreate(&stop)  );
      
      
      cutilSafeCall( cudaEventRecord( start, 0) );
      
	time_t t0 = time(NULL);
	waiter(1);
	time_t t1 = time(NULL);
	printf("\nTime Check! Time elapsed (measured by CPU): %.0f second(s)\n", float(t1-t0) );
	
      cutilSafeCall( cudaEventRecord( stop, 0) );
      cutilSafeCall( cudaEventSynchronize( stop ) );
      cutilSafeCall( cudaEventElapsedTime( &t_c, start, stop) );
      printf("Time elapsed(GPU): %3.3f ms\n\n",t_c);

      
   
    printf("Now calculation on device...\n\n");
  
    size_t size = N * sizeof(int);
  
    // Allocate randomnumbers h_A and h_B in host memory
    h_A = (int*)malloc(size);
    if (h_A == 0) free(h_A);
    h_B = (int*)malloc(size);
    if (h_B == 0) free(h_B);
    h_C = (int*)malloc(size);
    if (h_C == 0) free(h_C);
    h_D = (double*)malloc(N*sizeof(double));
    if (h_D == 0) free(h_D);
/*
    h_E = (int*)malloc(size);    
    if (h_E == 0) free(h_E);
    h_F = (int*)malloc(size);
    if (h_F == 0) free(h_F);
*/    

//     cutilSafeCall( cudaHostAlloc((void**)&h_result, N*sizeof(cuda::complex<double>)*4, cudaHostAllocDefault) );    
    h_result = (cuda::complex<double>*)malloc(1+N*sizeof(cuda::complex<double>));
    if (h_result == 0) free(h_result);

    gRandom->SetSeed(23434543);
    // Initialize randomnumbers etc.
//     RandomInit(h_A, N);
//     RandomInit(h_B, N);
/*    // cgc test
    cgcInitj(h_A, N);
    cgcInitm(h_B, N);
    cgcInitj(h_C, N);
    cgcInitm(h_D, N);    
    cgcInitJ(h_E, N);
    cgcInitM(h_F, N);  
*/
    // dfunc test
    cgcInitj(h_A, N);
    cgcInitm(h_B, N);
    cgcInitm(h_C, N);
    cgcInittheta(h_D, N);      
    
    // Allocate randomnumbers in device memory
    cutilSafeCall( cudaMalloc((void**)&d_A, size) );
    cutilSafeCall( cudaMalloc((void**)&d_B, size) );
    cutilSafeCall( cudaMalloc((void**)&d_C, size) );
    cutilSafeCall( cudaMalloc((void**)&d_D, N*sizeof(double)) );
/*
    cutilSafeCall( cudaMalloc((void**)&d_E, size) );
    cutilSafeCall( cudaMalloc((void**)&d_F, size) );
*/
    cutilSafeCall( cudaMalloc((void**)&d_result, N*sizeof(cuda::complex<double>)) );    
   
    // Copy randomnumbers from host memory to device memory
    cutilSafeCall( cudaEventRecord( start, 0) );

    cutilSafeCall( cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(d_C, h_C, size, cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(d_D, h_D, N*sizeof(double), cudaMemcpyHostToDevice) );
/*
    cutilSafeCall( cudaMemcpy(d_E, h_E, size, cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(d_F, h_F, size, cudaMemcpyHostToDevice) );
*/
    
    cutilSafeCall( cudaEventRecord( stop, 0) );
    cutilSafeCall( cudaEventSynchronize( stop ) );
    cutilSafeCall( cudaEventElapsedTime( &t_htd, start, stop));
    printf("%.3f ms required for uploading\n",t_htd);
    printf("%f GB/s Upload\n", ((7*size)/(t_htd/1e3))/(1024*1024*1024) );

    
    // Invoke kernel
    cutilSafeCall( cudaEventRecord( start, 0) );
    
    int threadsPerBlock = 512;
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
//     square<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, d_D, d_E, d_F, d_result);
    square<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, d_D, d_result);

    
//     cutilSafeCall( cudaThreadSynchronize() );
    cutilSafeCall( cudaEventRecord( stop, 0) );
    cutilSafeCall( cudaEventSynchronize( stop ) );
    cutilSafeCall( cudaEventElapsedTime( &t_cod, start, stop));

    printf("%.3f ms required for calculation on device\n", t_cod );
    printf("Bandwidth: %f GB/s\n", ((11*size)/(t_cod/1e3))/(1024*1024*1024) );
    
    
    // Copy result from device memory to host memory
    cutilSafeCall( cudaEventRecord( start, 0) );
    
    cutilSafeCall( cudaMemcpy(h_result, d_result, N*sizeof(cuda::complex<double>), cudaMemcpyDeviceToHost) );
    
    cutilSafeCall( cudaEventRecord( stop, 0) );
    cutilSafeCall( cudaEventSynchronize( stop ) );
    cutilSafeCall( cudaEventElapsedTime( &t_dth, start, stop));
    printf("%.3f ms required for downloading\n", t_dth ); 
    printf("%f GB/s Download\n\n", ((4*size)/(t_dth/1e3))/(1024*1024*1024) );

    
///     Human testing here 
       int r=0;
       for (r=0;r<N;r++) { 
	std::cout << "GPU: " << h_result[r] << std::endl;
// 	std::cout << "CPU: " << DFunction< std::complex<double> >(1,1,1,PI,PI/2,0) << std::endl;
	std::cout << "CPU: " << clebschGordanCoeff<double>(8,0,0,0,4,0) << std::endl;
// 	printf("output gpu: %f\n", h_result[r]);
// 	printf("output cpu: %f\n", clebschGordanCoeff<double>(8,0,0,0,4,0));
// 	printf("output cpu: %f\n", dFunction<double>(1,1,1,3.1415));
// 	printf("%i. -- GPU: < %2.1f %2.1f %2.1f %2.1f | %2.1f %2.1f > = %2.6f\n", r,  (float)h_A[r]/2, (float)h_B[r]/2, (float)h_C[r]/2, (float)h_D[r]/2, (float)h_E[r]/2, (float)h_F[r]/2, h_result[r]); 
// 	printf("%i. -- CPU: < %2.1f %2.1f %2.1f %2.1f | %2.1f %2.1f > = %2.6f\n", r,  (float)h_A[r]/2, (float)h_B[r]/2, (float)h_C[r]/2, (float)h_D[r]/2, (float)h_E[r]/2, (float)h_F[r]/2, clebschGordanCoeff<double>(h_A[r],h_B[r],h_C[r],h_D[r],h_E[r],h_F[r]) );
// 	printf("         %i,   %i,   %i,   %i,    %i,   %i\n", h_A[r], h_B[r], h_C[r], h_D[r], h_E[r], h_F[r]);
       }

    
    
    
    // Calculation on host
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    printf("\nNow calculation on host: Time: %s \n\n", asctime(timeinfo) );
    timespec h_cal_beg, h_cal_end; 
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &h_cal_beg); 
    
    int i=0;
    for (i=0;i<N;i++) {
//       CPUpow[i] = clebschGordanCoeff<double>(h_A[i],h_B[i],h_C[i],h_D[i],h_E[i],h_F[i]);
      CPUpow[i] = dFunction<double>(h_A[i],h_B[i],h_C[i],h_D[i]);
    }
	/// CPU function here

					      
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &h_cal_end); 
    printf("%.3f ms required for calculation on host\n\n", ( ( (h_cal_end.tv_sec - h_cal_beg.tv_sec) ) + (h_cal_end.tv_nsec - h_cal_beg.tv_nsec) / 1e9 ) * 1e3 ); 
    printf("Pure calculation speed up: GPU is %3.1f times faster than the CPU.\n",  ( ( ( (h_cal_end.tv_sec - h_cal_beg.tv_sec) ) + (h_cal_end.tv_nsec - h_cal_beg.tv_nsec) / 1e9 ) * 1e3 ) / (t_cod) );
    printf("Over-all speed up:         GPU is %3.1f times faster than the CPU.\n\n",  ( ( ( (h_cal_end.tv_sec - h_cal_beg.tv_sec) ) + (h_cal_end.tv_nsec - h_cal_beg.tv_nsec) / 1e9 ) * 1e3 ) / (t_cod+t_htd+t_dth) );  
     
    
/*
    
    
    // Verify result
    int k;
    for (k=0;k<N;k++) {

      if (fabs(h_result[k] - CPUpow[k]) > 1e-5) {
	printf("%i. paires of variates incorrect!\n",k);
// 	printf("GPU: < %2.1f %2.1f %2.1f %2.1f | %2.1f %2.1f > = %2.6f ----- CPU: < %2.1f %2.1f %2.1f %2.1f | %2.1f %2.1f > = %2.6f\n", (float)h_A[k]/2, (float)h_B[k]/2, (float)h_C[k]/2, (float)h_D[k]/2, (float)h_E[k]/2, (float)h_F[k]/2, h_result[k],
// 																	(float)h_A[k]/2, (float)h_B[k]/2, (float)h_C[k]/2, (float)h_D[k]/2, (float)h_E[k]/2, (float)h_F[k]/2, CPUpow[k]); 
	printf("GPU: dWigner(%f,%f,%f,%f) = %2.6f ----- CPU: dWigner(%f,%f,%f,%f) = %2.6f\n", (float)h_A[k]/2, (float)h_B[k]/2, (float)h_C[k]/2, h_D[k], h_result[k], (float)h_A[k]/2, (float)h_B[k]/2, (float)h_C[k]/2, h_D[k], CPUpow[k]); 
	printf("%20.10f - %20.10f = %.10f > 1e-9\n", h_result[k] , CPUpow[k], h_result[k] - CPUpow[k]);
	break;
      }
    }
    
    printf("TEST %s \n", (k == N) ? "PASSED" : "FAILED");

    

    
    
    // Visualize results via root   

    // create a root file and a TTree
    // the structure to hold the variables for the branch
    struct errors_t {
      Double_t errabs; 		// absolute errors
      Double_t errrel; 		// relative errors
      Double_t inputA;		// self-explaining
      Double_t inputB;
      Double_t inputC;
      Double_t inputD;
//       Double_t inputE;
//       Double_t inputF;
      Double_t output;
    };
    errors_t errors;


    // create a new ROOT file
    TFile *f = new TFile("errors.root","RECREATE");
    // create a TTree
    TTree *tree = new TTree("T","Absolute und relative Errors");
    // create one branch with all information from the stucture
    tree->Branch("Errors",&errors.errabs,"errabs/D:errrel:inputA:inputB:inputC:inputD:output");
    // fill the tree 
    for (int i=0;i<N;i++) {
      errors.errabs = h_result[i] - CPUpow[i];
      errors.errrel = 1 - (h_result[i] / CPUpow[i]);
      errors.inputA = h_A[i];
      errors.inputB = h_B[i];
      errors.inputC = h_C[i];
      errors.inputD = h_D[i];
//       errors.inputE = h_E[i];
//       errors.inputF = h_F[i];
      errors.output = h_result[i];
      tree->Fill();
    }
    // check what the tree looks like
    tree->Print();
    f->Write();


*/
    
    
    
    // Free host memory
        free(h_A);
        free(h_B);
        free(h_C);
        free(h_D);
//         free(h_E);
//         free(h_F);
	free(h_result);
//         cutilSafeCall( cudaFreeHost(h_result) ); 
    
    // Free device memory
        cudaFree(d_A);
        cudaFree(d_B);
        cudaFree(d_C);
        cudaFree(d_D);
//         cudaFree(d_E);
//         cudaFree(d_F);
	cudaFree(d_result);
        cutilSafeCall( cudaEventDestroy(start) );
	cutilSafeCall( cudaEventDestroy(stop) );

    
    
    
    
}
    
//    Cleanup
void Cleanup(void)
{
    exit(0);
}


// Allocates an array with random float entries.
void RandomInit(double* data, int n)
{
//     srand( time(NULL) );
    for (int i = 0; i < n; ++i) {
//       data[i] = (rand()%10000+1)/1239. ;
      data[i] = gRandom->Rndm();
    }
}


void cgcInitj(int* data, int n)
{
  for (int i = 0; i <= n; i++) {
    data[i] = 1 + abs((int)(gRandom->Rndm()*2));
  }
}

void cgcInitm(int* data, int n)
{
  for (int i = 0; i <= n; i++) {
    data[i] = -3 + (int)(gRandom->Rndm()*6);
  }
}
void cgcInittheta(double* data, int n)
{
  for (int i = 0; i <= n; i++) {
    data[i] = 0 + (gRandom->Rndm()*3.1415);
  }
}
void cgcInitJ(int* data, int n)
{
  for (int i = 0; i <= n; i++) {
    data[i] = abs((int)(gRandom->Rndm()*6));
  }
}
void cgcInitM(int* data, int n)
{
  for (int i = 0; i <= n; i++) {
    data[i] = -5 + (int)(gRandom->Rndm()*10);
  }
}


// Timestamp
void waiter ( int seconds ) // ja kellner
{
  clock_t endwait;
  endwait = clock () + seconds * CLOCKS_PER_SEC ;
  while (clock() < endwait) {}
}







