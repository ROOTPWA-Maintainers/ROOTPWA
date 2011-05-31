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


#define N 500000
  
// Functions
    void Cleanup(void);
    void RandomInit(double*, int);
    void waiter(int);


// Kernel that executes on the CUDA device

__global__ void square(const double *A, const double *B, double *C)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
//   if (idx<N) C[idx] = (2.0*7.5*B[idx]*B[idx])/(3*A[idx]*A[idx]);						/// GPU function here
//   if (idx<N) C[idx] = A[idx]+A[idx]*B[idx];
  if (idx<N) C[idx] = (A[idx]*A[idx]*A[idx]*A[idx]*A[idx]*A[idx]*A[idx]*A[idx]*A[idx]*A[idx])+B[idx] ;
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
    printf("Clock rate: %d MHz\n", prop.clockRate/1000 );
    printf("Multiprocessor count: %d\n", prop.multiProcessorCount );   
    printf("Total global memory: %d Byte\n", prop.totalGlobalMem );
    printf("Threads in warp: %d\n", prop.warpSize );
  }
  
  
  
  
  
  
  // Variables
    double* h_A;
    double* h_B;
    double* h_C;
    double* d_A;
    double* d_B;
    double* d_C;
    double* CPUpow;
    
    // Allocating memory for array on host to avoid segmentation fault
    CPUpow = (double *)malloc(N*sizeof(double));
      if(NULL == CPUpow) {
      printf("Error at malloc(gamma)...?\n");
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
  
    size_t size = N * sizeof(double);
  
    // Allocate randomnumbers h_A and h_B in host memory
    h_A = (double*)malloc(size);
    if (h_A == 0) free(h_A);
    h_B = (double*)malloc(size);
    if (h_B == 0) free(h_B);
    h_C = (double*)malloc(size);
    if (h_C == 0) free(h_C);

    gRandom->SetSeed(23434543);
    // Initialize randomnumbers
    RandomInit(h_A, N);
    RandomInit(h_B, N);
  
    // Allocate randomnumbers in device memory
    cutilSafeCall( cudaMalloc((void**)&d_A, size) );
    cutilSafeCall( cudaMalloc((void**)&d_B, size) );
    cutilSafeCall( cudaMalloc((void**)&d_C, size) );

   
    // Copy randomnumbers from host memory to device memory
    cutilSafeCall( cudaEventRecord( start, 0) );

    cutilSafeCall( cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice) );
    
    cutilSafeCall( cudaEventRecord( stop, 0) );
    cutilSafeCall( cudaEventSynchronize( stop ) );
    cutilSafeCall( cudaEventElapsedTime( &t_htd, start, stop));
    printf("%.3f ms required for uploading\n",t_htd);
    printf("%f GB/s Upload\n", ((2*size)/(t_htd/1e3))/(1024*1024*1024) );

    
    // Invoke kernel
    cutilSafeCall( cudaEventRecord( start, 0) );
    
    int threadsPerBlock = 512;
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
    square<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C); //, N);
    
//     cutilSafeCall( cudaThreadSynchronize() );
    cutilSafeCall( cudaEventRecord( stop, 0) );
    cutilSafeCall( cudaEventSynchronize( stop ) );
    cutilSafeCall( cudaEventElapsedTime( &t_cod, start, stop));

    printf("%.3f ms required for calculation on device\n", t_cod );
    printf("Bandwidth: %f GB/s\n", ((3*size)/(t_cod/1e3))/(1024*1024*1024) );
    
    
    // Copy result from device memory to host memory
    cutilSafeCall( cudaEventRecord( start, 0) );
    
    cutilSafeCall( cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost) );
    
    cutilSafeCall( cudaEventRecord( stop, 0) );
    cutilSafeCall( cudaEventSynchronize( stop ) );
    cutilSafeCall( cudaEventElapsedTime( &t_dth, start, stop));
    printf("%.3f ms required for downloading\n", t_dth ); 
    printf("%f GB/s Download\n", ((size)/(t_dth/1e3))/(1024*1024*1024) );

    
///     Human testing here 
//      int r;
//      for (r=0;r<N;r++) { printf("GPU: %6.3f to the power of %6.3f = %20.10f\n", h_A[r], h_B[r], h_C[r]); }
//      for (r=0;r<N;r++) { printf("CPU: %6.3f to the power of %6.3f = %20.10f\n", h_A[r], h_B[r], pow(h_A[r],h_B[r]) ); }
    
    
    
    // Calculation on host
    printf("\nNow calculation on host...\n\n");
    timespec h_cal_beg, h_cal_end; 
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &h_cal_beg); 
    
    int i;
    for (i = 0; i<N; i++) {      CPUpow[i] = 
					      // (2.0*7.5*h_B[i]*h_B[i])/(3*h_A[i]*h_A[i]) ;   }					/// CPU function here
					      // h_A[i]+h_A[i]*h_B[i] ; }
					      (h_A[i]*h_A[i]*h_A[i]*h_A[i]*h_A[i]*h_A[i]*h_A[i]*h_A[i]*h_A[i]*h_A[i])+h_B[i] ; }
					      
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &h_cal_end); 
    printf("%.3f ms required for calculation on host\n\n", ( ( (h_cal_end.tv_sec - h_cal_beg.tv_sec) ) + (h_cal_end.tv_nsec - h_cal_beg.tv_nsec) / 1e9 ) * 1e3 ); 
    printf("Pure calculation speed up: GPU is %3.1f times faster than the CPU.\n",  ( ( ( (h_cal_end.tv_sec - h_cal_beg.tv_sec) ) + (h_cal_end.tv_nsec - h_cal_beg.tv_nsec) / 1e9 ) * 1e3 ) / (t_cod) );
    printf("Over-all speed up:         GPU is %3.1f times faster than the CPU.\n\n",  ( ( ( (h_cal_end.tv_sec - h_cal_beg.tv_sec) ) + (h_cal_end.tv_nsec - h_cal_beg.tv_nsec) / 1e9 ) * 1e3 ) / (t_cod+t_htd+t_dth) );  
     
    // Verify result
    int k;
    for (k=0;k<N;k++) {

      if (fabs(h_C[k] - CPUpow[k]) > 1e-9) {
	printf("%i. paires of variates incorrect!\n",k);
	printf("GPU: %6.3f to the power of %6.3f = %20.10f. CPU: %6.3f to the power of %6.3f = %20.10f.\n", h_A[k], h_B[k], h_C[k], h_A[k], h_B[k], CPUpow[k] );
	printf("%20.10f - %20.10f = %.10f > 1e-9\n", h_C[k] , CPUpow[k], h_C[k] - CPUpow[k]);
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
      Double_t outputC;
    };
    errors_t errors;


    // create a new ROOT file
    TFile *f = new TFile("errors.root","RECREATE");
    // create a TTree
    TTree *tree = new TTree("T","Absolute und relative Errors");
    // create one branch with all information from the stucture
    tree->Branch("Errors",&errors.errabs,"errabs/D:errrel:inputA:inputB:outputC");
    // fill the tree 
    for (int i=0;i<N;i++) {
      errors.errabs = h_C[i] - CPUpow[i];
      errors.errrel = 1 - (h_C[i] / CPUpow[i]);
      errors.inputA = h_A[i];
      errors.inputB = h_B[i];
      errors.outputC = h_C[i];
      tree->Fill();
    }
    // check what the tree looks like
    tree->Print();
    f->Write();

    
    // Free host memory
        free(h_A);
        free(h_B);
        free(h_C);
    
    
    // Free device memory
        cudaFree(d_A);
        cudaFree(d_B);
        cudaFree(d_C);
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

// Timestamp
void waiter ( int seconds ) // ja kellner
{
  clock_t endwait;
  endwait = clock () + seconds * CLOCKS_PER_SEC ;
  while (clock() < endwait) {}
}


