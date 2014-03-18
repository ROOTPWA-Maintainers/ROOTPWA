#include <math.h>
#include <stdio.h>
#include <cuda.h>
#include <helper_cuda.h>
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
#include <typeinfo>
#include "templtest_kernel.cu"
#include "complex.cuh"
// #include "clebschGordanCoeff.hpp"



using namespace rpwa;
// #include "sharedmem.cuh" // later on for shared memory
// #define N 10



/////////////////////////////
// Declaration of functions//
/////////////////////////////
template <class T>
void tmpltest(int);

// Allocates an array with random float entries.
template<class T>
T RandomInit(T* data, int n)
{
    srand( time(NULL) );
    for (int i = 0; i < n; ++i) {
      data[i] = (rand()%10000+1)/1239. ;
    }
    return 0;
}
 




int main (int argc, char *argv[]) {
 
    int N = atoi( argv[1] );
  
    printf("\nTesting template for %lu doubles on GPU...\n", N);
    tmpltest<double>(N);

    printf("\nTesting template for %lu floats on GPU...\n", N);
    tmpltest<float>(N);
    
  
}

template<class T>
void
tmpltest(int N) 
{

    // Variables  
    T* h_A; 
    T* h_B; 
    T* h_C;
//    cuda::complex<T>* h_cmplxC; 
    T* d_A; 
    T* d_B; 
    T* d_C;
//    cuda::complex<T>* d_cmplxC; 

    T* CPUresult;
//    cuda::complex<T>* CPUcmplxresult;
    
    
    // Allocating memory for array on host to avoid segmentation fault
    CPUresult = (T *)malloc(N*sizeof(T));
//    CPUcmplxresult = (cuda::complex<T> *)malloc(N*sizeof(cuda::complex<double>));
    if(NULL == CPUresult) {
      printf("Error at malloc...?\n");
    }
    
    // Timestamps
    float t_htd;	// host to device
    float t_dth;	// device to host
    float t_cod;	// calculation on device
    
    // Set cuda event for measuring time and bandwidth
    cudaEvent_t start, stop;
    checkCudaErrors( cudaEventCreate(&start) );
    checkCudaErrors( cudaEventCreate(&stop)  );

    size_t size = N * sizeof(T); 	// used for allocating memory on host and device
    size_t sdoc = sizeof(T);            // used for distinguishing from doubles, floats, ...


    // Allocate masses and width in host memory
//     h_m = (T*)malloc(size);
//     if (h_m == 0) free(h_m);
//     h_q = (T*)malloc(size);
//     if (h_q == 0) free(h_q);
//     h_BW = (cuda::complex<T>*)malloc(size);
//     if (h_BW == 0) free(h_BW);    
    
    // Testing page-locked memory here to speed up down and uploading; allocating values in host memory
    checkCudaErrors( cudaHostAlloc((void**)&h_A, size, cudaHostAllocDefault) );
    checkCudaErrors( cudaHostAlloc((void**)&h_B, size, cudaHostAllocDefault) ); 
    checkCudaErrors( cudaHostAlloc((void**)&h_C, size, cudaHostAllocDefault) );    
    
    // Initialize randomnumbers
    gRandom->SetSeed(111990);
    RandomInit(h_A,N);
    RandomInit(h_B,N);
  
    // Allocate randomnumbers in device memory
    checkCudaErrors( cudaMalloc((void**)&d_A, size) );
    checkCudaErrors( cudaMalloc((void**)&d_B, size) );
    checkCudaErrors( cudaMalloc((void**)&d_C, size) );    
  
    // Copy randomnumbers from host memory to device memory
    checkCudaErrors( cudaEventRecord( start, 0) );    
    checkCudaErrors( cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaEventRecord( stop, 0) );
    checkCudaErrors( cudaEventSynchronize( stop ) );
    checkCudaErrors( cudaEventElapsedTime( &t_htd, start, stop));
    printf("%.3f ms required for uploading\n",t_htd);
    printf("%f GB/s Upload\n", ((2*size)/(t_htd/1e3))/(1024*1024*1024) );    
    
    // Invoke kernel 
    checkCudaErrors( cudaEventRecord( start, 0) );    
    int threadsPerBlock = 512;
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
    testPow<T><<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, N);
    checkCudaErrors( cudaThreadSynchronize() );
    checkCudaErrors( cudaEventRecord( stop, 0) );
    checkCudaErrors( cudaEventSynchronize( stop ) );
    checkCudaErrors( cudaEventElapsedTime( &t_cod, start, stop));
    printf("%.3f ms required for calculation on device\n", t_cod );
    printf("Bandwidth: %f GB/s\n", ((size*3)/(t_cod/1e3))/(1024*1024*1024) );    ///  bandwidth = (bytes_read + bytes_written) / time 
    
    // Copy result from device memory to host memory
    checkCudaErrors( cudaEventRecord( start, 0) );
    checkCudaErrors( cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaEventRecord( stop, 0) );
    checkCudaErrors( cudaEventSynchronize( stop ) );
    checkCudaErrors( cudaEventElapsedTime( &t_dth, start, stop));
    printf("%.3f ms required for download\n", t_dth ); 
    printf("%f GB/s Download\n", ((size)/(t_dth/1e3))/(1024*1024*1024) ); // because of complex numbers doubles the size of memory
    
    
///   Human Testing here ///   
//                long unsigned int r;
//               for (r=0;r<N;r++) { 
//         	printf("mass: %f\n", h_m[r]);
//         	printf("q-value: %f\n", h_q[r]);	
//      	printf("L: %i\n", L);	
//           	std::cout << "GPU: " <<  h_C[r] << std::endl; 
//            	std::cout << "CPU: " <<  (m0*Gamma0)/(m0*m0 - h_m[r]*h_m[r] - imag*m0*(Gamma0 * (m0 / h_m[r]) * (h_q[r] / q0) * ( BF(L,h_q[r])*BF(L,h_q[r]) / ( BF(L,q0)*BF(L,q0) ) ) ) ) << std::endl; 
//		std::cout << "CPU: " << h_A[r] / h_B[r] << std::endl;
//                 }
       
    // Calculation on host
    printf("\nNow calculation the same on host...\n");
    timespec h_cal_beg_bw, h_cal_end_bw; 
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &h_cal_beg_bw);     
//    for (long unsigned int j = 0; j<N; j++) { CPUcmplxresult[j] = (m0*Gamma0)/(m0*m0 - h_m[j]*h_m[j] - imag*m0*(Gamma0 * (m0 / h_m[j]) * (h_q[j] / q0) * ( BF(L,h_q[j])*BF(L,h_q[j]) / ( BF(L,q0)*BF(L,q0) ) ) ) ) ;   }
    for (long unsigned int j = 0; j<N; j++) { CPUresult[j] = pow(h_A[j],h_B[j]) ;   }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &h_cal_end_bw); 
    printf("%.3f ms required for calculation on host\n", ( ( (h_cal_end_bw.tv_sec - h_cal_beg_bw.tv_sec) ) + (h_cal_end_bw.tv_nsec - h_cal_beg_bw.tv_nsec) / 1e9 ) * 1e3 );     
    printf("Pure calculation speed up: GPU is %3.1f times faster than the CPU.\n",  ( ( ( (h_cal_end_bw.tv_sec - h_cal_beg_bw.tv_sec) ) + (h_cal_end_bw.tv_nsec - h_cal_beg_bw.tv_nsec) / 1e9 ) * 1e3 ) / (t_cod) );
    printf("Over-all speed up:         GPU is %3.1f times faster than the CPU.\n\n",  ( ( ( (h_cal_end_bw.tv_sec - h_cal_beg_bw.tv_sec) ) + (h_cal_end_bw.tv_nsec - h_cal_beg_bw.tv_nsec) / 1e9 ) * 1e3 ) / (t_cod+t_htd+t_dth) );
   
    // Verify result 
    printf("Verifying results for bw...\n");
    long unsigned int l = 0;
    for (l=0;l<N;l++) {
      if (abs(h_C[l] - CPUresult[l]) > 1e-9) {
	printf("%li. pair of variates incorrect!\n",l);
	std::cout << "GPU: " << h_C[l] << " CPU: " << CPUresult[l] << std::endl;
	std::cout << "absolute error (complex absolute...): " << abs( h_C[l] - CPUresult[l] )  << std::endl;
	break;
      }
    }
    printf("TEST %s \n\n", (l == N) ? "PASSED" : "FAILED");
    
    
    
    // Visualize results via root
    // Create a root file and a TTree
    // Structure to hold the variables for the branch
    struct errors_t {
      Double_t abs_d; // absolute errors for doubles
      Double_t rel_d; // relative errors for doubles
      Double_t A_d;
      Double_t B_d;
      Double_t C_d;
      Double_t abs_f; // absolute errors for floats
      Double_t rel_f; // relative errors for floats
      Double_t A_f;
      Double_t B_f;       
      Double_t C_f;
     
    };
    errors_t errors;
    

    
    // Distinction of cases for doubles(8) and floats(4)
    if (sdoc==8) {
          // Create a new ROOT file
	  TFile *f = new TFile("errorsdoubletest.root","RECREATE");
	  // Create a TTree
	  TTree *tree = new TTree("T","Absolute und relative errors");
	 // Create one branch with all information from the stucture
 	 tree->Branch("ErrorsD",&errors.abs_d,"abs_d/D:rel_d:A_d:B_d:C_d");
 	 // Fill the tree 
 	 for (long unsigned int i=0;i<N;i++) {
	   errors.A_d = h_A[i];
	   errors.B_d = h_B[i];
	   errors.C_d = h_C[i];	   
	   errors.abs_d = h_C[i] - CPUresult[i];
 	   errors.rel_d = 1 - abs(h_C[i] / CPUresult[i]);
	   tree->Fill();
 	 }     
       	 // Check what the tree looks like
	 tree->Print();
	 f->Write();
    }
    else if (sdoc==4) {
          // Create a new ROOT file
	  TFile *f = new TFile("errorsfloattest.root","RECREATE");
	  // Create a TTree
	  TTree *tree = new TTree("T","Absolute und relative errors");
   	 tree->Branch("ErrorsF",&errors.abs_f,"abs_f/D:rel_f:A_f:B_f:C_f");
 	 for (long unsigned int i=0;i<N;i++) {
	   errors.A_f = h_A[i];
	   errors.B_f = h_B[i];
	   errors.C_f = h_C[i];	   
	   errors.abs_f = h_C[i] - CPUresult[i];
 	   errors.rel_f = 1 - abs(h_C[i] / CPUresult[i]);
	   tree->Fill();
 	 }     
	 tree->Print();
	 f->Write();
    }
    else { printf("default bla bla"); }
		
    
    // Free host memory
//     free(h_m);
//     free(h_q);
//     free(h_BW);
    
    // Freeing allocated host memory from cudaHostAlloc
    checkCudaErrors( cudaFreeHost(h_A) );
    checkCudaErrors( cudaFreeHost(h_B) );    
    checkCudaErrors( cudaFreeHost(h_C) );    
    
    // Free device memory
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
    checkCudaErrors( cudaEventDestroy(start) );
    checkCudaErrors( cudaEventDestroy(stop) );
    
}
  
  


