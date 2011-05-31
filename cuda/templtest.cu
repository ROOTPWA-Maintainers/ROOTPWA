// function template
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
#include <typeinfo>
#include "templtest_kernel.cu"
using namespace std;
using namespace rpwa;

// #include "sharedmem.cuh" // later on for shared memory
// #define N 10



/////////////////////////////
// Declaration of functions//
/////////////////////////////
template <class T>
void tmpltest();

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

int main () {
  
  printf("\nTesting template for doubles: Exponentiating %lu randomnumbers on GPU...\n", N);
  tmpltest<double>();
  
//   tmpltest<int>(); // insert casts here for integers?

//   printf("\nTesting template for floats: Exponentiating %lu randomnumbers on GPU...\n", N);
//   tmpltest<float>();
  
}

template<class T>
void
tmpltest() 
{

    // Variables
    T* h_A;
    T* h_B;
    cuda::complex<T>* h_C;
    T* d_A;
    T* d_B;
    cuda::complex<T>* d_C;
    cuda::complex<double> CPUpow[N];   
    
    // Allocating memory for array on host to avoid segmentation fault
//     CPUpow = (double *)malloc(N*sizeof(double));
//     if(NULL == CPUpow) {
//       printf("Error at malloc(pow)...?\n");
//     }
    
    // Timestamps
    float t_htd;	// host to device
    float t_dth;	// device to host
    float t_cod;	// calculation on device
    
    // Set cuda event for measuring time and bandwidth
    cudaEvent_t start, stop;
    cutilSafeCall( cudaEventCreate(&start) );
    cutilSafeCall( cudaEventCreate(&stop)  );

    size_t size = N * sizeof(T); 	// used for allocating memory on host and device
//     size_t sdoc = sizeof(T); 		// used for distinguishing from doubles, floats, ...
  
    // Allocate randomnumbers h_A and h_B in host memory
    h_A = (T*)malloc(size);
    if (h_A == 0) free(h_A);
    h_B = (T*)malloc(size);
    if (h_B == 0) free(h_B);
    h_C = (cuda::complex<T>*)malloc(2*size);
    if (h_C == 0) free(h_C);
    
    // Initialize randomnumbers
    RandomInit(h_A, N); 
    RandomInit(h_B, N);
  
    // Allocate randomnumbers in device memory
    cutilSafeCall( cudaMalloc((void**)&d_A, size) );
    cutilSafeCall( cudaMalloc((void**)&d_B, size) );
    cutilSafeCall( cudaMalloc((void**)&d_C, 2*size) );
//     cutilSafeCall( cudaMalloc((void**)&d_D, size) );
  
    // Copy randomnumbers from host memory to device memory
    cutilSafeCall( cudaEventRecord( start, 0) );    
    cutilSafeCall( cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaEventRecord( stop, 0) );
    cutilSafeCall( cudaEventSynchronize( stop ) );
    cutilSafeCall( cudaEventElapsedTime( &t_htd, start, stop));
    printf("%.3f ms required for uploading\n",t_htd);
    printf("%f GB/s Upload\n", ((2*size)/(t_htd/1e3))/(1024*1024*1024) );       
    
    // Invoke kernel (pow-function)
    cutilSafeCall( cudaEventRecord( start, 0) );    
    int threadsPerBlock = 512;
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
    testPow<T><<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C);
    cutilSafeCall( cudaThreadSynchronize() );
    cutilSafeCall( cudaEventRecord( stop, 0) );
    cutilSafeCall( cudaEventSynchronize( stop ) );
    cutilSafeCall( cudaEventElapsedTime( &t_cod, start, stop));
    printf("%.3f ms required for calculation on device\n", t_cod );
    printf("???? insert knowledge here for bandwidth on GPU ????\n");  
      
    // Copy result from device memory to host memory
    cutilSafeCall( cudaEventRecord( start, 0) );
//     cutilSafeCall( cudaMemcpy(d_C, d_D, size, cudaMemcpyDeviceToDevice) );    
//     cutilSafeCall( cudaThreadSynchronize() );    
    cutilSafeCall( cudaMemcpy(h_C, d_C, 2*size, cudaMemcpyDeviceToHost) );
    cutilSafeCall( cudaEventRecord( stop, 0) );
    cutilSafeCall( cudaEventSynchronize( stop ) );
    cutilSafeCall( cudaEventElapsedTime( &t_dth, start, stop));
    printf("%.3f ms required for download\n", t_dth ); 
    printf("%f GB/s Download\n", ((size)/(t_dth/1e3))/(1024*1024*1024) );  
    
///      Human Testing here ///
	cuda::complex<T> imag(0.0,1.0);
        long unsigned int r;
        for (r=0;r<N;r++) { 
  	cout << "GPU: " <<  h_C[r] << endl; 
  	cout << "CPU: " <<  pow(h_A[r],h_B[r])*imag + 1.0 << endl; 
        }
       	
       	
//        	cuda::complex<T> bla[N];
// 	for (long unsigned int i=0;i<N;i++) {
//  	bla[i] = cuda::complex<T>(1.0,i);
// 	}
// 	cout << "imaginaries test: i = " << imag << endl;
// 	for (long unsigned int i=0;i<N;i++) {
// 	cout << "imaginaries test: 1 + " << i << "i = " << bla[i] << endl;
// 	}
          
    // Calculation on host
    printf("\nNow calculation the same on host...\n");
    timespec h_cal_beg, h_cal_end; 
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &h_cal_beg);     
    long unsigned int i;
    for (i = 0; i<N; i++) {      CPUpow[i] = (pow(h_A[i],h_B[i]))*imag+1.0 ;   }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &h_cal_end); 
    printf("%.3f ms required for calculation on host\n\n", ( ( (h_cal_end.tv_sec - h_cal_beg.tv_sec) ) + (h_cal_end.tv_nsec - h_cal_beg.tv_nsec) / 1e9 ) * 1e3 );     
    
    // Verify result pow
    printf("Verifying results for pow...\n");
    long unsigned int k;
    for (k=0;k<N;k++) {

      if ( abs( (h_C[k]) - (CPUpow[k]) )  > 1e-9) {
	cout << k << ". pair of variates incorrect!" << endl;
	cout << "GPU: " << h_C[k] << endl;
	cout << "CPU: " << pow(h_A[k],h_B[k])*imag + 1.0 << endl;
	cout << "absolute error (complex absolute...): " << abs( (h_C[k]) - (CPUpow[k]) )  << endl;
	break;
      }
    }
    printf("TEST %s \n\n", (k == N) ? "PASSED" : "FAILED");
   
    // Visualize results via root
    // Create a root file and a TTree
    // Structure to hold the variables for the branch
//     struct errors_t {
//       Double_t abs_d; // absolute errors for doubles
//       Double_t rel_d; // relative errors for doubles
//       Double_t abs_f; // absolute errors for floats
//       Double_t rel_f; // relative errors for floats
//     };
//     errors_t errors;
    

    
    // Distinction of cases for doubles(8) and floats(4)
//     if (sdoc==8) {
//           // Create a new ROOT file
// 	  TFile *f = new TFile("errorsdouble.root","UPDATE");
// 	  // Create a TTree
// 	  TTree *tree = new TTree("T","Absolute und relative errors");
// 	 // Create one branch with all information from the stucture
//  	 tree->Branch("ErrorsD",&errors.abs_d,"abs_d/D:rel_d");
//  	 // Fill the tree 
//  	 for (long unsigned int i=0;i<N;i++) {
// 	   errors.abs_d = h_C[i] - CPUpow[i];
//  	   errors.rel_d = 1 - (h_C[i] / CPUpow[i]);
// 	   tree->Fill();
//  	 }     
//        	 // Check what the tree looks like
// 	 tree->Print();
// 	 f->Write();
//     }
//     else if (sdoc==4) {
//           // Create a new ROOT file
// 	  TFile *f = new TFile("errorsfloat.root","UPDATE");
// 	  // Create a TTree
// 	  TTree *tree = new TTree("T","Absolute und relative errors");
//    	 tree->Branch("ErrorsF",&errors.abs_f,"abs_f/D:rel_f");
//  	 for (long unsigned int i=0;i<N;i++) {
// 	   errors.abs_f = h_C[i] - CPUpow[i];
//  	   errors.rel_f = 1 - (h_C[i] / CPUpow[i]);
// 	   tree->Fill();
//  	 }     
// 	 tree->Print();
// 	 f->Write();
//     }
//     else { printf("default bla bla"); }
		
    
    // Free host memory
    free(h_A);
    free(h_B);
    free(h_C);
    
    
    // Free device memory
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
//     cudaFree(d_D);
    cutilSafeCall( cudaEventDestroy(start) );
    cutilSafeCall( cudaEventDestroy(stop) );
    
}
  
  


