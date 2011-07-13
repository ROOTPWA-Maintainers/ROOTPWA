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
 

// Random mass and breakup-momenta entries for Breit-Wigner
template<class T>
T RanMassInit(T* mass, T* bum, int n)
{
  for (int i=0;i<n;++i) {
    mass[i] = gRandom->Rndm()+0.2; //    0.75549;				/// rho(770)
//     T scale = gRandom->Rndm()*0.00002+0.49999;
    bum[i] = sqrt((mass[i]*mass[i]-(0.18*mass[i]+0.18*mass[i])*(0.18*mass[i]+0.18*mass[i]))*(mass[i]*mass[i]-(0.18*mass[i]-0.18*mass[i])*(0.18*mass[i]-0.18*mass[i]))) / (2*mass[i]) ; /// pi+pi-, m(pi)/m(rho) = 0.17997...
  }
  return 0;
}


// Random q (mass dependend) entries for Breit-Wigner
template<class T>
T RanBumInit(T* bum, int n)
{
  for (int i=0;i<n;++i) {
    bum[i] = 0.3 ;
  }
  return 0;
}



// Blatt-Weisskopf Barrierfunction, i failed to implement... somehow.. so i added it here extra. templated too :D
template<class T>
T BF( int L, T q)
{
// #define Pr 0.1973 // Gev/c corresponds to 1 fermi
  double ret;
  double z = (q/Pr) * (q/Pr);
  int m = L/2;
  switch (m) {
  case 0:
    ret = 1.0;
    break;
  case 1:
    ret = sqrt( (2.0 * z)/(z + 1));
    break;
  case 2:
    ret = sqrt( (13.0 * z * z)/(pow(z - 3.0,2.0) + 9.0 * z));
    break;
  case 3:
    ret = sqrt( (277.0 * pow(z,3.0))/(z * pow(z - 15.0,2.0) + 9.0 * pow(2.0 * z - 5.0,2.0)));
    break;
  case 4:
    ret = sqrt( (12746.0 * pow(z,4.0))/( pow(z * z - 45.0 * z + 105.0,2.0) + 25.0 * z * pow(2.0 * z - 21.0,2.0)));
    break;
  case 5:
    ret =  sqrt( (998881.0 * z*z*z*z*z)/(893025.0 +99225.0*z +6300.0*z*z +315.0*z*z*z +15.0*z*z*z*z +z*z*z*z*z));
    break;
  case 6:
    ret =  sqrt( (118394977.0 * z*z*z*z*z*z)/(108056025.0 + z*(9823275.0 + z*(496125.0 + z*(18900.0 + z*(630.0 + z*(21.0 + z))))))) ;
    break;
  case 7:
    ret =  sqrt( (19727003738.0 * z*z*z*z*z*z*z)/((T)18261468225.0 + z*((T)1404728325.0 + z*((T)58939650.0 + z*((T)1819125.0 + z*((T)47250.0 + z*((T)1134.0 + z*((T)28.0 + z)))))))) ;
    break;
  default:
    std::cerr << "Blatt-Weisskopf called for undefined L = " << L/2 << std::endl;
    ret = 1.0;
    break;
  }
  return ret;
}



int main () {
  
//     double CGC = clebschGordanCoeff(0.5,0.5,0.5,0.5,1,0);
//     printf("bla %f\n", CGC);
    printf("\nTesting template for doubles: BW-Gamma for %lu randomentries on GPU...\n", N);
    tmpltest<double>();

//   tmpltest<int>(); // insert casts here for integers?

    printf("\nTesting template for floats: BW-Gamma for %lu randomnumbers on GPU...\n", N);
    tmpltest<float>();
    
  
}

template<class T>
void
tmpltest() 
{

    // Variables  
    T* h_m; 				// masses for BW
    T* h_q; 				// breakup momentum for BW
    cuda::complex<T>* h_BW; 		// complex results for BW
    T* d_m; 
    T* d_q; 
    cuda::complex<T>* d_BW; 
    cuda::complex<double>* CPUbw;
    
    
    // Allocating memory for array on host to avoid segmentation fault
    CPUbw = (cuda::complex<double> *)malloc(2*N*sizeof(double));
    if(NULL == CPUbw) {
      printf("Error at malloc(gamma)...?\n");
    }
    
    // Timestamps
    float t_htd;	// host to device
    float t_dth;	// device to host
    float t_cod;	// calculation on device
    
    // Set cuda event for measuring time and bandwidth
    cudaEvent_t start, stop;
    cutilSafeCall( cudaEventCreate(&start) );
    cutilSafeCall( cudaEventCreate(&stop)  );

    size_t size 	= N * sizeof(T); 	// used for allocating memory on host and device
    size_t scmplx 	= N * sizeof(cuda::complex<T>);
    size_t sdoc 	= sizeof(T); 		// used for distinguishing from doubles, floats, ...
//    size_t sint 	= sizeof(int);		// used for L
  
    // Allocate masses and width in host memory
//     h_m = (T*)malloc(size);
//     if (h_m == 0) free(h_m);
//     h_q = (T*)malloc(size);
//     if (h_q == 0) free(h_q);
//     h_BW = (cuda::complex<T>*)malloc(size);
//     if (h_BW == 0) free(h_BW);    
    
    // Testing page-locked memory here to speed up down and uploading; allocating values in host memory
    cutilSafeCall( cudaHostAlloc((void**)&h_m, size, cudaHostAllocDefault) );
    cutilSafeCall( cudaHostAlloc((void**)&h_q, size, cudaHostAllocDefault) ); 
    cutilSafeCall( cudaHostAlloc((void**)&h_BW, scmplx, cudaHostAllocDefault) );    
    
    // Initialize masses, breakup momentum
    gRandom->SetSeed(111990);
    RanMassInit(h_m, h_q, N);
//     RanBumInit(h_q, N);
  
    // Allocate masses and width in device memory
    cutilSafeCall( cudaMalloc((void**)&d_m, size) );
    cutilSafeCall( cudaMalloc((void**)&d_q, size) );
    cutilSafeCall( cudaMalloc((void**)&d_BW, scmplx) );    
  
    // Copy masses and width from host memory to device memory
    cutilSafeCall( cudaEventRecord( start, 0) );    
    cutilSafeCall( cudaMemcpy(d_m, h_m, size, cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(d_q, h_q, size, cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaEventRecord( stop, 0) );
    cutilSafeCall( cudaEventSynchronize( stop ) );
    cutilSafeCall( cudaEventElapsedTime( &t_htd, start, stop));
    printf("%.3f ms required for uploading\n",t_htd);
    printf("%f GB/s Upload\n", ((2*size)/(t_htd/1e3))/(1024*1024*1024) );    
    
    // Invoke kernel (BW-function)
    cutilSafeCall( cudaEventRecord( start, 0) );    
    int threadsPerBlock = 512;
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
    imagtest<T><<<blocksPerGrid, threadsPerBlock>>>(d_m, d_q, d_BW);
    cutilSafeCall( cudaThreadSynchronize() );
    cutilSafeCall( cudaEventRecord( stop, 0) );
    cutilSafeCall( cudaEventSynchronize( stop ) );
    cutilSafeCall( cudaEventElapsedTime( &t_cod, start, stop));
    printf("%.3f ms required for calculation on device\n", t_cod );
    printf("Bandwidth: %f GB/s\n", ((size*2+scmplx)/(t_cod/1e3))/(1024*1024*1024) );    ///  bandwidth = (bytes_read + bytes_written) / time 
    
    // Copy result from device memory to host memory
    cutilSafeCall( cudaEventRecord( start, 0) );
    cutilSafeCall( cudaMemcpy(h_BW, d_BW, scmplx, cudaMemcpyDeviceToHost) );
    cutilSafeCall( cudaEventRecord( stop, 0) );
    cutilSafeCall( cudaEventSynchronize( stop ) );
    cutilSafeCall( cudaEventElapsedTime( &t_dth, start, stop));
    printf("%.3f ms required for download\n", t_dth ); 
    printf("%f GB/s Download\n", ((scmplx)/(t_dth/1e3))/(1024*1024*1024) ); // because of complex numbers doubles the size of memory
    
    
///   Human Testing here ///   
//                 long unsigned int r;
//                  for (r=0;r<N;r++) { 
//         	printf("mass: %f\n", h_m[r]);
//         	printf("q-value: %f\n", h_q[r]);	
//    //      	printf("L: %i\n", L);	
//            	std::cout << "GPU: " <<  h_BW[r] << std::endl; 
//            	std::cout << "CPU: " <<  (m0*Gamma0)/(m0*m0 - h_m[r]*h_m[r] - imag*m0*(Gamma0 * (m0 / h_m[r]) * (h_q[r] / q0) * ( BF(L,h_q[r])*BF(L,h_q[r]) / ( BF(L,q0)*BF(L,q0) ) ) ) ) << std::endl; 
//                  }
       
    // Calculation on host
    printf("\nNow calculation the same on host...\n");
    timespec h_cal_beg_bw, h_cal_end_bw; 
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &h_cal_beg_bw);     
    for (long unsigned int j = 0; j<N; j++) { CPUbw[j] = (m0*Gamma0)/(m0*m0 - h_m[j]*h_m[j] - imag*m0*(Gamma0 * (m0 / h_m[j]) * (h_q[j] / q0) * ( BF(L,h_q[j])*BF(L,h_q[j]) / ( BF(L,q0)*BF(L,q0) ) ) ) ) ;   }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &h_cal_end_bw); 
    printf("%.3f ms required for calculation on host\n", ( ( (h_cal_end_bw.tv_sec - h_cal_beg_bw.tv_sec) ) + (h_cal_end_bw.tv_nsec - h_cal_beg_bw.tv_nsec) / 1e9 ) * 1e3 );     
    printf("Pure calculation speed up: GPU is %3.1f times faster than the CPU.\n",  ( ( ( (h_cal_end_bw.tv_sec - h_cal_beg_bw.tv_sec) ) + (h_cal_end_bw.tv_nsec - h_cal_beg_bw.tv_nsec) / 1e9 ) * 1e3 ) / (t_cod) );
    printf("Over-all speed up:         GPU is %3.1f times faster than the CPU.\n\n",  ( ( ( (h_cal_end_bw.tv_sec - h_cal_beg_bw.tv_sec) ) + (h_cal_end_bw.tv_nsec - h_cal_beg_bw.tv_nsec) / 1e9 ) * 1e3 ) / (t_cod+t_htd+t_dth) );
   
    // Verify result BW
    printf("Verifying results for bw...\n");
    long unsigned int l = 0;
    for (l=0;l<N;l++) {
      if (abs(h_BW[l] - (cuda::complex<T>)CPUbw[l]) > 1e-9) {
	printf("%li. pair of variates incorrect!\n",l);
	std::cout << "GPU: " << h_BW[l] << " CPU: " << CPUbw[l] << std::endl;
	std::cout << "absolute error (complex absolute...): " << abs( h_BW[l] - (cuda::complex<T>)CPUbw[l] )  << std::endl;
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
      Double_t q_d;
      Double_t masses_d;
      Double_t bwdata_d;
      Double_t abs_f; // absolute errors for floats
      Double_t rel_f; // relative errors for floats
      Double_t q_f;
      Double_t masses_f;       
      Double_t bwdata_f;
     
    };
    errors_t errors;
    

    
    // Distinction of cases for doubles(8) and floats(4)
    if (sdoc==8) {
          // Create a new ROOT file
	  TFile *f = new TFile("errorsdoublebw.root","UPDATE");
	  // Create a TTree
	  TTree *tree = new TTree("T","Absolute und relative errors");
	 // Create one branch with all information from the stucture
 	 tree->Branch("ErrorsD",&errors.abs_d,"abs_d/D:rel_d:q_d:masses_d:bwdata_d");
 	 // Fill the tree 
 	 for (long unsigned int i=0;i<N;i++) {
	   errors.q_d = h_q[i];
	   errors.masses_d = h_m[i];
	   errors.bwdata_d = abs(h_BW[i]);	   
	   errors.abs_d = abs(h_BW[i] - (cuda::complex<T>)CPUbw[i]);
 	   errors.rel_d = 1 - abs(h_BW[i] / (cuda::complex<T>)CPUbw[i]);
	   tree->Fill();
 	 }     
       	 // Check what the tree looks like
	 tree->Print();
	 f->Write();
    }
    else if (sdoc==4) {
          // Create a new ROOT file
	  TFile *f = new TFile("errorsfloatbw.root","UPDATE");
	  // Create a TTree
	  TTree *tree = new TTree("T","Absolute und relative errors");
   	 tree->Branch("ErrorsF",&errors.abs_f,"abs_f/D:rel_f:q_f:masses_f:bwdata_f");
 	 for (long unsigned int i=0;i<N;i++) {
	   errors.q_f = h_q[i];
	   errors.masses_f = h_m[i];
	   errors.bwdata_f = abs(h_BW[i]);	   
	   errors.abs_f = abs(h_BW[i] - (cuda::complex<T>)CPUbw[i]);
 	   errors.rel_f = 1 - abs(h_BW[i] / (cuda::complex<T>)CPUbw[i]);
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
    cutilSafeCall( cudaFreeHost(h_m) );
    cutilSafeCall( cudaFreeHost(h_q) );    
    cutilSafeCall( cudaFreeHost(h_BW) );    
    
    // Free device memory
    cudaFree(d_m);
    cudaFree(d_q);
    cudaFree(d_BW);
    cutilSafeCall( cudaEventDestroy(start) );
    cutilSafeCall( cudaEventDestroy(stop) );
    
}
  
  


