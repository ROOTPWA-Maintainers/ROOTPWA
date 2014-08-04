#include <stdio.h>
#include <complex>
#include "complex.cuh"
#include <cuda.h>
#include <helper_cuda.h>
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
#include <time.h>

#include "helampltest.cuh"

using namespace std;
using namespace rpwa;

// !NOTE! spins and projection quantum numbers are in units of hbar/2

namespace cupwa {
	
	struct vertexdata {
		int JX, Lambda, J1, L1, S1, J2, L2, S2;
		double theta1, phi1, theta2, phi2, wX, wf2, massf2, massX, wpi, qX, qf2;
	};
  
	//Kernel for cuPWA using struct
	__global__ void
	GPUAmp2(vertexdata* data, cuda::complex<double>* result_GPUAmp, int N)
	{

		const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
		// if(idx<N) { result_GPUAmp[idx] = 0;}
		if(idx<N) { 
			result_GPUAmp[idx] = cuHelAmplitude(data[idx].theta1,data[idx].phi1,data[idx].theta2,data[idx].phi2,data[idx].wX,data[idx].wf2,data[idx].qX,data[idx].qf2,data[idx].JX,data[idx].Lambda,data[idx].J1,data[idx].L1,data[idx].S1,data[idx].J2,data[idx].L2,data[idx].S2); 
		}
		
	}
	
	void DMAccess(vertexdata* data, int N)
	{
		printf("address0 before: %p\n", (void*) &data[0].JX);
		checkCudaErrors( cudaHostAlloc((void**)&data, N*sizeof(vertexdata), cudaHostAllocDefault) );
		printf("address0 after: %p\n", (void*) &data[0].JX);    
	}
	
	void GPUAmpcall2(vertexdata* data,int N, const vector<complex<double> >& cpudata)
	{
	
		printf("\nNow calculating on device...\n");
		timespec bartime_beg, bartime_end; 
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &bartime_beg);
		cudaEvent_t start, stop;
		checkCudaErrors( cudaEventCreate(&start) );
		checkCudaErrors( cudaEventCreate(&stop)  );
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &bartime_end); 
		printf("\n\n%.3f ms required for creating cuda-events (timespec)\n", ( ( (bartime_end.tv_sec - bartime_beg.tv_sec) ) + (bartime_end.tv_nsec - bartime_beg.tv_nsec) / 1e9 ) * 1e3 ); 
		
		float timer[10];
		size_t sdata = sizeof(vertexdata);
		size_t sres  = sizeof(cuda::complex<double>);
		cuda::complex<double>* result_GPUAmp;
		cuda::complex<double>* dev_result_GPUAmp;
		vertexdata* dev_data; 
		
		checkCudaErrors( cudaEventRecord( start, 0) );
		timespec footime_beg, footime_end; 
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &footime_beg); 
		checkCudaErrors( cudaHostAlloc((void**)&result_GPUAmp, N*sres, cudaHostAllocDefault) );
		// checkCudaErrors( cudaHostAlloc((void**)&data, N*sdata, cudaHostAllocDefault) ); //doesn't work properly
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &footime_end); 
		printf("\n\n%.3f ms required for allocating memory (timespec)\n", ( ( (footime_end.tv_sec - footime_beg.tv_sec) ) + (footime_end.tv_nsec - footime_beg.tv_nsec) / 1e9 ) * 1e3 ); 
		checkCudaErrors( cudaEventRecord( stop, 0) );
		checkCudaErrors( cudaEventSynchronize( stop ) );
		checkCudaErrors( cudaEventElapsedTime( &timer[0], start, stop));
		printf("%.3f ms required for allocating memory\n",timer[0]);
		
		checkCudaErrors( cudaMalloc((void**)&dev_result_GPUAmp, N*sres) );
		checkCudaErrors( cudaMalloc((void**)&dev_data, N*sdata) );
		
		checkCudaErrors( cudaEventRecord( start, 0) );
		checkCudaErrors( cudaMemcpy(dev_data, data, N*sdata, cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaEventRecord( stop, 0) );
		checkCudaErrors( cudaEventSynchronize( stop ) );
		checkCudaErrors( cudaEventElapsedTime( &timer[1], start, stop));
		printf("%.3f ms required for uploading\n",timer[1]);
		printf("%f GB/s bandwidth for upload\n",(N*(sdata/(timer[1]/1e3))/(1024*1024*1024)));
		
		int threadsPerBlock = 512;
		int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;    
		checkCudaErrors( cudaEventRecord( start, 0) );
		GPUAmp2<<<blocksPerGrid, threadsPerBlock>>>(dev_data,dev_result_GPUAmp,N);
		checkCudaErrors( cudaThreadSynchronize() );
		checkCudaErrors( cudaEventRecord( stop, 0) );
		checkCudaErrors( cudaEventSynchronize( stop ) );
		checkCudaErrors( cudaEventElapsedTime( &timer[2], start, stop));
		printf("%.3f ms required for calculating on device\n",timer[2]);
		printf("%f GB/s Bandwidth\n", (N*(sdata+sres)/(timer[2]/1e3))/(1024*1024*1024) );
		checkCudaErrors( cudaEventRecord( start, 0) );
		checkCudaErrors( cudaMemcpy(result_GPUAmp, dev_result_GPUAmp, N*sres, cudaMemcpyDeviceToHost) );
		checkCudaErrors( cudaEventRecord( stop, 0) );
		checkCudaErrors( cudaEventSynchronize( stop ) );
		checkCudaErrors( cudaEventElapsedTime( &timer[3], start, stop));
		printf("%.3f ms required for downloading\n",timer[3]);
		printf("%f GB/s bandwidth for download\n",(N*(sres/(timer[3]/1e3))/(1024*1024*1024)));
		printf("%.3f ms required for GPU call without malloc stuff, etc.\n",timer[1]+timer[2]+timer[3]);
		
		//for(int e=0;e<N;e++) {
		//	printf("\nGPUAmp[%d]: %2.20f + i*(%2.20f)",e, result_GPUAmp[e].real(), result_GPUAmp[e].imag() );
		//	printf("\nCPUAmp[%d]: %2.20f + i*(%2.20f)",e, cpudata[e].real(), cpudata[e].imag() );
		//}
		
		timespec roottime_beg, roottime_end; 
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &roottime_beg); 
		checkCudaErrors( cudaEventRecord( start, 0) );
		// Visualize results via root
		// Create a root file and a TTree
		// Structure to hold the variables for the branch
		struct output_t {
		Double_t mass;
		Double_t amplimag;
		Double_t amplreal;
		Double_t amplabs;
		Double_t errrelimag;
		Double_t errrelreal;
		Double_t errabsimag;
		Double_t errabsreal;
		Double_t cpuimag;
		Double_t cpureal;
		Double_t cpuabs;
		};
		output_t output;
		
		// Create a new ROOT file
		TFile *f = new TFile("Xtopipipi.root","RECREATE");
		// Create a TTree
		TTree *tree = new TTree("T","output");
		// Create one branch with all information from the stucture
		tree->Branch("Output",&output.mass,"mass/D:amplimag:amplreal:amplabs:errrelimag:errrelreal:errabsimag:errabsreal:cpuimag:cpureal:cpuabs");
		// Fill the tree 
		for (int i=0;i<N;i++) {
		output.mass 	= data[i].wX;
		output.amplimag 	= result_GPUAmp[i].imag();
		output.amplreal 	= result_GPUAmp[i].real();
		output.amplabs  	= abs(result_GPUAmp[i]);
		output.errrelimag = 1 - (result_GPUAmp[i].imag() / cpudata[i].imag());
		output.errrelreal = 1 - (result_GPUAmp[i].real() / cpudata[i].real());
		output.errabsimag = result_GPUAmp[i].imag() - cpudata[i].imag();
		output.errabsreal = result_GPUAmp[i].real() - cpudata[i].real();
		output.cpuimag	= cpudata[i].imag();
		output.cpureal	= cpudata[i].real();
		output.cpuabs 	= abs(cpudata[i]);
		tree->Fill();
		}     
		// Check what the tree looks like
		tree->Print();
		f->Write();
		checkCudaErrors( cudaEventRecord( stop, 0) );
		checkCudaErrors( cudaEventSynchronize( stop ) );
		checkCudaErrors( cudaEventElapsedTime( &timer[4], start, stop));
		printf("%.3f ms required for root stuff\n",timer[4]);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &roottime_end); 
		printf("\n\n%.3f ms required for root stuff (time spec)\n", ( ( (roottime_end.tv_sec - roottime_beg.tv_sec) ) + (roottime_end.tv_nsec - roottime_beg.tv_nsec) / 1e9 ) * 1e3 ); 
		
		checkCudaErrors( cudaFreeHost(result_GPUAmp) );
		// checkCudaErrors( cudaFreeHost(data) );
		cudaFree(dev_result_GPUAmp);
		cudaFree(dev_data);
		
		checkCudaErrors( cudaEventDestroy(start) );
		checkCudaErrors( cudaEventDestroy(stop) );
		
	}
	
}
