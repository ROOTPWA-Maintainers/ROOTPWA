#include <stdio.h>
#include <complex>
#include "complex.cuh"
#include <cuda.h>
#include <cutil_inline.h>
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
//#include "helampltestInterface.cuh"
//#include "helampltestKernel.cu"

using namespace std;
using namespace rpwa;
//using namespace cupwa;

#define FOURPI 4*3.1415926

  ///      !NOTE! spins and projection quantum numbers are in units of hbar/2

namespace cupwa {

  
    __device__ double
    cuBreakUpMomentum(const double M,   // mass of mother particle
		      const double m1,  // mass of daughter particle 1
		      const double m2)  // mass of daughter particle 2
    {
      if (M < m1 + m2)
	return 0;
      return sqrt((M - m1 - m2) * (M + m1 + m2) * (M - m1 + m2) * (M + m1 - m2)) / (2 * M);
    }


    // note that for speed up purposes the barrier factor is defined without the squareroot here!
    __device__ double
    BFactor( int L, double q )
    {
      double Pr = 0.1973;
      double z = (q * q)/(Pr * Pr);
      double ret;
      int m = L/2;
      switch (m) {
      case 0:
	ret =  1.0;
	break;
      case 1:
	ret =  ((double)2.0 * z)/(z + 1) ;
	break;
      case 2:
	ret =  ((double)13.0 * z * z)/( (double)9.0 + z*(z+(double)3.0) ) ;								//
	break;
      case 3:
	ret =  ((double)277.0 * z*z*z)/((double)225.0 + z*((double)45.0 + z*((double)6.0 + z))) ;						//
	break;
      case 4:
	ret =  ((double)12746.0 * z*z*z*z)/((double)11025.0 + z*((double)1575.0 + z*((double)135.0 + z*((double)10.0 + z))) ) ; 			//
	break;
      case 5:
	ret =  ((double)998881.0 * z*z*z*z*z)/((double)893025.0 + z*((double)99225.0 + z*((double)6300.0 + z*((double)315.0 + z*((double)15.0 + z))))) ; 			//  
	break;
      case 6:
	ret =  ((double)118394977.0 * z*z*z*z*z*z)/((double)108056025.0 + z*((double)9823275.0 + z*((double)496125.0 + z*((double)18900.0 + z*((double)630.0 + z*((double)21.0 + z)))))) ;
	break;
      case 7:
	ret =  ((double)19727003738.0 * z*z*z*z*z*z*z)/((double)18261468225.0 + z*((double)1404728325.0 + z*((double)58939650.0 + z*((double)1819125.0 + z*((double)47250.0 + z*((double)1134.0 + z*((double)28.0 + z))))))) ;
	break;
      default:
	ret =  1.0;
	break;
      }
      return ret;
    }
    
    
    __device__ int
    cufac( int n )
    {
      int cufac = 1;
      if (n == 0) { return 1; }
      else if (n > 0) {
	for (int i = 1; i<=n; i++) { cufac =  cufac * i; }
	return cufac;
      }
      else { return 0; }
    }
    
    
    __device__ double
    cuabs( double x )
    {
      if (x<0) { return -x; }
      else { return x; }
    }
    
    
    __device__ int
    cupmo( const int exponent )
    {
      if (exponent & 0x1) { return -1; }
      else { return +1; }
    }
    
    
    __device__ bool
    cuisOdd(const int val) 
    {
      return val & 0x1;
    }
    
    
    template<class T>
    __device__ void
    cuswap ( T a, T b )
    {
      T c = a;
      a=b;
      b=c;
    }
    
    
    __device__ int
    cuReflFactor(const int j,
		const int P,
		const int m,
		const int refl)
    {
      if (cuabs(P) != 1) {
	return 0;
      }
      if (m < 0) {
	return 0;
      }
      if (cuabs(refl) != 1) {
	return 0;
      }
      return refl * P * cupmo( (j - m) / 2 );
    }
    
    
    __device__ double
    cucgc(	 int j1,
	    int m1,
	    int j2,
	    int m2,
	    int J,
	    int M)
    {
      double clebschVal;
      if ( cuisOdd(j1-m1) or cuisOdd(j2-m2) or cuisOdd(J-M) ) {
	return 0;
      }
      if ( (cuabs(m1) > j1) or (cuabs(m2) > j2) or (cuabs(M) > J) ) {
	return 0;
      }
      if ( m1 + m2 != M ) {
	return 0;
      }
      if ((j1 < 0) or (j2 < 0) or (J < 0)) { 
	return 0;
      }
      if ((J < cuabs(j1 - j2)) or (J > j1 + j2)) { 
	return 0;
      }
      if (cuisOdd(j1 + j2 - J)) { 
	return 0;
      }
      else {
	int nu = 0;
	while ( ( (j1 - j2 - M) / 2 + nu < 0) or ( (j1-m1) / 2 + nu < 0 ) ) {
	  nu++;
	}
	double sum = 0;
	int d1, d2, n1;
	while ( ((d1 = (J - j1 + j2) / 2 - nu) >= 0) and ((d2 = (J + M) / 2 - nu) >= 0) and ((n1 = (j2 + J + m1) / 2 - nu) >= 0) ) {
	  const int d3 = (j1 - j2 - M) / 2 + nu;
	  const int n2 = (j1 - m1) / 2 + nu;
	  sum += ( (double)cupmo(nu + (j2 + m2) / 2) * cufac(n1) * cufac(n2) ) / ( cufac(nu) * cufac(d1) * cufac(d2) * cufac(d3));
	  nu++;
	}
	if ( sum == 0 ) {
	  return 0;
	}
	const double N1 = cufac((J  + j1 - j2) / 2);
	const double N2 = cufac((J  - j1 + j2) / 2);
	const double N3 = cufac((j1 + j2 - J ) / 2);
	const double N4 = cufac((J + M) / 2);
	const double N5 = cufac((J - M) / 2);
	    
	const double D0 = cufac((j1 + j2 + J) / 2 + 1);
	const double D1 = cufac((j1 - m1) / 2);
	const double D2 = cufac((j1 + m1) / 2);
	const double D3 = cufac((j2 - m2) / 2);
	const double D4 = cufac((j2 + m2) / 2);
    
	const double A  = ( (J + 1) * N1 * N2 * N3 * N4 * N5 ) / (D0 * D1 * D2 * D3 * D4);
	    
	clebschVal = ::sqrt(A) * sum;
    //     if (clebschVal > 1) { clebschVal = clebschVal / 2 ; }
      }
      return clebschVal;
    }
    
    
    __device__ double
    cudfunc( int j,
	    int m,
	    int n,
	    double theta)
    {
      if ( (j<0) or (cuabs(m) > j) or (cuabs(n) > j) ) {
	return 0;
      }
      if ( cuisOdd(j-m) or cuisOdd(j-n) ) {
	return 0;
      }
      
      if (j==0) {
	return 1;
      }
      //swap spinprojections for negative angle
      int _m = m;
      int _n = n;
      double thetaHalf = theta / 2;
      if (theta < 0) {
	thetaHalf = cuabs(thetaHalf);
	cuswap<int>( _m , _n );
      }
      
      const double cosThetaHalf = cos(thetaHalf);
      const double sinThetaHalf = sin(thetaHalf);
      
      double dFuncVal = 0;
      
      const int jpm       		= ( j + _m ) / 2;
      const int jpn       		= ( j + _n ) / 2;
      const int jmm       		= ( j - _m ) / 2;
      const int jmn       		= ( j - _n ) / 2;
      const double   kk        	= cufac(jpm) * cufac(jmm) * cufac(jpn) * cufac(jmn);
      const double   constTerm 	= cupmo(jpm) * ::sqrt(kk);
      
      double sumTerm = 0;
      const int mpn  = ( _m + _n ) / 2;
      const int kMin = max(0, mpn);
      const int kMax = min(jpm, jpn);
      for (int k = kMin; k <= kMax; ++k) {
	const int kmn1 = 2 * k - (_m + _n) / 2;
	const int jmnk = j + (_m + _n) / 2 - 2 * k;
	const int jmk  = (j + _m) / 2 - k;
	const int jnk  = (j + _n) / 2 - k;
	const int kmn2 = k - (_m + _n) / 2;
	const double factor = (  cufac(k) * cufac(jmk) * cufac(jnk) * cufac(kmn2) ) * cupmo(k);	/// why devide by? same effect as multiplying??
	sumTerm += ::pow(cosThetaHalf, kmn1) * ::pow(sinThetaHalf, jmnk) / factor;
      }
      dFuncVal = constTerm * sumTerm;
      return dFuncVal;
      
    }
    
    
    __device__ cuda::complex<double>
    cuDfunc(const int 			j,
	    const int 			m,
	    const int 			n,
	    const double		alpha,
	    const double		beta,
	    const double		gamma)
    {
      const double arg 			= - ( ((double)m / 2) * alpha + ((double)n / 2) * gamma ) ;
      const cuda::complex<double> argimag(0, arg);
      const cuda::complex<double> DFuncVal 	= cuda::exp(argimag) * cudfunc(j,m,n,beta);
      return DFuncVal;
    }
    
    
    __device__ cuda::complex<double>
    cuDfuncConj(const int 			j,
		const int 			m,
		const int 			n,
		const double			alpha,
		const double			beta,
		const double			gamma)
    {
      const cuda::complex<double> DFuncVal = cuda::conj(cuDfunc(j,m,n,alpha,beta,gamma));
      return DFuncVal;
    }
    
    
    __device__ cuda::complex<double>
    cuSpherHarm(const int 		l,
		const int 		m,
		const double 		theta,
		const double 		phi)
    {
      const cuda::complex<double> Yarg(0, ((double)m*phi) / 2 );
      const cuda::complex<double> YVal = ::sqrt((l+1)/(FOURPI)) * cuda::exp(Yarg) * cudfunc(l,m,0,theta);
      return YVal;
    }
    
    
    __device__ cuda::complex<double>
    cuDfuncRefl(const int 		j,
		const int 		m,
		const int 		n,
		const int 		P,
		const int 		refl,
		const double 		alpha,
		const double 		beta,
		const double	 	gamma)
    {
      cuda::complex<double> DFuncVal;
      const int reflFactor = cuReflFactor(j,P,m,refl);
      if (m == 0) {
	if (reflFactor == +1) {
	  return 0;
	}
	else {
	  DFuncVal = cuDfunc(j,0,n,alpha,beta,gamma);
	}
      }
      else {
	DFuncVal = ( cuDfunc(j,+m,n,alpha,beta,gamma) - (double)reflFactor * cuDfunc(j,-m,n,alpha,beta,gamma) ) / ::sqrt(2.)  ;
      }
      return DFuncVal;
    }
    
    
    __device__ cuda::complex<double>
    cuDfuncReflConj(const int j,
		    const int m,
		    const int n,
		    const int P,
		    const int refl,
		    const double alpha,
		    const double beta,
		    const double gamma)
    {
      const cuda::complex<double> DFuncVal = cuda::conj(cuDfuncRefl(j,m,n,P,refl,alpha,beta,gamma));
      return DFuncVal;
    }
    
    
    __device__ double
    GammaFunc( double m, double m0, double Gamma0, double q, double q0, int L) 
    {
	
      double BFValue 		= BFactor(L,q);
      double BF0	 	= BFactor(L, q0);
      
    double Gamma = ( Gamma0 * m0 * q  * BFValue ) / ( m * q0  * BF0 ) ; 
    return Gamma;
    }
    
    
    __device__ cuda::complex<double>
    bwig(double m, double m0, double Gamma0, double q, double q0, int L)
    {
      
	cuda::complex<double> imag(0,1);
	
      cuda::complex<double> BWVal;
      double GammaValue = GammaFunc(m,m0,Gamma0,q,q0,L);
      
      BWVal = ( Gamma0*m0) / ( (m0*m0 - m*m) - imag*m0*GammaValue) ;
      return BWVal;
    }
    
    
    __device__ cuda::complex<double>
    cuHelAmplitude( const double theta1,
		    const double phi1,
		    const double theta2,
		    const double phi2,
		    const double wX,
		    const double wf2,
		    const double qX,
		    const double qf2,
		    const int JX,
		    const int MX,
		    const int J1,
		    const int L1,
		    const int S1,
		    const int J2,
		    const int L2,
		    const int S2)
    {
 //     double wPIpm 	= 0.13957;
//      double qX 	= 0.50514288800993423; //cuBreakUpMomentum(wX,0.587483,wPIpm);
//      double qf2 	= 0.25846533761461093; //cuBreakUpMomentum(0.587483,wPIpm,wPIpm);
      cuda::complex<double> amplsum(0,0);
//      int lambda2 = 0;
      double Gamma0 = 0.1851;
      double q0 = 0.622239;
//      cuda::complex<double> bwigdat		(0.183925,0.007331);
      
	for(int lambda1 = -4; lambda1 <= +4; lambda1+=2) {
	  
	      amplsum +=	sqrt(L1+1.0) * cucgc(J1,lambda1,J2,0,S1,lambda1) * cucgc(L1,0,S1,lambda1,JX,lambda1) * cuDfuncReflConj(JX,MX,lambda1,-1,1,phi1,theta1,0) * sqrt(BFactor(L1,qX)) * 
				sqrt(L2+1.0) * cucgc(0,0,0,0,0,0) * cucgc(L2,0,S2,0,L2,0) * cuDfuncConj(4,lambda1,0,phi2,theta2,0) * sqrt(BFactor(L2,qf2)) * bwig(wf2,1.2754,Gamma0,qf2,q0,L2);

	}
     

      return amplsum;
      /*
      for(int lambda1 = -4; lambda1 <= 4; lambda1+=2) {
	amplsum += sqrt(L1+1.0) * cucgc(4,lambda1,0,0,4,lambda1)*cucgc(4,0,4,lambda1,4,lambda1)*cuDfuncReflConj(JX,MX,lambda1,-1,1,-0.564732,1.29395,0)*  sqrt(BFactor(L1,qX)) * 
		   sqrt(L2+1.0) * cucgc(0,0,0,0,0,0) * cucgc(4,0,0,0,4,0) * cuDfuncConj(4,lambda1,0,-1.6509,1.11757,0) * sqrt(BFactor(L2,qf2)) * bwig(0.58748315940941132,wf2,Gamma0,qf2,q0,L2);
			 
      }
      return   amplsum;
      */
    }
    
        struct vertexdata {
	    
	    int JX, Lambda, J1, L1, S1, J2, L2, S2;
	    double theta1, phi1, theta2, phi2, wX, wf2, massf2, massX, wpi, qX, qf2;
	    
	    };
  
      //Kernel for cuPWA using struct
    __global__ void
    GPUAmp2(vertexdata* data, cuda::complex<double>* result_GPUAmp, int N)
    {

      const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
//       if(idx<N) { result_GPUAmp[idx] = 0;}
     if(idx<N) { result_GPUAmp[idx] = cuHelAmplitude(data[idx].theta1,data[idx].phi1,data[idx].theta2,data[idx].phi2,data[idx].wX,data[idx].wf2,data[idx].qX,data[idx].qf2,data[idx].JX,data[idx].Lambda,data[idx].J1,data[idx].L1,data[idx].S1,data[idx].J2,data[idx].L2,data[idx].S2); }
	  
    }


    //Kernel for cuPWA
    __global__ void
    GPUAmp(double* theta1, double* phi1, double* theta2, double* phi2, double* wX, double* wf2, int* JX, int* Lambda, int* J1, int* L1, int* S1, int* J2, int* L2, int* S2, cuda::complex<double>* result_GPUAmp)
    {
      int N = 1;
      const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
       if(idx<N) { result_GPUAmp[idx] = cuHelAmplitude(theta1[idx],phi1[idx],theta2[idx],phi2[idx],wX[idx],wf2[idx],0.50514288800993423,0.25846533761461093,JX[idx],Lambda[idx],J1[idx],L1[idx],S1[idx],J2[idx],L2[idx],S2[idx]);  }
    }
  
  void DMAccess(vertexdata* data, int N)
  {
    printf("address0 before: %x\n",&data[0].JX);
    cutilSafeCall( cudaHostAlloc((void**)&data, N*sizeof(vertexdata), cudaHostAllocDefault) );
    printf("address0 after: %x\n",&data[0].JX);    
  }
    
  void GPUAmpcall2(vertexdata* data,int N, const vector<complex<double> >& cpudata)
  {
   
    printf("\nNow calculating on device...\n");
	  timespec bartime_beg, bartime_end; 
	  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &bartime_beg);
      cudaEvent_t start, stop;
      cutilSafeCall( cudaEventCreate(&start) );
      cutilSafeCall( cudaEventCreate(&stop)  );
	  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &bartime_end); 
	  printf("\n\n%.3f ms required for creating cuda-events (timespec)\n", ( ( (bartime_end.tv_sec - bartime_beg.tv_sec) ) + (bartime_end.tv_nsec - bartime_beg.tv_nsec) / 1e9 ) * 1e3 ); 
          
	  float timer[10];
      size_t sdata = sizeof(vertexdata);
      size_t sres  = sizeof(cuda::complex<double>);
    cuda::complex<double>* result_GPUAmp;
    cuda::complex<double>* dev_result_GPUAmp;
    vertexdata* dev_data; 
	      
	  cutilSafeCall( cudaEventRecord( start, 0) );
	  timespec footime_beg, footime_end; 
	  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &footime_beg); 
    cutilSafeCall( cudaHostAlloc((void**)&result_GPUAmp, N*sres, cudaHostAllocDefault) );
// //    cutilSafeCall( cudaHostAlloc((void**)&data, N*sdata, cudaHostAllocDefault) ); //doesn't work properly
	  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &footime_end); 
	  printf("\n\n%.3f ms required for allocating memory (timespec)\n", ( ( (footime_end.tv_sec - footime_beg.tv_sec) ) + (footime_end.tv_nsec - footime_beg.tv_nsec) / 1e9 ) * 1e3 ); 
          cutilSafeCall( cudaEventRecord( stop, 0) );
	  cutilSafeCall( cudaEventSynchronize( stop ) );
	  cutilSafeCall( cudaEventElapsedTime( &timer[0], start, stop));
	  printf("%.3f ms required for allocating memory\n",timer[0]);
     
    cutilSafeCall( cudaMalloc((void**)&dev_result_GPUAmp, N*sres) );
    cutilSafeCall( cudaMalloc((void**)&dev_data, N*sdata) );
	  
	  cutilSafeCall( cudaEventRecord( start, 0) );
    cutilSafeCall( cudaMemcpy(dev_data, data, N*sdata, cudaMemcpyHostToDevice) );
          cutilSafeCall( cudaEventRecord( stop, 0) );
	  cutilSafeCall( cudaEventSynchronize( stop ) );
	  cutilSafeCall( cudaEventElapsedTime( &timer[1], start, stop));
	  printf("%.3f ms required for uploading\n",timer[1]);
	  printf("%f GB/s bandwidth for upload\n",(N*(sdata/(timer[1]/1e3))/(1024*1024*1024)));
	  
    int threadsPerBlock = 512;
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;    
	  cutilSafeCall( cudaEventRecord( start, 0) );
    GPUAmp2<<<blocksPerGrid, threadsPerBlock>>>(dev_data,dev_result_GPUAmp,N);
	  cutilSafeCall( cudaThreadSynchronize() );
	  cutilSafeCall( cudaEventRecord( stop, 0) );
	  cutilSafeCall( cudaEventSynchronize( stop ) );
	  cutilSafeCall( cudaEventElapsedTime( &timer[2], start, stop));
	  printf("%.3f ms required for calculating on device\n",timer[2]);
	  printf("%f GB/s Bandwidth\n", (N*(sdata+sres)/(timer[2]/1e3))/(1024*1024*1024) );
    
          cutilSafeCall( cudaEventRecord( start, 0) );
    cutilSafeCall( cudaMemcpy(result_GPUAmp, dev_result_GPUAmp, N*sres, cudaMemcpyDeviceToHost) );
          cutilSafeCall( cudaEventRecord( stop, 0) );
	  cutilSafeCall( cudaEventSynchronize( stop ) );
	  cutilSafeCall( cudaEventElapsedTime( &timer[3], start, stop));
	  printf("%.3f ms required for downloading\n",timer[3]);
	  printf("%f GB/s bandwidth for download\n",(N*(sres/(timer[3]/1e3))/(1024*1024*1024)));
	  printf("%.3f ms required for GPU call without malloc stuff, etc.\n",timer[1]+timer[2]+timer[3]);

	  
  if(0) {
    for(int e=0;e<N;e++) {
    printf("\nGPUAmp[%d]: %2.20f + i*(%2.20f)",e, result_GPUAmp[e].real(), result_GPUAmp[e].imag() );
    printf("\nCPUAmp[%d]: %2.20f + i*(%2.20f)",e, cpudata[e].real(), cpudata[e].imag() );
    }
  }

  if (1) {
    timespec roottime_beg, roottime_end; 
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &roottime_beg); 
    cutilSafeCall( cudaEventRecord( start, 0) );
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
    cutilSafeCall( cudaEventRecord( stop, 0) );
    cutilSafeCall( cudaEventSynchronize( stop ) );
    cutilSafeCall( cudaEventElapsedTime( &timer[4], start, stop));
    printf("%.3f ms required for root stuff\n",timer[4]);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &roottime_end); 
    printf("\n\n%.3f ms required for root stuff (time spec)\n", ( ( (roottime_end.tv_sec - roottime_beg.tv_sec) ) + (roottime_end.tv_nsec - roottime_beg.tv_nsec) / 1e9 ) * 1e3 ); 
 
  }
	
    cutilSafeCall( cudaFreeHost(result_GPUAmp) );
//    cutilSafeCall( cudaFreeHost(data) );
    cudaFree(dev_result_GPUAmp);
    cudaFree(dev_data);
   
    cutilSafeCall( cudaEventDestroy(start) );
    cutilSafeCall( cudaEventDestroy(stop) );

  }




















  void GPUAmpcall(double theta1, double phi1, double theta2, double phi2, double wX, double wf2, int JX, int Lambda, int J1, int L1, int S1, int J2, int L2, int S2)
  {
     printf("\ntheta1 = %f", theta1);
     printf("\nphi1 = %f", phi1);
     printf("\ntheta2 = %f", theta2);
     printf("\nphi2 = %f", phi2);
     printf("\nwX = %f", wX);
     printf("\nJX = %i\n", JX);
     printf("\nMX = Lambda = %i\n", Lambda);
     printf("\nJ1 = %i\n", J1);
     printf("\nL1 = %i\n", L1);
     printf("\nS1 = %i\n", S1);
     printf("\nJ2 = %i\n", J2);
     printf("\nL2 = %i\n", L2);
     printf("\nS2 = %i\n", S2);
		
//    printf("B = %f", B);
  
    
    cuda::complex<double>* result_GPUAmp;
    cuda::complex<double>* dev_result_GPUAmp;
    double* dev_theta1;
    double* dev_theta2;
    double* dev_phi1;
    double* dev_phi2;
    double* dev_wX;
    double* dev_wf2;
    int* dev_JX;
    int* dev_Lambda;
    int* dev_J1;
    int* dev_L1;
    int* dev_S1;
    int* dev_J2;
    int* dev_L2;
    int* dev_S2;
    

    
  
    result_GPUAmp = (cuda::complex<double>*)malloc(sizeof(cuda::complex<double>));
    if (result_GPUAmp == 0) free(result_GPUAmp);
  
    cutilSafeCall( cudaMalloc((void**)&dev_result_GPUAmp, sizeof(cuda::complex<double>)) );
    
    cutilSafeCall( cudaMalloc((void**)&dev_theta1, sizeof(double)) );
    cutilSafeCall( cudaMalloc((void**)&dev_theta2, sizeof(double)) );
    cutilSafeCall( cudaMalloc((void**)&dev_phi1, sizeof(double)) );
    cutilSafeCall( cudaMalloc((void**)&dev_phi2, sizeof(double)) );
    cutilSafeCall( cudaMalloc((void**)&dev_wX, sizeof(double)) );
    cutilSafeCall( cudaMalloc((void**)&dev_wf2, sizeof(double)) );
    cutilSafeCall( cudaMalloc((void**)&dev_JX, sizeof(int)) );
    cutilSafeCall( cudaMalloc((void**)&dev_Lambda, sizeof(int)) );
    cutilSafeCall( cudaMalloc((void**)&dev_J1, sizeof(int)) );
    cutilSafeCall( cudaMalloc((void**)&dev_L1, sizeof(int)) );
    cutilSafeCall( cudaMalloc((void**)&dev_S1, sizeof(int)) );
    cutilSafeCall( cudaMalloc((void**)&dev_J2, sizeof(int)) );
    cutilSafeCall( cudaMalloc((void**)&dev_L2, sizeof(int)) );
    cutilSafeCall( cudaMalloc((void**)&dev_S2, sizeof(int)) );
    
    cutilSafeCall( cudaMemcpy(dev_theta1, &theta1, sizeof(double), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dev_theta2, &theta2, sizeof(double), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dev_phi1,	  &phi1,   sizeof(double), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dev_phi2,	  &phi2,   sizeof(double), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dev_wX,	  &wX,     sizeof(double), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dev_wf2,	  &wf2,    sizeof(double), cudaMemcpyHostToDevice) );
    
    cutilSafeCall( cudaMemcpy(dev_Lambda, &Lambda, sizeof(int), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dev_JX, &JX, sizeof(int), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dev_J1, &J1, sizeof(int), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dev_L1, &L1, sizeof(int), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dev_S1, &S1, sizeof(int), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dev_J2, &J2, sizeof(int), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dev_L2, &L2, sizeof(int), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dev_S2, &S2, sizeof(int), cudaMemcpyHostToDevice) );
    
    
    
    int threadsPerBlock = 512;
    int blocksPerGrid = (1 + threadsPerBlock - 1) / threadsPerBlock;
    
    GPUAmp<<<blocksPerGrid, threadsPerBlock>>>(dev_theta1,dev_phi1,dev_theta2,dev_phi2,dev_wX,dev_wf2,dev_JX,dev_Lambda,dev_J1,dev_L1,dev_S1,dev_J2,dev_L2,dev_S2,dev_result_GPUAmp);
    
    cutilSafeCall( cudaMemcpy(result_GPUAmp, dev_result_GPUAmp, sizeof(cuda::complex<double>), cudaMemcpyDeviceToHost) );
    
    
    printf("\nresult: %2.20f + i*(%2.20f)\n\n", result_GPUAmp[0].real(), result_GPUAmp[0].imag() );
  
    free(result_GPUAmp);
    cudaFree(dev_result_GPUAmp);
    cudaFree(dev_theta1);
    cudaFree(dev_theta2);
    cudaFree(dev_phi1);
    cudaFree(dev_phi2);
    cudaFree(dev_wX);
    cudaFree(dev_wf2);
    cudaFree(dev_JX);
    cudaFree(dev_Lambda);
    cudaFree(dev_J1);
    cudaFree(dev_L1);
    cudaFree(dev_S1);
    cudaFree(dev_J2);
    cudaFree(dev_L2);
    cudaFree(dev_S2);

  //nvidia-smi -a -l

  
  }

} //namespace cupwa