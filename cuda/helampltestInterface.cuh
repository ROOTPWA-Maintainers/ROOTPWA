#include <iostream>
#include <complex>
#include "complex.cuh"

using namespace rpwa;




namespace cupwa {
  
//  void GPUAmp(double*,double*,double*);

    struct vertexdata {
	    
	    int JX, Lambda, J1, L1, S1, J2, L2, S2;
	    double theta1, phi1, theta2, phi2, wX, wf2, massf2, massX, wpi, qX, qf2;
	    
	    };

  void GPUAmpcall(double,double,double,double,double,double,double,double,int,int,int,int,int,int,int,int);
  void GPUAmpcall2(vertexdata*, int, const std::vector<std::complex<double> >& );
  void DMAccess(vertexdata*, int);



}