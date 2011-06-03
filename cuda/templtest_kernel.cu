#include <stdio.h>
#include "complex.cuh"
using namespace rpwa;



const long unsigned int N = 500000;


// Variables... for the host
const double m0 = 0.77549; 			// Gev/c^2, rho(770)
const double q0 = 0.361755;  			// breakup momentum --> expected from kinematics, here pi+pi-
const double Gamma0 = 0.1462; 			// GeV, rho(770)
int L = 10;					// angular momentum
const double Pr = 0.1973; 			// Gev/c, hbarc
cuda::complex<double> imag(0.0,1.0);		// imaginary unit


///Variables for the device functions are inside of them, such as m0, q0, Gamma0, L

/// ///////////////////////////////////// ///
/// templated pow function (complex test) ///
/// ///////////////////////////////////// ///
template<class T>
__global__ void 
testPow( T* A, T* B, cuda::complex<T>* C) 
{
  cuda::complex<T> imag(0.0,1.0);
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx<N) C[idx] = pow(A[idx],B[idx])*imag + 1.0 ;
}

/// ////////////////////////////////////////////// ///
/// templated pow function (interger exponentiate) ///
/// ////////////////////////////////////////////// ///
template<class T>
__device__ T 
powint( T* z, int n) 
{
  return pow( *z , n) ;
}





/// ////////////////////////////// ///
/// Blatt-Weisskopf barrier factor ///
/// ////////////////////////////// ///
template<class T>
__device__ T
BFactor( int L, T* q )
{
  T z = (*q * *q)/(Pr * Pr);			// substitution
  T ret;					// return value
  int m = L/2;
  /// actually the return values are squarerooted, but in the GammaFunction they are squared again, so...
  /// this implies a speedup-factor of 1.74 for every calculation
  switch (m) {
  case 0:
    ret =  (T)1.0;
    break;
  case 1:
    ret =  ((T)2.0 * z)/(z + 1) ;
    break;
  case 2:
    ret =  ((T)13.0 * z * z)/( (T)9.0 + z*(z+(T)3.0) ) ;								//
    break;
  case 3:
    ret =  ((T)277.0 * z*z*z)/((T)225.0 + z*((T)45.0 + z*((T)6.0 + z))) ;						//
    break;
  case 4:
    ret =  ((T)12746.0 * z*z*z*z)/((T)11025.0 + z*((T)1575.0 + z*((T)135.0 + z*((T)10.0 + z))) ) ; 			//
    break;
  case 5:
    ret =  z*z*z*z*z/((T)893025.0 + z*((T)99225.0 + z*((T)6300.0 + z*((T)315.0 + z*((T)15.0 + z))))) ; 			//  
    break;
  default:
    ret =  1.0;
    break;
  }
  return ret;
}

/// /////////////////////////////////////////////////////////////////////////////////////////// ///
/// templated breit-wigner function (complex) with GammaFunc and Blatt-Weisskopf barrier factor ///
/// /////////////////////////////////////////////////////////////////////////////////////////// ///
template<class T>
__global__ void
imagtest(T* m, T* q, cuda::complex<T>* BW )
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
    T m0 = 0.77549;							/// insert expected peak for mass
    T Gamma0 = 0.1462;							/// self explaining
    cuda::complex<T> imag(0,1);

  T GammaValue = GammaFunc(&m[idx],&q[idx]);
   
  if (idx<N) BW[idx] = ( Gamma0*m0) / ( (m0*m0 - m[idx]*m[idx]) - imag*m0*GammaValue) ; /// will test to calculate the brackets before adding an complex number ... no notable speedup
}

/// //////// ///
/// Gamma(m) ///
/// //////// ///
template<class T>
__device__ T
GammaFunc( T* m, T* q ) 
{
    int L = 10;								/// insert angular momentum
    T q0 = 0.361755;							/// insert expected breakup momentum by calculating kinematics
    
  T BFValue = BFactor(L,q);
  T BF0 = BFactor(L, &q0);
  
//  T Gamma = Gamma0 * (m0 / *m) * (*q / q0) * (BFValue / BF0) * (BFValue / BF0) ;
 T Gamma = ( Gamma0 * m0 * *q  * BFValue ) / ( *m * q0  * BF0 ) ; 	/// nearly 2 times (1.88) faster than the formula above, means first calculating the nominator and denominator and then dividing
 return Gamma;
}










/// //////////////////////////// ///
/// breit-wigner function (real) ///
/// //////////////////////////// ///
/// ...well just for testing.... ///
template<class T>
__global__ void
finBW( T* m, T* q, cuda::complex<T>* BW)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
    cuda::complex<T> imag(0,1);
    T GammaValue = GammaFunc(m,q);
  if (idx<N) BW[idx] =  ( m0 * Gamma0) / pow(( pow((m0*m0 - m[idx]*m[idx]),2) + m0*m0*GammaValue*GammaValue),0.5)  ;
}

