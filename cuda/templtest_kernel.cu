#include <stdio.h>
#include "complex.cuh"
using namespace rpwa;



const long unsigned int N = 500000;


// Variables... for the host
const double m0 = 1000; 			// Mev/c^2
const double q0 = 0.13975;  			// breakup momentum --> expected from kinematics
const double Gamma0 = 100; 			// Mev
int L = 8;					// angular momentum
const double Pr = 0.1973; 			// Gev/c
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
    ret = 1.0;
    break;
  case 1:
    ret =  (2.0 * z)/(z + 1) ;
    break;
  case 2:
//     ret =  (13.0 * z * z)/( (z - 3.0)*(z - 3.0) + 9.0 * z ) ;
    ret =  (13.0 * z * z)/( 9.0 + z*(z+3.0) ) ;								// small difference but faster, again...
//     ret =  (13.0 * z * z)/( z*z + 3.0*z + 9.0); 							// nope
    break;
  case 3:
//     ret =  (277.0 * z*z*z)/(z * (z - 15.0)*(z - 15.0) + 9.0 * (2.0 * z - 5.0)*(2.0 * z - 5.0) ) ;
    ret =  (277.0 * z*z*z)/(225.0 + z*(45.0 + z*(6.0 + z))) ;						// fastest again
//     ret =  (277 * z*z*z)/(225 + z*(45 + z*(6 + z))) ;						// fastest again but iwth integers, only influece to floats...
//     ret =  (277.0 * z*z*z)/(z*z*z+6.0*z*z+48.0*z+225) ;						// nope
    break;
  case 4:
//     ret =  (12746.0 * z*z*z*z)/( (z * z - 45.0 * z + 105.0)*(z * z - 45.0 * z + 105.0) + 25.0 * z * (2.0 * z - 21.0)*(2.0 * z - 21.0) );
     ret =  ((T)12746.0 * z*z*z*z)/((T)11025.0 + z*((T)1575.0 + z*((T)135.0 + z*((T)10.0 + z))) ) ; 			// fastest!!
//    ret =  (12746 * z*z*z*z)/( 11025 + z*(1575 + z*(135 + z*(10 + z))) ) ; 				// fastest but with integers
//     ret =  (12746.0 * powint(&z , 4))/( 11025.0 + z*(1575.0 + z*(135.0 + z*(10.0 + z))) ) ;		// not as fast as i would expect, okay gpu is not fully available... 
//     ret =  (12746.0 * z*z*z*z)/( z*z*z*z + 10.0*z*z*z + 135.0*z*z + 1575.0*z + 11025.0 ) ; 		// seems to be, that this is not faster?!?!?!? less computations but higher numbers... hmmm....
    break;
  case 5:
//     ret =   z*z*z*z*z/(893025.0 +99225.0*z +6300.0*z*z +315.0*z*z*z +15.0*z*z*z*z +z*z*z*z*z) ;
    ret =   z*z*z*z*z/(893025.0 + z*(99225.0 + z*(6300.0 + z*(315.0 + z*(15.0 + z))))) ; 		// again faster, maybe fastest? true, true
//     ret =   z*z*z*z*z/(893025 + z*(99225 + z*(6300 + z*(315 + z*(15 + z))))) ; 			// whats the point if i do not use floats or doubles? speed up...    
//     ret =   powint(&z,5)/(893025.0 + z*(99225.0 + z*(6300.0 + z*(315.0 + z*(15.0 + z))))) ; 		// not as fast as expected
//     ret =   z*z*z*z*z/( (z+13.2987) * (z*z-14.6962*z+329.652) * (z*z+16.3975*z+203.704) ) ; 		// wow i think that makes a difference, not a huge, but a difference. actually it is wrong, not exact enough...
    break;
  default:
    ret = 1.0;
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
    T m0 = 1000;							/// insert expected peak for mass
    T Gamma0 = 100;							/// self explaining
    cuda::complex<T> imag(0,1);

  T GammaValue = GammaFunc(&m[idx],&q[idx]);
   
  if (idx<N) BW[idx] = ( Gamma0*m0) / ( (m0*m0 - m[idx]*m[idx]) - imag*m0*GammaValue) ; /// will test to calculate the brackets before adding an complex number
//   if (idx<N) BW[idx] = ( Gamma0*m0*m0*m0 - Gamma0*m0*m[idx]*m[idx] ) / ( m0*m0*m0*m0 - 2*m0*m0*m[idx]*m[idx] - m[idx]*m[idx]*m[idx]*m[idx] + m0*m0*GammaValue*GammaValue ) + imag*((m0*m0*Gamma0*GammaValue)/( m0*m0*m0*m0 - 2*m0*m0*m[idx]*m[idx] - m[idx]*m[idx]*m[idx]*m[idx] + m0*m0*GammaValue*GammaValue ))  ;  
/// well somehow this doesn't work, but even if it works, there would be bandwidth loss, too much computation
}

/// //////// ///
/// Gamma(m) ///
/// //////// ///
template<class T>
__device__ T
GammaFunc( T* m, T* q ) 
{
    int L = 8;								/// insert angular momentum
    T q0 = 0.13975;							/// insert expected breakup momentum by calculating kinematics
    
  T BFValue = BFactor(L,q);
  T BF0 = BFactor(L, &q0);
  
//  T Gamma = Gamma0 * (m0 / *m) * (*q / q0) * (BFValue / BF0) * (BFValue / BF0) ;
 T Gamma = ( Gamma0 * m0 * *q  * BFValue ) / ( *m * q0  * BF0 ) ; /// nearly 2 times (1.88) faster than the formula above, means first calculating the nominator and denominator and the dividing
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





// (Gamma0*Gamma0 * (m0*m0 / (m[idx]*m[idx])) * (q[idx]*q[idx] / (q0*q0)) * bf) /// === GammaValue
