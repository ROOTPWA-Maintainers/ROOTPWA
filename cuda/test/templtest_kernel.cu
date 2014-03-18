#include <stdio.h>
#include "complex.cuh"
using namespace rpwa;

#define FOURPI 4*3.1415926

//const long unsigned int N = 25000000;


// Variables... for the host
const double m0 = 0.77549; 			// Gev/c^2, rho(770)
const double q0 = 0.361755;  			// breakup momentum --> expected from kinematics, here pi+pi-
const double Gamma0 = 0.1462; 			// GeV, rho(770)
int L = 0;					// angular momentum
const double Pr = 0.1973; 			// Gev/c, hbarc
cuda::complex<double> imag(0.0,1.0);		// imaginary unit


///Variables for the device functions are inside of them, such as m0, q0, Gamma0, L




/// ///////////////////////////////////// ///
/// templated pow function (complex test) ///
/// ///////////////////////////////////// ///
template<class T>
__global__ void 
testPow( T* A, T* B, T* C, int N) 
{
//  cuda::complex<T> imag(0.0,1.0);
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx<N) C[idx] = pow(A[idx],B[idx]);
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

/// //////////////// ///
/// Breakup Momentum ///
/// //////////////// ///
__device__ double
cuBreakUpMomentum(const double M,   // mass of mother particle
                  const double m1,  // mass of daughter particle 1
                  const double m2)  // mass of daughter particle 2
{
  if (M < m1 + m2)
    return 0;
  return sqrt((M - m1 - m2) * (M + m1 + m2) * (M - m1 + m2) * (M + m1 - m2)) / (2 * M);
}


/// ////////////////////////////// ///
/// Blatt-Weisskopf barrier factor ///
/// ////////////////////////////// ///
// note that for speed up purposes the barrier factor is defined without the squareroot here!
template<class T>
__device__ T
BFactor( int L, T q )
{
  T z = (q * q)/(Pr * Pr);			// substitution
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
    ret =  ((T)998881.0 * z*z*z*z*z)/((T)893025.0 + z*((T)99225.0 + z*((T)6300.0 + z*((T)315.0 + z*((T)15.0 + z))))) ; 			//  
    break;
  case 6:
    ret =  ((T)118394977.0 * z*z*z*z*z*z)/((T)108056025.0 + z*((T)9823275.0 + z*((T)496125.0 + z*((T)18900.0 + z*((T)630.0 + z*((T)21.0 + z)))))) ;
    break;
  case 7:
    ret =  ((T)19727003738.0 * z*z*z*z*z*z*z)/((T)18261468225.0 + z*((T)1404728325.0 + z*((T)58939650.0 + z*((T)1819125.0 + z*((T)47250.0 + z*((T)1134.0 + z*((T)28.0 + z))))))) ;
    break;
  default:
    ret =  1.0;
    break;
  }
  return ret;
}

/// ////////////////////////////// ///
/// Blatt-Weisskopf barrier factor ///
/// ////////////////////////////// ///
// note that for speed up purposes the barrier factor is defined without the squareroot here!
template<class T>
__device__ T
BFactorTest( int L, T* q )
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
    ret =  ((T)998881.0 * z*z*z*z*z)/((T)893025.0 + z*((T)99225.0 + z*((T)6300.0 + z*((T)315.0 + z*((T)15.0 + z))))) ; 			//  
    break;
  case 6:
    ret =  ((T)118394977.0 * z*z*z*z*z*z)/((T)108056025.0 + z*((T)9823275.0 + z*((T)496125.0 + z*((T)18900.0 + z*((T)630.0 + z*((T)21.0 + z)))))) ;
    break;
  case 7:
    ret =  ((T)19727003738.0 * z*z*z*z*z*z*z)/((T)18261468225.0 + z*((T)1404728325.0 + z*((T)58939650.0 + z*((T)1819125.0 + z*((T)47250.0 + z*((T)1134.0 + z*((T)28.0 + z))))))) ;
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

  T GammaValue = GammaFuncTest(&m[idx],&q[idx]);
  int N = 1;
  if (idx<N) BW[idx] = ( Gamma0*m0) / ( (m0*m0 - m[idx]*m[idx]) - imag*m0*GammaValue) ; /// will test to calculate the brackets before adding an complex number ... no notable speedup
}


/// //////// ///
/// Gamma(m) ///
/// //////// ///
template<class T>
__device__ T
GammaFuncTest( T* m, T* q) 
{
  int L = 0;								/// insert angular momentum
  T q0 = 0.361755;
    
  T BFValue 	= BFactorTest(L,q);
  T BF0 	= BFactorTest(L,&q0);
  
//  T Gamma = Gamma0 * (m0 / *m) * (*q / q0) * (BFValue / BF0) * (BFValue / BF0) ;
 T Gamma = ( Gamma0 * m0 * *q  * BFValue ) / ( *m * q0  * BF0 ) ; 	/// nearly 2 times (1.88) faster than the formula above, means first calculating the nominator and denominator and then dividing
 return Gamma;
}

/// ///////// ///
/// factorial ///
/// ///////// ///
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

/// ///////////////// ///
/// absolute function ///
/// ///////////////// ///
__device__ double
cuabs( double x )
{
  if (x<0) { return -x; }
  else { return x; }
}

/// ////// ///
/// (-1)^n ///
/// ////// ///
__device__ int
cupmo( const int exponent )
{
  if (exponent & 0x1) { return -1; }
  else { return +1; }
}

/// ///// ///
/// isOdd ///
/// ///// ///
__device__ bool
cuisOdd(const int val)  ///< returns whether val is an odd number 
{
  return val & 0x1;
}

/// //// ///
/// swap ///
/// //// ///
template<class T>
__device__ void
cuswap ( T a, T b )
{
  T c = a;
  a=b;
  b=c;
}

/// /////////////////// ///
/// Reflectivity Factor ///
/// /////////////////// ///
__device__ int
cuReflFactor(const int j,
             const int P,
             const int m,
             const int refl)  ///< calculates prefactor for reflectibity symmetrization
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




///      !NOTE! spins and projection quantum numbers are in units of hbar/2

/// /////////////////////////////////////////// ///
/// Clebsch-Gordan-Coefficients (wiki formulas) ///
/// /////////////////////////////////////////// ///
// doens't work properly till now
__device__ double
cucgcwiki(	int j1,
		int m1,
		int j2,
		int m2,
		int j,
		int m)
{
  double clebschVal;
  if ( cuisOdd(j1-m1) or cuisOdd(j2-m2) or cuisOdd(j-m) ) {
    return 0;
  }
  if ( (cuabs(m1) > j1) or (cuabs(m2) > j2) or (cuabs(m) > j) ) {
    return 0;
  }
  if ( m1 + m2 != m ) {
    return 0;
  }
  else {
    int k = 0;
    while ( ((j - j2 + m1) / 2 + k < 0) or ((j - j1 - m2) / 2 + k < 0) ) {
      k++;
    }
    double sum = 0;
    int d1, d2, d3;
    while ( ((d1 = (j1 + j2 - j) / 2 - k) >= 0) and ((d2 = (j1 - m1) / 2 - k) >= 0) and ((d3 = (j2 + m2) / 2 - k) >= 0) ) {
      const int d4 = (j - j2 + m1) / 2 + k;
      const int d5 = (j - j1 - m2) / 2 + k;
      sum += ( cupmo(k) / ( cufac(k) * cufac(d1) * cufac(d2) * cufac(d3) * cufac(d4) * cufac(d5) ) );
      k++;
    }
    if ( sum == 0 ) {
      return 0;
    }
    const double N1 = cufac((j + j1 - j2) / 2);
    const double N2 = cufac((j - j1 + j2) / 2);
    const double N3 = cufac((j1 + j2 - j) / 2);
    const double N4 = cufac((j + m) / 2);
    const double N5 = cufac((j - m) / 2);
    const double N6 = cufac((j1 - m1) / 2);
    const double N7 = cufac((j1 + m1) / 2);
    const double N8 = cufac((j2 - m2) / 2);
    const double N9 = cufac((j2 + m2) / 2);
    
    const double D0 = cufac((j1 + j2 + j) / 2 + 1);
    
    double A = ( (j + 1) * N1 * N2 * N3 * N4 * N5 * N6 * N7 * N8 * N9 ) / D0 ;
    clebschVal = sqrtf(A) * sum;
  }
  return clebschVal;
}


/// /////////////////////////// ///
/// Clebsch-Gordan-Coefficients ///
/// /////////////////////////// ///
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
  if ((j1 < 0) or (j2 < 0) or (J < 0)) { // negative spins are not allowed
    return 0;
  }
  if ((J < cuabs(j1 - j2)) or (J > j1 + j2)) { // make sure J is in physical allowed range
    return 0;
  }
  if (cuisOdd(j1 + j2 - J)) { // check that J is in the half-integer or integer series, respectively
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


// Description:
//      optimized Wigner d-function d^j_{m n}(theta) with caching
//      used as basis for optimized spherical harmonics Y_l^m(theta, phi)
//      as well as for optimized D-function D^j_{m n}(alpha, beta,
//      gamma) and D-function in reflectivity basis
//
//     	based on PWA2000 function d_jmn_b() in pputil.cc
//
//      !NOTE! spin j and projections m and n are in units of hbar/2
//



/// ///////////////// ///
/// Wigner d-Function ///
/// ///////////////// ///
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


/// ////////// ///
/// D-Function ///
/// ////////// ///
__device__ cuda::complex<double>
cuDfunc(const int 			j,
	const int 			m,
	const int 			n,
	const double			alpha,
	const double			beta,
	const double			gamma)
{
  const double arg 			= - ( ((double)m / 2) * alpha + ((double)n / 2) * gamma ) ;
  const cuda::complex<double> argimag(0, arg);
  const cuda::complex<double> DFuncVal 	= cuda::exp(argimag) * cudfunc(j,m,n,beta);
  return DFuncVal;
}

/// ///////////////////// ///
/// D-Function conjugated ///
/// ///////////////////// ///
__device__ cuda::complex<double>
cuDfuncConj(const int 			j,
	    const int 			m,
	    const int 			n,
	    const double		alpha,
	    const double		beta,
	    const double		gamma)
{
  const cuda::complex<double> DFuncVal = cuda::conj(cuDfunc(j,m,n,alpha,beta,gamma));
  return DFuncVal;
}

/// /////////////////// ///
/// Spherical Harmonics ///
/// /////////////////// ///
__device__ cuda::complex<double>
cuSpherHarm(const int 		l,
	    const int 		m,
	    const double 	theta,
	    const double 	phi)
{
  const cuda::complex<double> Yarg(0, ((double)m*phi) / 2 );
  const cuda::complex<double> YVal = ::sqrt((l+1)/(FOURPI)) * cuda::exp(Yarg) * cudfunc(l,m,0,theta);
  return YVal;
}

/// //////////////// ///
/// D-Function Refl. ///
/// //////////////// ///
__device__ cuda::complex<double>
cuDfuncRefl(const int 		j,
	    const int 		m,
	    const int 		n,
	    const int 		P,
	    const int 		refl,
	    const double 	alpha,
	    const double 	beta,
	    const double 	gamma)
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

/// ////////////////////// ///
/// D-Function Refl. Conj. ///
/// ////////////////////// ///
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
  int N = 1;
  if (idx<N) BW[idx] =  ( m0 * Gamma0) / pow(( pow((m0*m0 - m[idx]*m[idx]),2) + m0*m0*GammaValue*GammaValue),0.5)  ;
}

/// //////// ///
/// Gamma(m) ///
/// //////// ///
template<class T>
__device__ T
GammaFunc( T m, T m0, T Gamma0, T q, T q0, int L) 
{
//    int L = 8;								/// insert angular momentum

    
  T BFValue 	= BFactor(L,q);
  T BF0 	= BFactor(L, q0);
  
//  T Gamma = Gamma0 * (m0 / *m) * (*q / q0) * (BFValue / BF0) * (BFValue / BF0) ;
 T Gamma = ( Gamma0 * m0 * q  * BFValue ) / ( m * q0  * BF0 ) ; 	/// nearly 2 times (1.88) faster than the formula above, means first calculating the nominator and denominator and then dividing
 return Gamma;
}

/// //////////// ///
/// Breit-Wigner ///
/// //////////// ///
__device__ cuda::complex<double>
bwig(double m, double m0, double Gamma0, double q, double q0, int L)
{
//    double m0 			= 0.77549;							/// insert expected peak for mass
//    double Gamma0 		= 0.1462;							/// self explaining
    cuda::complex<double> imag(0,1);
    
  cuda::complex<double> BWVal;
  double GammaValue = GammaFunc<double>(m,m0,Gamma0,q,q0,L);
   
  BWVal = ( Gamma0*m0) / ( (m0*m0 - m*m) - imag*m0*GammaValue) ;
  return BWVal;
}

/*
/// ///////////////// ///
/// Amplitude test... ///
/// ///////////////// ///
__device__ cuda::complex<double>
cuHelAmplitude(	const double theta1,
		const double phi1,
		const double theta2,
		const double phi2,
		const double wX,
		const double wf2,
		const int JX,
		const int MX,
		const int J1,
		const int L1,
		const int S1,
		const int J2,
		const int L2,
		const int S2)
{
  double wPIpm 	= 0.13957;
  double qX 	= cuBreakUpMomentum(wX,wf2,wPIpm);
  double qf2 	= cuBreakUpMomentum(wf2,wPIpm,wPIpm);
  cuda::complex<double> amplsum(0,0);
  for(int lambda1 = -J1; lambda1 <= +J1; lambda1+=2) {
    for(int lambda2 = -J2; lambda2 <= +J2; lambda2+=2) {
      amplsum += 	::sqrt(2.0*L1+1.0) * cucgc(J1,lambda1,J2,lambda2,J1,lambda1-lambda2) * cucgc(L1,0,S1,lambda1-lambda2,JX,lambda1-lambda2) * cuDfuncConj(JX,MX,lambda1,theta1,phi1,0) * ::sqrt((double)BFactor(L1,qX)) * 
			::sqrt(2.0*L2+1.0) * cucgc(0,0,0,0,0,0) * cucgc(L2,0,0,0,L2,0) * cuDfuncConj(4,lambda1,0,theta2,phi2,0) * ::sqrt((double)BFactor(L2,qf2)) * bwig(wf2,qf2,L2);
    }
  }
  return amplsum;
}
*/


