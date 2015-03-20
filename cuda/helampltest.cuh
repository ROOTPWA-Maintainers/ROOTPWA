///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      general isobar decay amplitude in helicity formalism
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef CUISOBARHELICITYAMPLITUDE_H
#define CUISOBARHELICITYAMPLITUDE_H

#include "complex.cuh"

using namespace rpwa;

#define FOURPI (4 * 3.1415926)

namespace cupwa { 

/// //////////////// ///
/// Breakup Momentum ///
/// //////////////// ///
template<class T>
__device__ T
cuBreakUpMomentum(const T M,   // mass of mother particle
                  const T m1,  // mass of daughter particle 1
                  const T m2)  // mass of daughter particle 2
{
	if (M < m1 + m2) return 0;
	return sqrt((M - m1 - m2) * (M + m1 + m2) * (M - m1 + m2) * (M + m1 - m2)) / (2 * M);
}


/// ////////////////////////////// ///
/// Blatt-Weisskopf barrier factor ///
/// ////////////////////////////// ///
// note that for speed up purposes the barrier factor is defined without the squareroot here!
template<class T>
__device__ T
BFactor(int L, T q) 
{
	
	T Pr = 0.1973;
	T z = (q * q) / (Pr * Pr); // substitution
	int m = L / 2;
	
	/// actually the return values are squarerooted, but in the GammaFunction they are squared again, so...
	/// this implies a speedup-factor of 1.74 for every calculation
	
	switch (m) {
		case 0:
			return (T)1.0;
		case 1:
			return ((T)2.0 * z)/(z + 1);
		case 2:
			return ((T)13.0 * z*z)/( (T)9.0 + z*(z + (T)3.0) );
		case 3:
			return ((T)277.0 * z*z*z)/((T)225.0 + z*((T)45.0 + z*((T)6.0 + z)));
		case 4:
			return ((T)12746.0 * z*z*z*z)/((T)11025.0 + z*((T)1575.0 + z*((T)135.0 + z*((T)10.0 + z))) );
		case 5:
			return ((T)998881.0 * z*z*z*z*z)/((T)893025.0 + z*((T)99225.0 + z*((T)6300.0 + z*((T)315.0 + z*((T)15.0 + z)))));
		case 6:
			return ((T)118394977.0 * z*z*z*z*z*z)/((T)108056025.0 + z*((T)9823275.0 + z*((T)496125.0 + z*((T)18900.0 + z*((T)630.0 + z*((T)21.0 + z))))));
		case 7:
			return ((T)19727003738.0 * z*z*z*z*z*z*z)/((T)18261468225.0 + z*((T)1404728325.0 + z*((T)58939650.0 + z*((T)1819125.0 + z*((T)47250.0 + z*((T)1134.0 + z*((T)28.0 + z)))))));
		default:
			return (T) 1.0;
	}
	
}

/// ///////// ///
/// factorial ///
/// ///////// ///
__device__ int
cufac(int n)
{
	if (n == 0) {
		return 1;
	}
	if (n < 0) {
		return 0;
	}
	int cufac = 1;
	for (int i = 2; i <= n; i++) { 
		cufac = cufac * i; 
	}
	return cufac;
}

/// ///////////////// ///
/// absolute function ///
/// ///////////////// ///
template<class T>
__device__ double
cuabs(T x)
{
	return ((x < 0) ? -x : x);
}

/// ////// ///
/// (-1)^n ///
/// ////// ///
__device__ int
cupmo(int exponent)
{
	return ((exponent & 0x1) ? -1 : 1);
}

/// ///// ///
/// isOdd ///
/// ///// ///
__device__ bool
cuisOdd(int val)  ///< returns whether val is an odd number 
{
	return val & 0x1;
}

/// //// ///
/// swap ///
/// //// ///
template<class T>
__device__ void
cuswap (T& a, T& b)
{
	T temp = a;
	a = b;
	b = temp;
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
	return refl * P * cupmo((j - m) / 2);
}


/// !NOTE! spins and projection quantum numbers are in units of hbar/2

/// /////////////////////////// ///
/// Clebsch-Gordan-Coefficients ///
/// /////////////////////////// ///
template<class T>
__device__ T
cucgc(const int j1,
	  const int m1,
	  const int j2,
	  const int m2,
	  const int J,
	  const int M)
{
	if (cuisOdd(j1 - m1) || cuisOdd(j2 - m2) || cuisOdd(J - M)) {
		return 0;
	}
	if (cuabs(m1) > j1 || cuabs(m2) > j2 || cuabs(M) > J) {
		return 0;
	}
	if (m1 + m2 != M) {
		return 0;
	}
	if (j1 < 0 || j2 < 0 || J < 0) { // negative spins are not allowed
		return 0;
	}
	if (J < cuabs(j1 - j2) || J > j1 + j2) { // make sure J is in physical allowed range
		return 0;
	}
	if (cuisOdd(j1 + j2 - J)) { // check that J is in the half-integer or integer series, respectively
		return 0;
	}
	else {
		T clebschVal;
		int nu = 0;
		while ((j1 - j2 - M) / 2 + nu < 0 || (j1 - m1) / 2 + nu < 0) {
			nu++;
		}
		T sum = 0;
		int d1, d2, n1;
		while ( ((d1 = (J - j1 + j2) / 2 - nu) >= 0) and ((d2 = (J + M) / 2 - nu) >= 0) and ((n1 = (j2 + J + m1) / 2 - nu) >= 0) ) {
			const int d3 = (j1 - j2 - M) / 2 + nu;
			const int n2 = (j1 - m1) / 2 + nu;
			sum += ( (T)cupmo(nu + (j2 + m2) / 2) * cufac(n1) * cufac(n2) ) / ( cufac(nu) * cufac(d1) * cufac(d2) * cufac(d3));
			nu++;
		}
		if (sum == 0) {
			return 0;
		}
		const T N1 = cufac((J  + j1 - j2) / 2);
		const T N2 = cufac((J  - j1 + j2) / 2);
		const T N3 = cufac((j1 + j2 - J ) / 2);
		const T N4 = cufac((J + M) / 2);
		const T N5 = cufac((J - M) / 2);
		
		const T D0 = cufac((j1 + j2 + J) / 2 + 1);
		const T D1 = cufac((j1 - m1) / 2);
		const T D2 = cufac((j1 + m1) / 2);
		const T D3 = cufac((j2 - m2) / 2);
		const T D4 = cufac((j2 + m2) / 2);
		
		const T A  = ((J + 1) * N1 * N2 * N3 * N4 * N5) / (D0 * D1 * D2 * D3 * D4);
		
		clebschVal = ::sqrt(A) * sum;
		// if (clebschVal > 1) { clebschVal = clebschVal / 2 ; }
		return clebschVal;
	}
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
template<class T>
__device__ T
cudfunc(const int j,
	    const int m,
	    const int n,
	    const T theta)
{
	if (j < 0 || cuabs(m) > j || cuabs(n) > j) {
		return 0;
	}
	if (cuisOdd(j-m) || cuisOdd(j-n)) {
		return 0;
	}
	if (j == 0) {
		return 1;
	}
	
	//swap spinprojections for negative angle
	int _m = m;
	int _n = n;
	T thetaHalf = theta / 2;
	if (theta < 0) {
		thetaHalf = cuabs(thetaHalf);
		cuswap<int>(_m , _n);
	}
	
	const T cosThetaHalf = cos(thetaHalf);
	const T sinThetaHalf = sin(thetaHalf);
	
	const int jpm = (j + _m) / 2;
	const int jpn = (j + _n) / 2;
	const int jmm = (j - _m) / 2;
	const int jmn = (j - _n) / 2;
	const T kk = cufac(jpm) * cufac(jmm) * cufac(jpn) * cufac(jmn);
	const T constTerm = cupmo(jpm) * ::sqrt(kk);
	
	const int mpn  = (_m + _n) / 2;
	const int kMin = max(0, mpn);
	const int kMax = min(jpm, jpn);
	
	T sumTerm = 0;
	for (int k = kMin; k <= kMax; ++k) {
		const int kmn1 = 2 * k - (_m + _n) / 2;
		const int jmnk = j + (_m + _n) / 2 - 2 * k;
		const int jmk  = (j + _m) / 2 - k;
		const int jnk  = (j + _n) / 2 - k;
		const int kmn2 = k - (_m + _n) / 2;
		const T factor = cufac(k) * cufac(jmk) * cufac(jnk) * cufac(kmn2) * cupmo(k); /// why devide by? same effect as multiplying??
		sumTerm += ::pow(cosThetaHalf, kmn1) * ::pow(sinThetaHalf, jmnk) / factor;
	}
	
	return constTerm * sumTerm;
	
}


/// ////////// ///
/// D-Function ///
/// ////////// ///
template<class T>
__device__ cuda::complex<T>
cuDfunc(const int j,
		const int m,
		const int n,
		const T alpha,
		const T beta,
		const T gamma)
{
	const T arg = - (((T)m / 2) * alpha + ((T)n / 2) * gamma);
	return cuda::exp(cuda::complex<T>(0, arg)) * cudfunc(j, m, n, beta);
}

/// ///////////////////// ///
/// D-Function conjugated ///
/// ///////////////////// ///
template<class T>
__device__ cuda::complex<T>
cuDfuncConj(const int j,
			const int m,
			const int n,
			const T alpha,
			const T beta,
			const T gamma)
{
	return cuda::conj(cuDfunc(j, m, n, alpha, beta, gamma));
}

/// /////////////////// ///
/// Spherical Harmonics ///
/// /////////////////// ///
template<class T>
__device__ cuda::complex<T>
cuSpherHarm(const int l,
			const int m,
			const T theta,
			const T phi)
{
	const cuda::complex<T> Yarg(0, m * phi / 2);
	return ::sqrt((l + 1) / FOURPI) * cuda::exp(Yarg) * cudfunc(l, m, 0, theta);
}

/// //////////////// ///
/// D-Function Refl. ///
/// //////////////// ///
template<class T>
__device__ cuda::complex<T>
cuDfuncRefl(const int j,
			const int m,
			const int n,
			const int P,
			const int refl,
			const T alpha,
			const T beta,
			const T gamma)
{
	const int reflFactor = cuReflFactor(j, P, m, refl);
	if(m != 0) {
		return (cuDfunc(j, +m, n, alpha, beta, gamma) - (T)reflFactor * cuDfunc(j, -m, n, alpha, beta, gamma)) / ((T) ::sqrt(2.0));
	}
	else if (reflFactor == +1) {
		return 0;
	}
	else {
		return cuDfunc(j, 0, n, alpha, beta, gamma);
	}
}

/// ////////////////////// ///
/// D-Function Refl. Conj. ///
/// ////////////////////// ///
template<class T>
__device__ cuda::complex<T>
cuDfuncReflConj(const int j,
				const int m,
				const int n,
				const int P,
				const int refl,
				const T alpha,
				const T beta,
				const T gamma)
{
	return cuda::conj(cuDfuncRefl(j, m, n, P, refl, alpha, beta, gamma));
}

/// //////// ///
/// Gamma(m) ///
/// //////// ///
template<class T>
__device__ T
GammaFunc(T m, T m0, T Gamma0, T q, T q0, int L) 
{
	
	// int L = 8; /// insert angular momentum
	T BFValue = BFactor(L, q);
	T BF0 = BFactor(L, q0);
	
	// return Gamma0 * (m0 / *m) * (*q / q0) * (BFValue / BF0) * (BFValue / BF0);
	return ( Gamma0 * m0 * q  * BFValue ) / ( m * q0  * BF0 ); /// nearly 2 times (1.88) faster than the formula above, means first calculating the nominator and denominator and then dividing
	
}

/// //////////// ///
/// Breit-Wigner ///
/// //////////// ///
template<class T>
__device__ cuda::complex<T>
bwig(T m, T m0, T Gamma0, T q, T q0, int L)
{
	// double m0 = 0.77549; /// insert expected peak for mass
	// double Gamma0 = 0.1462; /// self explaining (?)
	cuda::complex<T> imag(0, 1);
	
	T GammaValue = GammaFunc<T>(m, m0, Gamma0, q, q0, L);
	
	return (Gamma0 * m0) / (m0 * m0 - m * m - imag * m0 * GammaValue);
}

/// ///////////////// ///
/// Amplitude test... ///
/// ///////////////// ///
/*
template<class T>
__device__ cuda::complex<T>
cuHelAmplitude(const T theta1,
			   const T phi1,
			   const T theta2,
			   const T phi2,
			   const T wX,
			   const T wf2,
			   const int JX,
			   const int MX,
			   const int J1,
			   const int L1,
			   const int S1,
			   const int J2,
			   const int L2,
			   const int S2)
{
	T wPIpm 	= 0.13957;
	T qX 	= cuBreakUpMomentum(wX,wf2,wPIpm);
	T qf2 	= cuBreakUpMomentum(wf2,wPIpm,wPIpm);
	T Gamma0 = 1;
	T q0 = 1;
	cuda::complex<T> amplsum(0,0);
	for(int lambda1 = -J1; lambda1 <= +J1; lambda1+=2) {
		for(int lambda2 = -J2; lambda2 <= +J2; lambda2+=2) {
			amplsum += 
				::sqrt(L1+1.0) 
				* cucgc(J1,lambda1,J2,lambda2,J1,lambda1-lambda2) 
				* cucgc(L1,0,S1,lambda1-lambda2,JX,lambda1-lambda2) 
				* cuDfuncConj(JX,MX,lambda1,theta1,phi1,0) 
				* ::sqrt((double)BFactor(L1,qX)) 
				* ::sqrt(L2+1.0) 
				* cucgc(0,0,0,0,0,0) 
				* cucgc(L2,0,0,0,L2,0) 
				* cuDfuncConj(4,lambda1,0,theta2,phi2,0) 
				* ::sqrt((double)BFactor(L2,qf2)) 
				* bwig(wf2,1.2754,Gamma0,qf2,q0,L2);
		}
	}
	return amplsum;
}
*/

template<class T>
__device__ cuda::complex<T>
cuHelAmplitude(const T theta1,
			   const T phi1,
			   const T theta2,
			   const T phi2,
			   const T wX,
			   const T wf2,
			   const T qX,
			   const T qf2,
			   const int JX,
			   const int MX,
			   const int J1,
			   const int L1,
			   const int S1,
			   const int J2,
			   const int L2,
			   const int S2)
{
	// double wPIpm = 0.13957;
	// double qX = 0.50514288800993423; //cuBreakUpMomentum(wX,0.587483,wPIpm);
	// double qf2 = 0.25846533761461093; //cuBreakUpMomentum(0.587483,wPIpm,wPIpm);
	// int lambda2 = 0;
	double Gamma0 = 0.1851;
	double q0 = 0.622239;
	// cuda::complex<double> bwigdat (0.183925,0.007331);
	
	cuda::complex<double> amplsum(0, 0);
	
	for(int lambda1 = -4; lambda1 <= +4; lambda1 += 2) {
		
		amplsum += 
			sqrt(L1 + (T)1) 
				* cucgc<T>(J1, lambda1, J2, 0, S1, lambda1) 
				* cucgc<T>(L1, 0, S1, lambda1, JX, lambda1) 
				* cuDfuncReflConj<T>(JX, MX, lambda1, -1, 1, phi1, theta1, 0) 
				* sqrt(BFactor(L1, qX)) 
				* sqrt(L2 + (T)1) 
				* cucgc<T>(0, 0, 0, 0, 0, 0) 
				* cucgc<T>(L2, 0, S2, 0, L2, 0) 
				* cuDfuncConj<T>(4, lambda1, 0, phi2, theta2, 0) 
				* sqrt(BFactor(L2, qf2)) 
				* bwig(wf2, (T)1.2754, Gamma0, qf2, q0, L2);
				
	}
	
	/*
	for(int lambda1 = -4; lambda1 <= 4; lambda1+=2) {
		
		amplsum += 
			sqrt(L1+1.0) 
			* cucgc(4,lambda1,0,0,4,lambda1)
			* cucgc(4,0,4,lambda1,4,lambda1)
			* cuDfuncReflConj(JX,MX,lambda1,-1,1,-0.564732,1.29395,0) 
			* sqrt(BFactor(L1,qX)) 
			* sqrt(L2+1.0) 
			* cucgc(0,0,0,0,0,0) 
			* cucgc(4,0,0,0,4,0) 
			* cuDfuncConj(4,lambda1,0,-1.6509,1.11757,0) 
			* sqrt(BFactor(L2,qf2)) 
			* bwig(0.58748315940941132,wf2,Gamma0,qf2,q0,L2);
			
	}
	*/
	
	return amplsum;
	
}

}

#endif
