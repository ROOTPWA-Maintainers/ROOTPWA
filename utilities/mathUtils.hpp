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
//      collection of useful mathematical functions
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef MATHUTILS_H
#define MATHUTILS_H


#include <cmath>
#include <algorithm>
#include <complex>

#include "pputil.h"
#include "utilities.h"
#include "factorial.hpp"


namespace rpwa {


	// redefine some standard functions to make switching of complex data types easier
	template<typename T> inline T exp(const T& x) { return std::exp(x); }


  // some wrappers for libpp functions
  // !NOTE! all angular momenta and spin projections are in units of hbar/2
  inline
  double
  normFactor(const int  L,
             const bool debug = false)  ///< standard normalization factor in amplitudes
  {
	  const double norm = std::sqrt(L + 1);
    if (debug)
      printInfo << "normalization factor sqrt(2 * L = " << 0.5 * L << " + 1) = "
                << maxPrecision(norm) << std::endl;
    return norm;
  }
  

	inline
	int
	powMinusOne(const int exponent)  ///< optimized function for pow(-1, n)
	{
		if (exponent & 0x1)  // exponent is odd
			return -1;
		else                 // exponent is even
			return +1;
	}


	// based on PWA2000 implementation in pputil.cc
	template<typename T>
	T
	dFunc(const int  j,
	      const int  m,
	      const int  n,
	      const T&   theta,
	      const bool debug = false)  ///< Wigner d-function d^j_{m n}(theta)
	{
		if ((j < 0) or (std::abs(m) > j) or (std::abs(n) > j)) {
			printErr << "illegal argument for Wigner d^{J = " << 0.5 * j << "}"
			         << "_{M = " << 0.5 * m << ", " << "M' = " << 0.5 * n << "}"
			         << "(theta = " << theta << "). aborting." << std::endl;
			throw;
		}

		// swap spin projections for negative angle
		int _m        = m;
		int _n        = n;
		T   thetaHalf = theta / 2;
		if (theta < 0) {
			thetaHalf = std::abs(thetaHalf);
			std::swap(_m, _n);
		}

		const int j_p_m     = (j + _m) / 2;
		const int j_m_m     = (j - _m) / 2;
		const int j_p_n     = (j + _n) / 2;
		const int j_m_n     = (j - _n) / 2;
		const T   kk        =   rpwa::factorial<T>::instance()(j_p_m)
			                    * rpwa::factorial<T>::instance()(j_m_m)
			                    * rpwa::factorial<T>::instance()(j_p_n)
			                    * rpwa::factorial<T>::instance()(j_m_n);
		const T   constTerm = powMinusOne(j_p_m) * std::sqrt(kk);	
 
		const T   cosThetaHalf = cos(thetaHalf);
		const T   sinThetaHalf = sin(thetaHalf);
		const int m_p_n        = (_m + _n) / 2;
		const int kMin         = std::max(0,     m_p_n);
		const int kMax         = std::min(j_p_m, j_p_n);
		T         sumTerm      = 0;
		for (int k = kMin; k <= kMax; ++k) {
			const int kmn1 = 2 * k - (_m + _n) / 2;
			const int jmnk = j + (_m + _n) / 2 - 2 * k;
			const int jmk  = (j + _m) / 2 - k;
			const int jnk  = (j + _n) / 2 - k;
			const int kmn2 = k - (_m + _n) / 2;
			sumTerm += powMinusOne(k) * std::pow(cosThetaHalf, kmn1) * std::pow(sinThetaHalf, jmnk)
				         / (  rpwa::factorial<T>::instance()(k)   * rpwa::factorial<T>::instance()(jmk)
				            * rpwa::factorial<T>::instance()(jnk) * rpwa::factorial<T>::instance()(kmn2));
		}
		const T dFuncVal = constTerm * sumTerm;

		if (debug)
			printInfo << "Wigner d^{J = " << 0.5 * j << "}" << "_{M = " << 0.5 * m << ", "
			          << "M' = " << 0.5 * n << "}" << "(theta = " << theta << ") = "
			          << maxPrecision(dFuncVal) << std::endl;
		return dFuncVal;
	}


	template<typename complexT, typename scalarT>
	inline
	complexT
	DFunc(const int      j,
	      const int      m,
	      const int      n,
	      const scalarT& alpha,
	      const scalarT& beta,
	      const scalarT& gamma = 0,
	      const bool     debug = false)  ///< Wigner D-function D^j_{m n}(alpha, beta, gamma)
	{
		const scalarT  arg      = ((scalarT)m / 2) * alpha + ((scalarT)n / 2) * gamma;
		const complexT DFuncVal = rpwa::exp(complexT(0, -arg)) * dFunc(j, m, n, beta);
		if (debug)
			printInfo << "Wigner D^{J = " << 0.5 * j << "}" << "_{M = " << 0.5 * m << ", "
			          << "M' = " << 0.5 * n << "}" << "(alpha = " << alpha << ", beta = " << beta << ", "
			          << "gamma = " << gamma << ") = " << maxPrecisionDouble(DFuncVal) << std::endl;
		return DFuncVal;
	}


  inline
  std::complex<double>
  DFuncConj(const int    j,
            const int    m,
            const int    n,
            const double phi,
            const double theta,
            const bool   debug = false)  ///< conjugated Wigner D-function D^{J *}_{M lambda}(phi, theta, 0)
  {
    const std::complex<double> DFunc = conj(D(phi, theta, 0, j, m, n));
    if (debug)
      printInfo << "Wigner D^{J = " << 0.5 * j << " *}" << "_{M = " << 0.5 * m << ", "
                << "M' = " << 0.5 * n << "}" << "(alpha = " << phi << ", beta = " << theta << ", "
                << "gamma = 0) = " << maxPrecisionDouble(DFunc) << std::endl;
    return DFunc;
  }
  
  
  inline
  std::complex<double>
  DFuncConjRefl(const int    j,
                const int    m,
                const int    n,
                const int    P,
                const int    refl,
                const double phi,
                const double theta,
                const bool   debug = false)  ///< conjugated Wigner D-function {^epsilon}D^{J P *}_{M lambda}(phi, theta, 0) in reflectivity basis
  {
    if (m < 0) {
      printWarn << "in reflectivity basis M = " << 0.5 * m << " < 0 is not allowed. "
                << "returning 0." << std::endl;
      return 0;
    }
    if (abs(refl) != 1) {
      printWarn << "reflectivity value epsilon = " << refl << " != +-1 is not allowed. "
                << "returning 0." << std::endl;
      return 0;
    }
    const double               preFactor  = (m == 0 ? 0.5 : 1 / sqrt(2));
    const double               reflFactor = (double)refl * (double)P * pow(-1, 0.5 * (j - m));
    const std::complex<double> DFunc     
      =  preFactor * (               DFuncConj(j,  m, n, phi, theta)
                      - reflFactor * DFuncConj(j, -m, n, phi, theta));
    if (debug)
      printInfo << "Wigner D^{J = " << 0.5 * j << ", P = " << sign(P) << ", "
                << "refl = " << sign(refl) << " *}" << "_{M = " << 0.5 * m << ", "
                << "M' = " << 0.5 * n << "}(alpha = " << phi << ", "
                << "beta = " << theta << ", gamma = 0) = " << maxPrecisionDouble(DFunc) << std::endl;
    return DFunc;
  }

  
  inline
  double
  cgCoeff(const int  J1,
          const int  M1,
          const int  J2,
          const int  M2,
          const int  J,
          const int  M,
          const bool debug = false)  ///< Clebsch-Gordan coefficient (J1 M1, J2 M2 | J M)
  {
    const double cg = clebsch(J1, J2, J, M1, M2, M);
    if (debug)
      printInfo << "Clebsch-Gordan (J_1 = " << 0.5 * J1 << ", M_1 = " << 0.5 * M1 << "; "
                << "J_2 = " << 0.5 * J2 << ", M_2 = " << 0.5 * M2
                << " | J = " << 0.5 * J << ", M = " << 0.5 * M << ") = "
                << maxPrecision(cg) << std::endl;
    return cg;
  }
  
  
  inline
  double
  barrierFactor(const int    L,
                const double breakupMom,
                const bool   debug = false)  ///< Blatt-Weisskopf barrier factor
  {
    const double bf = F(L, breakupMom);
    if (debug)
      printInfo << "Blatt-Weisskopf barrier factor(L = " << 0.5 * L << ", "
                << "q = " << breakupMom << " GeV) = " << maxPrecision(bf) << std::endl;
    return bf;
  }

  
  inline
  std::complex<double>
  breitWigner(const double m,
              const double m0,
              const double Gamma0,
              const int    L,
              const double q,
              const double q0)  ///< relativistic Breit-Wigner with mass-dependent width
  {
    const double Gamma  = Gamma0 * (m0 / m) * (q / q0) * (pow(F(L, q), 2) / pow(F(L, q0), 2));
    return (m0 * Gamma0) / (m0 * m0 - m * m - imag * m0 * Gamma);
  }
  
  
} // namespace rpwa


#endif  // MATHUTILS_H
