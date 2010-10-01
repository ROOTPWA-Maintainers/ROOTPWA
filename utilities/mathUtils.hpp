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

#include <boost/math/tools/promotion.hpp>

#include "pputil.h"
#include "utilities.h"


namespace rpwa {


	//////////////////////////////////////////////////////////////////////////////
	// define aliases for some math functions so implementations may be switched easliy
	template<typename T> inline T abs(const T& x) { return std::abs (x); }
	template<typename T>
	inline typename boost::math::tools::promote_args<T>::type sqrt(const T& x) { return std::sqrt(x); }
	template<typename T>
	inline typename boost::math::tools::promote_args<T>::type exp (const T& x) { return std::exp (x); }
	template<typename T1, typename T2>
	inline typename boost::math::tools::promote_args<T1, T2>::type pow(const T1& base,
	                                                                   const T2& exponent)
	{ return std::pow(base, exponent); }


	//////////////////////////////////////////////////////////////////////////////
	// various small helper functions
  inline
  double
  normFactor(const int  l,
             const bool debug = false)  ///< standard normalization factor in amplitudes
  {
	  const double norm = rpwa::sqrt(l + 1);
    if (debug)
      printInfo << "normalization factor sqrt(2 * L = " << 0.5 * l << " + 1) = "
                << maxPrecision(norm) << std::endl;
    return norm;
  }
  

	inline
	int
	powMinusOne(const int exponent)  ///< optimized function for (-1)^n
	{
		if (exponent & 0x1)  // exponent is odd
			return -1;
		else                 // exponent is even
			return +1;
	}

	
	inline
	int
	reflectivityFactor(const int j,
	                   const int P,
	                   const int m,
	                   const int refl)  ///< calculates prefactor for reflectibity symmetrization
	{
    if (rpwa::abs(P) != 1) {
      printWarn << "parity value P = " << P << " != +-1 is not allowed. "
                << "returning 0." << std::endl;
      return 0;
    }
    if (m < 0) {
      printWarn << "in reflectivity basis M = " << 0.5 * m << " < 0 is not allowed. "
                << "returning 0." << std::endl;
      return 0;
    }
    if (rpwa::abs(refl) != 1) {
      printWarn << "reflectivity value epsilon = " << refl << " != +-1 is not allowed. "
                << "returning 0." << std::endl;
      return 0;
    }
    return refl * P * powMinusOne((j - m) / 2);
	}


	//////////////////////////////////////////////////////////////////////////////
  // some wrappers for libpp functions
  // !NOTE! all angular momenta and spin projections are in units of hbar/2
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
	  const double Gamma  =   Gamma0 * (m0 / m) * (q / q0)
		                      * (rpwa::pow(F(L, q), 2) / rpwa::pow(F(L, q0), 2));
    return (m0 * Gamma0) / (m0 * m0 - m * m - imag * m0 * Gamma);
  }
  
  
} // namespace rpwa


#endif  // MATHUTILS_H
