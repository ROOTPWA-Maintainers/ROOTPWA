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
//      collection of useful physics functions and constants
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef PHYSUTILS_H
#define PHYSUTILS_H


#include "../pwa2000/libpp/pputil.h"
#include "mathUtils.hpp"
#include "reportingUtils.hpp"
#include "conversionUtils.hpp"


namespace rpwa {


	// computes breakup momentum of 2-body decay
	inline
	double
	breakupMomentum(const double M,   // mass of mother particle
	                const double m1,  // mass of daughter particle 1
	                const double m2)  // mass of daughter particle 2
	{
		if (M < m1 + m2)
			return 0;
		return rpwa::sqrt((M - m1 - m2) * (M + m1 + m2) * (M - m1 + m2) * (M + m1 - m2)) / (2 * M);
	}


	// kinematic border in Dalitz plot; PDG 2008 eq. 38.22a, b
	// for decay M -> m0 m1 m2
	inline
	double
	dalitzKinematicBorder(const double  mass_2,      // 2-body mass squared on x-axis
	                      const double  M,           // 3-body mass
	                      const double* m,           // array with the 3 daughter masses
	                      const bool    min = true)  // switches between curves for minimum and maximum mass squared on y-axis
	{
		if (mass_2 < 0)
			return 0;
		const double  mass   = rpwa::sqrt(mass_2);
		const double  M_2    = M * M;                                    // 3-body mass squared
		const double  m_2[3] = {m[0] * m[0], m[1] * m[1], m[2] * m[2]};  // daughter masses squared

		// calculate energies of particles 1 and 2 in m01 RF
		const double E1 = (mass_2 - m_2[0] + m_2[1]) / (2 * mass);
		const double E2 = (M_2    - mass_2 - m_2[2]) / (2 * mass);
		const double E1_2  = E1 * E1;
		const double E2_2  = E2 * E2;
		if ((E1_2 < m_2[1]) or (E2_2 < m_2[2]))
			return 0;

		// calculate m12^2
		const double p1     = rpwa::sqrt(E1_2 - m_2[1]);
		const double p2     = rpwa::sqrt(E2_2 - m_2[2]);
		const double Esum_2 = (E1 + E2) * (E1 + E2);
		if (min)
			return Esum_2 - (p1 + p2) * (p1 + p2);
		else
			return Esum_2 - (p1 - p2) * (p1 - p2);
	}


	inline
	double
	angMomNormFactor(const int  L,
	                 const bool debug = false)  ///< angular momentum normalization factor in amplitudes
	{
		const double norm = rpwa::sqrt(L + 1);
		if (debug)
			printDebug << "normalization factor sqrt(2 * L = " << spinQn(L) << " + 1) = "
			           << maxPrecision(norm) << std::endl;
		return norm;
	}


	inline
	int
	reflectivityFactor(const int J,
	                   const int P,
	                   const int M,
	                   const int refl)  ///< calculates prefactor for reflectibity symmetrization
	{
		if (rpwa::abs(P) != 1) {
			printWarn << "parity value P = " << P << " != +-1 is not allowed. "
			          << "returning 0." << std::endl;
			return 0;
		}
		if (M < 0) {
			printWarn << "in reflectivity basis M = " << spinQn(M) << " < 0 is not allowed. "
			          << "returning 0." << std::endl;
			return 0;
		}
		if (rpwa::abs(refl) != 1) {
			printWarn << "reflectivity value epsilon = " << refl << " != +-1 is not allowed. "
			          << "returning 0." << std::endl;
			return 0;
		}
		return refl * P * powMinusOne((J - M) / 2);
	}


	//////////////////////////////////////////////////////////////////////////////
	// some wrappers for libpp functions
	// !NOTE! all angular momenta and spin projections are in units of hbar/2
	inline
	double
	barrierFactor(const int    L,
	              const double breakupMom,
	              const bool   debug = false)  ///< Blatt-Weisskopf barrier factor
	{
		const double bf = F(L, breakupMom);
		if (debug)
			printDebug << "Blatt-Weisskopf barrier factor(L = " << spinQn(L) << ", "
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


#endif  // PHYSUTILS_H
