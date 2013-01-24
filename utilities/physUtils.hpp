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


#include "mathUtils.hpp"
#include "reportingUtils.hpp"
#include "conversionUtils.hpp"


namespace rpwa {


	// computes squared breakup momentum of 2-body decay
	inline
	double
	breakupMomentumSquared(const double M,   // mass of mother particle
	                       const double m1,  // mass of daughter particle 1
	                       const double m2,  // mass of daughter particle 2
	                       const bool   allowSubThr = false)  // if set sub-threshold decays with negative return values are allowed
	{
		if (not allowSubThr and (M < m1 + m2)) {
			printErr << "mother mass " << M << " GeV/c^2 is smaller than sum of daughter masses "
			         << m1 << " + " << m2 << " GeV/c^2. this should never happen. aborting." << std::endl;
			throw;
		}
		return (M - m1 - m2) * (M + m1 + m2) * (M - m1 + m2) * (M + m1 - m2) / (4 * M * M);
	}


	// computes breakup momentum of 2-body decay
	inline
	double
	breakupMomentum(const double M,   // mass of mother particle
	                const double m1,  // mass of daughter particle 1
	                const double m2)  // mass of daughter particle 2
	{
		return rpwa::sqrt(breakupMomentumSquared(M, m1, m2, false));
	}


	// computes breakup momentum of 2-body decay
	// complex version with analytic continuation below threshold as used in K-matrix formalism
	inline
	std::complex<double>
	breakupMomentumComplex(const double M,   // mass of mother particle
	                       const double m1,  // mass of daughter particle 1
	                       const double m2)  // mass of daughter particle 2
	{
		const double q2 = rpwa::breakupMomentumSquared(M, m1, m2, true);
		const double q  = sqrt(fabs(q2));
		if (q2 < 0)
			return std::complex<double>(0, q);
		return std::complex<double>(q, 0);
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


	// computes square of Blatt-Weisskopf barrier factor for 2-body decay
	// !NOTE! L is units of hbar/2
	inline
	double
	barrierFactorSquared(const int    L,               // relative orbital angular momentum
	                     const double breakupMom,      // breakup momentum of 2-body decay [GeV/c]
	                     const bool   debug = false,
	                     const double Pr    = 0.1973)  // momentum scale 0.1973 GeV/c corresponds to 1 fm interaction radius
	{
		const double z   = (breakupMom * breakupMom) / (Pr * Pr);
		double       bf2 = 0;
		switch (L) {
		case 0:  // L = 0
			bf2 = 1;
			break;
		case 2:  // L = 1
			bf2 = (2 * z) / (z + 1);
			break;
		case 4:  // L = 2
			bf2 = (13 * z * z) / (z * (z + 3) + 9);
			break;
		case 6:  // L = 3
			bf2 = (277 * z * z * z) / (z * (z * (z + 6) + 45) + 225);
			break;
		case 8:  // L = 4
			{
				const double z2 = z * z;
				bf2 = (12746 * z2 * z2) / (z * (z * (z * (z + 10) + 135) + 1575) + 11025);
			}
			break;
		case 10:  // L = 5
			{
				const double z2 = z * z;
				bf2 = (998881 * z2 * z2 * z)
					/ (z * (z * (z * (z * (z + 15) + 315) + 6300) + 99225) + 893025);
			}
			break;
		case 12:  // L = 6
			{
				const double z3 = z * z * z;
				bf2 = (118394977 * z3 * z3)
					/ (z * (z * (z * (z * (z * (z + 21) + 630) + 18900) + 496125) + 9823275) + 108056025);
			}
			break;
		case 14:  // L = 7
			{
				const double z3 = z * z * z;
				bf2 = (19727003738LL * z3 * z3 * z)
					/ (z * (z * (z * (z * (z * (z * (z + 28) + 1134) + 47250) + 1819125) + 58939650)
					        + 1404728325L) + 18261468225LL);
			}
			break;
		default:
			printDebug << "calculation of Blatt-Weisskopf barrier factor is not (yet) implemented for L = "
			           << spinQn(L) << ". returning 0." << std::endl;
			return 0;
		}
		if (debug)
			printDebug << "squared Blatt-Weisskopf barrier factor(L = " << spinQn(L) << ", "
			           << "q = " << maxPrecision(breakupMom) << " GeV/c; P_r = " << Pr << " GeV/c) = "
			           << maxPrecision(bf2) << std::endl;
		return bf2;
	}


	// computes Blatt-Weisskopf barrier factor for 2-body decay
	// !NOTE! L is units of hbar/2
	inline
	double
	barrierFactor(const int    L,               // relative orbital angular momentum
	              const double breakupMom,      // breakup momentum of 2-body decay [GeV/c]
	              const bool   debug = false,
	              const double Pr    = 0.1973)  // momentum scale 0.1973 GeV/c corresponds to 1 fm interaction radius
	{
		const double bf = rpwa::sqrt(barrierFactorSquared(L, breakupMom, false, Pr));
		if (debug)
			printDebug << "Blatt-Weisskopf barrier factor(L = " << spinQn(L) << ", "
			           << "q = " << maxPrecision(breakupMom) << " GeV/c; P_r = " << Pr << " GeV/c) = "
			           << maxPrecision(bf) << std::endl;
		return bf;
	}


	// computes relativistic Breit-Wigner amplitude with mass-dependent width for 2-body decay
	// !NOTE! L is units of hbar/2
	inline
	std::complex<double>
	breitWigner(const double M,       // mass
	            const double M0,      // peak position
	            const double Gamma0,  // total width
	            const int    L,       // relative orbital angular momentum
	            const double q,       // 2-body breakup momentum
	            const double q0)      // 2-body breakup momentum at peak position
	{
		if (q0 == 0)
			return 0;
		const double Gamma  = Gamma0 * (M0 / M) * (q / q0)
			                    * (barrierFactorSquared(L, q) / barrierFactorSquared(L, q0));
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = M0 * Gamma0;
		const double B = M0 * M0 - M * M;
		const double C = M0 * Gamma;
		return (A / (B * B + C * C)) * std::complex<double>(B, C);
		// return (M0 * Gamma0) / (M0 * M0 - M * M - imag * M0 * Gamma);
	}


} // namespace rpwa


#endif  // PHYSUTILS_H
