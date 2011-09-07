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


#include <cmath>
#include <algorithm>

#include <boost/tuple/tuple.hpp>

#include "mathUtils.hpp"
#include "clebschGordanCoeff.hpp"


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
		return sqrt((M - m1 - m2) * (M + m1 + m2) * (M - m1 + m2) * (M + m1 - m2)) / (2 * M);
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
		const double  mass   = sqrt(mass_2);
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
		const double p1     = sqrt(E1_2 - m_2[1]);
		const double p2     = sqrt(E2_2 - m_2[2]);
		const double Esum_2 = (E1 + E2) * (E1 + E2);
		if (min)
			return Esum_2 - (p1 + p2) * (p1 + p2);
		else
			return Esum_2 - (p1 - p2) * (p1 - p2);
	}


	//////////////////////////////////////////////////////////////////////////////
	// functions related to spin algebra
	// !NOTE! all spin quantum numbers are in units of hbar/2

	// computes minimum and maximum possible spins that can be made by
	// coupling spin A and spin B
	inline
	boost::tuples::tuple<int, int>
	getSpinRange(const int spinA,
	             const int spinB, 
	             bool*     valid = 0)
	{
		// make sure that allowedRange can always be directly used in for loops
		boost::tuples::tuple<int, int> allowedRange(0, -1);
		if ((spinA < 0) or (spinB < 0)) {
			if (valid)
				*valid = false;
			return allowedRange;
		}
		allowedRange = boost::tuples::make_tuple(std::abs(spinA - spinB), spinA + spinB);
		if (valid)
			*valid = true;
		return allowedRange;
	}

	// computes spin quantum number larger than or equal to lowerBound
	// in spin series defined by spinQn
	inline
	int
	spinQnLargerEqual(const int spinQn,
	                  const int lowerBound)
	{
		if (isOdd(spinQn))
			if (isOdd(lowerBound))
				return std::max(spinQn, lowerBound);
			else
				return std::max(spinQn, lowerBound + 1);
		else
			if (isEven(lowerBound))
				return std::max(spinQn, lowerBound);
			else
				return std::max(spinQn, lowerBound + 1);
	}

	// computes spin quantum number smaller than or equal to upperBound in
	// spin series defined by spinQn
	inline
	int
	spinQnSmallerEqual(const int spinQn,
	                   const int upperBound)
	{
		if (isOdd(spinQn))
			if (isOdd(upperBound))
				return std::min(spinQn, upperBound);
			else
				return std::min(spinQn, upperBound - 1);
		else
			if (isEven(upperBound))
				return std::min(spinQn, upperBound);
			else
				return std::min(spinQn, upperBound - 1);
	}

	// computes minimum and maximum possible spins that can be made by
	// coupling spin A and spin B taking into account externally defined
	// spin range
	inline
	boost::tuples::tuple<int, int>
	getSpinRange(const int                             spinA,
	             const int                             spinB,
	             const boost::tuples::tuple<int, int>& demandedRange,
	             bool*                                 valid = 0)
	{
		bool validRange;
		boost::tuples::tuple<int, int> allowedRange = getSpinRange(spinA, spinB, &validRange);
		if (not validRange) {
			if (valid)
				*valid = false;
			return allowedRange;
		}
		boost::tuples::get<0>(allowedRange) = spinQnLargerEqual (boost::tuples::get<0>(allowedRange),
		                                                         boost::tuples::get<0>(demandedRange));
		boost::tuples::get<1>(allowedRange) = spinQnSmallerEqual(boost::tuples::get<1>(allowedRange),
		                                                         boost::tuples::get<1>(demandedRange));
		if (valid) {
			if (boost::tuples::get<0>(allowedRange) <= boost::tuples::get<1>(allowedRange))
				*valid = true;
			else
				*valid = false;
		}
		return allowedRange;
	}


	// checks that spin and its projection quantum number are consistent
	inline
	bool
	spinAndProjAreCompatible(const int spin,
	                         const int spinProj)
	{ return (abs(spinProj) <= spin) and isEven(spin - spinProj); }


	// checks whether to daughter spin states can couple to given parent spin state
	inline
	bool
	spinStatesCanCouple(const int daughterSpins    [2],
	                    const int daughterSpinProjs[2],
	                    const int parentSpin,
	                    const int parentSpinProj)
	{
		//!!! the first two checks should probably be part of clebschGordanCoeff
		for (unsigned i = 0; i < 2; ++i)
			if (not spinAndProjAreCompatible(daughterSpins[i], daughterSpinProjs[i]))
				return false;
		if (not spinAndProjAreCompatible(parentSpin, parentSpinProj))
			return false;
		if (clebschGordanCoeff<double>(daughterSpins[0], daughterSpinProjs[0],
		                               daughterSpins[1], daughterSpinProjs[1],
		                               parentSpin,       parentSpinProj) == 0)
			return false;
		return true;		
	}


	// checks whether JPC combination is exotic
	inline
	bool
	jpcIsExotic(const int J,
	            const int P,
	            const int C)
	{
		// for baryons all JPCs are allowed
		if (isOdd(J))
			return false;
		// quark model restrictions for mesons: P == C is always allowed
		// check that P = (-1)^(J + 1)
		if (    (P != C)
		    and (   (C != ((J % 4     == 0) ? 1 : -1))
		         or (P != ((J + 2 % 4 == 0) ? 1 : -1))))
			return true;
		return false;
	}


	// checks whether JPG combination is exotic
	inline
	bool
	jpgIsExotic(const int J,
	            const int P,
	            const int G)
	{
		// for baryons all JPGs are allowed
		if (isOdd(J))
			return false;
		return false;
	}


} // namespace rpwa


#endif  // PHYSUTILS_H
