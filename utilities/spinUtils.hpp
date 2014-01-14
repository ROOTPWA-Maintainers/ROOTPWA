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
//      functions related to spin algebra and functor that calculates
//      Clebsch-Gordan coefficients and caches results
//
//      !NOTE! spins and projection quantum numbers are in units of hbar/2
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef CLEBSCHGORDANCOEFF_HPP
#define CLEBSCHGORDANCOEFF_HPP


#include <vector>
#include <limits>

#include <boost/tuple/tuple.hpp>

#include "mathUtils.hpp"
#include "factorial.hpp"
#include "conversionUtils.hpp"


namespace rpwa {

	inline
	boost::tuples::tuple<int, int>
	getSpinRange(const int spinA,
	             const int spinB,
	             bool*     valid = 0)  ///< computes minimum and maximum possible spins that can be made by coupling spin A and spin B
	{
		// make sure that allowedRange can always be directly used in for loops
		boost::tuples::tuple<int, int> allowedRange(0, -1);
		if ((spinA < 0) or (spinB < 0)) {
			if (valid)
				*valid = false;
			return allowedRange;
		}
		allowedRange = boost::tuples::make_tuple(rpwa::abs(spinA - spinB), spinA + spinB);
		if (valid)
			*valid = true;
		return allowedRange;
	}


	inline
	int
	spinQnLargerEqual(const int spinQn,
	                  const int lowerBound)  ///< computes spin quantum number larger than or equal to lowerBound in spin series defined by spinQn
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


	inline
	int
	spinQnSmallerEqual(const int spinQn,
	                   const int upperBound)  ///< computes spin quantum number smaller than or equal to upperBound in spin series defined by spinQn
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


	inline
	boost::tuples::tuple<int, int>
	getSpinRange(const int                             spinA,
	             const int                             spinB,
	             const boost::tuples::tuple<int, int>& demandedRange,
	             bool*                                 valid = 0)  ///< computes minimum and maximum possible spins that can be made by coupling spin A and spin B taking into account externally defined spin range
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


	inline
	bool
	spinAndProjAreCompatible(const int spin,
	                         const int spinProj)  ///< checks that spin and its projection quantum number are consistent
	{
		return (spin >= 0) and (rpwa::abs(spinProj) <= spin) and isEven(spin - spinProj);
	}


	inline
	bool
	spinStatesCanCouple(const int j1,
	                    const int j2,
	                    const int J)  ///< returns, whether j1 and j2 can couple to J
	{
		bool validRange;
		boost::tuples::tuple<int, int> allowedRange = getSpinRange(j1, j2, &validRange);
		if (not validRange)
			return false;
		// make sure J is in physical allowed range
		if (   (J < boost::tuples::get<0>(allowedRange))
		    or (J > boost::tuples::get<1>(allowedRange)))
			return false;
		if (isOdd(j1 + j2 - J))  // check that J is in the half-integer or integer series, respectively
			return false;
		return true;
	}


	bool
	spinStatesCanCouple(const int j1,
	                    const int m1,
	                    const int j2,
	                    const int m2,
	                    const int J,
	                    const int M);  ///< checks whether two daughter spin states can couple to given parent spin state


	inline
	bool
	spinStatesCanCouple(const int daughterSpins    [2],
	                    const int daughterSpinProjs[2],
	                    const int parentSpin,
	                    const int parentSpinProj)  ///< checks whether two daughter spin states can couple to given parent spin state
	{
		return spinStatesCanCouple(daughterSpins[0], daughterSpinProjs[0],
		                           daughterSpins[1], daughterSpinProjs[1],
		                           parentSpin,       parentSpinProj);
	}


	inline
	bool
	jpcIsExotic(const int _J,
	            const int P,
	            const int C)  ///< checks whether JPC combination is exotic
	{
		// for baryons and flavoured mesons C is not defined
		if ((P == 0) or (C == 0) or (_J < 0) or isOdd(_J))
			return false;
		const int J = _J / 2;  // spin in units of hbar
		// quark model restrictions for mesons: P = (-1)^(L + 1), C = (-1)^(L + S)
		// for S = 1 is P == C; each J can be made with odd and even L,
		// except for J = 0, which requires L = 1 -> 0-- is forbidden
		// for S = 0 is J = L so that P = (-1)^(J + 1) and C = -P for all J
		if (J == 0) {
			if ((P == -1) and (C == -1))
				return true;
		}
		if ((C == -P) and (P != powMinusOne(J + 1)))
			return true;
		return false;
	}


	inline
	bool
	igjpIsExotic(const int _I,
	             const int G,
	             const int _J,
	             const int P)  ///< checks whether IGJP combination is exotic
	{
		// for baryons and flavoured mesons G is not defined
		if ((P == 0) or (G == 0) or (_I < 0) or isOdd(_I) or (_J < 0) or isOdd(_J))
			return false;
		const int I = _I / 2;  // isospin in units of hbar
		const int J = _J / 2;  // spin in units of hbar
		// quark model restrictions for mesons: P = (-1)^(L + 1), G = (-1)^(I + L + S)
		// for S = 1 is P == (-1)^I * G; each J can be made with odd and even L,
		// except for J = 0, which requires L = 1 -> 0-0- and 1+0- are forbidden
		// for S = 0 is J = L so that P = (-1)^(J + 1) and G = (-1)^(I + 1) * P for all J
		if (J == 0) {
			if ((P == -1) and (G == powMinusOne(I) * P))
				return true;
		}
		if (G == powMinusOne(I + 1) * P and P != powMinusOne(J + 1))
			return true;
		return false;
	}


	//////////////////////////////////////////////////////////////////////////////
	// Clebsch-Gordan coefficient functor
	template<typename T>
	class clebschGordanCoeffCached {

	public:

		static clebschGordanCoeffCached& instance() { return _instance; }  ///< get singleton instance
		T operator ()(const int j1,
		              const int m1,
		              const int j2,
		              const int m2,
		              const int J,
		              const int M)  ///< returns Clebsch-Gordan coefficient (j1 m1 j2 m2 | J M)
		{
			// check input parameters
			if (_useCache and ((j1 >= _maxJ) or (j2 >= _maxJ) or (J >= _maxJ))) {
				printErr << "spins are too large. maximum allowed spin is "
				         << spinQn(_maxJ - 1) << ". aborting." << std::endl;
				throw;
			}
			if (   not spinAndProjAreCompatible(j1, m1)
			    or not spinAndProjAreCompatible(j2, m2)
			    or not spinAndProjAreCompatible(J,  M )) {
				if (_debug)
					printWarn << "spins and spin projections are inconsistent: "
					          << "(j1 = " << spinQn(j1) << ", m1 = " << spinQn(m1) << ", "
					          << "j2 = "  << spinQn(j2) << ", m2 = " << spinQn(m2) << ", "
					          << "J = "   << spinQn(J)  << ", M = "  << spinQn(M)  << ")" << std::endl;
				return 0;
			}
			if (not spinStatesCanCouple(j1, j2, J)) {
				if (_debug)
					printWarn << "spins j1 = " << spinQn(j1) << " and j2 = " << spinQn(j2)
					          << " cannot couple to J = "  << spinQn(J) << std::endl;
				return 0;
			}
			if (m1 + m2 != M) {
				if (_debug)
					printWarn << "spin projections m1 = " << spinQn(m1) << " and m2 = " << spinQn(m2)
					          << " cannot couple to M = " << spinQn(M) << std::endl;
				return 0;
			}

			T        clebschVal = 0;
			double*& cacheEntry = _cache[j1][m1][j2][m2][J];
			if (_useCache and cacheEntry) {
				// calculate function value using cache
				clebschVal = *cacheEntry;
			} else {
				// calculate function value and put intermediate values into cache
				int nu = 0;
				while (    ((j1 - j2 - M) / 2 + nu < 0)
				        or ((j1 - m1)     / 2 + nu < 0))
					nu++;

				T   sum = 0;
				int d1, d2, n1;
				while (     ((d1 = (J - j1 + j2) / 2 - nu) >= 0)
				        and ((d2 = (J + M)       / 2 - nu) >= 0)
				        and ((n1 = (j2 + J + m1) / 2 - nu) >= 0)) {
					const int d3 = (j1 - j2 - M) / 2 + nu;
					const int n2 = (j1 - m1)     / 2 + nu;
					sum +=   powMinusOne(nu + (j2 + m2) / 2) * rpwa::factorial<T>(n1) * rpwa::factorial<T>(n2)
						     / (  rpwa::factorial<T>(nu) * rpwa::factorial<T>(d1)
						        * rpwa::factorial<T>(d2) * rpwa::factorial<T>(d3));
					nu++;
				}

				if (sum == 0)
					return 0;

				const T N1 = rpwa::factorial<T>((J  + j1 - j2) / 2);
				const T N2 = rpwa::factorial<T>((J  - j1 + j2) / 2);
				const T N3 = rpwa::factorial<T>((j1 + j2 - J ) / 2);
				const T N4 = rpwa::factorial<T>((J + M) / 2);
				const T N5 = rpwa::factorial<T>((J - M) / 2);

				const T D0 = rpwa::factorial<T>((j1 + j2 + J) / 2 + 1);
				const T D1 = rpwa::factorial<T>((j1 - m1) / 2);
				const T D2 = rpwa::factorial<T>((j1 + m1) / 2);
				const T D3 = rpwa::factorial<T>((j2 - m2) / 2);
				const T D4 = rpwa::factorial<T>((j2 + m2) / 2);

				const T A  = (J + 1) * N1 * N2 * N3 * N4 * N5 / (D0 * D1 * D2 * D3 * D4);

				clebschVal = rpwa::sqrt(A) * sum;
				if (_useCache) {
					cacheEntry  = new double;
					*cacheEntry = clebschVal;
				}
			}

			return clebschVal;
		}

		static bool useCache()                              { return _useCache;     }  ///< returns caching flag
		static void setUseCache(const bool useCache = true) { _useCache = useCache; }  ///< sets caching flag
		static unsigned int cacheSize()  ///< returns cache size in bytes
		{
			unsigned int size = _maxJ * (2 * _maxJ - 1) * _maxJ * (2 * _maxJ - 1) * 2 * _maxJ * sizeof(T*);
			for (int j1 = 0; j1 < _maxJ; ++j1)
				for (int m1 = -j1; m1 <= j1; ++m1)
					for (int j2 = 0; j2 < _maxJ; ++j2)
						for (int m2 = -j2; m2 <= j2; ++m2)
							for (int J = 0; J < 2 * _maxJ; ++J)
								if (_cache[j1][m1][j2][m2][J])
									size += sizeof(T);
			return size;
		}

		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	private:

		clebschGordanCoeffCached () { }
		~clebschGordanCoeffCached() { }
		clebschGordanCoeffCached (const clebschGordanCoeffCached&);
		clebschGordanCoeffCached& operator =(const clebschGordanCoeffCached&);

		static clebschGordanCoeffCached _instance;  ///< singleton instance
		static bool                     _useCache;  ///< if set to true, cache is used
		static bool                     _debug;     ///< if set to true, debug messages are printed

		static const int _maxJ = 18;  ///< maximum allowed angular momentum * 2 + 1
		static T*        _cache[_maxJ][2 * _maxJ - 1][_maxJ][2 * _maxJ - 1][2 * _maxJ];  ///< cache for intermediate terms [j1][m1][j2][m2][J]
	};


	template<typename T> clebschGordanCoeffCached<T> clebschGordanCoeffCached<T>::_instance;
	template<typename T> bool                        clebschGordanCoeffCached<T>::_useCache = true;
	template<typename T> bool                        clebschGordanCoeffCached<T>::_debug    = false;

	template<typename T> T*
	clebschGordanCoeffCached<T>::_cache[_maxJ][2 * _maxJ - 1][_maxJ][2 * _maxJ - 1][2 * _maxJ];


	template<typename T>
	inline
	T
	clebschGordanCoeff(const int  j1,
	                   const int  m1,
	                   const int  j2,
	                   const int  m2,
	                   const int  J,
	                   const int  M,
	                   const bool debug = false)  ///< returns Clebsch-Gordan coefficient (j1 m1 j2 m2 | J M)
	{
		const T cgCoeff = clebschGordanCoeffCached<T>::instance()(j1, m1, j2, m2, J, M);
		if (debug)
			printDebug << "Clebsch-Gordan (j1 = " << spinQn(j1) << ", m1 = " << spinQn(m1) << "; "
			           << "j2 = " << spinQn(j2) << ", m2 = " << spinQn(m2) << " | "
			           << "J = " << spinQn(J) << ", M = " << spinQn(M) << ") = "
		             << maxPrecisionDouble(cgCoeff) << std::endl;
		return cgCoeff;
	}
	//
	//////////////////////////////////////////////////////////////////////////////


	inline
	bool
	spinStatesCanCouple(const int j1,
	                    const int m1,
	                    const int j2,
	                    const int m2,
	                    const int J,
	                    const int M)
	{
		if (clebschGordanCoeff<double>(j1, m1, j2, m2, J, M) == 0)
			return false;
		return true;
	}


}  // namespace rpwa


#endif  // CLEBSCHGORDANCOEFF_HPP
