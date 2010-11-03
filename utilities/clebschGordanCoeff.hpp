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
//      functor that calculates Clebsch-Gordan coefficients and caches results
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

#include "utilities.h"
#include "mathUtils.hpp"
#include "factorial.hpp"


namespace rpwa {

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
			if ((j1 < 0) || (j2 < 0) || (J < 0)) {
				printErr << "negative spins are not allowed (j1 = " << 0.5 * j1 << ", "
				         << "j2 = " << 0.5 * j2 << ", J = " << 0.5 * J << "). aborting." << std::endl;
				throw;
			}
			if (_useCache && ((j1 >= _maxJ) || (j2 >= _maxJ) || (J >= _maxJ))) {
				printErr << "spins are too large. maximum allowed spin is "
				         << (_maxJ - 1) * 0.5 << ". aborting." << std::endl;
				throw;
			}
			if (isOdd(j1 - m1) || isOdd(j2 - m2) || isOdd(J - M)) {
				printErr << "integer spin projection for half integer spin or vice versa "
				         << "(j1 = " << 0.5 * j1 << ", m1 = " << 0.5 * m1 << ", j2 = " << 0.5 * j2 << ", "
				         << "m2 = " << 0.5 * m2 << ", J = " << 0.5 * J << ", M = " << 0.5 * M << "). "
				         << "aborting." << std::endl;
				throw;
			}
			if ((rpwa::abs(m1) > j1) || (rpwa::abs(m2) > j2) || (rpwa::abs(M) > J))
				return 0;
			if (m1 + m2 != M)
				return 0;
			if ((J < rpwa::abs(j1 - j2)) || J > j1 + j2)
				return 0;


			T        clebschVal = 0;
			double*& cacheEntry = _cache[j1][m1][j2][m2][J];
			if (_useCache and cacheEntry) {
				// calculate function value using cache
				clebschVal = *cacheEntry;
			} else {
				// calculate function value and put intermediate values into cache
				int nu = 0;
				while (    ((j1 - j2 - M) / 2 + nu < 0)
				        || ((j1 - m1)     / 2 + nu < 0))
					nu++;

				T   sum = 0;
				int d1, d2, n1;
				while (    ((d1 = (J - j1 + j2) / 2 - nu) >= 0)
				        && ((d2 = (J + M)       / 2 - nu) >= 0)
				        && ((n1 = (j2 + J + m1) / 2 - nu) >= 0)) {
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


	private:

		clebschGordanCoeffCached () { }
		~clebschGordanCoeffCached() { }
		clebschGordanCoeffCached (const clebschGordanCoeffCached&);
		clebschGordanCoeffCached& operator =(const clebschGordanCoeffCached&);

		static clebschGordanCoeffCached _instance;  ///< singleton instance
		static bool                     _useCache;  ///< if set to true, cache is used

		static const int _maxJ = 18;  ///< maximum allowed angular momentum * 2 + 1
		static T*        _cache[_maxJ][2 * _maxJ - 1][_maxJ][2 * _maxJ - 1][2 * _maxJ];  ///< cache for intermediate terms [j1][m1][j2][m2][J]
	};


	template<typename T> clebschGordanCoeffCached<T> clebschGordanCoeffCached<T>::_instance;
	template<typename T> bool                        clebschGordanCoeffCached<T>::_useCache = true;

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
		  printInfo << "Clebsch-Gordan (j1 = " << 0.5 * j1 << ", m1 = " << 0.5 * m1 << "; "
		            << "j2 = " << 0.5 * j2 << ", m2 = " << 0.5 * m2 << " | "
		            << "J = " << 0.5 * J << ", M = " << 0.5 * M << ") = "
		            << maxPrecisionDouble(cgCoeff) << std::endl;
		return cgCoeff;
	}
	

}  // namespace rpwa


#endif  // CLEBSCHGORDANCOEFF_HPP
