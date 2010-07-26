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
//      optimized Wigner d-function d^j_{m n} with caching
//     	based on PWA2000 function d_jmn_b() in pputil.cc
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef DFUNCTION_HPP
#define DFUNCTION_HPP


#include <cmath>
#include <algorithm>
#include <vector>

#include "utilities.h"
#include "mathUtils.hpp"
#include "factorial.hpp"


namespace rpwa {

	template<typename T>
	class dFunction {
	
	public:
			
		static dFunction& instance() { return _instance; }  ///< get singleton instance
		
		struct cacheEntryType {

			cacheEntryType()
				: valid    (false),
				  constTerm(0)
			{ }

			bool             valid;
			T                constTerm;
			std::vector<int> kmn1;
			std::vector<int> jmnk;
			std::vector<T>   factor;

		};

		T operator ()(const int j,
		              const int m,
		              const int n,
		              const T&  theta)  ///< returns d^j_{m n}(theta)
		{
			if (j >= (int)_maxJ) {
				printErr << "J = " << 0.5 * j << " is too large. maximum allowed J is "
				         << (_maxJ - 1) * 0.5 << ". aborting." << std::endl;
				throw;
			}
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

			const T cosThetaHalf = cos(thetaHalf);
			const T sinThetaHalf = sin(thetaHalf);

			T               dFuncVal   = 0;
			cacheEntryType& cacheEntry = _cache[j][j + _m][j + _n];
			if (cacheEntry.valid and _useCache) {
				// calculate function value using cache
				T sumTerm = 0;
				for (unsigned int i = 0; i < cacheEntry.factor.size(); ++i) {
					sumTerm +=   std::pow(cosThetaHalf, cacheEntry.kmn1[i])
						         * std::pow(sinThetaHalf, cacheEntry.jmnk[i]) / cacheEntry.factor[i];
				}
				dFuncVal = cacheEntry.constTerm * sumTerm;
			} else {
				// calculate function value and put values into cache
				const int jpm       = (j + _m) / 2;
				const int jpn       = (j + _n) / 2;
				const int jmm       = (j - _m) / 2;
				const int jmn       = (j - _n) / 2;
				const T   kk        =   rpwa::factorial<T>::instance()(jpm)
					                    * rpwa::factorial<T>::instance()(jmm)
					                    * rpwa::factorial<T>::instance()(jpn)
					                    * rpwa::factorial<T>::instance()(jmn);
				const T   constTerm = powMinusOne(jpm) * std::sqrt(kk);
				if (_useCache)
					cacheEntry.constTerm = constTerm;
				
				T         sumTerm = 0;
				const int mpn     = (_m + _n) / 2;
				const int kMin    = std::max(0,   mpn);
				const int kMax    = std::min(jpm, jpn);
				for (int k = kMin; k <= kMax; ++k) {
					const int kmn1   = 2 * k - (_m + _n) / 2;
					const int jmnk   = j + (_m + _n) / 2 - 2 * k;
					const int jmk    = (j + _m) / 2 - k;
					const int jnk    = (j + _n) / 2 - k;
					const int kmn2   = k - (_m + _n) / 2;
					const T   factor = (  rpwa::factorial<T>::instance()(k)
						                  * rpwa::factorial<T>::instance()(jmk)
						                  * rpwa::factorial<T>::instance()(jnk)
					                    * rpwa::factorial<T>::instance()(kmn2)) / powMinusOne(k);
					if (_useCache) {
						cacheEntry.kmn1.push_back(kmn1);
						cacheEntry.jmnk.push_back(jmnk);
						cacheEntry.factor.push_back(factor);
					}
					// using the 1 / factor here so that function value is the same as in PWA2000
					sumTerm += std::pow(cosThetaHalf, kmn1) * std::pow(sinThetaHalf, jmnk) / factor;
				}
				dFuncVal = constTerm * sumTerm;
				if (_useCache)
					cacheEntry.valid = true;
			}

			if (_debug)
				printInfo << "Wigner d^{J = " << 0.5 * j << "}" << "_{M = " << 0.5 * m << ", "
				          << "M' = " << 0.5 * n << "}" << "(theta = " << theta << ") = "
				          << maxPrecision(dFuncVal) << std::endl;
			return dFuncVal;
		}

		static bool useCache() { return _useCache; }                                   ///< returns caching flag
		static void setUseCache(const bool useCache = true) { _useCache = useCache; }  ///< sets caching flag

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag
    

	private:

		dFunction () { }
		~dFunction() { }
		dFunction (const dFunction&);
		dFunction& operator =(const dFunction&);

		static dFunction _instance;  ///< singleton instance
		static bool      _debug;     ///< if set to true, debug messages are printed
		static bool      _useCache;  ///< if set to true, cache is used

		static const unsigned int _maxJ = 41;                           ///< 2 * maximum allowed angular momentum + 1
		static cacheEntryType     _cache[_maxJ][2 * _maxJ][2 * _maxJ];  ///< cache for already calculated values [j][m][m']
		
	};


	template<typename T> dFunction<T> dFunction<T>::_instance;
	template<typename T> bool         dFunction<T>::_debug    = false;
	template<typename T> bool         dFunction<T>::_useCache = true;

	template<typename T> typename dFunction<T>::cacheEntryType
	  dFunction<T>::_cache[_maxJ][2 * _maxJ][2 * _maxJ];


	template<typename complexT, typename T>
	inline
	complexT
	DFunction(const int  j,
	          const int  m,
	          const int  n,
	          const T&   alpha,
	          const T&   beta,
	          const T&   gamma = T(0),
	          const bool debug = false)  ///< Wigner D-function D^j_{m n}(alpha, beta, gamma)
	{
		const T        arg      = ((T)m / 2) * alpha + ((T)n / 2) * gamma;
		const complexT DFuncVal =   rpwa::exp(complexT(0, -arg))
			                        * dFunction<T>::instance()(j, m, n, beta);
		if (debug)
			printInfo << "Wigner D^{J = " << 0.5 * j << "}" << "_{M = " << 0.5 * m << ", "
			          << "M' = " << 0.5 * n << "}" << "(alpha = " << alpha << ", beta = " << beta << ", "
			          << "gamma = " << gamma << ") = " << maxPrecisionDouble(DFuncVal) << std::endl;
		return DFuncVal;
	}


	template<typename complexT, typename T>
	inline
	complexT
	DFunctionConj(const int  j,
	              const int  m,
	              const int  n,
	              const T&   alpha,
	              const T&   beta,
	              const T&   gamma = T(0),
	              const bool debug = false)  ///< complex conjugate of Wigner D-function D^j_{m n}(alpha, beta, gamma)
	{
		const complexT DFuncVal = conj(DFunction<complexT>(j, m, n, alpha, beta, gamma, false));
    if (debug)
      printInfo << "Wigner D^{J = " << 0.5 * j << " *}" << "_{M = " << 0.5 * m << ", "
                << "M' = " << 0.5 * n << "}" << "(alpha = " << alpha << ", beta = " << beta << ", "
                << "gamma = " << gamma << ") = " << maxPrecisionDouble(DFuncVal) << std::endl;
		return DFuncVal;
	}


	template<typename complexT, typename T>
  inline
  complexT
  DFunctionRefl(const int  j,
                const int  m,
                const int  n,
                const int  P,
                const int  refl,
	              const T&   alpha,
	              const T&   beta,
	              const T&   gamma = T(0),
                const bool debug = false)  ///< Wigner D-function D^{j P refl}_{m n}(alpha, beta, gamma) in reflectivity basis
  {
    if (m < 0) {
      printWarn << "in reflectivity basis M = " << 0.5 * m << " < 0 is not allowed. "
                << "returning 0." << std::endl;
      return 0;
    }
    if (std::abs(refl) != 1) {
      printWarn << "reflectivity value epsilon = " << refl << " != +-1 is not allowed. "
                << "returning 0." << std::endl;
      return 0;
    }
    if (std::abs(P) != 1) {
      printWarn << "parity value P = " << P << " != +-1 is not allowed. "
                << "returning 0." << std::endl;
      return 0;
    }
    complexT  DFuncVal;
    const int reflFactor = refl * P * powMinusOne((j - m) / 2);
    if (m == 0)
	    DFuncVal = (reflFactor == +1) ? complexT(0) : DFunction<complexT>(j, 0, n, alpha, beta, gamma);
    else {
	    DFuncVal = 1 / std::sqrt((T)2)
		             * (                   DFunction<complexT>(j,  m, n, alpha, beta, gamma)
		                 - (T)reflFactor * DFunction<complexT>(j, -m, n, alpha, beta, gamma));
    }
    if (debug)
      printInfo << "Wigner D^{J = " << 0.5 * j << ", P = " << sign(P) << ", "
                << "refl = " << sign(refl) << "}" << "_{M = " << 0.5 * m << ", "
                << "M' = " << 0.5 * n << "}(alpha = " << alpha << ", " << "beta = " << beta << ", "
                << "gamma = " << gamma << ") = " << maxPrecisionDouble(DFuncVal) << std::endl;
    return DFuncVal;
  }

  
	template<typename complexT, typename T>
  inline
  complexT
  DFunctionReflConj(const int  j,
                    const int  m,
                    const int  n,
                    const int  P,
                    const int  refl,
                    const T&   alpha,
                    const T&   beta,
                    const T&   gamma = T(0),
                    const bool debug = false)  ///< complex conjugate of Wigner D-function D^{j P refl}_{m n}(alpha, beta, gamma) in reflectivity basis
  {
	  const complexT DFuncVal = conj(DFunctionRefl<complexT>(j, m, n, P, refl,
	                                                         alpha, beta, gamma, false));
    if (debug)
      printInfo << "Wigner D^{J = " << 0.5 * j << ", P = " << sign(P) << ", "
                << "refl = " << sign(refl) << " *}" << "_{M = " << 0.5 * m << ", "
                << "M' = " << 0.5 * n << "}(alpha = " << alpha << ", " << "beta = " << beta << ", "
                << "gamma = " << gamma << ") = " << maxPrecisionDouble(DFuncVal) << std::endl;
    return DFuncVal;
  }

  
}  // namespace rpwa


#endif  // DFUNCTION_HPP
