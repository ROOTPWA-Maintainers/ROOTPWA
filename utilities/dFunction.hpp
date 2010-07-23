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
//      optimized Wigner d-function d^j_{m m'} with caching
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
#include <map>

#include "boost/tuple/tuple.hpp"

#include "utilities.h"
#include "factorial.hpp"


namespace rpwa {

	template<typename T>
	class dFunction {
	
	public:
			
		static dFunction& instance() { return _instance; }  ///< get singleton instance
		

		// define cache types
		typedef boost::tuple<int, int, int> qnTupleType;

		struct compareQnTuple {
			bool operator() (const qnTupleType& lhs,
			                 const qnTupleType& rhs) const
			{
				const int lhsJ = boost::get<0>(lhs);
				const int rhsJ = boost::get<0>(rhs);
				if (lhsJ != rhsJ)
					return lhsJ < rhsJ;
				const int lhsM = boost::get<1>(lhs);
				const int rhsM = boost::get<1>(rhs);
				if (lhsM != rhsM)
					return lhsM < rhsM;
				const int lhsN = boost::get<2>(lhs);
				const int rhsN = boost::get<2>(rhs);
				return lhsN < rhsN;
			}
		};
		
		struct cacheEntryType {
			T                constTerm;
			std::vector<int> kmn1;
			std::vector<int> jmnk;
			std::vector<T>   factor;
		};


		T operator ()(const int  j,
		              const int  m,
		              const int  n,
		              const T&   theta)  ///< returns d^j_{m n}(theta)
		{
			if ((j < 0) or (std::abs(m) > j) or (std::abs(n) > j)) {
				printErr << "illegal argument for Wigner d^{J = " << 0.5 * j << "}"
				         << "_{M = " << 0.5 * m << ", " << "M' = " << 0.5 * n << "}"
				         << "(theta = " << theta << "). aborting." << std::endl;
				throw;
			}
			
			// swap spin projections for negative angle
			qnTupleType qnTuple(j, m, n);
			T           thetaHalf = theta / 2;
			if (theta < 0) {
				thetaHalf = std::abs(thetaHalf);
				std::swap(boost::get<1>(qnTuple), boost::get<2>(qnTuple));
			}

			const T cosThetaHalf = cos(thetaHalf);
			const T sinThetaHalf = sin(thetaHalf);

			T             dFuncVal = 0;
			cacheIterator entry    = _cache.find(qnTuple);
			if (entry == _cache.end()) {
				// calculate function value and put values into cache
				cacheEntryType& cacheEntry = _cache[qnTuple];
				const int&      _m         = boost::get<1>(qnTuple);
				const int&      _n         = boost::get<2>(qnTuple);
			
				const int jpm       = (j + _m) / 2;
				const int jpn       = (j + _n) / 2;
				const int jmm       = (j - _m) / 2;
				const int jmn       = (j - _n) / 2;
				const T   kk        =   rpwa::factorial<T>::instance()(jpm)
					                    * rpwa::factorial<T>::instance()(jmm)
					                    * rpwa::factorial<T>::instance()(jpn)
					                    * rpwa::factorial<T>::instance()(jmn);
				const T   constTerm = powMinusOne(jpm) * std::sqrt(kk);
				cacheEntry.constTerm = constTerm;
				
				T         sumTerm = 0;
				const int mpn     = (_m + _n) / 2;
				const int kMin    = std::max(0,   mpn);
				const int kMax    = std::min(jpm, jpn);
				for (int k = kMin; k <= kMax; ++k) {
					const int kmn1 = 2 * k - (_m + _n) / 2;
					const int jmnk = j + (_m + _n) / 2 - 2 * k;
					cacheEntry.kmn1.push_back(kmn1);
					cacheEntry.jmnk.push_back(jmnk);
					const int jmk    = (j + _m) / 2 - k;
					const int jnk    = (j + _n) / 2 - k;
					const int kmn2   = k - (_m + _n) / 2;
					const T   factor = (  rpwa::factorial<T>::instance()(k)
						                  * rpwa::factorial<T>::instance()(jmk)
						                  * rpwa::factorial<T>::instance()(jnk)
					                    * rpwa::factorial<T>::instance()(kmn2)) / powMinusOne(k);
					cacheEntry.factor.push_back(factor);
					// using the 1 / factor here so that function value is the same as in PWA2000
					sumTerm += std::pow(cosThetaHalf, kmn1) * std::pow(sinThetaHalf, jmnk) / factor;
				}
				dFuncVal = constTerm * sumTerm;
			} else {
				// calculate function value using cache
				const cacheEntryType& cacheEntry = entry->second;
				T                     sumTerm    = 0;
				for (unsigned int i = 0; i < cacheEntry.factor.size(); ++i) {
					sumTerm +=   std::pow(cosThetaHalf, cacheEntry.kmn1[i])
						         * std::pow(sinThetaHalf, cacheEntry.jmnk[i]) / cacheEntry.factor[i];
				}
				dFuncVal = cacheEntry.constTerm * sumTerm;
			}

			if (_debug)
				printInfo << "Wigner d^{J = " << 0.5 * j << "}" << "_{M = " << 0.5 * m << ", "
				          << "M' = " << 0.5 * n << "}" << "(theta = " << theta << ") = "
				          << maxPrecision(dFuncVal) << std::endl;
			return dFuncVal;
		}

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag
    

	private:

		dFunction () { }
		~dFunction() { }
		dFunction (const dFunction&);
		dFunction& operator =(const dFunction&);

		static dFunction _instance;  ///< singleton instance
		static bool      _debug;     ///< if set to true, debug messages are printed
		
		static std::map<qnTupleType, cacheEntryType, compareQnTuple> _cache;  ///< cache for already calculated values
		typedef typename std::map<qnTupleType, cacheEntryType,
				compareQnTuple>::const_iterator cacheIterator;
		
	};


	template<typename T> dFunction<T>   dFunction<T>::_instance;
	template<typename T> bool           dFunction<T>::_debug = false;
	
	template<typename T> std::map<typename dFunction<T>::qnTupleType,
	                              typename dFunction<T>::cacheEntryType,
	                              typename dFunction<T>::compareQnTuple> dFunction<T>::_cache;


}  // namespace rpwa


#endif  // DFUNCTION_HPP
