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
//      optimized Wigner d-function d^j_{m n}(theta) with caching
//      used as basis for optimized spherical harmonics Y_l^m(theta, phi)
//      as well as for optimized D-function D^j_{m n}(alpha, beta,
//      gamma) and D-function in reflectivity basis
//
//     	based on PWA2000 function d_jmn_b() in pputil.cc
//
//      !NOTE! spin j and projections m and n are in units of hbar/2
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

#include "physUtils.hpp"
#include "factorial.hpp"


namespace rpwa {

	template<typename T>
	class dFunctionCached {

	public:

		struct cacheEntryType {

			cacheEntryType()
				: constTerm(0)
			{ }

			T                constTerm;
			std::vector<int> kmn1;
			std::vector<int> jmnk;
			std::vector<T>   factor;

		};

		typedef std::vector<std::vector<std::vector<cacheEntryType> > > cacheType;


		static dFunctionCached& instance() { return _instance; }  ///< get singleton instance

		T operator ()(const int j,
		              const int m,
		              const int n,
		              const T&  theta)  ///< returns d^j_{m n}(theta)
		{
			// check input parameters
			if (j >= (int)_maxJ) {
				printErr << "J = " << 0.5 * j << " is too large. maximum allowed J is "
				         << (_maxJ - 1) * 0.5 << ". aborting." << std::endl;
				throw;
			}
			if ((j < 0) or (rpwa::abs(m) > j) or (rpwa::abs(n) > j)) {
				printErr << "illegal argument for Wigner d^{J = " << 0.5 * j << "}"
				         << "_{M = " << 0.5 * m << ", M' = " << 0.5 * n << "}"
				         << "(theta = " << maxPrecision(theta) << "). aborting." << std::endl;
				throw;
			}

			// trivial case
			if (j == 0)
				return 1;

			// swap spin projections for negative angle
			int _m        = m;
			int _n        = n;
			T   thetaHalf = theta / 2;
			if (theta < 0) {
				thetaHalf = rpwa::abs(thetaHalf);
				std::swap(_m, _n);
			}

			const T cosThetaHalf = cos(thetaHalf);
			const T sinThetaHalf = sin(thetaHalf);

			T                dFuncVal   = 0;
			cacheEntryType*& cacheEntry = _cache[j][(j + _m) / 2][(j + _n) / 2];
			if (_useCache and cacheEntry) {
				// calculate function value using cache
				T sumTerm = 0;
				for (unsigned int i = 0; i < cacheEntry->factor.size(); ++i) {
					sumTerm +=   rpwa::pow(cosThetaHalf, cacheEntry->kmn1[i])
						         * rpwa::pow(sinThetaHalf, cacheEntry->jmnk[i]) / cacheEntry->factor[i];
				}
				dFuncVal = cacheEntry->constTerm * sumTerm;
			} else {
				// calculate function value and put intermediate values into cache
				if (_useCache)
					cacheEntry = new cacheEntryType();
				const int jpm       = (j + _m) / 2;
				const int jpn       = (j + _n) / 2;
				const int jmm       = (j - _m) / 2;
				const int jmn       = (j - _n) / 2;
				const T   kk        =   rpwa::factorial<T>(jpm)
					                    * rpwa::factorial<T>(jmm)
					                    * rpwa::factorial<T>(jpn)
					                    * rpwa::factorial<T>(jmn);
				const T   constTerm = powMinusOne(jpm) * rpwa::sqrt(kk);
				if (cacheEntry)
					cacheEntry->constTerm = constTerm;

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
					const T   factor = (  rpwa::factorial<T>(k   )
						                  * rpwa::factorial<T>(jmk )
						                  * rpwa::factorial<T>(jnk )
					                    * rpwa::factorial<T>(kmn2)) / powMinusOne(k);
					if (cacheEntry) {
						cacheEntry->kmn1.push_back  (kmn1);
						cacheEntry->jmnk.push_back  (jmnk);
						cacheEntry->factor.push_back(factor);
					}
					// using the 1 / factor here so that function value is the same as in PWA2000
					sumTerm += rpwa::pow(cosThetaHalf, kmn1) * rpwa::pow(sinThetaHalf, jmnk) / factor;
				}
				dFuncVal = constTerm * sumTerm;
			}

			return dFuncVal;
		}

		static bool useCache()                              { return _useCache;     }  ///< returns caching flag
		static void setUseCache(const bool useCache = true) { _useCache = useCache; }  ///< sets caching flag
		static unsigned int cacheSize()  ///< returns cache size in bytes
		{
			unsigned int size = _maxJ * (_maxJ + 1) * (_maxJ + 1) * sizeof(_cache[0][0][0]);
			for (int j = 0; j < (int)_maxJ; ++j)
				for (int m = -j; m <= j; m += 2)
					for (int n = -j; n <= j; n += 2) {
						cacheEntryType*& cacheEntry = _cache[j][(j + m) / 2][(j + n) / 2];
						if (cacheEntry) {
							size += sizeof(*cacheEntry);
							size += cacheEntry->kmn1.capacity  () * sizeof(int);
							size += cacheEntry->jmnk.capacity  () * sizeof(int);
							size += cacheEntry->factor.capacity() * sizeof(T);
						}
					}
			return size;
		}


	private:

		dFunctionCached () { }
		~dFunctionCached()
		{
			for (int j = 0; j < (int)_maxJ; ++j)
				for (int m = -j; m <= j; m += 2)
					for (int n = -j; n <= j; n += 2)
						delete _cache[j][(j + m) / 2][(j + n) / 2];
		}
		dFunctionCached (const dFunctionCached&);
		dFunctionCached& operator =(const dFunctionCached&);

		static dFunctionCached _instance;  ///< singleton instance
		static bool            _useCache;  ///< if set to true, cache is used

		static const unsigned int _maxJ = 41;                           ///< maximum allowed angular momentum * 2 + 1
		static cacheEntryType*    _cache[_maxJ][_maxJ + 1][_maxJ + 1];  ///< cache for intermediate terms [j][m][n]
	};


	template<typename T> dFunctionCached<T> dFunctionCached<T>::_instance;
	template<typename T> bool               dFunctionCached<T>::_useCache = true;

	template<typename T> typename dFunctionCached<T>::cacheEntryType*
	  dFunctionCached<T>::_cache[_maxJ][_maxJ + 1][_maxJ + 1];


	template<typename T>
	inline
	T
	dFunction(const int  j,
	          const int  m,
	          const int  n,
	          const T&   theta,
	          const bool debug = false)  ///< Wigner d-function d^j_{m n}(theta)
	{
		const T dFuncVal = dFunctionCached<T>::instance()(j, m, n, theta);
		if (debug)
			printDebug << "Wigner d^{J = " << 0.5 * j << "}_{M = " << 0.5 * m << ", "
			           << "M' = " << 0.5 * n << "}(theta = " << maxPrecision(theta) << ") = "
			           << maxPrecision(dFuncVal) << std::endl;
		return dFuncVal;
	}


	template<typename complexT>
  inline
  complexT
  sphericalHarmonic
	(const int                            l,
	 const int                            m,
	 const typename complexT::value_type& theta,
	 const typename complexT::value_type& phi,
	 const bool                           debug = false)  ///< spherical harmonics Y_l^{m}(theta, phi)
	{
		typedef typename complexT::value_type T;
		// crude implementation using Wigner d-function
		const complexT YVal =   rpwa::sqrt((l + 1) / fourPi)
			                    * rpwa::exp(complexT(0, ((T)m / 2) * phi)) * dFunction(l, m, 0, theta);
	  if (debug)
		  printDebug << "spherical harmonic Y_{l = " << 0.5 * l << "}^{m = " << 0.5 * m << "}"
		             << "(phi = " << maxPrecision(phi) << ", theta = " << maxPrecision(theta) << ") = "
		             << maxPrecisionDouble(YVal) << std::endl;
    return YVal;
  }


	template<typename complexT>
	inline
	complexT
	DFunction(const int                            j,
	          const int                            m,
	          const int                            n,
	          const typename complexT::value_type& alpha,
	          const typename complexT::value_type& beta,
	          const typename complexT::value_type& gamma = typename complexT::value_type(),
	          const bool                           debug = false)  ///< Wigner D-function D^j_{m n}(alpha, beta, gamma)
	{
		typedef typename complexT::value_type T;
		const T        arg      = ((T)m / 2) * alpha + ((T)n / 2) * gamma;
		const complexT DFuncVal = rpwa::exp(complexT(0, -arg)) * dFunction(j, m, n, beta);
		if (debug)
			printDebug << "Wigner D^{J = " << 0.5 * j << "}_{M = " << 0.5 * m << ", "
			           << "M' = " << 0.5 * n << "}(alpha = " << maxPrecision(alpha)
			           << ", beta = " << maxPrecision(beta) << ", gamma = " << maxPrecision(gamma)
			           << ") = " << maxPrecisionDouble(DFuncVal) << std::endl;
		return DFuncVal;
	}


	template<typename complexT>
	inline
	complexT
	DFunctionConj(const int                            j,
	              const int                            m,
	              const int                            n,
	              const typename complexT::value_type& alpha,
	              const typename complexT::value_type& beta,
	              const typename complexT::value_type& gamma = typename complexT::value_type(),
	              const bool                           debug = false)  ///< complex conjugate of Wigner D-function D^j_{m n}(alpha, beta, gamma)
	{
		const complexT DFuncVal = conj(DFunction<complexT>(j, m, n, alpha, beta, gamma, false));
		if (debug)
			printDebug << "Wigner D^{J = " << 0.5 * j << "}*_{M = " << 0.5 * m << ", "
			           << "M' = " << 0.5 * n << "}(alpha = " << maxPrecision(alpha)
			           << ", beta = " << maxPrecision(beta) << ", gamma = " << maxPrecision(gamma)
			           << ") = " << maxPrecisionDouble(DFuncVal) << std::endl;
		return DFuncVal;
	}


	template<typename complexT>
  inline
  complexT
  DFunctionRefl(const int                            j,
                const int                            m,
                const int                            n,
                const int                            P,
                const int                            refl,
	              const typename complexT::value_type& alpha,
	              const typename complexT::value_type& beta,
                const typename complexT::value_type& gamma = typename complexT::value_type(),
                const bool                           debug = false)  ///< Wigner D-function D^{j P refl}_{m n}(alpha, beta, gamma) in reflectivity basis
  {
		typedef typename complexT::value_type T;
    complexT  DFuncVal;
    const int reflFactor = reflectivityFactor(j, P, m, refl);
    if (m == 0) {
	    if (reflFactor == +1) {
		    DFuncVal = complexT(0);
		    // printWarn << "Wigner D^{J = " << 0.5 * j << ", P = " << sign(P) << ", "
		    //           << "refl = " << sign(refl) << "}_{M = " << 0.5 * m << ", "
		    //           << "M' = " << 0.5 * n << "}(alpha, beta, gamma) is always zero." << endl;
	    } else
		    DFuncVal = DFunction<complexT>(j, 0, n, alpha, beta, gamma, false);
    } else {
	    DFuncVal = 1 / rpwa::sqrt(2)
		    * (                  DFunction<complexT>(j, +m, n, alpha, beta, gamma, false)
		       - (T)reflFactor * DFunction<complexT>(j, -m, n, alpha, beta, gamma, false));
    }
    if (debug)
	    printDebug << "Wigner D^{J = " << 0.5 * j << ", P = " << sign(P) << ", "
	               << "refl = " << sign(refl) << "}_{M = " << 0.5 * m << ", "
	               << "M' = " << 0.5 * n << "}(alpha = " << maxPrecision(alpha) << ", "
	               << "beta = " << maxPrecision(beta) << ", gamma = " << maxPrecision(gamma)
	               << ") = " << maxPrecisionDouble(DFuncVal) << std::endl;
    return DFuncVal;
  }


	template<typename complexT>
  inline
  complexT
  DFunctionReflConj(const int                            j,
                    const int                            m,
                    const int                            n,
                    const int                            P,
                    const int                            refl,
                    const typename complexT::value_type& alpha,
                    const typename complexT::value_type& beta,
                    const typename complexT::value_type& gamma = typename complexT::value_type(),
                    const bool                           debug = false)  ///< complex conjugate of Wigner D-function D^{j P refl}_{m n}(alpha, beta, gamma) in reflectivity basis
  {
	  const complexT DFuncVal = conj(DFunctionRefl<complexT>(j, m, n, P, refl,
	                                                         alpha, beta, gamma, false));
    if (debug)
	    printDebug << "Wigner D^{J = " << 0.5 * j << ", P = " << sign(P) << ", "
	               << "refl = " << sign(refl) << "}*_{M = " << 0.5 * m << ", "
	               << "M' = " << 0.5 * n << "}(alpha = " << maxPrecision(alpha) << ", "
	               << "beta = " << maxPrecision(beta) << ", gamma = " << maxPrecision(gamma) << ") = "
	               << maxPrecisionDouble(DFuncVal) << std::endl;
    return DFuncVal;
  }


}  // namespace rpwa


#endif  // DFUNCTION_HPP
