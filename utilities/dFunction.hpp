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

#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime_api.h>
#endif

#include "cudaUtils.hpp"
#include "physUtils.hpp"
#include "factorial.hpp"


namespace rpwa {

	template<typename T>
	class dFunctionCached {

	public:

		static const unsigned int _maxJ = 41; ///< maximum allowed angular momentum * 2 + 1

		// this type has to be a POD type because it is copied with memcpy
		struct cacheEntryType {
			T   constTerm;
			int numSumTerms;  // The number of valid entries in the following 3 arrays. (has to be <= _maxJ + 1)
			int kmn1  [_maxJ + 1];
			int jmnk  [_maxJ + 1];
			T   factor[_maxJ + 1];
		};

		struct cacheType {
			const typename factorialCached<T>::cacheType* cudaFactorialCache;
			cacheEntryType entries[_maxJ][_maxJ + 1][_maxJ + 1];  ///< cache for intermediate terms [j][m][n]
		};

		static T calculate(const int j,
				           const int m,
				           const int n,
				           const T&  theta)  ///< returns d^j_{m n}(theta)
		{

			if (_useCache and not _cacheIsInitialized) {
				printErr << "dFunction cache has not been initialized by calling initDFunction()" << std::endl;
				throw;
			}

			const cacheType* cacheData = _useCache ? (& _cache) : NULL;
			return calculate(j, m, n, theta, cacheData);

		}

		HOST_DEVICE
		static T calculate(const int j,
		                   const int m,
		                   const int n,
		                   const T&  theta,
		                   const cacheType* cacheData)  ///< returns d^j_{m n}(theta)
		{

			// check input parameters
			if (j >= (int)_maxJ) {
#ifndef __CUDA_ARCH__
				printErr << "J = " << 0.5 * j << " is too large. maximum allowed J is "
						 << (_maxJ - 1) * 0.5 << ". aborting." << std::endl;
				throw;
#else
				return 0;
#endif
			}

			if ((j < 0) or (rpwa::abs(m) > j) or (rpwa::abs(n) > j)) {
#ifndef __CUDA_ARCH__
				printErr << "illegal argument for Wigner d^{J = " << 0.5 * j << "}"
						 << "_{M = " << 0.5 * m << ", M' = " << 0.5 * n << "}"
						 << "(theta = " << maxPrecision(theta) << "). aborting." << std::endl;
				throw;
#else
				return 0;
#endif
			}

			const typename factorialCached<T>::cacheType* cudaFactorialCache = NULL;
#ifdef __CUDA_ARCH__
			cudaFactorialCache = cacheData->cudaFactorialCache;
#endif

			// trivial case
			if (j == 0)
				return 1;

			// swap spin projections for negative angle
			int _m        = m;
			int _n        = n;
			T   thetaHalf = theta / 2;
			if (theta < 0) {
				thetaHalf = rpwa::abs(thetaHalf);
				int temp = _m;
				_m = _n;
				_n = temp;
			}

			const T cosThetaHalf = cos(thetaHalf);
			const T sinThetaHalf = sin(thetaHalf);

			T               dFuncVal   = 0;
			if (cacheData) {
				// calculate function value using cache
				const cacheEntryType& cacheEntry = cacheData->entries[j][(j + _m) / 2][(j + _n) / 2];
				T sumTerm = 0;
				for (unsigned int i = 0; i < cacheEntry.numSumTerms; ++i) {
					sumTerm +=   rpwa::pow(cosThetaHalf, cacheEntry.kmn1[i])
						         * rpwa::pow(sinThetaHalf, cacheEntry.jmnk[i]) / cacheEntry.factor[i];
				}
				dFuncVal = cacheEntry.constTerm * sumTerm;
			} else {
				// calculate function value
				const int jpm       = (j + _m) / 2;
				const int jpn       = (j + _n) / 2;
				const int jmm       = (j - _m) / 2;
				const int jmn       = (j - _n) / 2;
				const T   kk        =   rpwa::factorial<T>(jpm, cudaFactorialCache)
					                    * rpwa::factorial<T>(jmm, cudaFactorialCache)
					                    * rpwa::factorial<T>(jpn, cudaFactorialCache)
					                    * rpwa::factorial<T>(jmn, cudaFactorialCache);
				const T   constTerm = powMinusOne(jpm) * rpwa::sqrt(kk);

				T         sumTerm = 0;
				const int mpn     = (_m + _n) / 2;
				const int kMin    = rpwa::max(0,   mpn);
				const int kMax    = rpwa::min(jpm, jpn);
				for (int k = kMin; k <= kMax; ++k) {
					const int kmn1   = 2 * k - (_m + _n) / 2;
					const int jmnk   = j + (_m + _n) / 2 - 2 * k;
					const int jmk    = (j + _m) / 2 - k;
					const int jnk    = (j + _n) / 2 - k;
					const int kmn2   = k - (_m + _n) / 2;
					const T   factor = (  rpwa::factorial<T>(k, cudaFactorialCache)
						                  * rpwa::factorial<T>(jmk, cudaFactorialCache)
						                  * rpwa::factorial<T>(jnk, cudaFactorialCache)
					                    * rpwa::factorial<T>(kmn2, cudaFactorialCache)) / powMinusOne(k);
					// using the 1 / factor here so that function value is the same as in PWA2000
					sumTerm += rpwa::pow(cosThetaHalf, kmn1) * rpwa::pow(sinThetaHalf, jmnk) / factor;
				}
				dFuncVal = constTerm * sumTerm;
			}

			return dFuncVal;
		}

		static bool         useCache   ()                           { return _useCache;      }  ///< returns caching flag
		static void         setUseCache(const bool useCache = true) { _useCache = useCache;  }  ///< sets caching flag
		static unsigned int cacheSize  ()                           { return sizeof(_cache); }  ///< returns cache size in bytes

		static void initCacheEntry(cacheEntryType& cacheEntry, int j, int m, int n)
		{
			const int jpm        = (j + m) / 2;
			const int jpn        = (j + n) / 2;
			const int jmm        = (j - m) / 2;
			const int jmn        = (j - n) / 2;
			const T   kk         = rpwa::factorial<T>(jpm)
									* rpwa::factorial<T>(jmm)
									* rpwa::factorial<T>(jpn)
									* rpwa::factorial<T>(jmn);
			cacheEntry.constTerm = powMinusOne(jpm) * rpwa::sqrt(kk);

			const int mpn     = (m + n) / 2;
			const int kMin    = std::max(0,   mpn);
			const int kMax    = std::min(jpm, jpn);
			cacheEntry.numSumTerms = kMax - kMin + 1;
			for (int k = kMin; k <= kMax; ++k) {
				const int kmn1   = 2 * k - (m + n) / 2;
				const int jmnk   = j + (m + n) / 2 - 2 * k;
				const int jmk    = (j + m) / 2 - k;
				const int jnk    = (j + n) / 2 - k;
				const int kmn2   = k - (m + n) / 2;
				const T   factor = (  rpwa::factorial<T>(k   )
									* rpwa::factorial<T>(jmk )
									* rpwa::factorial<T>(jnk )
									* rpwa::factorial<T>(kmn2)) / powMinusOne(k);
				int kRel = k - kMin;
				if (kRel >= j + 1) {
					printErr << "there are more than (j + 1) possible values for k, "
						  << "this should be impossible" << std::endl;
					throw;
				}
				cacheEntry.kmn1  [kRel] = kmn1;
				cacheEntry.jmnk  [kRel] = jmnk;
				cacheEntry.factor[kRel] = factor;
			}
		}

		static void initCache()
		{
			if (_cacheIsInitialized)
				return;

			if (_debug)
				printDebug << "filling dFunction cache" << std::endl;

			_cache.cudaFactorialCache = NULL;// only set in cuda version of this cache

			for (int j = 0; j < (int)_maxJ; ++j) {
				for (int m = -j; m <= j; m += 2) {
					for (int n = -j; n <= j; n += 2) {
						cacheEntryType& cacheEntry = _cache.entries[j][(j + m) / 2][(j + n) / 2];
						initCacheEntry(cacheEntry, j, m, n);
					}
				}
			}

#ifdef USE_CUDA
			_cache.cudaFactorialCache = getFactorialCudaCache<T>();
			cudaMalloc((void**)&_cudaCache, sizeof(cacheType));
			cudaMemcpy(_cudaCache, &_cache, sizeof(cacheType), cudaMemcpyHostToDevice);
			_cache.cudaFactorialCache = NULL;
#endif

			_cacheIsInitialized = true;
		}

#ifdef USE_CUDA
		static const cacheType* getCudaCache() { return _cudaCache; }
#endif

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	private:

		dFunctionCached () { }
		~dFunctionCached() { }
		dFunctionCached (const dFunctionCached&);
		dFunctionCached& operator =(const dFunctionCached&);

		static bool      _useCache;            ///< if set to true, cache is used
		static bool      _cacheIsInitialized;  ///< indicates, whether cache was initialized or not
		static bool      _debug;               ///< if set to true, debug messages are printed

		static cacheType _cache;

#ifdef USE_CUDA
		static cacheType* _cudaCache; ///< pointer to device memory containing the cache
#endif

	};

	template<typename T> bool dFunctionCached<T>::_useCache           = true;
	template<typename T> bool dFunctionCached<T>::_cacheIsInitialized = false;
	template<typename T> bool dFunctionCached<T>::_debug              = false;

	template<typename T> typename dFunctionCached<T>::cacheType dFunctionCached<T>::_cache;
#ifdef USE_CUDA
	template<typename T> typename dFunctionCached<T>::cacheType* dFunctionCached<T>::_cudaCache;
#endif

	template<typename T>
	inline
	void
	initDFunction() ///< Initializes the cache used for dFunction calculation
	{ dFunctionCached<T>::initCache(); }

#ifdef USE_CUDA
	template<typename T>
	inline
	const typename dFunctionCached<T>::cacheType*
	getDFunctionCudaCache() ///< Returns the device version of the cache
	{ return dFunctionCached<T>::getCudaCache(); }
#endif


	template<typename T>
	HOST_DEVICE
	inline
	T
	dFunction(const int  j,
	          const int  m,
	          const int  n,
	          const T&   theta,
	          const bool debug = false, ///< Wigner d-function d^j_{m n}(theta)
	          const typename dFunctionCached<T>::cacheType* cudaCacheData = NULL) ///< cuda cache memory, ignored when not in device code
	{
#ifndef __CUDA_ARCH__
		const T dFuncVal = dFunctionCached<T>::calculate(j, m, n, theta);
		if (debug)
			printDebug << "Wigner d^{J = " << 0.5 * j << "}_{M = " << 0.5 * m << ", "
			           << "M' = " << 0.5 * n << "}(theta = " << maxPrecision(theta) << ") = "
			           << maxPrecision(dFuncVal) << std::endl;
#else
		const T dFuncVal = dFunctionCached<T>::calculate(j, m, n, theta, cudaCacheData);
#endif
		return dFuncVal;
	}

	template<typename complexT>
	HOST_DEVICE
	inline
	complexT
	sphericalHarmonic(const int                            l,
			          const int                            m,
			          const typename complexT::value_type& theta,
			          const typename complexT::value_type& phi,
			          const bool                           debug = false,  ///< spherical harmonics Y_l^{m}(theta, phi)
    				  const typename dFunctionCached<typename complexT::value_type>::cacheType* cudaCacheData = NULL) // cuda cache memory, ignored when not in device code

	{
		typedef typename complexT::value_type T;

		double fourPi     = 4 *  2 * asin((double)1);
		// crude implementation using Wigner d-function
		const complexT YVal = rpwa::sqrt((l + 1) / fourPi)
		                      * rpwa::exp(complexT(0, ((T)m / 2) * phi)) * dFunction(l, m, 0, theta, false, cudaCacheData);
#ifndef __CUDA_ARCH__
		if (debug)
			printDebug << "spherical harmonic Y_{l = " << 0.5 * l << "}^{m = " << 0.5 * m << "}"
			           << "(phi = " << maxPrecision(phi) << ", theta = " << maxPrecision(theta) << ") = "
		               << maxPrecisionDouble(YVal) << std::endl;
#endif
		return YVal;
	}


	template<typename complexT>
	HOST_DEVICE
	inline
	complexT
	DFunction(const int                            j,
	          const int                            m,
	          const int                            n,
	          const typename complexT::value_type& alpha,
	          const typename complexT::value_type& beta,
	          const typename complexT::value_type& gamma = typename complexT::value_type(),
	          const bool                           debug = false,  ///< Wigner D-function D^j_{m n}(alpha, beta, gamma)
    		  const typename dFunctionCached<typename complexT::value_type>::cacheType* cudaCacheData = NULL) // cuda cache memory, ignored when not in device code
	{
		typedef typename complexT::value_type T;
		const T        arg      = ((T)m / 2) * alpha + ((T)n / 2) * gamma;
		const complexT DFuncVal = rpwa::exp(complexT(0, -arg)) * dFunction(j, m, n, beta, false, cudaCacheData);
#ifndef __CUDA_ARCH__
		if (debug)
			printDebug << "Wigner D^{J = " << 0.5 * j << "}_{M = " << 0.5 * m << ", "
			           << "M' = " << 0.5 * n << "}(alpha = " << maxPrecision(alpha)
			           << ", beta = " << maxPrecision(beta) << ", gamma = " << maxPrecision(gamma)
			           << ") = " << maxPrecisionDouble(DFuncVal) << std::endl;
#endif
		return DFuncVal;
	}


	template<typename complexT>
	HOST_DEVICE
	inline
	complexT
	DFunctionConj(const int                            j,
	              const int                            m,
	              const int                            n,
	              const typename complexT::value_type& alpha,
	              const typename complexT::value_type& beta,
	              const typename complexT::value_type& gamma = typename complexT::value_type(),
	              const bool                           debug = false,  ///< complex conjugate of Wigner D-function D^j_{m n}(alpha, beta, gamma)
    			  const typename dFunctionCached<typename complexT::value_type>::cacheType* cudaCacheData = NULL) // cuda cache memory, ignored when not in device code
	{
		const complexT DFuncVal = conj(DFunction<complexT>(j, m, n, alpha, beta, gamma, false, cudaCacheData));
#ifndef __CUDA_ARCH__
		if (debug)
			printDebug << "Wigner D^{J = " << 0.5 * j << "}*_{M = " << 0.5 * m << ", "
			           << "M' = " << 0.5 * n << "}(alpha = " << maxPrecision(alpha)
			           << ", beta = " << maxPrecision(beta) << ", gamma = " << maxPrecision(gamma)
			           << ") = " << maxPrecisionDouble(DFuncVal) << std::endl;
#endif
		return DFuncVal;
	}


	template<typename complexT>
	HOST_DEVICE
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
                  const bool                           debug = false,  ///< Wigner D-function D^{j P refl}_{m n}(alpha, beta, gamma) in reflectivity basis
    			  const typename dFunctionCached<typename complexT::value_type>::cacheType* cudaCacheData = NULL) // cuda cache memory, ignored when not in device code
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
				DFuncVal = DFunction<complexT>(j, 0, n, alpha, beta, gamma, false, cudaCacheData);
		} else {
			DFuncVal = 1 / rpwa::sqrt(2)
		    	* (                  DFunction<complexT>(j, +m, n, alpha, beta, gamma, false, cudaCacheData)
		    	   - (T)reflFactor * DFunction<complexT>(j, -m, n, alpha, beta, gamma, false, cudaCacheData));
		}
#ifndef __CUDA_ARCH__
		if (debug)
			printDebug << "Wigner D^{J = " << 0.5 * j << ", P = " << sign(P) << ", "
	               	   << "refl = " << sign(refl) << "}_{M = " << 0.5 * m << ", "
	               	   << "M' = " << 0.5 * n << "}(alpha = " << maxPrecision(alpha) << ", "
	               	   << "beta = " << maxPrecision(beta) << ", gamma = " << maxPrecision(gamma)
	               	   << ") = " << maxPrecisionDouble(DFuncVal) << std::endl;
#endif
		return DFuncVal;
	}


	template<typename complexT>
	HOST_DEVICE
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
                      const bool                           debug = false,  ///< complex conjugate of Wigner D-function D^{j P refl}_{m n}(alpha, beta, gamma) in reflectivity basis
    				  const typename dFunctionCached<typename complexT::value_type>::cacheType* cudaCacheData = NULL) // cuda cache memory, ignored when not in device code
	{
		const complexT DFuncVal = conj(DFunctionRefl<complexT>(j, m, n, P, refl,
	                                                         alpha, beta, gamma, false, cudaCacheData));
#ifndef __CUDA_ARCH__
		if (debug)
			printDebug << "Wigner D^{J = " << 0.5 * j << ", P = " << sign(P) << ", "
	               	   << "refl = " << sign(refl) << "}*_{M = " << 0.5 * m << ", "
	               	   << "M' = " << 0.5 * n << "}(alpha = " << maxPrecision(alpha) << ", "
	               	   << "beta = " << maxPrecision(beta) << ", gamma = " << maxPrecision(gamma) << ") = "
	               	   << maxPrecisionDouble(DFuncVal) << std::endl;
#endif
    	return DFuncVal;
	}


}  // namespace rpwa


#endif  // DFUNCTION_HPP
