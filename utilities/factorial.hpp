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
//      simple n! singleton functor with caching for small n
//      checks for overflows
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef FACTORIAL_HPP
#define FACTORIAL_HPP


#include <vector>
#include <limits>

#include "cudaUtils.hpp"
#include "reportingUtils.hpp"

#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime_api.h>
#endif

namespace rpwa {

	template<typename T>
	class factorialCached {

	public:

		struct cacheType {
			int numEntries;
			T* entries;
		};

		HOST_DEVICE
		static T calcFactorial(unsigned int n)
		{
			T result = 1;
			for(int i = 1; i <= n; ++i)
			{
				result *= i;
			}
			return result;
		}

		static void initCache()
		{
			if (_debug)
				printDebug << "adding ";
			if (_cache.size() > 0)
				return;
			_cache.clear();
			_cache.push_back(1);
			unsigned int i = 1;
			while((std::numeric_limits<T>::max() / (T)i) > _cache[i - 1]) {
				T newValue = ((T)i) * _cache[i - 1];
				if (_debug)
					std::cout << i << "! = " << newValue << "  ";
				_cache.push_back(newValue);
				++i;
			}
			if (_debug)
				std::cout << " to factorial cache" << std::endl;
#ifdef USE_CUDA

			// copy entries to device
			T* entries;
			cudaMalloc((void**)&entries, sizeof(T) * _cache.size());
			cudaMemcpy(entries, &(_cache[0]), sizeof(T) * _cache.size(), cudaMemcpyHostToDevice);

			// fill cache struct
			cacheType cudaCache;
			cudaCache.numEntries = _cache.size();
			cudaCache.entries = entries;

			// copy cache struct to device
			cudaMalloc((void**)&_cudaCache, sizeof(cacheType));
			cudaMemcpy(&_cudaCache, &cudaCache, sizeof(cacheType), cudaMemcpyHostToDevice);

#endif
		}

		HOST_DEVICE
		static T calculate(const unsigned int n, const cacheType* cudaCacheData = NULL)  ///< returns n!
		{
#ifdef __CUDA_ARCH__
			if (cudaCacheData == NULL || n >= cudaCacheData->numEntries)
				return T();
			return cudaCacheData->entries[n];
#else
			if (_cache.size() == 0)
				initCache();
			if (n >= _cache.size()) {
				printErr << "cannot calculate factorial of " << n << ": value out of range" << std::endl;
				throw;
			}
			return _cache[n];
#endif
		}

#ifdef USE_CUDA
		static const cacheType* getCudaCache() { return _cudaCache; }
#endif

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	private:

		factorialCached () { }
		~factorialCached() { }
		factorialCached (const factorialCached&);
		factorialCached& operator =(const factorialCached&);

 		static std::vector<T> _cache;     ///< cache for already calculated values
#ifdef USE_CUDA
 		static cacheType*     _cudaCache; ///< cache for already calculated values in device memory
#endif
		static bool           _debug;     ///< if set to true, debug messages are printed

	};

	template<typename T> std::vector<T>     factorialCached<T>::_cache;
#ifdef USE_CUDA
	template<typename T> typename factorialCached<T>::cacheType* factorialCached<T>::_cudaCache = NULL;
#endif
	template<typename T> bool               factorialCached<T>::_debug = false;


	template<typename T>
	inline
	void
	initFactorial()  ///< Initializes the cache containing all needed factorial values
	{ return factorialCached<T>::initCache();	}

#ifdef USE_CUDA
	template<typename T>
	inline
	const typename factorialCached<T>::cacheType*
	getFactorialCudaCache() ///< Returns the device version of the cache
	{ return factorialCached<T>::getCudaCache(); }
#endif


	template<typename T>
	HOST_DEVICE
	inline
	T
	factorial(const unsigned int n, const typename factorialCached<T>::cacheType* cudaCacheData = NULL)  ///< returns factorial of n
	{ return factorialCached<T>::calculate(n, cudaCacheData);	}


}  // namespace rpwa


#endif  // FACTORIAL_HPP
