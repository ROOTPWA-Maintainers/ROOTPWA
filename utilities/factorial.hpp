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

#include "reportingUtils.hpp"


namespace rpwa {

	template<typename T>
	class factorialCached {
	
	public:
			
		static factorialCached& instance() { return _instance; }  ///< get singleton instance
		T operator ()(const unsigned int n)                       ///< returns n!
		{
			const unsigned int cacheSize = _cache.size();
			if (n >= cacheSize) {
				_cache.resize(n + 1);
				for (unsigned int i = cacheSize; i <= n; ++i) {
					// check for overflow
					if ((std::numeric_limits<T>::max() / (T)i) < _cache[i - 1]) {
						printErr << "target data type too small to hold " << i << "! "
						         << "(" << i - 1 << "! = " << _cache[i - 1] << " * " << i << " > "
						         << std::numeric_limits<T>::max() << "). aborting." << std::endl;
						throw;
					}
					_cache[i] = ((T)i) * _cache[i - 1];
				}
			}
			return _cache[n];
		}
    

	private:

		factorialCached () { }
		~factorialCached() { }
		factorialCached (const factorialCached&);
		factorialCached& operator =(const factorialCached&);

		static factorialCached _instance;  ///< singleton instance
 		static std::vector<T>  _cache;     ///< cache for already calculated values

	};


	template<typename T> factorialCached<T> factorialCached<T>::_instance;
	template<typename T> std::vector<T>     factorialCached<T>::_cache(1, 1);


	template<typename T>
	inline
	T
	factorial(const unsigned int n)  ///< returns factorial of n
	{ return factorialCached<T>::instance()(n);	}
	

}  // namespace rpwa


#endif  // FACTORIAL_HPP
