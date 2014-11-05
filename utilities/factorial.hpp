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
		}

		T operator ()(const unsigned int n)  ///< returns n!
		{
			if (_cache.size() == 0)
				initCache();
			if (n >= _cache.size()) {
				printErr << "cannot calculate factorial of " << n << ": value out of range" << std::endl;
				throw;
			}
			return _cache[n];
		}

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	private:

		factorialCached () { }
		~factorialCached() { }
		factorialCached (const factorialCached&);
		factorialCached& operator =(const factorialCached&);

		static factorialCached _instance;  ///< singleton instance
 		static std::vector<T>  _cache;     ///< cache for already calculated values
		static bool            _debug;     ///< if set to true, debug messages are printed

	};


	template<typename T> factorialCached<T> factorialCached<T>::_instance;
	template<typename T> std::vector<T>     factorialCached<T>::_cache;
	template<typename T> bool               factorialCached<T>::_debug = false;


	template<typename T>
	inline
	void
	initFactorial()  ///< Initializes the cache containing all needed factorial values
	{ return factorialCached<T>::initCache();	}


	template<typename T>
	inline
	T
	factorial(const unsigned int n)  ///< returns factorial of n
	{ return factorialCached<T>::instance()(n);	}


}  // namespace rpwa


#endif  // FACTORIAL_HPP
