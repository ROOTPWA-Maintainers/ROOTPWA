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
//      collection of useful mathematical functions
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef MATHUTILS_H
#define MATHUTILS_H


#include <cmath>
#include <algorithm>
#include <complex>

#include <boost/math/tools/promotion.hpp>


namespace rpwa {


	// mathematical constants
	const double pi     = 2 * asin((double)1);
	const double piHalf = pi / 2;
	const double twoPi  = 2 * pi;
	const double fourPi = 4 * pi;

	const std::complex<double> imag(0, 1);
  
  
	//////////////////////////////////////////////////////////////////////////////
	// define aliases for some math functions so implementations may be switched easliy
	template<typename T> inline T abs(const T& x) { return std::abs (x); }
	template<typename T>
	inline typename boost::math::tools::promote_args<T>::type sqrt(const T& x) { return std::sqrt(x); }
	template<typename T>
	inline typename boost::math::tools::promote_args<T>::type exp (const T& x) { return std::exp (x); }
	template<typename T1, typename T2>
	inline typename boost::math::tools::promote_args<T1, T2>::type pow(const T1& base,
	                                                                   const T2& exponent)
	{ return std::pow(base, exponent); }


	//////////////////////////////////////////////////////////////////////////////
	// various small helper functions
	inline
	int
	powMinusOne(const int exponent)  ///< optimized function for (-1)^n
	{
		if (exponent & 0x1)  // exponent is odd
			return -1;
		else                 // exponent is even
			return +1;
	}

	
	template <typename T>
	inline
	bool
	isOdd(const T val)  ///< returns whether val is an odd number (assuming T is integer type)
	{
		return val & 0x1;
	}
  

	template <typename T>
	inline
	bool
	isEven(const T val)  ///< returns whether val is an even number (assuming T is integer type)
	{
		return not isOdd(val);
	}
  

} // namespace rpwa


#endif  // MATHUTILS_H
