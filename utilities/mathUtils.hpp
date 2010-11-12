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
//      collection of useful mathematical functions and constants
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


namespace rpwa {


	// mathematical constants
	const double pi     = 2 * asin((double)1);
	const double piHalf = pi / 2;
	const double twoPi  = 2 * pi;
	const double fourPi = 4 * pi;
  
  
	// computes n!
	inline
	unsigned int
	factorial(const unsigned int n)
	{
		unsigned int fac = 1;
		for (unsigned int i = 1; i <= n; ++i)
			fac *= i;
		return fac;
	}


} // namespace rpwa


#endif  // MATHUTILS_H
