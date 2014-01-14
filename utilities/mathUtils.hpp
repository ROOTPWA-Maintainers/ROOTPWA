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
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>


namespace rpwa {

	namespace ublas = boost::numeric::ublas;


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


	template<typename T>
	inline
	T signum(const T& val)  ///< extracts sign from value
	{
		if (val < 0)
			return -1;
		if (val > 0)
			return +1;
		return 0;
	}


	//////////////////////////////////////////////////////////////////////////////
	// matrix inversion routine using lu_factorize and lu_substitute
	// see http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?LU_Matrix_Inversion
	template<typename T>
	bool
	invertMatrix(const ublas::matrix<T>& A,
	             ublas::matrix<T>&       inverseA)
	{
		// create working copy of input
		ublas::matrix<T> M(A);
		// create permutation matrix for LU-factorization
		ublas::permutation_matrix<std::size_t> pM(M.size1());
		// perform LU-factorization
		if (ublas::lu_factorize(M, pM) != 0)
			return false;
		// create identity matrix of "inverse"
		inverseA.assign(ublas::identity_matrix<T>(M.size1()));
		// backsubstitute to get the inverse
		ublas::lu_substitute(M, pM, inverseA);
		return true;
	}


	template<typename T>
	ublas::matrix<T>
	invertMatrix(const ublas::matrix<T>& A,
	             bool&                   isSingular)
	{
		ublas::matrix<T> inverseA(A.size1(), A.size2());
		isSingular = !invert(A, inverseA);
		return inverseA;
	}


} // namespace rpwa


#endif  // MATHUTILS_H
