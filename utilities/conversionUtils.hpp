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
//      some simple conversion routines
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef CONVERSIONUTILS_HPP
#define CONVERSIONUTILS_HPP


#include <string>
#include <sstream>
#include <vector>


namespace rpwa {


	inline
	int
	sign(const char s)  ///< converts sign character into int
	{
		if (s == '+')
			return +1;
		if (s == '-')
			return -1;
		return 0;
	}


	template<typename T>
	inline
	std::string
	sign(const T& val)  ///< extracts sign string from a value
	{
		if (val < 0)
			return "-";
		if (val > 0)
			return "+";
		return "0";
	}


	template<typename T>
	inline
	const T*
	toArray(const std::vector<T>& vec)  ///< returns underlying array of vector
	{
		return &(*(vec.begin()));
	}

	template<typename T>
	inline
	T*
	toArray(std::vector<T>& vec)  ///< returns underlying array of vector
	{
		return &(*(vec.begin()));
	}

	template<typename T>
	inline
	const T*
	toArray(std::size_t&          size,
	        const std::vector<T>& vec)  ///< returns underlying array of vector and its size
	{
		size = vec.size();
		return &(*(vec.begin()));
	}

	template<typename T>
	inline
	T*
	toArray(std::size_t&    size,
	        std::vector<T>& vec)  ///< returns underlying array of vector and its size
	{
		size = vec.size();
		return &(*(vec.begin()));
	}


	template<typename T>
	inline
	std::string
	toString(const T& fromValue)  ///< converts any type into string that has output stream operator defined
	{
		std::ostringstream to;
		to << fromValue;
		return to.str();
	}


	inline
	std::string trueFalse(const bool val)  ///< converts bool to "true"/"false" string
	{
		if (val)
			return "true";
		else
			return "false";
	}

	inline
	std::string yesNo(const bool val)  ///< converts bool to "yes"/"no" string
	{
		if (val)
			return "yes";
		else
			return "no";
	}

	inline
	std::string onOff(const bool val)  ///< converts bool to "on"/"off" string
	{
		if (val)
			return "on";
		else
			return "off";
	}

	inline
	std::string enDisabled(const bool val)  ///< converts bool to "enabled"/"disabled" string
	{
		if (val)
			return "enabled";
		else
			return "disabled";
	}


	inline
	double
	spinQn(const int qn)  ///< convert internal spin quantum number into real quantum number in units of hbar
	{
		return 0.5 * qn;
	}


	inline
	std::string
	parityQn(const int qn)  ///< convert parity quantum number into string (+, -, or empty string)
	{
		return (qn != 0) ? sign(qn) : "";
	}
	

}  // namespace rpwa


#endif  // CONVERSIONUTILS_HPP
