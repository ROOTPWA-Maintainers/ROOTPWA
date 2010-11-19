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
#include <vector>


namespace rpwa {


	// converts sign character into int
	inline
	int
	sign(const char s)
	{
		if (s == '+')
			return +1;
		if (s == '-')
			return -1;
		return 0;
	}


	// extracts sign string from a value
	template<typename T>
	inline
	std::string
	sign(const T& val)
	{
		if (val < 0)
			return "-";
		if (val > 0)
			return "+";
		return "0";
	}


	// extracts sign from value
	template<typename T>
	inline
	T signum(const T& val)
	{
		if (val < 0)
			return -1;
		if (val > 0)
			return +1;
		return 0;
	}


	// conversion functions from std::vector to C array
	template<typename T>
	inline
	const T*
	toArray(const std::vector<T>& vec)
	{
		return &(*(vec.begin()));
	}

	template<typename T>
	inline
	T*
	toArray(std::vector<T>& vec)
	{
		return &(*(vec.begin()));
	}

	template<typename T>
	inline
	const T*
	toArray(std::size_t&          size,
	        const std::vector<T>& vec)
	{
		size = vec.size();
		return &(*(vec.begin()));
	}

	template<typename T>
	inline
	T*
	toArray(std::size_t&    size,
	        std::vector<T>& vec)
	{
		size = vec.size();
		return &(*(vec.begin()));
	}


	// converts bool to "true"/"false" string
	inline
	std::string trueFalse(const bool val)
	{
		if (val)
			return "true";
		else
			return "false";
	}

	// converts bool to "yes"/"no" string
	inline
	std::string yesNo(const bool val)
	{
		if (val)
			return "yes";
		else
			return "no";
	}

	// converts bool to "on"/"off" string
	inline
	std::string onOff(const bool val)
	{
		if (val)
			return "on";
		else
			return "off";
	}

	// converts bool to "enabled"/"disabled" string
	inline
	std::string enDisabled(const bool val)
	{
		if (val)
			return "enabled";
		else
			return "disabled";
	}


}  // namespace rpwa


#endif  // CONVERSIONUTILS_HPP
