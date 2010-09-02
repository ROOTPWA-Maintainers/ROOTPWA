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
//      some primitive standardized streams for reporting
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef REPORTINGUTILS_HPP
#define REPORTINGUTILS_HPP


#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <limits>
#include <vector>


namespace rpwa {


	//////////////////////////////////////////////////////////////////////////////
  // macros for printinf errors, warnings, and infos

  // cuts out block "className::methodName" from __PRETTY_FUNCTION__ output
  inline
  std::string
  getClassMethod__(std::string prettyFunction)
  {
    size_t pos = prettyFunction.find("(");
    if (pos == std::string::npos)
      return prettyFunction;           // something is not right
    prettyFunction.erase(pos);         // cut away signature
    pos = prettyFunction.rfind(" ");
    if (pos == std::string::npos)
      return prettyFunction;           // something is not right
    prettyFunction.erase(0, pos + 1);  // cut away return type
    return prettyFunction;
  }

#define printErr  std::cerr << "!!! " << __PRETTY_FUNCTION__ << " [" << __FILE__ << ":" << __LINE__ << "]: error: "   << std::flush
#define printWarn std::cerr << "??? " << __PRETTY_FUNCTION__ << " [" << __FILE__ << ":" << __LINE__ << "]: warning: " << std::flush
#define printInfo std::cout << ">>> " << getClassMethod__(__PRETTY_FUNCTION__) << "(): info: "  << std::flush


	//////////////////////////////////////////////////////////////////////////////
	// output stream manipulator that prints a value with its maximum precision

	template<typename T> class maxPrecisionValue__;

	template<typename T>
	inline
	maxPrecisionValue__<T>
	maxPrecision(const T& value)
	{ return maxPrecisionValue__<T>(value); }

	// output stream manipulator that prints a value with its maximum precision
	// in addition manipulator reserves space so that values will align
	template<typename T>
	inline
	maxPrecisionValue__<T>
	maxPrecisionAlign(const T& value)
	{ return maxPrecisionValue__<T>(value, maxPrecisionValue__<T>::ALIGN); }

	// output stream manipulator that prints a value with maximum precision for double
	template<typename T>
	inline
	maxPrecisionValue__<T>
	maxPrecisionDouble(const T& value)
	{ return maxPrecisionValue__<T>(value, maxPrecisionValue__<T>::DOUBLE); }

	// general helper class that encapsulates a value of type T
	template<typename T>
	class maxPrecisionValue__ {
	public:
		enum modeEnum { PLAIN,
		                ALIGN,
		                DOUBLE};  // forces precision for double
		maxPrecisionValue__(const T&       value,
		                    const modeEnum mode = PLAIN)
			: _value(value),
			  _mode (mode)
		{ }
		std::ostream& print(std::ostream& out) const
		{
			const int nmbDigits = (_mode != DOUBLE) ? std::numeric_limits<T>::digits10 + 1
				: std::numeric_limits<double>::digits10 + 1;
			std::ostringstream s;
			s.precision(nmbDigits);
			s.setf(std::ios_base::scientific, std::ios_base::floatfield);
			s << _value;
			switch (_mode) {
			case ALIGN:
				return out << std::setw(nmbDigits + 7) << s.str();  // make space for sign, dot, and exponent
			case PLAIN: case DOUBLE: default:
				return out << s.str();
			}
		}
	private:
		const T& _value;
		modeEnum _mode;
	};

	template<typename T>
	inline
	std::ostream& operator << (std::ostream&                 out,
	                           const maxPrecisionValue__<T>& value)
	{ return value.print(out); }


	//////////////////////////////////////////////////////////////////////////////
	// simple stream operators for some STL classes

	template<typename T1, typename T2>
	inline
	std::ostream&
	operator << (std::ostream&            out,
	             const std::pair<T1, T2>& pair)
	{
		return out << "(" << pair.first << ", " << pair.second << ")";
	}
	
	
	template<typename T>
	inline
	std::ostream&
	operator << (std::ostream&         out,
	             const std::vector<T>& vec)
	{
		out << "{";
		for (unsigned int i = 0; i < (vec.size() - 1); ++i)
			out << "[" << i << "] = " << vec[i] << ", ";
		return out << "[" << vec.size() - 1 << "] = " << vec[vec.size() - 1] << "}";
	}


}  // namespace rpwa


#endif  // REPORTINGUTILS_HPP
