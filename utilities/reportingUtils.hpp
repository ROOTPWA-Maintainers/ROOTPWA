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


namespace rpwa {


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

  // macros for printinf errors, warnings, and infos
#define printErr  std::cerr << "!!! " << __PRETTY_FUNCTION__ << " [" << __FILE__ << ":" << __LINE__ << "]: error: "   << std::flush
#define printWarn std::cerr << "??? " << __PRETTY_FUNCTION__ << " [" << __FILE__ << ":" << __LINE__ << "]: warning: " << std::flush
#define printInfo std::cout << ">>> " << getClassMethod__(__PRETTY_FUNCTION__) << "(): info: "  << std::flush


}  // namespace rpwa


#endif  // REPORTINGUTILS_HPP
