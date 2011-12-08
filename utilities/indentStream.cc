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
//      Boost.Iostream filter for output indentation
//      taken from Boost mailing list:
//      http://lists.boost.org/Archives/boost/2008/02/133679.php
//      copyright Roland Schwarz
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <iostream>

#include "indentStream.h"


#ifndef ROOT_CINT
rpwa::indentStream rpwa::cout(std::cout);
rpwa::indentStream rpwa::cerr(std::cerr);
#else
std::ostream rpwa::cout = std::cout;
std::ostream rpwa::cerr = std::cerr;
#endif
