//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
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
// $Rev::                        $: revision of last commit
// $Author::                     $: author of last commit
// $Date::                       $: date of last commit
//
// Description:
//      Angular momentum arithmetics class
//
//
// Author List:
//      Sebastian Neubert          TUM            (original author)
//
//
//-------------------------------------------------------------------------


// This Class' Header ------------------
#include "angMomCoupl.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------


// Class Member definitions -----------

using namespace rpwa;

std::vector<int> 
angMomCoupl::getCombinations(){
  std::vector<int> combo;
  for(int j=abs(_j1-_j2);j<=(_j1+_j2);j+=2)combo.push_back(j);
  return combo;
}
