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

#include <cstdlib>
#include <vector>

namespace rpwa {

  class angMomCoupl {
  public:
    angMomCoupl(int j1, int j2) : _j1(j1), _j2(j2) {}
    
    bool inRange(int j) const {return (j>=abs(_j1-_j2) && j<=_j1+_j2);}
    std::vector<int> getCombinations();
    
  private:
    int _j1;   
    int _j2;
    
  };
  
}
