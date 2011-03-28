///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
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
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////

/** @brief compound system with mass dependence
 */

#ifndef COMPOUND_HH
#define COMPOUND_HH

#include <vector>


// Base Class Headers ----------------


// Collaborating Class Declarations --


class absMassDep;

namespace rpwa { 

  
class compound {
public:

  compound();

  // Accessors -----------------------
  absMassDep* getMassDep() const {return _massDep;}
  int spin() const {return _j;}

  // Modifiers -----------------------

  // Operations ----------------------

private:

  // Private Data Members ------------
  absMassDep* _massDep;
  int _j;


  // Private Methods -----------------

};


}

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
