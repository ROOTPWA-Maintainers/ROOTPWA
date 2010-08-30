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

#ifndef ABSMASSDEP_HH
#define ABSMASSDEP_HH

#include <vector>
#include <complex>

// Base Class Headers ----------------


// Collaborating Class Declarations --


namespace rpwa { 

typedef std::complex<double> cd;

class absMassDep {
public:
  
  virtual ~absMassDep(){}

  // Accessors -----------------------
  virtual cd val(double m)=0;
  virtual double get_rho(double m, unsigned int i)const =0;
  virtual int l(unsigned int i)const  =0;
  virtual int l()const  =0;

  // Modifiers -----------------------

  // Operations ----------------------


};

}

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
