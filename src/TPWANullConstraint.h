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

//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Contraints an amplitude to 0
//
//
// Environment:
//      rootpwa
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef TPWANULLCONSTRAINT_HH
#define TPWANULLCONSTRAINT_HH

// Base Class Headers ----------------
#include "TPWAConstraint.h"

// Collaborating Class Headers -------
#include <complex>

// Collaborating Class Declarations --

class TPWANullConstraint : public TPWAConstraint {
public:

  // Constructors/Destructors ---------
  TPWANullConstraint(){}
  virtual ~TPWANullConstraint(){}

  virtual TPWAConstraint* clone(){return new TPWANullConstraint();}

   // Accessors -----------------------
  virtual int npar()const {return 0;} // returns number of free parameters 0,1 or even 2
  virtual std::string type()const {return "NullConstraint";}
  virtual std::string parname(unsigned int i) const {return "_NONE";}

   // Operations ----------------------
  virtual std::complex<double> cAmp(const std::complex<double>& amp)
  {return std::complex<double>(0,0);}
  virtual std::complex<double> dampdpar(unsigned int i)
  {return std::complex<double>(0,0);}

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
