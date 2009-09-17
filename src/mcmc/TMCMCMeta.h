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
//      Meta info for MCMC pwa
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef TMCMCMETA_HH
#define TMCMCMETA_HH

// Base Class Headers ----------------
#include "TObject.h"

// Collaborating Class Headers -------
#include <ostream> // remove if you do not need streaming op
#include <vector>
using std::vector;

#include "TString.h"
#include "TCMatrix.h"

// Collaborating Class Declarations --


class TMCMCMeta : public TObject {
public:

  // Constructors/Destructors ---------
  TMCMCMeta(){}
  virtual ~TMCMCMeta(){}

  // Operators
  //TMCMCMeta& operator=(const TMCMCMeta&);
  //friend bool operator== (const TMCMCMeta& lhs, const TMCMCMeta& rhs);
  //friend bool operator< (const TMCMCMeta& lhs, const TMCMCMeta& rhs);
  //friend std::ostream& operator<< (std::ostream& s, const TMCMCMeta& me);

  // Accessors -----------------------
  unsigned int npar() const {return parnames.size();}

  // Modifiers -----------------------


  // Operations ----------------------

  vector<TString> parnames;
  double stepsize;
  double virtmass;
  double low_mass;
  double high_mass;
  unsigned int Nleap;
  unsigned int NEvents; // number of events used in Likelihood
  TCMatrix Norm;

public:
  ClassDef(TMCMCMeta,2)
};



#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
