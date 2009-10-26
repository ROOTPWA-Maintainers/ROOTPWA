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
//      Implementation of class FSParticle
//      see FSParticle.hh for details
//
// Environment:
//      Software developed for the COMPASS at CERN.
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

// This Class' Header ------------------
#include "FSParticle.h"

// C/C++ Headers ----------------------
#include <iostream>

// Collaborating Class Headers --------


// Class Member definitions -----------

bool operator== (const FSParticle& lhs, const FSParticle& rhs){
  return lhs.q()==rhs.q() && lhs.p()==rhs.p();
}


FSParticle::FSParticle()
  : _q(99)
{}
  
FSParticle::FSParticle(const TLorentzVector& p,
		       const TVector3& v,
		       double q)
: _p(p), _v(v), _q(q)
{}

