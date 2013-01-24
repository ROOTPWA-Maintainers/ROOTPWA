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
//
// Description:
//      Final state particle
//
//
// Environment:
//      Software developed for the COMPASS Experiment at CERN
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef FSPARTICLE_HH
#define FSPARTICLE_HH

// Collaborating Class Headers -------
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TVector3.h"
#include <assert.h>

// Collaborating Class Declarations --

class FSParticle {
public:

  // Constructors/Destructors ---------
  FSParticle();
  FSParticle(const TLorentzVector& p,
	     const TVector3& v,
	     int q);

  virtual ~FSParticle(){}

  // Operators
  friend bool operator== (const FSParticle& lhs, const FSParticle& rhs);

  // Accessors -----------------------
  const TLorentzVector& p() const {return _p;}
  const TVector3& v() const {return _v;}
  int q() const {return _q;}

  bool IsEqual(FSParticle* pi) {return (*this)==(*pi);}
  
  void Transform(const TLorentzRotation& L){_p*=L;}
 

private:

  // Private Data Members ------------
  TLorentzVector _p; // 4-momentum
  TVector3 _v; // vertex
  int _q; // charge


  // Private Methods -----------------
};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
