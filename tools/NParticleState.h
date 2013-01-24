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
//      N-Particle state
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

#ifndef NPARTICLESTATE_HH
#define NPARTICLESTATE_HH


// Collaborating Class Headers ------
#include "TLorentzVector.h"
#include "TVector3.h"
#include <vector>

class FSParticle;


class NParticleState {
public:

  // Constructors/Destructors ---------
  NParticleState();
  virtual ~NParticleState();

  // Operators

  // Accessors -----------------------
  unsigned int n() const {return _n;}
  int q() const {return _q;} // total charge
  int qabs() const; // summed abs charge 
  TLorentzVector p() const; // Momentum of NParticleState
  TLorentzVector pfs(unsigned int i) const; // Momentum of ith daughter particle;
  FSParticle* getParticle(unsigned int i) const {return _fspart.at(i);}
  double Q2() const; // momentum transfer to target
  double t() const {return -Q2();}  // momentum transfer to target
  TVector3 vertex() const;  // Vertex
  const TLorentzVector& beam() const {return _beam;}
  double rapidity() const;

  bool Exclusive(double d=10);

  // Modifiers -----------------------
  bool addParticle(FSParticle* part); // returns false if double counting! 
  void setBeam(const TLorentzVector& beam);

  // Operations ----------------------
  bool isSubstate(const NParticleState* motherstate) const;
  bool isDisjunctFrom(const NParticleState* isobar) const;

private:

  // Private Data Members ------------
  unsigned int _n; // number of particles
  int _q; // total charge
  TLorentzVector _p; // total momentum
  TLorentzVector _beam; 
  std::vector<FSParticle*> _fspart; // 4-vectors final state particles

  // Private Methods -----------------
  
};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
