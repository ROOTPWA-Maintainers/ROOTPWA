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
//      class that describes decay vertex of isobar into two particles
//      the isobar -> particle1 + particle 2 vertex has exactly one
//      incoming mother and two outgoing daughter particle
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "utilities.h"
#include "isobarDecayVertex.h"
#include "angMomCoupl.h"

	
using namespace std;
using namespace rpwa;


isobarDecayVertex::isobarDecayVertex(particle&          mother,
				     particle&          daughter1,
				     particle&          daughter2,
				     const unsigned int L,
				     const unsigned int S)
  : interactionVertex(),
    _L               (L),
    _S               (S)
{
  interactionVertex::addInParticle (mother);
  interactionVertex::addOutParticle(daughter1);
  interactionVertex::addOutParticle(daughter2);
}


isobarDecayVertex::isobarDecayVertex(const isobarDecayVertex& vert)
{
  *this = vert;
}


isobarDecayVertex::~isobarDecayVertex()
{ }


isobarDecayVertex&
isobarDecayVertex::operator = (const isobarDecayVertex& vert)
{
  if (this != &vert) {
    interactionVertex::operator = (vert);
    _L = vert._L;
    _S = vert._S;
  }
  return *this;
}


bool 
isobarDecayVertex::checkConsistency(){
  bool result=true;
  // check multiplicative QN: G,P, C is superflous if we use G??
  if(mother().G()!=daughter1().G()*daughter2().G()){
    if(debug())
      printWarn << "G-Parity mismatch in "<< endl;
    result=false;
  }
  int lparity= (_L % 4 == 0) ? 1 : -1; // modulo 4 because L in units of Hbar/2
  if(mother().P()!=daughter1().P()*daughter2().P() * lparity){
    if(debug())
      printWarn << "Parity mismatch in "<< endl;
    result=false;
  }
  // if(mother().C()!=daughter1().C()*daughter2().C()){
//     if(debug())
//       printWarn << "C-Parity mismatch in isobarDecay "<< *this << flush;
//     result=false;
//   }
  
  // check charge
  if(mother().charge()!=daughter1().charge()+daughter2().charge()){
    if(debug())
      printWarn << "Charge mismatch in "<< endl;
    result=false;
  }


  // check flavour (Isospin, Strangeness, Charm, Beauty) coupling

  


  // check spin coupling: s in {|s1-s2| .. s1+s2}
  
  if(!angMomCoupl(daughter1().J(),daughter2().J()).inRange(_S)){
    if(debug())
      printWarn << "Spin-Spin coupling out of range in "<< endl;
    result=false;
  }
  
  // check l-s coupling: J in {|l-s| .. l+s}
  
  if(!angMomCoupl(_L,_S).inRange(mother().J())){
    if(debug())
      printWarn << "Spin-Orbit coupling out of range in "<< endl;
    result=false;
  }

  if(!result && debug())printWarn <<  *this << flush;
    
  return result;

}



const TLorentzVector&
isobarDecayVertex::updateMotherLzVec()
{
  if (debug())
    printInfo << "updating Lorentz-vector of " << mother().name()
	      << " p_before = " << mother().lzVec() << " GeV, " << flush;
  mother().setLzVec(daughter1().lzVec() + daughter2().lzVec());
  cout << "p_after = " << mother().lzVec() << " GeV" << endl;
  return mother().lzVec();
}


ostream&
isobarDecayVertex::print(ostream& out) const
{
  out << "isobar decay vertex "
      << "(data are " << ((!dataAreValid()) ? "not " : "") << "valid):" << endl
      << "    mother "     << *(inParticles()[0])  << endl
      << "    daughter 1 " << *(outParticles()[0]) << endl
      << "    daughter 2 " << *(outParticles()[1]) << endl
      << "    L = " << _L << ", 2S = " << _S << endl;
  return out;
}
