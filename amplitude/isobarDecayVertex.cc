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


bool isobarDecayVertex::_debug = false;


isobarDecayVertex::isobarDecayVertex(const particlePtr&       mother,
				     const particlePtr&       daughter1,
				     const particlePtr&       daughter2,
				     const unsigned int       L,
				     const unsigned int       S,
				     const massDependencePtr& massDep)
  : interactionVertex(),
    _L                (L),
    _S                (S),
    _massDep          (massDep)
{
  if (!mother) {
    printErr << "null pointer to mother particle. aborting." << endl;
    throw;
  }
  if (!daughter1) {
    printErr << "null pointer to daughter 1 particle. aborting." << endl;
    throw;
  }
  if (!daughter2) {
    printErr << "null pointer to daughter 2 particle. aborting." << endl;
    throw;
  }
  interactionVertex::addInParticle (mother);
  interactionVertex::addOutParticle(daughter1);
  interactionVertex::addOutParticle(daughter2);
  if (!_massDep) {
    _massDep = createFlatMassDependence();
    if (_debug)
      printWarn << "null pointer to mass dependence. setting " << *_massDep << endl;
  }
  if (_debug)
    printInfo << "constructed " << *this << endl;
}


isobarDecayVertex::isobarDecayVertex(const isobarDecayVertex& vert)
{
  *this = vert;
}


isobarDecayVertex::~isobarDecayVertex()
{ }


isobarDecayVertex&
isobarDecayVertex::operator =(const isobarDecayVertex& vert)
{
  if (this != &vert) {
    interactionVertex::operator =(vert);
    _L       = vert._L;
    _S       = vert._S;
    _massDep = vert._massDep;
  }
  return *this;
}

    
isobarDecayVertex*
isobarDecayVertex::clone(const bool cloneInParticles,
			 const bool cloneOutParticles) const
{
  isobarDecayVertex* vertexClone = new isobarDecayVertex(*this);
  if (cloneInParticles)
    vertexClone->cloneInParticles();
  if (cloneOutParticles)
    vertexClone->cloneOutParticles();
  return vertexClone;
}


const TLorentzVector&
isobarDecayVertex::calcMotherLzVec()
{
  if (_debug)
    printInfo << "calculating Lorentz-vector of particle " << mother()->name()
	      << " before = " << mother()->lzVec() << " GeV, " << flush;
  mother()->setLzVec(daughter1()->lzVec() + daughter2()->lzVec());
  if (_debug)
    cout << "after = " << mother()->lzVec() << " GeV" << endl;
  return mother()->lzVec();
}


bool
isobarDecayVertex::checkMultiplicativeQn(const int     mQn,
					 const int     d1Qn,
					 const int     d2Qn,
					 const string& qnName)
{
  if (mQn == d1Qn * d2Qn) {
    if (_debug)
      printInfo << *this << ": " << ((qnName == "") ? "multiplicative quantum number" : qnName)
		<< " is consistent" << endl;
    return true;
  } else {
    if (_debug)
      printWarn << ((qnName == "") ? "multiplicative quantum number" : qnName) << " mismatch: "
		<< mother()->name() << " = " << mQn << " != " << d1Qn * d2Qn  << " = "
		<< "(" << daughter1()->name() << " = " << d1Qn << ") * "
		<< "(" << daughter2()->name() << " = " << d2Qn << ")" << endl;
    return false;
  }
}


bool
isobarDecayVertex::checkAdditiveQn(const int     mQn,
				   const int     d1Qn,
				   const int     d2Qn,
				   const string& qnName)
{
  if (mQn == d1Qn + d2Qn) {
    if (_debug)
      printInfo << *this << ": " << ((qnName == "") ? "additive quantum number" : qnName)
		<< " is consistent" << endl;
    return true;
  } else {
    if (_debug)
      printWarn << ((qnName == "") ? "additive quantum number" : qnName) << " mismatch: "
		<< mother()->name() << " = " << mQn << " != " << d1Qn + d2Qn  << " = "
		<< "(" << daughter1()->name() << " = " << d1Qn << ") + "
		<< "(" << daughter2()->name() << " = " << d2Qn << ")" << endl;
    return false;
  }
}


bool 
isobarDecayVertex::checkConsistency()
{
  bool vertexConsistent = true;
  // check multiplicative quantum numbers
  // G-parity
  if (!checkMultiplicativeQn(mother()->G(), daughter1()->G(), daughter2()->G(), "G-parity"))
    vertexConsistent = false;
  // C-parity
  const int cParity = mother()->G() * (mother()->isospin() % 4 == 0 ? 1 : -1);
  if (cParity != mother()->C()) {
    vertexConsistent = false;
    if (_debug)
      printWarn << "C-parity mismatch: " << mother()->name() << " = " << mother()->C()
		<< " != " << cParity  << " = "
		<< "(" << mother()->G() << ")^" << mother()->isospin() << " = G^2I"<< endl;
  } else if (_debug)
    printInfo << *this << ": C-parity is consistent" << endl;
  // parity
  const int angMomParity = (_L % 4 == 0) ? 1 : -1;  // modulo 4 because L is in units of hbar / 2
  if (mother()->P() != daughter1()->P() * daughter2()->P() * angMomParity) {
    if (_debug)
      printWarn << "parity mismatch: "
		<< mother()->name() << " = " << mother()->P() << " != "
		<< daughter1()->P() * daughter2()->P() * angMomParity  << " = "
		<< "(" << daughter1()->name() << " = " << daughter1()->P() << ") * "
		<< "(" << daughter2()->name() << " = " << daughter2()->P() << ") * "
		<< "(ang. momentum = " << angMomParity << ")" << endl;
    vertexConsistent = false;
  } else if (_debug)
    printInfo << *this << ": parity is consistent" << endl;
  // check additive quantum numbers
  // charge
  if (!checkAdditiveQn(mother()->charge(),      daughter1()->charge(),      daughter2()->charge(),      "charge"))
    vertexConsistent = false;
  // baryon number
  if (!checkAdditiveQn(mother()->baryonNmb(),   daughter1()->baryonNmb(),   daughter2()->baryonNmb(),   "baryonNmb"))
    vertexConsistent = false;
  // strangeness
  if (!checkAdditiveQn(mother()->strangeness(), daughter1()->strangeness(), daughter2()->strangeness(), "strangeness"))
    vertexConsistent = false;
  // charm
  if (!checkAdditiveQn(mother()->charm(),       daughter1()->charm(),       daughter2()->charm(),       "charm"))
    vertexConsistent = false;
  // beautty
  if (!checkAdditiveQn(mother()->beauty(),      daughter1()->beauty(),      daughter2()->beauty(),      "beauty"))
    vertexConsistent = false;
  // check angular momentum like quantum numbers
  // spin coupling: S in {|s1 - s2|, ..., s1 + s2}
  if (!angMomCoupl(daughter1()->J(), daughter2()->J()).inRange(_S)) {
    if(_debug)
      printWarn << "spins "
		<< "(" << daughter1()->name() << " J = " << daughter1()->J() * 0.5 << ") and "
		<< "(" << daughter2()->name() << " J = " << daughter2()->J() * 0.5 << ") "
		<< "cannot couple to total spin S = " << _S * 0.5 << endl;
    vertexConsistent = false;
  } else if (_debug)
    printInfo << *this << ": spin-spin coupling is consistent" << endl;
  // L-S coupling: J in {|L - S|, ..., L + S}
  if (!angMomCoupl(_L, _S).inRange(mother()->J())) {
    if (_debug)
      printWarn << "orbital angular momentum L = " << _L * 0.5 << " and spin S = " << _S * 0.5
		<< " cannot couple to angular momentum J = " << mother()->J() * 0.5 << endl;
    vertexConsistent = false;
  } else if (_debug)
    printInfo << *this << ": L-S coupling is consistent" << endl;
  // isospin coupling: I in {|I_1 - I_2|, ..., I_1 + I_2}
  if (!angMomCoupl(daughter1()->isospin(), daughter2()->isospin()).inRange(mother()->isospin())) {
    if (_debug)
      printWarn << "isospins "
		<< "(" << daughter1()->name() << " I = " << daughter1()->isospin() * 0.5 << ") and "
		<< "(" << daughter2()->name() << " I = " << daughter2()->isospin() * 0.5 << ") "
		<< "cannot couple to total isospin I = " << mother()->isospin() * 0.5 << endl;
    vertexConsistent = false;
  } else if (_debug)
    printInfo << *this << ": isospin coupling is consistent" << endl;
  //!!! missing: spin projections
  if (_debug) {
    if (vertexConsistent)
      printInfo << "vertex data are consistent: " << *this << endl;
    else
      printWarn << "vertex data are inconsistent (see warnings above): "
		<< *this << endl;
  }
  return vertexConsistent;
}


ostream&
isobarDecayVertex::print(ostream& out) const
{
  out << "isobar decay vertex: "
      << mother()->qnSummary() << "  --->  "
      << daughter1()->qnSummary()
      << "  [L = " << _L * 0.5 << ", S = " << _S * 0.5 << "]  "
      << daughter2()->qnSummary()
      << "; " << *_massDep;
  return out;
}


ostream&
isobarDecayVertex::dump(ostream& out) const
{
  out << "isobar decay vertex:" << endl
      << "    mother: "     << *mother()    << endl
      << "    daughter 1: " << *daughter1() << endl
      << "    daughter 2: " << *daughter2() << endl
      << "    L = " << _L * 0.5 << ", S = " << _S * 0.5 << endl
      << "    " << *_massDep << endl;
  return out;
}


ostream&
isobarDecayVertex::printPointers(ostream& out) const
{
  out << "isobar decay vertex "  << this        << ": "
      << "mother particle: "     << mother()    << "; "
      << "daughter 1 particle: " << daughter1() << "; "
      << "daughter 2 particle: " << daughter2() << endl;
  return out;
}
