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
//      incoming parent and two outgoing daughter particle
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <sstream>

#include "reportingUtilsRoot.hpp"
#include "spinUtils.hpp"
#include "isobarDecayVertex.h"

	
using namespace std;
using namespace rpwa;


bool isobarDecayVertex::_debug = false;


isobarDecayVertex::isobarDecayVertex(const particlePtr&       parent,
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
	if (not parent) {
		printErr << "null pointer to parent particle. aborting." << endl;
		throw;
	}
	if (not daughter1) {
		printErr << "null pointer to daughter 1 particle. aborting." << endl;
		throw;
	}
	if (not daughter2) {
		printErr << "null pointer to daughter 2 particle. aborting." << endl;
		throw;
	}
	interactionVertex::addInParticle (parent);
	interactionVertex::addOutParticle(daughter1);
	interactionVertex::addOutParticle(daughter2);
	if (not _massDep) {
		_massDep = createFlatMassDependence();
		if (_debug)
			printWarn << "null pointer to mass dependence. setting " << *_massDep << endl;
	}
	if (_debug)
		printDebug << "constructed " << *this << endl;
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
isobarDecayVertex::doClone(const bool cloneInParticles,
                           const bool cloneOutParticles) const
{
	isobarDecayVertex* vertexClone = new isobarDecayVertex(*this);
	if (cloneInParticles)
		vertexClone->cloneInParticles();
	if (cloneOutParticles)
		vertexClone->cloneOutParticles();
	if (_debug)
		printDebug << "cloned " << *this << "; " << this << " -> " << vertexClone << " "
		           << ((cloneInParticles ) ? "in" : "ex") << "cluding incoming particles, "
		           << ((cloneOutParticles) ? "in" : "ex") << "cluding outgoing particles" << std::endl;
	return vertexClone;
}


bool
isobarDecayVertex::isEqualTo(const interactionVertex& vert) const
{
	const isobarDecayVertex* isoVert = dynamic_cast<const isobarDecayVertex*>(&vert);
	if (not isoVert)
		return false;
	if (not interactionVertex::isEqualTo(vert))
		return false;
	if (   (L()                 != isoVert->L())
	    or (S()                 != isoVert->S())
	    or (*(massDependence()) != *(isoVert->massDependence())))
		return false;
	return true;
}


bool
isobarDecayVertex::addInParticle(const particlePtr&)
{
	if (_debug)
		printWarn << "cannot add incoming particle to " << *this << endl;
	return false;
}


bool
isobarDecayVertex::addOutParticle(const particlePtr&)
{
	if (_debug)
		printWarn << "cannot add outgoing particle to " << *this << endl;
	return false;
}


const TLorentzVector&
isobarDecayVertex::calcParentLzVec()
{
	if (_debug)
		printDebug << "calculating Lorentz-vector of parent particle " << parent()->name()
		           << " before = " << parent()->lzVec() << " GeV, " << flush;
	parent()->setLzVec(daughter1()->lzVec() + daughter2()->lzVec());
	if (_debug)
		cout << "after = " << parent()->lzVec() << " GeV" << endl;
	return parent()->lzVec();
}


int
isobarDecayVertex::calcParentCharge()
{
	if (_debug)
		printDebug << "calculating charge of parent particle " << parent()->name()
		           << " before = " << parent()->charge() << ", " << flush;
	parent()->setCharge(daughter1()->charge() + daughter2()->charge());
	if (_debug)
		cout << "after = " << parent()->charge() << endl;
	return parent()->charge();
}


int
isobarDecayVertex::calcParentBaryonNmb()
{
	if (_debug)
		printDebug << "calculating baryon number of parent particle " << parent()->name()
		           << " before = " << parent()->baryonNmb() << ", " << flush;
	parent()->setBaryonNmb(daughter1()->baryonNmb() + daughter2()->baryonNmb());
	if (_debug)
		cout << "after = " << parent()->baryonNmb() << endl;
	return parent()->baryonNmb();
}


bool
isobarDecayVertex::checkMultiplicativeQn(const int     mQn,
                                         const int     d1Qn,
                                         const int     d2Qn,
                                         const string& qnName)
{
	if (mQn == d1Qn * d2Qn) {
		if (_debug)
			printDebug << "success: " << *this << ": "
			           << ((qnName == "") ? "multiplicative quantum number" : qnName)
			           << " is consistent" << endl;
		return true;
	} else {
		printWarn << ((qnName == "") ? "multiplicative quantum number" : qnName) << " mismatch: "
		          << parent()->name() << " = " << mQn << " != " << d1Qn * d2Qn  << " = "
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
			printDebug << "success: " << *this << ": "
			           << ((qnName == "") ? "additive quantum number" : qnName)
			           << " is consistent" << endl;
		return true;
	} else {
		printWarn << ((qnName == "") ? "additive quantum number" : qnName) << " mismatch: "
		          << parent()->name() << " = " << mQn << " != " << d1Qn + d2Qn  << " = "
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
	if (not checkMultiplicativeQn(parent()->G(), daughter1()->G(), daughter2()->G(), "G-parity"))
		vertexConsistent = false;
	// parity
	const int angMomParity = (_L % 4 == 0) ? 1 : -1;  // modulo 4 because L is in units of hbar / 2
	if (parent()->P() != daughter1()->P() * daughter2()->P() * angMomParity) {
		printWarn << "parity mismatch: "
		          << parent()->name() << " = " << parent()->P() << " != "
		          << daughter1()->P() * daughter2()->P() * angMomParity  << " = "
		          << "(" << daughter1()->name() << " = " << daughter1()->P() << ") * "
		          << "(" << daughter2()->name() << " = " << daughter2()->P() << ") * "
		          << "(ang. momentum = " << angMomParity << ")" << endl;
		vertexConsistent = false;
	} else if (_debug)
		printDebug << "success: " << *this << ": parity is consistent" << endl;
	// check additive quantum numbers
	// charge
	if (not checkAdditiveQn(parent()->charge(), daughter1()->charge(),
	                        daughter2()->charge(), "charge"))
		vertexConsistent = false;
	// baryon number
	if (not checkAdditiveQn(parent()->baryonNmb(), daughter1()->baryonNmb(),
	                        daughter2()->baryonNmb(), "baryonNmb"))
		vertexConsistent = false;
	// strangeness
	if (not checkAdditiveQn(parent()->strangeness(), daughter1()->strangeness(),
	                        daughter2()->strangeness(), "strangeness"))
		vertexConsistent = false;
	// charm
	if (not checkAdditiveQn(parent()->charm(), daughter1()->charm(),
	                        daughter2()->charm(), "charm"))
		vertexConsistent = false;
	// beauty
	if (not checkAdditiveQn(parent()->beauty(), daughter1()->beauty(),
	                        daughter2()->beauty(), "beauty"))
		vertexConsistent = false;
	// check orbital angular momentum
	if (not isEven(_L)) {
		printWarn << "relative orbital angular momentum L = " << spinQn(_L)
		          << " is not an integer" << endl;
		vertexConsistent = false;
	} else if (_debug)
		printDebug << "success: " << *this
		           << ": relative orbital angular momentum L = " << spinQn(_L) << " is integer" << endl;
	// spin coupling: S in {|s1 - s2|, ..., s1 + s2}
	if (not spinStatesCanCouple(daughter1()->J(), daughter2()->J(), _S)) {
		printWarn << "spins "
		          << "(" << daughter1()->name() << " J = " << spinQn(daughter1()->J()) << ") and "
		          << "(" << daughter2()->name() << " J = " << spinQn(daughter2()->J()) << ") "
		          << "cannot couple to total spin S = " << spinQn(_S) << endl;
		vertexConsistent = false;
	} else if (_debug)
		printDebug << "success: " << *this << ": spin-spin coupling is consistent" << endl;
	// L-S coupling: J in {|L - S|, ..., L + S}
	if (not spinStatesCanCouple(_L, _S, parent()->J())) {
		printWarn << "orbital angular momentum L = " << spinQn(_L) << " and spin S = " << spinQn(_S)
		          << " cannot couple to angular momentum J = " << spinQn(parent()->J()) << endl;
		vertexConsistent = false;
	} else if (_debug)
		printDebug << "success: " << *this << ": L-S coupling is consistent" << endl;
	// check, whether isospins can couple
	if (not spinStatesCanCouple(daughter1()->isospin(), daughter1()->isospinProj(),
	                            daughter2()->isospin(), daughter2()->isospinProj(),
	                            parent   ()->isospin(), parent   ()->isospinProj())) {
		printWarn << "isospin Clebsch-Gordan "
		          << "("   << spinQn(daughter1()->isospin()) << ", " << spinQn(daughter1()->isospinProj())
		          << "; "  << spinQn(daughter2()->isospin()) << ", " << spinQn(daughter2()->isospinProj())
		          << " | " << spinQn(parent   ()->isospin()) << ", " << spinQn(parent   ()->isospinProj())
		          << ") == 0" << endl;
		vertexConsistent = false;
	} else if (_debug)
		printDebug << *this << ": isospin Clebsch-Gordan is not zero" << endl;
	// check consistency of C-, G-parity, and isospin
	if (   (    (daughter1()->charge() == 0)
	        and (daughter1()->C() != daughter1()->G() * (daughter1()->isospin() % 4 == 0 ? 1 : -1)))
	    or (    (daughter1()->charge() != 0)
	        and (daughter1()->C() != 0))
	    or (    (daughter2()->charge() == 0)
          and (daughter2()->C() != daughter2()->G() * (daughter2()->isospin() % 4 == 0 ? 1 : -1)))
	    or (    (daughter2()->charge() != 0)
	        and (daughter2()->C() != 0))
	    or (    (parent   ()->charge() == 0)
          and (parent   ()->C() != parent   ()->G() * (parent   ()->isospin() % 4 == 0 ? 1 : -1)))
	    or (    (parent   ()->charge() != 0)
	        and (parent   ()->C() != 0))) {
		printWarn << "C-, G-parity, and isospin are inconsistent "
		          << daughter1()->qnSummary() << ", " << daughter2()->qnSummary() << ", "
		          << parent   ()->qnSummary() << endl;
		vertexConsistent = false;
	} else if (_debug)
		printDebug << *this << ": C-, G-parity, and isospin are consistent" << endl;

	//!!! missing: spin projections
	if (vertexConsistent) {
		if (_debug)
			printDebug << "data for " << *this  << " are consistent: " << endl;
	} else
		printWarn << "data for " << *this << " are inconsistent (see warnings above)" << endl;
	return vertexConsistent;
}


ostream&
isobarDecayVertex::print(ostream& out) const
{
	out << name() << ": "
	    << parent()->qnSummary() << "  --{" << *_massDep << "}-->  "
	    << daughter1()->qnSummary()
	    << "  [L = " << spinQn(_L) << ", S = " << spinQn(_S) << "]  "
	    << daughter2()->qnSummary();
	return out;
}


ostream&
isobarDecayVertex::dump(ostream& out) const
{
	out << name() << ":" << endl
	    << "    parent: "     << *parent()    << endl
	    << "    daughter 1: " << *daughter1() << endl
	    << "    daughter 2: " << *daughter2() << endl
	    << "    L = " << spinQn(_L) << ", S = " << spinQn(_S) << endl
	    << "    " << *_massDep << endl;
	return out;
}


ostream&
isobarDecayVertex::printPointers(ostream& out) const
{
	out << name() << " "  << this << ": "
	    << "parent particle: "     << parent()    << "; "
	    << "daughter 1 particle: " << daughter1() << "; "
	    << "daughter 2 particle: " << daughter2() << endl;
	return out;
}
