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
//      container class for all external particle related information
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "reportingUtils.hpp"
#include "reportingUtilsRoot.hpp"
#include "conversionUtils.hpp"
#include "particle.h"

	
using namespace std;
using namespace rpwa;


bool particle::_debug = false;


particle::particle()
	: particleProperties(),
	  _spinProj         (0),
	  _lzVec            (),
	  _index            (-1),
	  _refl             (0)
{ }


particle::particle(const particle& part)
{
	*this = part;
}


particle::particle(const particleProperties& partProp,
                   const int                 index,
                   const int                 spinProj,
                   const int                 refl,
                   const TVector3&           momentum)
	: particleProperties(partProp),
	  _spinProj         (spinProj),
	  _lzVec            (TLorentzVector(momentum, sqrt(momentum.Mag2() + mass() * mass()))),
	  _index            (index),
	  _refl             (refl)
{ }

	
particle::particle(const string&   partName,
                   const bool      requirePartInTable,
                   const int       index,
                   const int       spinProj,
                   const int       refl,
                   const TVector3& momentum)
	: particleProperties(),
	  _spinProj(spinProj),
	  _index   (index),
	  _refl    (refl)
{
	if (not fillFromDataTable(partName, requirePartInTable))
		// set at least name
		setName(partName);
	_lzVec = TLorentzVector(momentum, sqrt(momentum.Mag2() + mass() * mass()));
}

	
particle::particle(const string& partName,
                   const int     isospin,
                   const int     G,
                   const int     J,
                   const int     P,
                   const int     C,
                   const int     spinProj,
                   const int     refl,
                   const int     index)
	: particleProperties(partName, isospin, G, J, P, C),
	  _spinProj(spinProj),
	  _index   (index),
	  _refl    (refl)
{ }


particle::~particle()
{ }


particle&
particle::operator =(const particle& part)
{
	if (this != &part) {
		particleProperties::operator =(part);
		_spinProj = part._spinProj;
		_lzVec    = part._lzVec;
		_index    = part._index;
		_refl     = part._refl;
	}
	return *this;
}


particle*
particle::doClone() const
{
	particle* particleClone(new particle(*this));
	if (_debug)
		printDebug << "cloned " << *this << "; " << this << " -> " << particleClone << endl;
	return particleClone;
}


bool
particle::isEqualTo(const particleProperties& partProp) const
{
	// printErr << "!!!particle ==" << endl;
	const particle* part = dynamic_cast<const particle*>(&partProp);
	if (not part)
		return false;
	if (not particleProperties::isEqualTo(partProp))
		return false;
	return (    (spinProj    () == part->spinProj    ())
	        and (index       () == part->index       ())
	        and (reflectivity() == part->reflectivity()));
}


void
particle::setProperties(const particleProperties& prop)
{
	if (this != &prop)
		particleProperties::operator =(prop);
}


string
particle::qnSummary() const
{
	ostringstream out;
	out << name() << "[" << spinQn(isospin()) << parityQn(G())
	    << "(" << spinQn(J()) << parityQn(P()) << parityQn(C()) << ")"
	    << spinQn(spinProj()) << parityQn(reflectivity()) << "]";
	return out.str();
}


ostream&
particle::print(ostream& out) const
{
	particleProperties::print(out);
	out << ", "
	    << "spin proj. = "     << spinQn(_spinProj) << ", "
	    << "reflectivity = "   << _refl             << ", "
	    << "Lorentz-vector = " << _lzVec            << " GeV, "
	    << "index = "          << _index;
	return out;
}


string
particle::label() const
{
	ostringstream out;
	out << name() << "[" << spinQn(isospin()) << parityQn(G())
	    << "(" << spinQn(J()) << parityQn(P()) << parityQn(C()) << ")]";
	return out.str();
}
