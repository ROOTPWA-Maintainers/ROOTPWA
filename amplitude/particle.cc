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


#include "utilities.h"
#include "particle.h"

	
using namespace std;
using namespace rpwa;


bool particle::_debug = false;


particle::particle()
  : particleProperties(),
    _charge           (0),
    _spinProj         (0),
    _lzVec            ()
{ }


particle::particle(const particle& part)
{
  *this = part;
}


particle::particle(const particleProperties& partProp,
		   const int                 charge,
		   const int                 spinProj,
		   const TVector3&           momentum)
  : particleProperties(partProp),
    _charge           (charge),
    _spinProj         (spinProj),
    _lzVec            (TLorentzVector(momentum, sqrt(momentum.Mag2() + mass() * mass())))
{ }

	
particle::particle(const string&   partName,
		   const int       charge,
		   const int       spinProj,
		   const TVector3& momentum)
  : _charge  (charge),
    _spinProj(spinProj)
{
  fillFromDataTable(partName);
  _lzVec = TLorentzVector(momentum, sqrt(momentum.Mag2() + mass() * mass()));
}

	
particle::~particle()
{ }


particle&
particle::operator = (const particle& part)
{
  if (this != &part) {
    particleProperties::operator = (part);
    _charge   = part._charge;
    _spinProj = part._spinProj;
    _lzVec    = part._lzVec;
  }
  return *this;
}


ostream&
particle::print(ostream& out) const
{
  particleProperties::print(out);
  out << ", "
      << "charge = "          << sign(_charge) << ", "
      << "spin projection = " << _spinProj     << ", "
      << "Lorentz-vector = "  << _lzVec;
  return out;
}
