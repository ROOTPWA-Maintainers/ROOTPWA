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
//      container class for particle properties
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "utilities.h"
#include "particleDataTable.h"
#include "particleProperties.h"

	
using namespace std;
using namespace rpwa;


bool particleProperties::_debug = false;

	
particleProperties::particleProperties()
  : _name       (""),
    _mass       (0),
    _width      (0),
    _baryonNmb  (0),
    _isospin    (0),
    _strangeness(0),
    _charm      (0),
    _beauty     (0),
    _G          (0),
    _J          (0),
    _P          (0),
    _C          (0)
{ }


particleProperties::particleProperties(const particleProperties& partProp)
{
  *this = partProp;
}


particleProperties::particleProperties(const std::string& partName,
				       const int          isospin,
				       const int          G,
				       const int          J,
				       const int          P,
				       const int          C)
  : _name       (partName),
    _mass       (0),
    _width      (0),
    _baryonNmb  (0),
    _isospin    (isospin),
    _strangeness(0),
    _charm      (0),
    _beauty     (0),
    _G          (G),
    _J          (J),
    _P          (P),
    _C          (C)
{
}


particleProperties::~particleProperties()
{ }


particleProperties&
particleProperties::operator = (const particleProperties& partProp)
{
  if (this != &partProp) {
    _name        = partProp._name;
    _mass        = partProp._mass;
    _width       = partProp._width;
    _baryonNmb   = partProp._baryonNmb;
    _isospin     = partProp._isospin;
    _strangeness = partProp._strangeness;
    _charm       = partProp._charm;
    _beauty      = partProp._beauty;
    _G           = partProp._G;
    _J           = partProp._J;
    _P           = partProp._P;
    _C           = partProp._C;
  }
  return *this;
}


bool
operator == (const particleProperties& lhsProp,
	     const particleProperties& rhsProp)
{
  return (   (lhsProp.name()        == rhsProp.name()       )
          && (lhsProp.mass()        == rhsProp.mass()       )
          && (lhsProp.width()       == rhsProp.width()      )
          && (lhsProp.baryonNmb()   == rhsProp.baryonNmb()  )
          && (lhsProp.isospin()     == rhsProp.isospin()    )
          && (lhsProp.strangeness() == rhsProp.strangeness())
          && (lhsProp.charm()       == rhsProp.charm()      )
          && (lhsProp.beauty()      == rhsProp.beauty()     )
          && (lhsProp.G()           == rhsProp.G()          )
          && (lhsProp.J()           == rhsProp.J()          )
          && (lhsProp.P()           == rhsProp.P()          )
	  && (lhsProp.C()           == rhsProp.C()          ));
}


bool
particleProperties::fillFromDataTable(const string& name)
{
  const particleProperties* partProp = particleDataTable::instance().entry(name);
  if (!partProp) {
    printWarn << "trying to fill particle properties for '" << name << "' from non-existing table entry" << endl;
    return false;
  } else {
    *this = *partProp;
    if (_debug)
      printInfo << "succesfully filled particle properties for '" << name << "': "
		<< *this << endl;
    return true;
  }
}


ostream&
particleProperties::print(ostream& out) const
{
  out << "particle '"  << _name         << "': "
      << "mass = "     << _mass         << " GeV/c^2, "
      << "width = "    << _width        << " GeV/c^2, "
      << "baryon # = " << _baryonNmb    << ", "
      << "(2I)^G(2J)^PC = ("  << _isospin << ")^" << sign(_G)
      << "(" << _J << ")^" << sign(_P) << sign(_C) << ", "
      << "strangeness = " << _strangeness << ", "
      << "charm = "       << _charm       << ", "
      << "beauty = "      << _beauty;
  return out;
}


ostream&
particleProperties::dump(ostream& out) const
{
  out << _name        << "\t"
      << _mass        << "\t"
      << _width       << "\t"
      << _baryonNmb   << "\t"
      << _isospin     << "\t"
      << _strangeness << "\t" 
      << _charm       << "\t" 
      << _beauty      << "\t" 
      << _G           << "\t" 
      << _J           << "\t" 
      << _P           << "\t" 
      << _C;
  return out;
}


bool
particleProperties::read(istringstream& line)
{
  if (_debug)
    printInfo << "trying to read particle properties from line '" << line.str() << "' ... " << flush;
  if (line >> _name
           >> _mass
           >> _width
           >> _baryonNmb
           >> _isospin
           >> _strangeness
           >> _charm
           >> _beauty
           >> _G
           >> _J
           >> _P
           >> _C) {
    if (_debug) {
      cout << "success" << endl;
      printInfo << "read "<< *this << endl;
    }
    return true;
  } else {
    printWarn << "problems reading particle data from line '" << line.str() << "'" << endl;
    return false;
  }
}
