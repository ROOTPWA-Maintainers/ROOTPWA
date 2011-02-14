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


#include <cstdlib>

#include "conversionUtils.hpp"
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
particleProperties::operator =(const particleProperties& partProp)
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
rpwa::operator ==(const particleProperties& lhsProp,
                  const particleProperties& rhsProp)
{
	return (    (lhsProp.name()        == rhsProp.name()       )
	        and (lhsProp.mass()        == rhsProp.mass()       )
	        and (lhsProp.width()       == rhsProp.width()      )
	        and (lhsProp.baryonNmb()   == rhsProp.baryonNmb()  )
	        and (lhsProp.isospin()     == rhsProp.isospin()    )
	        and (lhsProp.strangeness() == rhsProp.strangeness())
	        and (lhsProp.charm()       == rhsProp.charm()      )
	        and (lhsProp.beauty()      == rhsProp.beauty()     )
	        and (lhsProp.G()           == rhsProp.G()          )
	        and (lhsProp.J()           == rhsProp.J()          )
	        and (lhsProp.P()           == rhsProp.P()          )
	        and (lhsProp.C()           == rhsProp.C()          ));
}


// selector string can contain any of the following:
// I, G, J, P, C, strangeness, charm, beauty, baryonNmb, or allQn
bool 
rpwa::operator ==(const particleProperties&               lhsProp,
                  const pair<particleProperties, string>& rhs)
{
	const particleProperties& rhsProp    = rhs.first;
	const string&             selector   = rhs.second;
	const bool                checkAllQn = (selector.find("allQn") != string::npos);
	return (    (   ((selector.find("baryonNmb") == string::npos) and not checkAllQn)
	             or (lhsProp.baryonNmb()         == rhsProp.baryonNmb()))
	        and (   ((selector.find("I") == string::npos) and not checkAllQn)
	             or (lhsProp.isospin()   == rhsProp.isospin()))
	        and (   ((selector.find("strangeness") == string::npos) and not checkAllQn)
	             or (lhsProp.strangeness()         == rhsProp.strangeness()))
	        and (   ((selector.find("charm") == string::npos) and not checkAllQn)
	             or (lhsProp.charm()         == rhsProp.charm()))
	        and (   ((selector.find("beauty") == string::npos) and not checkAllQn)
	             or (lhsProp.beauty()         == rhsProp.beauty()))
	        and (   ((selector.find("G") == string::npos) and not checkAllQn)
	             or (lhsProp.G()         == rhsProp.G()))
	        and (   ((selector.find("J") == string::npos) and not checkAllQn)
	             or (lhsProp.J()         == rhsProp.J()))
	        and (   ((selector.find("P") == string::npos) and not checkAllQn)
	             or (lhsProp.P()         == rhsProp.P()))
	        and (   ((selector.find("C") == string::npos) and not checkAllQn)
	             or (lhsProp.C()         == rhsProp.C())));
}


bool
particleProperties::isXParticle() const
{
	return (_name == "X") or (_name == "X-") or (_name == "X0") or (_name == "X+");
}


bool
particleProperties::fillFromDataTable(const string& partName)
{
	string       name         = partName;
	const string strippedName = stripChargeFromName(partName);
	if (not particleDataTable::instance().isInTable(partName))
		// try with charge stripped from name
		name = strippedName;
	const particleProperties* partProp = particleDataTable::instance().entry(name);
	if (not partProp) {
		printWarn << "trying to fill particle properties for '" << partName << "' from non-existing table entry" << endl;
		return false;
	} else {
		*this = *partProp;
		if (_debug)
			printInfo << "succesfully filled particle properties for '" << partName << "': "
			          << *this << endl;
		return true;
	}
}


void
particleProperties::setSCB(const int strangeness,
                           const int charm,
                           const int beauty)
{
	setStrangeness(strangeness);
	setCharm      (charm);
	setBeauty     (beauty);
}


void
particleProperties::setIGJPC(const int isospin,
                             const int G,
                             const int J,
                             const int P,
                             const int C)
{
	setIsospin(isospin);
	setG      (G);
	setJ      (J);
	setP      (P);
	setC      (C);
}


string
particleProperties::qnSummary() const
{
	ostringstream out;
	out << name() << "[" << 0.5 * isospin() << ((G() != 0) ? sign(G()) : "")
	    << "(" << 0.5 * J() << ((P() != 0) ? sign(P()) : "") << ((C() != 0) ? sign(C()) : "") << ")]";
	return out.str();
}


ostream&
particleProperties::print(ostream& out) const
{
	out << "particle '"  << name()         << "': "
	    << "mass = "     << mass()         << " GeV/c^2, "
	    << "width = "    << width()        << " GeV/c^2, "
	    << "baryon # = " << baryonNmb()    << ", "
	    << "I" << ((G() != 0) ? "^G" : "") << " J";
	if (P() != 0)
		out << "^P" << ((C() != 0) ? "C" : "") << " = ";
  else
		out << ((C() != 0) ? "^C" : "") << " = ";
	out << 0.5 * isospin();
	if (G() != 0)
		out << "^" << sign(G());
	out << " " << 0.5 * J();
	if (P() != 0)
		out << "^" << sign(P()) << ((C() != 0) ? sign(C()) : "");
	else
		if (C() != 0)
			out << "^" << sign(C());
	out << ", "
	    << "strangeness = " << strangeness() << ", "
	    << "charm = "       << charm()       << ", "
	    << "beauty = "      << beauty();
	return out;
}


ostream&
particleProperties::dump(ostream& out) const
{
	out << name()        << "\t"
	    << mass()        << "\t"
	    << width()       << "\t"
	    << baryonNmb()   << "\t"
	    << isospin()     << "\t"
	    << strangeness() << "\t" 
	    << charm()       << "\t" 
	    << beauty()      << "\t" 
	    << G()           << "\t" 
	    << J()           << "\t" 
	    << P()           << "\t" 
	    << C();
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


string
particleProperties::chargeFromName(const string& name,
                                   int&          charge)
{
	// assumes that for singly charged particles the charge is
	// designated by the last character in the name (e.g. 'pi+')
	// if no charge character is found charge 0 is assumed
	// for multiply charged particles the charge state is the
	// second-to-last character (e.g. 'Delta(1232)2+')
	charge = 0;
	string     strippedName = name;
	const char last         = name[name.length() - 1];
	if ((last == '+') or (last == '0') or (last == '-')) {
		charge = sign(last);
		strippedName.erase(strippedName.length() - 1);
		const char secondLast[] = {name[name.length() - 2], '\0'};
		const int  chargeVal    = atoi(secondLast);
		if (chargeVal != 0) {
			charge *= chargeVal;
			strippedName.erase(strippedName.length() - 1);
		}
	}
	return strippedName;
}


string
particleProperties::stripChargeFromName(const string& name)
{
	int dummy;
	return chargeFromName(name, dummy);
}
