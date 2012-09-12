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
#include <iterator>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "spinUtils.hpp"
#include "conversionUtils.hpp"
#include "particleDataTable.h"
#include "particleProperties.h"

	
using namespace std;
using namespace rpwa;
using namespace boost;


bool particleProperties::_debug = false;

	
particleProperties::particleProperties()
	: _name        (""),
	  _antiPartName(""),
	  _charge      (0),
	  _mass        (0),
	  _mass2       (0),
	  _width       (0),
	  _baryonNmb   (0),
	  _isospin     (0),
	  _strangeness (0),
	  _charm       (0),
	  _beauty      (0),
	  _G           (0),
	  _J           (0),
	  _P           (0),
	  _C           (0)
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
	: _antiPartName(""),
	  _mass        (0),
	  _mass2       (0),
	  _width       (0),
	  _baryonNmb   (0),
	  _strangeness (0),
	  _charm       (0),
	  _beauty      (0)
{
	setName   (partName);
	setIsospin(isospin);
	setG      (G);
	setJ      (J);
	setP      (P);
	setC      (C);
}


particleProperties::~particleProperties()
{ }


particleProperties&
particleProperties::operator =(const particleProperties& partProp)
{
	if (this != &partProp) {
		_name         = partProp._name;
		_antiPartName = partProp._antiPartName;
		_charge       = partProp._charge;
		_mass         = partProp._mass;
		_mass2        = partProp._mass2;
		_width        = partProp._width;
		_baryonNmb    = partProp._baryonNmb;
		_isospin      = partProp._isospin;
		_strangeness  = partProp._strangeness;
		_charm        = partProp._charm;
		_beauty       = partProp._beauty;
		_G            = partProp._G;
		_J            = partProp._J;
		_P            = partProp._P;
		_C            = partProp._C;
	}
	return *this;
}


bool
rpwa::operator ==(const particleProperties& lhsProp,
                  const particleProperties& rhsProp)
{
	return (    (lhsProp.name()         == rhsProp.name()        )
	        and (lhsProp.antiPartName() == rhsProp.antiPartName())
	        and (lhsProp.charge()       == rhsProp.charge()      )
	        and (lhsProp.mass()         == rhsProp.mass()        )
	        and (lhsProp.width()        == rhsProp.width()       )
	        and (lhsProp.baryonNmb()    == rhsProp.baryonNmb()   )
	        and (lhsProp.isospin()      == rhsProp.isospin()     )
	        and (lhsProp.strangeness()  == rhsProp.strangeness() )
	        and (lhsProp.charm()        == rhsProp.charm()       )
	        and (lhsProp.beauty()       == rhsProp.beauty()      )
	        and (lhsProp.G()            == rhsProp.G()           )
	        and (lhsProp.J()            == rhsProp.J()           )
	        and (lhsProp.P()            == rhsProp.P()           )
	        and (lhsProp.C()            == rhsProp.C()           ));
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
	return (    (   ((selector.find("charge") == string::npos) and not checkAllQn)
	             or (lhsProp.charge()         == rhsProp.charge()))
	        and (   ((selector.find("baryonNmb") == string::npos) and not checkAllQn)
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
particleProperties::isSpinExotic() const
{
	return (isMeson() and igjpIsExotic(isospin(), G(), J(), P()));
}


bool
particleProperties::fillFromDataTable(const string& partName,
                                      const bool    warnIfNotExistent)
{
	string name = partName;
	if (not particleDataTable::instance().isInTable(partName)
	    and (partName == stripChargeFromName(partName)))
		// if no charge is given in the name assume charge 0
		name += '0';
	const particleProperties* partProp = particleDataTable::instance().entry(name, warnIfNotExistent);
	if (not partProp) {
		if (warnIfNotExistent)
			printWarn << "trying to fill particle properties for '"
			          << partName << "' from non-existing table entry" << endl;
		return false;
	} else {
		*this = *partProp;
		if (_debug)
			printDebug << "succesfully filled particle properties for '" << partName << "': "
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


particleProperties
particleProperties::antiPartProperties() const
{
	particleProperties antiPartProp;
	antiPartProp.setName        (_antiPartName);
	antiPartProp.setAntiPartName(_name        );
	antiPartProp.setCharge      (-_charge     );
	antiPartProp.setMass        (_mass        );
	antiPartProp.setWidth       (_width       );
	antiPartProp.setBaryonNmb   (-_baryonNmb  );
	antiPartProp.setIsospin     (_isospin     );
	antiPartProp.setStrangeness (-_strangeness);
	antiPartProp.setCharm       (-_charm      );
	antiPartProp.setBeauty      (-_beauty     );
	antiPartProp.setG           (_G           );
	antiPartProp.setJ           (_J           );
	antiPartProp.setP           (_P           );
	antiPartProp.setC           (_C           );
	return antiPartProp;
}


string
particleProperties::bareNameLaTeX() const
{
	// split name into mass and symbol parts
	// const string name = bareName();
	// typename boost::iterator_range<std::string> startPos;
	// find_first(name, "(");
	// iterator_range<string> endPos   = find_last (name, ")");
	// iterator_range<string> range    = make_iterator_range(startPos, endPos);
	// const string mass = copy_range(range);
	// printDebug << "!!!HERE " << name << ": " << mass << "; " << range << endl;
	// // handle antiparticle
	// if ()
	// // handle * particles
	
	// // setup particle-name dictionary
	// map<string, string> partNameDict;
	// isobars["gamma"  ] = "\\gamma";
	// isobars["mu"     ] = "\\mu";
	// isobars["pi"     ] = "\\pi";
	// isobars["eta"    ] = "\\eta";
	// isobars["sigma"  ] = "\\sigma";
	// isobars["rho"    ] = "\\rho";
	// isobars["omega"  ] = "\\omega";
	// isobars["phi"    ] = "\\phi";
	// isobars["kappa"  ] = "\\kappa";
	// isobars["nucleon"] = "N";
	// isobars["Delta"  ] = "\\Delta";
	// isobars["Lambda" ] = "\\Lambda";


	// vector<iterator_range<string::iterator> > foundPos;
	// find_all(foundPos, name, "(");
	// if (foundPos.size() > 1) {
	// 	printErr << "particle name '" << name << "' contains more than one '('. "
	// 	         << "cannot construct LaTeX name" << endl;
	// 	return "";
	// }
	return "";
}



string
particleProperties::qnSummary() const
{
	ostringstream out;
	out << name() << "[" << spinQn(isospin()) << parityQn(G())
	    << "(" << spinQn(J()) << parityQn(P()) << parityQn(C()) << ")]";
	return out.str();
}


ostream&
particleProperties::print(ostream& out) const
{
	out << "particle '"  << name()      << "': "
	    << "mass = "     << mass()      << " GeV/c^2, "
	    << "width = "    << width()     << " GeV/c^2, "
	    << "baryon # = " << baryonNmb() << ", "
	    << "I" << ((G() != 0) ? "^G" : "") << " J";
	if (P() != 0)
		out << "^P" << ((C() != 0) ? "C" : "") << " = ";
	else
		out << ((C() != 0) ? "^C" : "") << " = ";
	out << spinQn(isospin());
	if (G() != 0)
		out << "^" << sign(G());
	out << " " << spinQn(J());
	if (P() != 0)
		out << "^" << sign(P()) << parityQn(C());
	else
		if (C() != 0)
			out << "^" << sign(C());
	out << ", "
	    << "strangeness = "             << strangeness()     << ", "
	    << "charm = "                   << charm()           << ", "
	    << "beauty = "                  << beauty()          << ", "
	    << "is meson = "                << yesNo(isMeson())  << ", ";
	if (isMeson())
		out << "is spin-exotic = " << yesNo(isSpinExotic()) << ", ";
	out << "is baryon = "               << yesNo(isBaryon()) << ", "
	    << "is lepton = "               << yesNo(isLepton()) << ", "
	    << "is photon = "               << yesNo(isPhoton()) << ", "
	    << "antiparticle '"             << antiPartName()    << "', "
	    << "is its own antiparticle = " << yesNo(isItsOwnAntiPart());

	// decay products
	unsigned int ndec=nDecays();
	if(ndec > 0){
		out << "\n Known decay modes: " << endl;
		for(unsigned int idec=0;idec<ndec;++idec){
			std::copy(_decaymodes[idec].begin(), _decaymodes[idec].end(), std::ostream_iterator<string>(std::cout, " "));
			out << endl;
		}
	}

	return out;
}


ostream&
particleProperties::dump(ostream& out) const
{
	out << name        () << "\t"
	    << antiPartName() << "\t"
	    << mass        () << "\t"
	    << mass2       () << "\t"
	    << width       () << "\t"
	    << baryonNmb   () << "\t"
	    << isospin     () << "\t"
	    << strangeness () << "\t" 
	    << charm       () << "\t" 
	    << beauty      () << "\t" 
	    << G           () << "\t" 
	    << J           () << "\t" 
	    << P           () << "\t" 
	    << C           ();
	return out;
}


bool
particleProperties::read(istringstream& line)
{
	if (_debug)
		printDebug << "trying to read particle properties from line '" << line.str() << "' ... " << flush;
	std::string name = "", antiPartName = "";
	double      mass = 0, width = 0;
	int         baryonNmb = 0, isospin = 0, strangeness = 0, charm = 0, beauty = 0;
	int         G = 0, J = 0, P = 0, C = 0;
	if (line >> name
	         >> antiPartName
	         >> mass
	         >> width
	         >> baryonNmb
	         >> isospin
	         >> strangeness
	         >> charm
	         >> beauty
	         >> G
	         >> J
	         >> P
	         >> C) {
		setName        (name        );
		setAntiPartName(antiPartName);
		setMass        (mass        );
		setWidth       (width       );
		setBaryonNmb   (baryonNmb   );
		setIsospin     (isospin     );
		setStrangeness (strangeness );
		setCharm       (charm       );
		setBeauty      (beauty      );
		setG           (G           );
		setJ           (J           );
		setP           (P           );
		setC           (C           );
		if (_debug) {
			cout << "success" << endl;
			printDebug << "read "<< *this << endl;
		}
		return true;
	} else {
		printWarn << "problems reading particle data from line '" << line.str() << "'" << endl;
		return false;
	}
}


string
particleProperties::nameWithCharge(const string& bareName,
                                   const int     charge)
{
	string name = bareName;
	if (bareName != "") {
		if (abs(charge) > 1) {
			string q = "";
			try {
				q = lexical_cast<char>(abs(charge));
			} catch (bad_lexical_cast&) { }
			name += q;
		}
		name += sign(charge);
	}
	return name;
}


string
particleProperties::chargeFromName(const string& partName,
                                   int&          charge)
{
	// naming scheme: <bare name>XY
	//
	// for singly charged particles the charge is designated by the last
	// character Y in the particle name (e.g. 'pi0', 'rho-', 'p+'). If Y
	// is neither of '-', '0', or '+', the particle charge is set to 0.
	//
	// for multiply charged particles the absolute value of the charge
	// is defined by the second-to-last character X, the sign by the
	// lats character Y (e.g. 'Delta(1232)2+').
	//
	// !NOTE! <bare name> should not end on a digit, as this would be
	// interpreted as the absolute value of the charge, so instead of
	// e.g. 'pi2+' one should use 'pi2(1670)+'
	charge = 0;
	string     strippedName = partName;
	const char last         = partName[partName.length() - 1];
	if ((last == '+') or (last == '0') or (last == '-')) {
		charge = sign(last);
		strippedName.erase(strippedName.length() - 1);
		int chargeVal = 0;
		try {
			chargeVal = lexical_cast<int>(partName[partName.length() - 2]);
		} catch (bad_lexical_cast&) { }
		if (chargeVal != 0) {
			charge *= chargeVal;
			strippedName.erase(strippedName.length() - 1);
		}
	}
	return strippedName;
}


string
particleProperties::stripChargeFromName(const string& partName)
{
	int dummy;
	return chargeFromName(partName, dummy);
}


bool 
particleProperties::hasDecay(const set<string>& daughters) const {
  return (find(_decaymodes.begin(),_decaymodes.end(),daughters)!=_decaymodes.end());
}
