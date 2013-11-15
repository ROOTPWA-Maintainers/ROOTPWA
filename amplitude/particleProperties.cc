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


particleProperties::decayMode::decayMode(const multiset<string>& daughters,
                                         const int               L,
                                         const int               S)
	: _daughters(daughters),
	  _L        (L),
	  _S        (S)
{ }


particleProperties::decayMode::~decayMode()
{ }


bool
particleProperties::decayMode::operator ==(const decayMode& rhsDecay) const
{
	if (_daughters != rhsDecay._daughters)
		return false;
	// compare L and S only if they are defined in left- and right-hand side
	if ((_L != -1) and (rhsDecay._L != -1) and (_L != rhsDecay._L))
		return false;
	if ((_S != -1) and (rhsDecay._S != -1) and (_S != rhsDecay._S))
		return false;
	return true;
}


ostream&
particleProperties::decayMode::print(ostream& out) const
{
	copy(_daughters.begin(), _daughters.end(), ostream_iterator<string>(out, "  "));
	if (_L != -1)
		out << "[L = " << spinQn(_L);
	if (_S != -1) {
		if (_L != -1)
			out << ", ";
		else
			out << "[";
		out << "S = " << spinQn(_S) << "]";
	} else if (_L != -1)
		out << "]";
	return out;
}


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
	  _C           (0),
	  _decayModes  ()
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
	  _beauty      (0),
	  _decayModes  ()
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
		_decayModes   = partProp._decayModes;
	}
	return *this;
}


namespace rpwa {

	// selector string can contain any of the following:
	// I, G, J, P, C, strangeness, charm, beauty, baryonNmb, or allQn
	bool
	operator ==(const particleProperties&               lhsProp,
	            const pair<particleProperties, string>& rhsPropSel)
	{
		const particleProperties& rhsProp    = rhsPropSel.first;
		const string&             selector   = rhsPropSel.second;
		const bool                checkAllQn = (selector.find("allQn") != string::npos);
		return (    (   ((selector.find("charge")      == string::npos) and not checkAllQn)
		             or (lhsProp.charge()              == rhsProp.charge()))
		        and (   ((selector.find("baryonNmb")   == string::npos) and not checkAllQn)
		             or (lhsProp.baryonNmb()           == rhsProp.baryonNmb()))
		        and (   ((selector.find("I")           == string::npos) and not checkAllQn)
		             or (lhsProp.isospin()             == rhsProp.isospin()))
		        and (   ((selector.find("strangeness") == string::npos) and not checkAllQn)
		             or (lhsProp.strangeness()         == rhsProp.strangeness()))
		        and (   ((selector.find("charm")       == string::npos) and not checkAllQn)
		             or (lhsProp.charm()               == rhsProp.charm()))
		        and (   ((selector.find("beauty")      == string::npos) and not checkAllQn)
		             or (lhsProp.beauty()              == rhsProp.beauty()))
		        and (   ((selector.find("G")           == string::npos) and not checkAllQn)
		             or (lhsProp.G()                   == rhsProp.G()))
		        and (   ((selector.find("J")           == string::npos) and not checkAllQn)
		             or (lhsProp.J()                   == rhsProp.J()))
		        and (   ((selector.find("P")           == string::npos) and not checkAllQn)
		             or (lhsProp.P()                   == rhsProp.P()))
		        and (   ((selector.find("C")           == string::npos) and not checkAllQn)
		             or (lhsProp.C()                   == rhsProp.C())));
	}

}


unsigned int
particleProperties::geantId() const
{
	return particleDataTable::geantIdFromParticleName(this->name());
}


bool
particleProperties::isSpinExotic() const
{
	return (isMeson() and igjpIsExotic(isospin(), G(), J(), P()));
}


bool
particleProperties::hasDecay(const decayMode& decay) const
{
  return find(_decayModes.begin(), _decayModes.end(), decay) != _decayModes.end();
}

bool
particleProperties::isStable() const
{
	return _name == "pi";
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
particleProperties::antiPartProperties(const bool convertDecaysModes) const
{
	particleProperties antiPartProp;
	if (isItsOwnAntiPart()) {
		if (convertDecaysModes)
			return *this;
		else {
			antiPartProp = *this;
			antiPartProp.deleteDecayModes();
		}
	} else {
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
		// convert decay modes to antiparticles
		if (convertDecaysModes)
			for (size_t i = 0 ; i < nmbDecays(); ++i) {
				multiset<string> antiDaughters;
				decayMode decay = _decayModes[i];
				// !NOTE! for associative containers with value type == key type iterator is const_iterator
				for (multiset<string>::const_iterator it = decay._daughters.begin();
				     it != decay._daughters.end(); ++it) {
					particleProperties daughterProp;
					daughterProp.fillFromDataTable(*it);
					antiDaughters.insert(daughterProp.antiPartName());
				}
				decay._daughters = antiDaughters;
				antiPartProp.addDecayMode(decay);
			}
	}
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
	    << "strangeness = "             << strangeness()         << ", "
	    << "charm = "                   << charm()               << ", "
	    << "beauty = "                  << beauty()              << ", "
	    << "is meson = "                << yesNo(isMeson())      << ", ";
	if (isMeson())
		out << "is spin-exotic = "        << yesNo(isSpinExotic()) << ", ";
	out << "is baryon = "               << yesNo(isBaryon())     << ", "
	    << "is lepton = "               << yesNo(isLepton())     << ", "
	    << "is photon = "               << yesNo(isPhoton())     << ", "
	    << "antiparticle '"             << antiPartName()        << "', "
	    << "is its own antiparticle = " << yesNo(isItsOwnAntiPart());
	// decay products
	const unsigned int nmbDecays = this->nmbDecays();
	if (nmbDecays > 0) {
		out << endl << "    decay modes:" << endl;
		for (unsigned int i = 0; i < nmbDecays; ++i) {
			out << "        -> " << _decayModes[i];
			if (i < nmbDecays - 1)
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
particleProperties::isEqualTo(const particleProperties& rhsProp) const
{
	return (    (name        () == rhsProp.name        ())
	        and (antiPartName() == rhsProp.antiPartName())
	        and (charge      () == rhsProp.charge      ())
	        and (mass        () == rhsProp.mass        ())
	        and (width       () == rhsProp.width       ())
	        and (baryonNmb   () == rhsProp.baryonNmb   ())
	        and (isospin     () == rhsProp.isospin     ())
	        and (strangeness () == rhsProp.strangeness ())
	        and (charm       () == rhsProp.charm       ())
	        and (beauty      () == rhsProp.beauty      ())
	        and (G           () == rhsProp.G           ())
	        and (J           () == rhsProp.J           ())
	        and (P           () == rhsProp.P           ())
	        and (C           () == rhsProp.C           ()));
}
