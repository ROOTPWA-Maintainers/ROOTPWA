///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2016 Sebastian Uhl (TUM)
//
//    This file is part of ROOTPWA
//
//    ROOTPWA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ROOTPWA is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ROOTPWA.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      implementation information for a single parameter of the resonance
//      fit
//
//-------------------------------------------------------------------------


#include "parameter.h"


rpwa::resonanceFit::parameter::parameter()
	: _name("INVALID"),
	  _startValue(0.0),
	  _startError(0.0),
	  _fixed(false),
	  _limitLower(0.0),
	  _limitedLower(false),
	  _limitUpper(0.0),
	  _limitedUpper(false),
	  _step(0.0)
{
}


std::ostream&
rpwa::resonanceFit::parameter::print(std::ostream& out, const bool newLine) const
{
	out << "parameter '" << _name << "' "
	    << "start value: " << _startValue << " +/- " << _startError << " ";
	if(_limitedLower and _limitedUpper) {
		out << "limits: " << _limitLower << "-" << _limitUpper;
	} else if(_limitedLower) {
		out << "lower limit: " << _limitLower;
	} else if(_limitedUpper) {
		out << "upper limit: " << _limitUpper;
	} else {
		out << "unlimited";
	}
	out << " " << (_fixed ? "(FIXED) " : "") << "step size: " << _step;

	if(newLine) {
		out << std::endl;
	}

	return out;
}
