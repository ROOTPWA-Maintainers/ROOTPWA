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
//      production vertex virtual base class
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "reportingUtils.hpp"
#include "productionVertex.h"


using namespace std;
using namespace rpwa;


bool productionVertex::_debug = false;


productionVertex::productionVertex()
	: interactionVertex()
{
	if (_debug)
		printDebug << "constructed " << *this << endl;
}


productionVertex::~productionVertex()
{ }


// default implementation
std::vector<std::complex<double> >
productionVertex::productionAmp() const
{
	int numEvents = XParticle()->lzVecs().size();
	return std::vector<std::complex<double> >(numEvents, 1);
}
