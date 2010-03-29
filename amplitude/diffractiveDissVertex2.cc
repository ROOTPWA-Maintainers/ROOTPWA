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
//      class that describes production vertex in diffractive dissociation
//      beam-Reggeon-(X-system) vertex has exactly one incoming beam and
//      one outgoing X particle, which unambiguously defines the Reggeon
//      kinematics
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "utilities.h"
#include "diffractiveDissVertex2.h"

	
using namespace std;
using namespace rpwa;


bool diffractiveDissVertex2::_debug = false;


diffractiveDissVertex2::diffractiveDissVertex2(const particlePtr& beam,
					       const particlePtr& XSystem)
  : interactionVertex2()
{
  if (!beam) {
    printErr << "null pointer to beam particle. aborting." << endl;
    throw;
  }
  if (!XSystem) {
    printErr << "null pointer to particle representing X system. aborting." << endl;
    throw;
  }
  interactionVertex2::addInParticle (beam);
  interactionVertex2::addOutParticle(XSystem);
  if (_debug)
    printInfo << "contructed " << *this << endl;
}


diffractiveDissVertex2::diffractiveDissVertex2(const diffractiveDissVertex2& vert)
{
  *this = vert;
}


diffractiveDissVertex2::~diffractiveDissVertex2()
{ }


ostream&
diffractiveDissVertex2::print(ostream& out) const
{
  out << "diffractive dissociation vertex: "
      << "beam " << beam()->summary() << "  --->  "
      << XSystem()->summary();
  return out;
}


ostream&
diffractiveDissVertex2::dump(ostream& out) const
{
  out << "diffractive dissociation vertex: " << endl
      << "    beam: "     << *beam()    << endl
      << "    X system: " << *XSystem() << endl;
  return out;
}
