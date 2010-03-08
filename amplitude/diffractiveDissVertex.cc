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
#include "diffractiveDissVertex.h"

	
using namespace std;
using namespace rpwa;


diffractiveDissVertex::diffractiveDissVertex(particle& beam,
					     particle& XSystem)
  : vertex()
{
  vertex::addInParticle (beam);
  vertex::addOutParticle(XSystem);
}


diffractiveDissVertex::diffractiveDissVertex(const diffractiveDissVertex& vert)
{
  *this = vert;
}


diffractiveDissVertex::~diffractiveDissVertex()
{ }


diffractiveDissVertex&
diffractiveDissVertex::operator = (const diffractiveDissVertex& vert)
{
  if (this != &vert)
    vertex::operator = (vert);
  return *this;
}


bool
diffractiveDissVertex::dataValid() const
{
  if ((nmbInParticles() == 1) && (nmbOutParticles() == 1))
    return true;
  else
    return false;
}


ostream&
diffractiveDissVertex::print(ostream& out) const
{
  if (dataValid()) {
    out << "diffractive dissociation vertex beam particle:" << endl
	<< *(inParticles()[0]) << endl
	<< "diffractive dissociation vertex X-system:" << endl
	<< *(outParticles()[0]) << endl;
  } else {
    printWarn << "diffractive dissociation vertex data are incomplete or not valid" << endl;
    vertex::print(out);
  }
  return out;
}
