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
//      class that describes final state vertex decay topology
//      class is just used for internal book keeping
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "utilities.h"
#include "fsVertex.h"

	
using namespace std;
using namespace rpwa;


bool fsVertex::_debug = false;


fsVertex::fsVertex(particle& fsParticle)
  : interactionVertex()
{
  if (_debug)
    printInfo << "contructing final state vertex: "
	      << "final state particle = " << fsParticle.name() << endl;
  interactionVertex::addInParticle(fsParticle);
}


fsVertex::fsVertex(const fsVertex& vert)
{
  *this = vert;
}


fsVertex::~fsVertex()
{ }


ostream&
fsVertex::print(ostream& out) const
{
  out << "final state vertex "
      // << "(data are " << ((!dataAreValid()) ? "not " : "") << "valid):" << endl
      << "    final state particle " << fsParticle()  << endl;
  return out;
}
