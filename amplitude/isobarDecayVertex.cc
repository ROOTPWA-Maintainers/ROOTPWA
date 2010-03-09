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
//      class that describes decay vertex of isobar into two particles
//      the isobar -> particle1 + particle 2 vertex has exactly one
//      incoming mother and two outgoing daughter particle
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "utilities.h"
#include "isobarDecayVertex.h"

	
using namespace std;
using namespace rpwa;


isobarDecayVertex::isobarDecayVertex(particle&          mother,
				     particle&          daughter1,
				     particle&          daughter2,
				     const unsigned int L,
				     const unsigned int S)
  : interactionVertex(),
    _L    (L),
    _S    (S)
{
  interactionVertex::addInParticle (mother);
  interactionVertex::addOutParticle(daughter1);
  interactionVertex::addOutParticle(daughter2);
}


isobarDecayVertex::isobarDecayVertex(const isobarDecayVertex& vert)
{
  *this = vert;
}


isobarDecayVertex::~isobarDecayVertex()
{ }


isobarDecayVertex&
isobarDecayVertex::operator = (const isobarDecayVertex& vert)
{
  if (this != &vert) {
    interactionVertex::operator = (vert);
    _L = vert._L;
    _S = vert._S;
  }
  return *this;
}


ostream&
isobarDecayVertex::print(ostream& out) const
{
  out << "isobar decay vertex data are "
      << ((!dataAreValid()) ? "not " : "") << "valid:" << endl
      << "    mother "     << *(inParticles()[0])  << endl
      << "    daughter 1 " << *(outParticles()[0]) << endl
      << "    daughter 2 " << *(outParticles()[1]) << endl
      << "    L = " << _L << ", 2S = " << _S << endl;
  return out;
}
