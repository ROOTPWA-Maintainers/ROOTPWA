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
//      base class that desbribes general interaction vertex between particles
//      !NOTE! class stores pointers to particles, which are inserted via
//             references in order to ensure existence of objects; the calling
//             code has to ensure that lifetime of the particle instances is
//             longer than life time of the vertex instances they are assigned to
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "utilities.h"
#include "interactionVertex.h"

	
using namespace std;
using namespace rpwa;


bool interactionVertex::_debug = false;


interactionVertex::interactionVertex()
  : _inParticles (),
    _outParticles(),
    _dataValid   (false)
{ }


interactionVertex::interactionVertex(const interactionVertex& vert)
{
  *this = vert;
}


interactionVertex::~interactionVertex()
{ }


interactionVertex&
interactionVertex::operator = (const interactionVertex& vert)
{
  if (this != &vert) {
    _inParticles  = vert._inParticles;
    _outParticles = vert._outParticles;
    _dataValid    = vert._dataValid;
  }
  return *this;
}


// interactionVertex&
// interactionVertex::operator *= (const lorentzTransform& L)
// {
//   for (list<particle>::iterator i = _children.begin(); i != _children.end(); ++i)
//     *i *= L;
//   return *this;
// }


interactionVertex&
interactionVertex::clone() const
{
  interactionVertex* newVertex = new interactionVertex(*this);
  return *newVertex;
}


bool
interactionVertex::addInParticle (particle& part)
{
  if (_debug)
    printInfo << "adding incoming " << part << endl;
  _inParticles.push_back(&part);
  return true;
}


bool
interactionVertex::addOutParticle(particle& part)
{
  if (_debug)
    printInfo << "adding outgoing " << part << endl;
  _outParticles.push_back(&part);
  return true;
}


ostream&
interactionVertex::print(ostream& out) const
{
  out << "interaction vertex "
      << "(data are " << ((!dataAreValid()) ? "not " : "") << "valid):" << endl
      << "    incoming particles:" << endl;
  for (unsigned int i = 0; i < _inParticles.size(); ++i)
    out << "        " << *(_inParticles[i]) << endl;
  out << "    outgoing particles:" << endl;
  for (unsigned int i = 0; i < _outParticles.size(); ++i)
    out << "        " << *(_outParticles[i]) << endl;
  return out;
}
