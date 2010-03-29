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
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "utilities.h"
#include "interactionVertex2.h"

	
using namespace std;
using namespace rpwa;


bool interactionVertex2::_debug = false;


interactionVertex2::interactionVertex2()
  : _inParticles (),
    _outParticles()
{ }


interactionVertex2::interactionVertex2(const interactionVertex2& vert)
{
  *this = vert;
}


interactionVertex2::~interactionVertex2()
{ }


interactionVertex2&
interactionVertex2::operator =(const interactionVertex2& vert)
{
  if (this != &vert) {
    _inParticles  = vert._inParticles;
    _outParticles = vert._outParticles;
  }
  return *this;
}


interactionVertex2*
interactionVertex2::clone(const bool cloneInParticles,
			  const bool cloneOutParticles) const
{
  interactionVertex2* newVertex = new interactionVertex2(*this);
  if (cloneInParticles)
    for (unsigned int i = 0; i < newVertex->nmbInParticles(); ++i) {
      particlePtr newPart(&(newVertex->inParticles()[i]->clone()));
      newVertex->inParticles()[i] = newPart;
    }
  if (cloneOutParticles)
    for (unsigned int i = 0; i < newVertex->nmbOutParticles(); ++i) {
      particlePtr newPart(&(newVertex->outParticles()[i]->clone()));
      newVertex->outParticles()[i] = newPart;
    }
  return newVertex;
}


bool
interactionVertex2::addInParticle(const particlePtr& part)
{
  if (!part) {
    printErr << "trying to add null pointer to incoming particles. aborting." << endl;
    throw;
  }
  if (_debug)
    printInfo << "adding incoming " << *part << endl;
  _inParticles.push_back(part);
  return true;
}


bool
interactionVertex2::addOutParticle(const particlePtr& part)
{
  if (!part) {
    printErr << "trying to add null pointer to outgoing particles. aborting." << endl;
    throw;
  }
  if (_debug)
    printInfo << "adding outgoing " << *part << endl;
  _outParticles.push_back(part);
  return true;
}


void
interactionVertex2::transformOutParticles(const TLorentzRotation& L)
{
  for (unsigned int i = 0; i < nmbOutParticles(); ++i)
    _outParticles[i]->transform(L);
}


ostream&
interactionVertex2::print(ostream& out) const
{
  out << "interaction vertex: ";
  for (unsigned int i = 0; i < _inParticles.size(); ++i) {
    out << _inParticles[i]->summary();
    if (i < _inParticles.size() - 1)
      out << "  +  ";
  }
  out << "  --->  ";
  for (unsigned int i = 0; i < _outParticles.size(); ++i) {
    out << _outParticles[i]->summary();
    if (i < _outParticles.size() - 1)
      out << "  +  ";
  }
  return out;
}


ostream&
interactionVertex2::dump(ostream& out) const
{
  out << "interaction vertex:" << endl;
  for (unsigned int i = 0; i < _inParticles.size(); ++i)
    out << "    incoming[" << i << "]: " << *_inParticles[i] << endl;
  for (unsigned int i = 0; i < _outParticles.size(); ++i)
    out << "    outgoing[" << i << "]: " << *_outParticles[i] << endl;
  return out;
}
