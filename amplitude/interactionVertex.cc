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


#include "TClonesArray.h"

#include "utilities.h"
#include "interactionVertex.h"

	
using namespace std;
using namespace rpwa;


bool interactionVertex::_debug = false;


interactionVertex::interactionVertex()
  : _inParticles (),
    _outParticles()
{ }


interactionVertex::interactionVertex(const interactionVertex& vert)
{
  *this = vert;
}


interactionVertex::~interactionVertex()
{ }


interactionVertex&
interactionVertex::operator =(const interactionVertex& vert)
{
  if (this != &vert) {
    _inParticles  = vert._inParticles;
    _outParticles = vert._outParticles;
  }
  return *this;
}


interactionVertex*
interactionVertex::clone(const bool cloneInParticles,
			 const bool cloneOutParticles) const
{
  interactionVertex* vertexClone = new interactionVertex(*this);
  if (cloneInParticles)
    vertexClone->cloneInParticles();
  if (cloneOutParticles)
    vertexClone->cloneOutParticles();
  return vertexClone;
}


void
interactionVertex::clear()
{
  _inParticles.clear();
  _outParticles.clear();
}


bool
interactionVertex::addInParticle(const particlePtr& part)
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
interactionVertex::addOutParticle(const particlePtr& part)
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
interactionVertex::transformOutParticles(const TLorentzRotation& L)
{
  for (unsigned int i = 0; i < nmbOutParticles(); ++i)
    _outParticles[i]->transform(L);
}


bool
interactionVertex::readData(const TClonesArray&,
			    const TClonesArray&)
{
  printWarn << "function not implemented" << endl;
  return false;
}


ostream&
interactionVertex::print(ostream& out) const
{
  out << "interaction vertex: ";
  for (unsigned int i = 0; i < _inParticles.size(); ++i) {
    out << _inParticles[i]->qnSummary();
    if (i < _inParticles.size() - 1)
      out << "  +  ";
  }
  out << "  --->  ";
  for (unsigned int i = 0; i < _outParticles.size(); ++i) {
    out << _outParticles[i]->qnSummary();
    if (i < _outParticles.size() - 1)
      out << "  +  ";
  }
  return out;
}


ostream&
interactionVertex::dump(ostream& out) const
{
  out << "interaction vertex:" << endl;
  for (unsigned int i = 0; i < _inParticles.size(); ++i)
    out << "    incoming[" << i << "]: " << *_inParticles[i] << endl;
  for (unsigned int i = 0; i < _outParticles.size(); ++i)
    out << "    outgoing[" << i << "]: " << *_outParticles[i] << endl;
  return out;
}


ostream&
interactionVertex::printPointers(ostream& out) const
{
  out << "interaction vertex " << this << ": incoming particles: ";
  for (unsigned int i = 0; i < _inParticles.size(); ++i) {
    out << "[" << i << "] = " << _inParticles[i];
    if (i < _inParticles.size() - 1)
      out << ", ";
  }
  out << "; outgoing particles: ";
  for (unsigned int i = 0; i < _outParticles.size(); ++i) {
    out << "[" << i << "] = " << _outParticles[i];
    if (i < _outParticles.size() - 1)
      out << ", ";
  }
  out << endl;
  return out;
}


void
interactionVertex::cloneInParticles()
{
  for (unsigned int i = 0; i < nmbInParticles(); ++i) {
    particlePtr newPart(inParticles()[i]->clone());
    inParticles()[i] = newPart;
  }
}


void
interactionVertex::cloneOutParticles()
{
  for (unsigned int i = 0; i < nmbOutParticles(); ++i) {
    particlePtr newPart(outParticles()[i]->clone());
    outParticles()[i] = newPart;
  }
}
