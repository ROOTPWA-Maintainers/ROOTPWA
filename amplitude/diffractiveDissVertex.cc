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


#include "TClonesArray.h"
#include "TClass.h"
#include "TObjString.h"
#include "TVector3.h"

#include "utilities.h"
#include "diffractiveDissVertex.h"

	
using namespace std;
using namespace rpwa;


bool diffractiveDissVertex::_debug = false;


diffractiveDissVertex::diffractiveDissVertex(const particlePtr& beam,
					     const particlePtr& XSystem)
  : interactionVertex(),
    _beamMomCache    ()
{
  if (!beam) {
    printErr << "null pointer to beam particle. aborting." << endl;
    throw;
  }
  if (!XSystem) {
    printErr << "null pointer to particle representing X system. aborting." << endl;
    throw;
  }
  interactionVertex::addInParticle (beam);
  interactionVertex::addOutParticle(XSystem);
  if (_debug)
    printInfo << "constructed " << *this << endl;
}


diffractiveDissVertex::diffractiveDissVertex(const diffractiveDissVertex& vert)
{
  *this = vert;
}


diffractiveDissVertex::~diffractiveDissVertex()
{ }


diffractiveDissVertex&
diffractiveDissVertex::operator =(const diffractiveDissVertex& vert)
{
  if (this != &vert) {
    interactionVertex::operator =(vert);
    _beamMomCache = vert._beamMomCache;
  }
  return *this;
}


diffractiveDissVertex*
diffractiveDissVertex::clone(const bool cloneInParticles,
			     const bool cloneOutParticles) const
{
  diffractiveDissVertex* vertexClone = new diffractiveDissVertex(*this);
  if (cloneInParticles)
    vertexClone->cloneInParticles();
  if (cloneOutParticles)
    vertexClone->cloneOutParticles();
  if (_debug)
    printInfo << "cloned " << *this << "; " << this << " -> " << vertexClone << " "
	      << ((cloneInParticles ) ? "in" : "ex") << "cluding incoming particles, "
	      << ((cloneOutParticles) ? "in" : "ex") << "cluding outgoing particles" << std::endl;
  return vertexClone;
}


bool
diffractiveDissVertex::readData(const TClonesArray& prodKinParticles,
				const TClonesArray& prodKinMomenta)
{
  _beamMomCache = TVector3();
  // check inital state data
  bool success = true;
  const string partClassName = prodKinParticles.GetClass()->GetName();
  if (partClassName != "TObjString") {
    printWarn << "production kinematics particle names are of type " << partClassName
	      << " and not TObjString. cannot read production kinematics." << endl;
    success = false;
  }
  const string momClassName = prodKinMomenta.GetClass()->GetName();
  if (momClassName != "TVector3") {
    printWarn << "production kinematics momenta are of type " << momClassName
	      << " and not TVector3. cannot read production kinematics." << endl;
    success = false;
  }
  if (prodKinParticles.GetEntriesFast() != prodKinMomenta.GetEntriesFast()) {
    printWarn << "arrays of production kinematics particles and momenta have different sizes: "
	      << prodKinParticles.GetEntriesFast() << " vs. "
	      << prodKinMomenta.GetEntriesFast  () << endl;
    success = false;
  }
  if (!success)
    return false;
  // set inital state
  const string partName = ((TObjString*)prodKinParticles[0])->GetString().Data();
  if (partName != beam()->name()) {
    printWarn << "cannot find entry for beam particle '" << beam()->name() << "' in data." << endl;
    return false;
  }
  if (_debug)
    printInfo << "setting momentum of beam particle " << partName
	      << " to " << *((TVector3*)prodKinMomenta[0]) << " GeV" << endl;
  beam()->setMomentum(*((TVector3*)prodKinMomenta[0]));
  _beamMomCache = beam()->lzVec().Vect();
  return true;
}


bool
diffractiveDissVertex::revertMomenta()
{
  if (_debug)
    printInfo << "resetting beam momentum to " << _beamMomCache << " GeV" << endl;
  beam()->setMomentum(_beamMomCache);
  return true;
}


ostream&
diffractiveDissVertex::print(ostream& out) const
{
  out << "diffractive dissociation vertex: "
      << "beam " << beam()->qnSummary() << "  --->  "
      << XSystem()->qnSummary();
  return out;
}


ostream&
diffractiveDissVertex::dump(ostream& out) const
{
  out << "diffractive dissociation vertex: " << endl
      << "    beam: "     << *beam()    << endl
      << "    X system: " << *XSystem() << endl;
  return out;
}


ostream&
diffractiveDissVertex::printPointers(ostream& out) const
{
  out << "diffractive dissociation vertex " << this << ": "
      << "beam particle: "     << beam()    << "; "
      << "X system particle: " << XSystem() << endl;
  return out;
}
