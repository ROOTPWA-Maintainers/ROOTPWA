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
//      the beam-Reggeon-X vertex has exactly one incoming beam and
//      one outgoing particle (X), which unambiguously defines the Reggeon
//      kinematics; in addition the target has to be specified; if the recoil
//      particle is not specified, elastic scattering is assumed
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
                                             const particlePtr& target,
                                             const particlePtr& XParticle,
                                             const particlePtr& recoil)
	: productionVertex(),
	  _beamMomCache   (),
	  _targetMomCache (),
	  _recoilMomCache ()
{
	if (not beam) {
		printErr << "null pointer to beam particle. aborting." << endl;
		throw;
	}
	if (not target) {
		printErr << "null pointer to target particle. aborting." << endl;
		throw;
	}
	if (not XParticle) {
		printErr << "null pointer to particle representing X system. aborting." << endl;
		throw;
	}
	interactionVertex::addInParticle (beam);
	interactionVertex::addInParticle (target);
	interactionVertex::addOutParticle(XParticle);
	if (not recoil) {
		if (_debug)
			printWarn << "recoil not specified. assuming elastic scattering." << endl;
		interactionVertex::addOutParticle(createParticle(*target));
	}
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
		_beamMomCache   = vert._beamMomCache;
		_targetMomCache = vert._targetMomCache;
		_recoilMomCache = vert._recoilMomCache;
	}
	return *this;
}


diffractiveDissVertex*
diffractiveDissVertex::doClone(const bool cloneInParticles,
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
diffractiveDissVertex::addInParticle(const particlePtr&)
{
	if (_debug)
		printWarn << "cannot add incoming particle to " << *this << endl;
	return false;
}


bool
diffractiveDissVertex::addOutParticle(const particlePtr&)
{
	if (_debug)
		printWarn << "cannot add outgoing particle to " << *this << endl;
	return false;
}


complex<double>
diffractiveDissVertex::productionAmp() const
{
	return 1;
}


bool
diffractiveDissVertex::readData(const TClonesArray& prodKinParticles,
                                const TClonesArray& prodKinMomenta)
{
	_beamMomCache   = TVector3();
	_targetMomCache = TVector3();
	_recoilMomCache = TVector3();
	// check production vertex data
	bool         success       = true;
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
	const int nmbEntries = prodKinParticles.GetEntriesFast();
	if (nmbEntries != prodKinMomenta.GetEntriesFast()) {
		printWarn << "arrays of production kinematics particles and momenta have different sizes: "
		          << nmbEntries << " vs. " << prodKinMomenta.GetEntriesFast  () << endl;
		success = false;
	}
	if (nmbEntries < 2) {
		printWarn << "arrays of production kinematics particles and momenta have " << nmbEntries
		          << " entry/ies. need at least beam (index 0) and target (index 1); "
		          << "recoil (index 2) is optional." << endl;
		success = false;
	}
	if (not success)
		return false;
	// set beam
	const string beamName = ((TObjString*)prodKinParticles[0])->GetString().Data();
	if (beamName != beam()->name()) {
		printWarn << "cannot find entry for beam particle '" << beam()->name() << "' "
		          << "at index 0 in data." << endl;
		success = false;
	} else {
		if (_debug)
			printInfo << "setting momentum of beam particle " << beamName
			          << " to " << *((TVector3*)prodKinMomenta[0]) << " GeV" << endl;
		beam()->setMomentum(*((TVector3*)prodKinMomenta[0]));
		_beamMomCache = beam()->lzVec().Vect();
	}
	// set target
	const string targetName = ((TObjString*)prodKinParticles[1])->GetString().Data();
	if (targetName != target()->name()) {
		printWarn << "cannot find entry for target particle '" << target()->name() << "' "
		          << "at index 1 in data." << endl;
		success = false;
	} else {
		if (_debug)
			printInfo << "setting momentum of target particle " << targetName
			          << " to " << *((TVector3*)prodKinMomenta[1]) << " GeV" << endl;
		target()->setMomentum(*((TVector3*)prodKinMomenta[1]));
		_targetMomCache = target()->lzVec().Vect();
	}
	// set recoil (optional)
	if (nmbEntries >= 3) {
		const string recoilName = ((TObjString*)prodKinParticles[2])->GetString().Data();
		if (recoilName != recoil()->name()) {
			printWarn << "cannot find entry for recoil particle '" << recoil()->name() << "' "
			          << "at index 2 in data." << endl;
			success = false;
		} else {
			if (_debug)
				printInfo << "setting momentum of recoil particle " << recoilName
				          << " to " << *((TVector3*)prodKinMomenta[2]) << " GeV" << endl;
			recoil()->setMomentum(*((TVector3*)prodKinMomenta[2]));
			_recoilMomCache = recoil()->lzVec().Vect();
		}
	}
	return success;
}


bool
diffractiveDissVertex::revertMomenta()
{
	if (_debug)
		printInfo << "resetting beam momentum to " << _beamMomCache << " GeV" << endl;
	beam()->setMomentum(_beamMomCache);
	if (_debug)
		printInfo << "resetting target momentum to " << _targetMomCache << " GeV" << endl;
	target()->setMomentum(_targetMomCache);
	if (_debug)
		printInfo << "resetting recoil momentum to " << _recoilMomCache << " GeV" << endl;
	recoil()->setMomentum(_recoilMomCache);
	return true;
}


ostream&
diffractiveDissVertex::print(ostream& out) const
{
	out << name() << ": "
	    << "beam " << beam()->qnSummary() << "  +  target " << target()->qnSummary()
	    << "  --->  " << XParticle()->qnSummary() << "  +  recoil " << recoil()->qnSummary();
	return out;
}


ostream&
diffractiveDissVertex::dump(ostream& out) const
{
	out << name() << ": " << endl
	    << "    beam ..... " << *beam()      << endl
	    << "    target ... " << *target()    << endl
	    << "    X ........ " << *XParticle() << endl
	    << "    recoil ... " << *recoil()    << endl;
	return out;
}


ostream&
diffractiveDissVertex::printPointers(ostream& out) const
{
	out << name() << " " << this << ": "
	    << "beam particle = "   << beam()      << "; "
	    << "target particle = " << target()    << "; "
	    << "X particle = "      << XParticle() << "; "
	    << "recoil particle = " << recoil()    << endl;
	return out;
}
