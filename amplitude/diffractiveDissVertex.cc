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
//      class that describes production vertex in diffractive
//      dissociation the beam-Reggeon-X vertex has exactly one
//      incoming beam and one outgoing particle (X), which
//      unambiguously defines the Reggeon kinematics; if the target
//      momentum is not specified, a fixed target is assumed; if the
//      recoil particle is not specified, elastic scattering is
//      assumed
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <boost/numeric/conversion/cast.hpp>

#include "TClonesArray.h"
#include "TClass.h"
#include "TObjString.h"
#include "TVector3.h"

#include "reportingUtils.hpp"
#include "reportingUtilsRoot.hpp"
#include "diffractiveDissVertex.h"


using namespace std;
using namespace rpwa;

using boost::numeric_cast;


bool diffractiveDissVertex::_debug = false;


diffractiveDissVertex::diffractiveDissVertex(const particlePtr& beam,
                                             const particlePtr& target,
                                             const particlePtr& XParticle,
                                             const particlePtr& recoil)
	: productionVertex(),
	  _beamMomCache   (),
	  _recoilMomCache (),
	  _targetMomCache ()
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
	if (recoil)
		interactionVertex::addOutParticle(recoil);
	else {
		if (_debug)
			printWarn << "recoil not specified. assuming elastic scattering." << endl;
		interactionVertex::addOutParticle(createParticle(*target));
	}
	if (_debug)
		printDebug << "constructed " << *this << endl;
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
		_recoilMomCache = vert._recoilMomCache;
		_targetMomCache = vert._targetMomCache;
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
		printDebug << "cloned " << *this << "; " << this << " -> " << vertexClone << " "
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


std::vector<std::complex<double> >
diffractiveDissVertex::productionAmps() const
{
	size_t numEvents = _beamMomCache.size();
	if (numEvents == 0) {
		printErr << "no data to calculate production amplitude. aborting." << endl;
		throw;
	}
	return vector<complex<double> >(numEvents, 1);
}


void
diffractiveDissVertex::setXFlavorQN()
{
	particle& X    = *XParticle();
	particle& beam = *(this->beam());
	X.setBaryonNmb  (beam.baryonNmb());
	X.setStrangeness(beam.strangeness());
	X.setCharm      (beam.charm());
	X.setBeauty     (beam.beauty());
}


bool
diffractiveDissVertex::initKinematicsData(const TClonesArray& prodKinPartNames)
{
	_nmbProdKinPart = 0;

	// check production vertex data
	const string partClassName = prodKinPartNames.GetClass()->GetName();
	if (partClassName != "TObjString") {
		printWarn << "production kinematics particle names are of type '" << partClassName
		          << "' and not TObjString." << endl;
		return false;
	}
	_nmbProdKinPart = numeric_cast<size_t>(prodKinPartNames.GetEntriesFast());
	if (_nmbProdKinPart > 3) {
		printWarn << "array of production kinematics particle names has wrong size: "
		          << _nmbProdKinPart << ". need at least beam (index 0); recoil (index 1) and "
		          << "target (index 2) are optional." << endl;
		return false;
	}

	// beam at index 0 (mandatory)
	bool success = true;
	const string beamName = ((TObjString*)prodKinPartNames[0])->GetString().Data();
	if (beamName != beam()->name()) {
		printWarn << "wrong particle at index 0 in production kinematics input data: "
		          << "read '" << beamName << "', "
		          << "expected beam particle '" << beam()->name() << "'" << endl;
		success = false;
	}

	// recoil at index 1 (optional)
	if (_nmbProdKinPart >= 2) {
		const string recoilName = ((TObjString*)prodKinPartNames[1])->GetString().Data();
		if (recoilName != recoil()->name()) {
			printWarn << "wrong particle at index 1 in production kinematics input data: "
			          << "read '" << recoilName << "', "
			          << "expected recoil particle '" << recoil()->name() << "'" << endl;
			success = false;
		}
	}

	// target at index 2 (optional)
	if (_nmbProdKinPart >= 3) {
		const string targetName = ((TObjString*)prodKinPartNames[2])->GetString().Data();
		if (targetName != target()->name()) {
			printWarn << "wrong particle at index 2 in production kinematics input data: "
			          << "read '" << targetName << "', "
			          << "expected target particle '" << target()->name() << "'" << endl;
			success = false;
		}
	}

	return success;
}


bool
diffractiveDissVertex::readKinematicsData(const vector<vector<TVector3> >& prodKinMomenta)
{
	_beamMomCache.clear  ();
	_recoilMomCache.clear();
	_targetMomCache.clear();

	// check production vertex data
	const size_t nmbProdKinMom = prodKinMomenta.size();
	if (nmbProdKinMom != _nmbProdKinPart) {
		printWarn << "array of production kinematics particle momenta has wrong size: "
		          << nmbProdKinMom << " (expected " << _nmbProdKinPart << "). "
		          << "cannot read production kinematics." << endl;
		return false;
	}

	bool success = true;
	for (size_t i = 0; i < prodKinMomenta[0].size(); ++i) {

		// set beam
		const TVector3& beamMom = prodKinMomenta[0][i];
		if (_debug)
			printDebug << "setting momentum of beam particle '" << beam()->name()
					   << "' to " << beamMom << " GeV" << endl;
		_beamMomCache.push_back(beamMom);

		// set recoil (optional)
		if (_nmbProdKinPart >= 2) {
			const TVector3& recoilMom = prodKinMomenta[1][i];
			if (_debug)
				printDebug << "setting momentum of recoil particle '" << recoil()->name()
						   << "' to " << recoilMom << " GeV" << endl;
			_recoilMomCache.push_back(recoilMom);
		}

		// set target (optional); if not defined fixed target is assumed
		if (_nmbProdKinPart >= 3) {
			const TVector3& targetMom = prodKinMomenta[2][i];
			if (_debug)
				printDebug << "setting momentum of target particle '" << target()->name()
						   << "' to " << targetMom << " GeV" << endl;
			_targetMomCache.push_back(targetMom);
		}
	}

	return success;
}


bool
diffractiveDissVertex::revertMomenta()
{
	if (_debug) {
		printDebug << "resetting beam momentum to "       << _beamMomCache   << " GeV" << endl
		           << "    resetting recoil momentum to " << _recoilMomCache << " GeV" << endl
		           << "    resetting target momentum to " << _targetMomCache << " GeV" << endl;
	}
	beam  ()->setMomenta(_beamMomCache  );
	recoil()->setMomenta(_recoilMomCache);
	target()->setMomenta(_targetMomCache);
	return true;
}


ostream&
diffractiveDissVertex::print(ostream& out) const
{
	out << name() << ": beam " << beam()->qnSummary() << "  +  target " << target()->qnSummary()
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
