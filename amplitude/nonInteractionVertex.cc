
#include "nonInteractionVertex.h"

#include<TClass.h>
#include<TClonesArray.h>

#include "reportingUtils.hpp"
#include "reportingUtilsRoot.hpp"

using namespace rpwa;
using namespace std;


bool nonInteractionVertex::_debug = false;


nonInteractionVertex::nonInteractionVertex(const particlePtr& XParticle)
	: productionVertex(),
	  _XParticleCache()
{
	if (not XParticle) {
		printErr << "null pointer to XParticle. aborting." << endl;
		throw;
	}
	interactionVertex::addOutParticle(XParticle);
	if (_debug) {
		printDebug << "constructed " << *this << endl;
	}
}


bool
nonInteractionVertex::addInParticle(const particlePtr&)
{
	if (_debug) {
		printWarn << "cannot add incoming particle to " << *this << endl;
	}
	return false;
}


bool
nonInteractionVertex::addOutParticle(const particlePtr&)
{
	if (_debug) {
		printWarn << "cannot add outgoing particle to " << *this << endl;
	}
	return false;
}


bool
nonInteractionVertex::initKinematicsData(const TClonesArray& prodKinPartNames)
{
	// check production vertex data
	const string partClassName = prodKinPartNames.GetClass()->GetName();
	if (partClassName != "TObjString") {
		printWarn << "production kinematics particle names are of type '" << partClassName
		          << "' and not TObjString." << endl;
		return false;
	}
	if (prodKinPartNames.GetEntriesFast() != 1) {
		printWarn << "array of production kinematics particle names has wrong size: "
		          << prodKinPartNames.GetEntriesFast() << ", has to be 1." << endl;
		return false;
	}

	const string particleName = ((TObjString*)prodKinPartNames[0])->GetString().Data();
	if (particleName != XParticle()->name()) {
		printWarn << "wrong particle in production kinematics input data: "
		          << "read '" << particleName << "', "
		          << "expected particle '" << XParticle()->name() << "'" << endl;
		return false;
	}
	return true;
}


bool
nonInteractionVertex::readKinematicsData(const TClonesArray& prodKinMomenta)
{
	_XParticleCache   = TVector3();

	// check production vertex data
	const int nmbProdKinMom = prodKinMomenta.GetEntriesFast();
	if (nmbProdKinMom != 1) {
		printWarn << "array of production kinematics particle momenta has wrong size: "
		          << nmbProdKinMom << " (expected 1). "
		          << "cannot read production kinematics." << endl;
		return false;
	}

	TVector3* XParticleMom = dynamic_cast<TVector3*>(prodKinMomenta[0]);
	if (XParticleMom) {
		if (_debug) {
			printDebug << "setting momentum of beam particle '" << XParticle()->name()
			           << "' to " << *XParticleMom << " GeV" << endl;
		}
		XParticle()->setMomentum(*XParticleMom);
		_XParticleCache = XParticle()->momentum();
	} else {
		printWarn << "production kinematics data entry [0] is not of type TVector3. "
		          << "cannot read beam particle momentum." << endl;
		return false;
	}

	return true;
}


bool
nonInteractionVertex::revertMomenta()
{
	if (_debug) {
		printDebug << "resetting XParticle momentum to " << _XParticleCache << " GeV" << endl;
	}
	XParticle()->setMomentum(_XParticleCache);
	return true;
}


ostream&
nonInteractionVertex::print(ostream& out) const
{
	out << name() << ": XParticle " << XParticle()->qnSummary();
	return out;
}
