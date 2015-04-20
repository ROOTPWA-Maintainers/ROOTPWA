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
//      virtual base class for isobar decay amplitude independent of formalism
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <algorithm>
#include <cassert>

#include "TLorentzRotation.h"
#include "TMath.h"

#include "conversionUtils.hpp"
#include "factorial.hpp"
#include "isobarAmplitude.h"


using namespace std;
using namespace boost;
using namespace rpwa;


bool isobarAmplitude::_debug = false;


isobarAmplitude::isobarAmplitude()
	: _decay               (),
	  _useReflectivityBasis(true),
	  _boseSymmetrize      (true),
	  _isospinSymmetrize   (true),
	  _doSpaceInversion    (false),
	  _doReflection        (false)
{ }


isobarAmplitude::isobarAmplitude(const isobarDecayTopologyPtr& decay)
	: _useReflectivityBasis(true),
	  _boseSymmetrize      (true),
	  _isospinSymmetrize   (true),
	  _doSpaceInversion    (false),
	  _doReflection        (false)
{
	setDecayTopology(decay);
}


isobarAmplitude::~isobarAmplitude()
{ }


void
isobarAmplitude::setDecayTopology(const isobarDecayTopologyPtr& decay)
{
	if (not decay) {
		printErr << "null pointer to decay topology. Aborting..." << endl;
		throw;
	}
	if (not decay->checkTopology()) {
		printErr << "decay does not have the correct topology. Aborting..." << endl;
		throw;
	}
	if (not decay->checkConsistency()) {
		printErr << "decay is not consistent. Aborting..." << endl;
		throw;
	}
	_decay = decay;
	_decay->saveDecayToVertices(_decay);
}


void
isobarAmplitude::init()
{
	_symTermMaps.clear();
	// create first symmetrization entry with identity permutation map
	vector<unsigned int> identityPermMap;
	for (unsigned int i = 0; i < _decay->nmbFsParticles() ; ++i)
		identityPermMap.push_back(i);
	_symTermMaps.push_back(symTermMap(1, identityPermMap));
	// create final-state symmetriztion terms
	if(not initSymTermMaps()) {
		printErr << "Could not initialize amplitude symetrization maps." << endl;
		throw;
	}
}


// at the moment it is not possible to fix the helicity state of
// final state particles for cases where the amplitudes for the
// different helicity states have to be added incoherently
complex<double>
isobarAmplitude::amplitude() const
{
	const unsigned int nmbSymTerms = _symTermMaps.size();
	if (nmbSymTerms < 1) {
		printErr << "array of symmetrization terms is empty. make sure isobarAmplitude::init() "
		         << "was called. cannot calculate amplitude. returning 0." << endl;
		return 0;
	}
	// loop over all symmetrization terms; assumes that init() was called before
	complex<double> amp = 0;
	for (unsigned int i = 0; i < nmbSymTerms; ++i)
		amp += _symTermMaps[i].factor * symTermAmp(_symTermMaps[i].fsPartPermMap);
	return amp;
}


TLorentzRotation
isobarAmplitude::gjTransform(const TLorentzVector& beamLv,  // beam Lorentz-vector
                             const TLorentzVector& XLv)     // X  Lorentz-vector
{
	TLorentzVector beam    = beamLv;
	TLorentzVector X       = XLv;
	const TVector3 yGjAxis = beam.Vect().Cross(X.Vect());  // y-axis of Gottfried-Jackson frame
	// rotate so that yGjAxis becomes parallel to y-axis and beam momentum ends up in (x, z)-plane
	TRotation rot1;
	rot1.RotateZ(piHalf - yGjAxis.Phi());
	rot1.RotateX(yGjAxis.Theta() - piHalf);
	beam *= rot1;
	X    *= rot1;
	// boost to X RF
	TLorentzRotation boost;
	boost.Boost(-X.BoostVector());
	beam *= boost;
	// rotate about yGjAxis so that beam momentum is along z-axis
	TRotation rot2;
	rot2.RotateY(-signum(beam.X()) * beam.Theta());
	// construct total transformation
	TLorentzRotation gjTransform(rot1);
	gjTransform.Transform(boost);
	gjTransform.Transform(rot2);
	return gjTransform;
}


void
isobarAmplitude::spaceInvertDecay() const
{
	if (_debug)
		printDebug << "space inverting final state momenta." << endl;
	// transform final state particles into X rest frame
	const TLorentzVector&  beamLv  = _decay->productionVertex()->referenceLzVec();
	const TLorentzVector&  XLv     = _decay->XParticle()->lzVec();
	const TLorentzRotation gjTrans = gjTransform(beamLv, XLv);
	_decay->transformFsParticles(gjTrans);
	// perform parity transformation on final state particles in X rest frame
	for (unsigned int i = 0; i < _decay->nmbFsParticles(); ++i) {
		const particlePtr&    part     = _decay->fsParticles()[i];
		const TLorentzVector& fsPartLv = part->lzVec();
		part->setLzVec(TLorentzVector(-fsPartLv.Vect(), fsPartLv.E()));
	}
	// transform final state particles back to lab frame
	_decay->transformFsParticles(gjTrans.Inverse());
}


void
isobarAmplitude::reflectDecay() const
{
	if (_debug)
		printDebug << "reflecting final state momenta through production plane." << endl;
	// transform final state particles into X rest frame
	const TLorentzVector&  beamLv  = _decay->productionVertex()->referenceLzVec();
	const TLorentzVector&  XLv     = _decay->XParticle()->lzVec();
	const TLorentzRotation gjTrans = gjTransform(beamLv, XLv);
	_decay->transformFsParticles(gjTrans);
	// reflect final state particles through production plane
	for (unsigned int i = 0; i < _decay->nmbFsParticles(); ++i) {
		const particlePtr&    part     = _decay->fsParticles()[i];
		const TLorentzVector& fsPartLv = part->lzVec();
		part->setLzVec(TLorentzVector(fsPartLv.X(), -fsPartLv.Y(), fsPartLv.Z(), fsPartLv.E()));
	}
	// transform final state particles back to lab frame
	_decay->transformFsParticles(gjTrans.Inverse());
}


// assumes that all particles in the decay are mesons
// !!! in this primitive recursion scheme daughter amplitudes might be
//     called multiple times; for amplitudes with high-spin isobars a
//     large speed gain might be achievable by precalculating the
//     amplitudes at each vertex for all possible helicity values;
//     additional speed can be gained by checking that daughter
//     amplitudes are not zero before calculating the remaining parts
//     of the amplitude
complex<double>
isobarAmplitude::twoBodyDecayAmplitudeSum(const isobarDecayVertexPtr& vertex,           // current vertex
                                          const bool                  topVertex) const  // switches special treatment of X decay vertex; needed for reflectivity basis
{
	const particlePtr& parent    = vertex->parent();
	const particlePtr& daughter1 = vertex->daughter1();
	const particlePtr& daughter2 = vertex->daughter2();
	complex<double>    ampSum    = 0;
	for (int lambda1 = -daughter1->J(); lambda1 <= +daughter1->J(); lambda1 += 2) {
		// calculate decay amplitude for daughter 1
		daughter1->setSpinProj(lambda1);
		const isobarDecayVertexPtr& daughter1Vertex =
			dynamic_pointer_cast<isobarDecayVertex>(_decay->toVertex(daughter1));
		complex<double> daughter1Amp = 0;
		if (daughter1Vertex) {
			daughter1Amp = twoBodyDecayAmplitudeSum(daughter1Vertex, false);
		} else
			daughter1Amp = 1;
		for (int lambda2 = -daughter2->J(); lambda2 <= +daughter2->J(); lambda2 += 2) {
			// calculate decay amplitude for daughter 2
			daughter2->setSpinProj(lambda2);
			const isobarDecayVertexPtr& daughter2Vertex =
				dynamic_pointer_cast<isobarDecayVertex>(_decay->toVertex(daughter2));
			complex<double> daughter2Amp = 0;
			if (daughter2Vertex) {
				daughter2Amp = twoBodyDecayAmplitudeSum(daughter2Vertex, false);
			} else
				daughter2Amp = 1;
			complex<double> parentAmp = twoBodyDecayAmplitude(vertex, topVertex);
			complex<double> amp       = parentAmp * daughter1Amp * daughter2Amp;
			if (_debug)
				printDebug << "amplitude term for : "
				           << parent->name()    << " [lambda = " << spinQn(parent->spinProj()) << "] -> "
				           << daughter1->name() << " [lambda = " << spinQn(lambda1           ) << "] + "
				           << daughter2->name() << " [lambda = " << spinQn(lambda2           ) << "] = "
				           << "parent amp. = " << maxPrecisionDouble(parentAmp)
				           << " * daughter_1 amp = " << maxPrecisionDouble(daughter1Amp)
				           << " * daughter_2 amp = " << maxPrecisionDouble(daughter2Amp) << " = "
				           << maxPrecisionDouble(amp) << endl;
			ampSum += amp;
		}
	}
	if (_debug)
		printDebug << "decay amplitude for " << *vertex << " = " << maxPrecisionDouble(ampSum) << endl;
	return ampSum;
}


complex<double>
isobarAmplitude::symTermAmp(const vector<unsigned int>& fsPartPermMap) const
{
	// (re)set final state momenta
	if (not _decay->revertMomenta(fsPartPermMap)) {
		printErr << "problems reverting momenta in decay topology. cannot calculate amplitude. "
		         << "returning 0." << endl;
		return 0;
	}
	// transform daughters into their respective RFs
	transformDaughters();
	// calculate amplitude
	return twoBodyDecayAmplitudeSum(_decay->XIsobarDecayVertex(), true);
}


bool
isobarAmplitude::initSymTermMaps()
{

	vector<symTermMap> isoSymTermMaps;
	vector<symTermMap> boseSymTermMaps;
	unsigned int nmbIsoSymTerms = 0;
	unsigned int nmbBoseSymTerms = 0;
	if(_isospinSymmetrize) {
		// get isospin permutation maps
		isoSymTermMaps = _decay->getIsospinSymmetrization();
		nmbIsoSymTerms = isoSymTermMaps.size();
		if (nmbIsoSymTerms < 1) {
			printErr << "array of isospin-symmetrization terms is empty. "
			         << "cannot isospin-symmetrize amplitude." << endl;
			return false;
		}
		if (nmbIsoSymTerms == 1) {
			printInfo << "no isospin symmetrization needed for this amplitude." << endl;
			isoSymTermMaps[0].factor = 1;
		} else {
			// check isospin permutation maps
			bool isoPermMapsOkay = true;
			for (unsigned int i = 0; i < nmbIsoSymTerms; ++i)
				if (isoSymTermMaps[i].fsPartPermMap.size() != _decay->nmbFsParticles()) {
					printErr << "final-state permutation map for isospin symmetrization has wrong size = "
							 << isoSymTermMaps[i].fsPartPermMap.size() << " expected "
							 << _decay->nmbFsParticles() << " . cannot isospin-symmetrize amplitude." << endl;
					isoPermMapsOkay = false;
				}
			if (not isoPermMapsOkay)
				return false;
			printInfo << "Found " << nmbIsoSymTerms << " isospin symmetrization terms." << endl;
		}
	}
	if(_boseSymmetrize) {
		boseSymTermMaps = _decay->getBoseSymmetrization();
		nmbBoseSymTerms = boseSymTermMaps.size();
		if (nmbBoseSymTerms < 1) {
			printErr << "array of Bose-symmetrization terms is empty. "
			         << "cannot Bose-symmetrize amplitude." << endl;
			return false;
		}
		if (nmbBoseSymTerms == 1) {
			printInfo << "no Bose symmetrization needed for this amplitude." << endl;
		} else {
			printInfo << "Found " << nmbBoseSymTerms << " Bose symmetrization terms." << endl;
		}
	}
	if(nmbIsoSymTerms + nmbBoseSymTerms > 0) {
		_symTermMaps.clear();
	}
	for(unsigned int iso_i = 0; iso_i < nmbIsoSymTerms; ++iso_i) {
			complex<double> isoFactor = isoSymTermMaps[iso_i].factor;
			vector<unsigned int> isoSymTermMap = isoSymTermMaps[iso_i].fsPartPermMap;
		for(unsigned int bose_i = 0; bose_i < nmbBoseSymTerms; ++bose_i) {
			complex<double> boseFactor = boseSymTermMaps[bose_i].factor;
			vector<unsigned int> boseSymTermMap = boseSymTermMaps[bose_i].fsPartPermMap;
			vector<unsigned int> newSymTermMap(isoSymTermMap.size(), 0);
			assert(isoSymTermMap.size() == boseSymTermMap.size());
			for(unsigned int i = 0; i < boseSymTermMap.size(); ++i) {
				newSymTermMap[i] = boseSymTermMap[isoSymTermMap[i]];
			}
			complex<double> newFactor = isoFactor * boseFactor;
			_symTermMaps.push_back(symTermMap(newFactor, newSymTermMap));
		}
	}

	if(_debug) {
		printDebug << "Found the folloing symmetrization parametrization:" << endl;
		cout << "isospin terms: " << endl;
		for(unsigned int i = 0; i < nmbIsoSymTerms; ++i) {
			vector<unsigned int> map = isoSymTermMaps[i].fsPartPermMap;
			cout << isoSymTermMaps[i].factor << ": [" << map[0];
			for(unsigned int j = 1; j < map.size(); ++j) {
				cout << ", " << map[j];
			}
			cout << "]" << endl;
		}
		cout << "Bose terms: " << endl;
		for(unsigned int i = 0; i < nmbBoseSymTerms; ++i) {
			vector<unsigned int> map = boseSymTermMaps[i].fsPartPermMap;
			cout << boseSymTermMaps[i].factor << ": [" << map[0];
			for(unsigned int j = 1; j < map.size(); ++j) {
				cout << ", " << map[j];
			}
			cout << "]" << endl;
		}
		cout << "combined terms: " << endl;
		for(unsigned int i = 0; i < _symTermMaps.size(); ++i) {
			vector<unsigned int> map = _symTermMaps[i].fsPartPermMap;
			cout << _symTermMaps[i].factor << ": [" << map[0];
			for(unsigned int j = 1; j < map.size(); ++j) {
				cout << ", " << map[j];
			}
			cout << "]" << endl;
		}
	}

	return true;

};


ostream&
isobarAmplitude::printParameters(ostream& out) const
{
	out << name() << ": " << endl
	    << "    reflectivity basis ............... " << enDisabled(_useReflectivityBasis) << endl
	    << "    Bose symmetrization .............. " << enDisabled(_boseSymmetrize      ) << endl
	    << "    isospin symmetrization ........... " << enDisabled(_isospinSymmetrize   ) << endl
	    << "    space inversion of FS momenta .... " << enDisabled(_doSpaceInversion    ) << endl
	    << "    reflection through prod. plane ... " << enDisabled(_doReflection        ) << endl;
	return out;
}


ostream&
isobarAmplitude::print(ostream& out) const
{
	printParameters(out);
	out << *_decay;
	return out;
}
