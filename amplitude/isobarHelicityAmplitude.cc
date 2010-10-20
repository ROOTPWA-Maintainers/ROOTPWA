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
//      general isobar decay amplitude in helicity formalism
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

#include "utilities.h"
#include "clebschGordanCoeff.hpp"
#include "dFunction.hpp"
#include "isobarHelicityAmplitude.h"

  
using namespace std;
using namespace boost;
using namespace rpwa;


bool isobarHelicityAmplitude::_debug = false;


isobarHelicityAmplitude::isobarHelicityAmplitude()
	: isobarAmplitude()
{ }


isobarHelicityAmplitude::isobarHelicityAmplitude(const isobarDecayTopologyPtr& decay)
	: isobarAmplitude(decay)
{ }


isobarHelicityAmplitude::~isobarHelicityAmplitude()
{ }


TLorentzRotation
isobarHelicityAmplitude::hfTransform(const TLorentzVector& daughterLv)
{
	TLorentzVector daughter = daughterLv;
	const TVector3 zAxisParent(0, 0, 1);  // take z-axis as defined in parent frame
	const TVector3 yHfAxis = zAxisParent.Cross(daughter.Vect());  // y-axis of helicity frame
	// rotate so that yHfAxis becomes parallel to y-axis and zHfAxis ends up in (x, z)-plane
	TRotation rot1;
	rot1.RotateZ(-yHfAxis.Phi());
	rot1.RotateY(piHalf - yHfAxis.Theta());
	rot1.RotateZ(piHalf);
	daughter *= rot1;
	// rotate about yHfAxis so that daughter momentum is along z-axis
	TRotation rot2;
	rot2.RotateY(-signum(daughter.X()) * daughter.Theta());
	daughter *= rot2;
	// boost to daughter RF
	rot1.Transform(rot2);
	TLorentzRotation hfTransform(rot1);
	hfTransform.Boost(-daughter.BoostVector());
	return hfTransform;
}


void
isobarHelicityAmplitude::transformDaughters() const
{
	// calculate Lorentz-vectors of all isobars
	_decay->calcIsobarLzVec();
	// modify event for testing purposes
	if (_doSpaceInversion) {
		spaceInvertDecay();
		// recalculate Lorentz-vectors of all isobars
		_decay->calcIsobarLzVec();
	}
	if (_doReflection) {
		reflectDecay();
		// recalculate Lorentz-vectors of all isobars
		_decay->calcIsobarLzVec();
	}
	// calculate Lorentz-transformations into the correct frames for the
	// daughters in the decay vertices
	// 1) transform daughters of all decay vertices into Gottfried-Jackson frame
	const TLorentzVector&  beamLv  = _decay->productionVertex()->referenceLzVec();
	const TLorentzVector&  XLv     = _decay->XParticle()->lzVec();
	const TLorentzRotation gjTrans = gjTransform(beamLv, XLv);
	for (unsigned int i = 0; i < _decay->nmbDecayVertices(); ++i) {
		const isobarDecayVertexPtr& vertex = _decay->isobarDecayVertices()[i];
		if (_debug)
			printInfo << "transforming outgoing particles of vertex " << *vertex
			          << " into " << vertex->parent()->name() << " Gottfried-Jackson RF" << endl;
		vertex->transformOutParticles(gjTrans);
	}
	// 2) transform daughters of isobar decay vertices to the respective helicity frames
	for (unsigned int i = 1; i < _decay->nmbDecayVertices(); ++i) {  // exclude X-decay vertex
		const isobarDecayVertexPtr& vertex = _decay->isobarDecayVertices()[i];
		if (_debug)
			printInfo << "transforming all child particles of vertex " << *vertex
			          << " into " << vertex->parent()->name() << " helicity RF" << endl;
		const TLorentzRotation hfTrans = hfTransform(vertex->parent()->lzVec());
		// get all particles downstream of this vertex
		decayTopologyGraphType subGraph = _decay->dfsSubGraph(vertex);
		decayTopologyGraphType::edgeIterator iEd, iEdEnd;
		for (tie(iEd, iEdEnd) = subGraph.edges(); iEd != iEdEnd; ++iEd) {
			const particlePtr& part = subGraph.particle(*iEd);
			if (_debug)
				cout << "    transforming " << part->name() << " into "
				     << vertex->parent()->name() << " helicity RF" << endl;
			part->transform(hfTrans);
		}
	}
}


// assumes that daughters were transformed into parent RF
complex<double>
isobarHelicityAmplitude::twoBodyDecayAmplitude(const isobarDecayVertexPtr& vertex,
                                               const bool                  topVertex) const
{
	if (_debug)
		printInfo << "calculating two-body decay amplitude in helicity formalism for " << *vertex << endl;

	const particlePtr& parent    = vertex->parent();
	const particlePtr& daughter1 = vertex->daughter1();
	const particlePtr& daughter2 = vertex->daughter2();

	// calculate Clebsch-Gordan coefficient for L-S coupling
	const int    L         = vertex->L();
	const int    S         = vertex->S();
	const int    J         = parent->J();
	const int    lambda1   = daughter1->spinProj();
	const int    lambda2   = daughter2->spinProj();
	const int    lambda    = lambda1 - lambda2;
	const double lsClebsch = clebschGordanCoeff<double>(L, 0, S, lambda, J, lambda, _debug);
	if (lsClebsch == 0)
		return 0;

	// calculate Clebsch-Gordan coefficient for S-S coupling
	const int    s1        = daughter1->J();
	const int    s2        = daughter2->J();
	const double ssClebsch = clebschGordanCoeff<double>(s1, lambda1, s2, -lambda2, S, lambda, _debug);
	if (ssClebsch == 0)
		return 0;

	// calculate D-function
	const int       Lambda = parent->spinProj();
	const int       P      = parent->P();
	const int       refl   = parent->reflectivity();
	const double    phi    = daughter1->lzVec().Phi();  // use daughter1 as analyzer
	const double    theta  = daughter1->lzVec().Theta();
	complex<double> DFunc;
	if (topVertex and _useReflectivityBasis)
		DFunc = DFunctionReflConj<complex<double> >(J, Lambda, lambda, P, refl, phi, theta, 0, _debug);
	else
		DFunc = DFunctionConj<complex<double> >(J, Lambda, lambda, phi, theta, 0, _debug);

	// calulate barrier factor
	const double q  = daughter1->lzVec().Vect().Mag();
	const double bf = barrierFactor(L, q, _debug);

	// calculate Breit-Wigner
	const complex<double> bw = vertex->massDepAmplitude();

	// calculate normalization factor
	const double norm = normFactor(L, _debug);

	// calculate decay amplitude
	complex<double> amp = norm * DFunc * lsClebsch * ssClebsch * bf * bw;
  
	if (_debug)
		printInfo << "two-body decay amplitude = " << maxPrecisionDouble(amp) << endl;
	return amp;
}
