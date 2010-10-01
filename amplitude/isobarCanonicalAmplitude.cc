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
//      general isobar decay amplitude in caninical formalism
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
#include "dFunction.hpp"
#include "isobarCanonicalAmplitude.h"

  
using namespace std;
using namespace boost;
using namespace rpwa;


bool isobarCanonicalAmplitude::_debug = false;


isobarCanonicalAmplitude::isobarCanonicalAmplitude()
	: isobarAmplitude()
{ }


isobarCanonicalAmplitude::isobarCanonicalAmplitude(const isobarDecayTopologyPtr& decay)
	: isobarAmplitude(decay)
{ }


isobarCanonicalAmplitude::~isobarCanonicalAmplitude()
{ }


void
isobarCanonicalAmplitude::transformDaughters() const
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
	// 2) transform daughters of isobar decay vertices to the respective rest frames
	for (unsigned int i = 1; i < _decay->nmbDecayVertices(); ++i) {  // exclude X-decay vertex
		const isobarDecayVertexPtr& vertex = _decay->isobarDecayVertices()[i];
		if (_debug)
			printInfo << "transforming all child particles of vertex " << *vertex
			          << " into " << vertex->parent()->name() << " daughter RF" << endl;
		// coordinate system does not change so this is just a simple Lorentz=boost
		const TVector3 rfBoost = -vertex->parent()->lzVec().BoostVector();
		// get all particles downstream of this vertex
		decayTopologyGraphType subGraph = _decay->dfsSubGraph(vertex);
		decayTopologyGraphType::edgeIterator iEd, iEdEnd;
		for (tie(iEd, iEdEnd) = subGraph.edges(); iEd != iEdEnd; ++iEd) {
			const particlePtr& part = subGraph.particle(*iEd);
			if (_debug)
				cout << "    transforming " << part->name() << " into "
				     << vertex->parent()->name() << " RF" << endl;
			part->transform(rfBoost);
		}
	}
}


// assumes that daughters were transformed into parent RF
complex<double>
isobarCanonicalAmplitude::twoBodyDecayAmplitude(const isobarDecayVertexPtr& vertex,
                                                const bool                  topVertex) const
{
	if (_debug)
		printInfo << "calculating two-body decay amplitude in canonical formalism "
		          << "for " << *vertex << endl;

	const particlePtr& parent    = vertex->parent();
	const particlePtr& daughter1 = vertex->daughter1();
	const particlePtr& daughter2 = vertex->daughter2();

	// calculate Clebsch-Gordan coefficient for S-S coupling
	const int    s1        = daughter1->J();
	const int    m1        = daughter1->spinProj();
	const int    s2        = daughter2->J();
	const int    m2        = daughter2->spinProj();
	const int    S         = vertex->S();
	const int    mS        = m1 + m2;
	const double ssClebsch = cgCoeff(s1, m1, s2, m2, S, mS, _debug);
	if (ssClebsch == 0)
		return 0;

	// calulate barrier factor
	const int    L  = vertex->L();
	const double q  = daughter1->lzVec().Vect().Mag();
	const double bf = barrierFactor(L, q, _debug);

	// calculate Breit-Wigner
	const complex<double> bw = vertex->massDepAmplitude();

	// calculate normalization factor
	const int    J    = parent->J();
	const double norm = normFactor(J, _debug);

	// sum over all possible spin projections of L
	const int       M     = parent->spinProj();
	const int       P     = parent->P();
	const int       refl  = parent->reflectivity();
	const double    phi   = daughter1->lzVec().Phi();  // use daughter1 as analyzer
	const double    theta = daughter1->lzVec().Theta();
	complex<double> amp	  = 0;
	for (int mL = -L; mL <= L; mL += 2) {
		// calculate Clebsch-Gordan coefficient for L-S coupling
		double LSClebsch;
		if (_useReflectivityBasis and topVertex) {
			// symmetrize L-S coupling term
			const int reflFactor = reflectivityFactor(J, P, M, refl);
			if (M == 0) {
				if (reflFactor == +1)
					LSClebsch = 0;
				else
					LSClebsch = cgCoeff(L, mL, S, mS, J, 0, _debug);
			} else {
				LSClebsch = 1 / rpwa::sqrt(2)
					* (               cgCoeff(L, mL, S, mS, J, +M, _debug)
					   - reflFactor * cgCoeff(L, mL, S, mS, J, -M, _debug));
			}
		} else
			LSClebsch = cgCoeff(L, mL, S, mS, J, M, _debug);
		if (LSClebsch == 0)
			continue;
		// multiply spherical harmonic
		amp += LSClebsch * sphericalHarmonicNoNorm<complex<double> >(L, mL, theta, phi, _debug);
	}

	// calculate decay amplitude
	amp *= norm * ssClebsch * bf * bw;
  
	if (_debug)
		printInfo << "two-body decay amplitude = " << maxPrecisionDouble(amp) << endl;
	return amp;
}
