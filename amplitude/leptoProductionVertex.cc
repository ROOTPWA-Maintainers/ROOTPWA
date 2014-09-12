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
//      class that describes leptoproduction vertex
//      the kinematics is defined by the incoming beam lepton, the
//      scattered lepton and the target; if the recoil particle is not
//      specified, elastic scattering is assumed
//
//      kinematic quantities are calculated based on COMPASS note
//      2009-6 using only Lorentz-scalars and without any
//      approximations concerning lepton mass or large Q^2
//
//      calculation of photon spin-density matrix is based on
//      Schilling and Wolf, NPB 61, 381 (1973)
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <cmath>

#include "TClonesArray.h"
#include "TClass.h"
#include "TObjString.h"
#include "TVector3.h"

#include "reportingUtils.hpp"
#include "reportingUtilsRoot.hpp"
#include "leptoProductionVertex.h"


using namespace std;
using namespace rpwa;


bool leptoProductionVertex::_debug = false;


leptoProductionVertex::leptoProductionVertex(const particlePtr& beamLepton,
                                             const particlePtr& target,
                                             const particlePtr& XParticle,
                                             const particlePtr& recoil)
	: productionVertex        (),
	  _longPol                (0),
	  _beamLeptonMomCache     (),
	  _scatteredLeptonMomCache(),
	  _recoilMomCache         (),
	  _targetMomCache         ()
{
	if (not beamLepton) {
		printErr << "null pointer to beam lepton particle. aborting." << endl;
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
	interactionVertex::addInParticle (beamLepton);
	interactionVertex::addInParticle (target);
	interactionVertex::addInParticle (createParticle("gamma"));  // virtual photon
	interactionVertex::addOutParticle(XParticle);
	if (not recoil) {
		if (_debug)
			printWarn << "recoil not specified. assuming elastic scattering." << endl;
		interactionVertex::addOutParticle(createParticle(*target));
	}
	interactionVertex::addOutParticle(createParticle(*beamLepton));  // lepton scatters elastically
	if (_debug)
		printDebug << "constructed " << *this << endl;
}


leptoProductionVertex::leptoProductionVertex(const leptoProductionVertex& vert)
{
	*this = vert;
}


leptoProductionVertex::~leptoProductionVertex()
{ }


leptoProductionVertex&
leptoProductionVertex::operator =(const leptoProductionVertex& vert)
{
	if (this != &vert) {
		interactionVertex::operator =(vert);
		_beamLeptonMomCache      = vert._beamLeptonMomCache;
		_scatteredLeptonMomCache = vert._scatteredLeptonMomCache;
		_recoilMomCache          = vert._recoilMomCache;
		_targetMomCache          = vert._targetMomCache;
	}
	return *this;
}


leptoProductionVertex*
leptoProductionVertex::doClone(const bool cloneInParticles,
                               const bool cloneOutParticles) const
{
	leptoProductionVertex* vertexClone = new leptoProductionVertex(*this);
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
leptoProductionVertex::addInParticle(const particlePtr&)
{
	if (_debug)
		printWarn << "cannot add incoming particle to " << *this << endl;
	return false;
}


bool
leptoProductionVertex::addOutParticle(const particlePtr&)
{
	if (_debug)
		printWarn << "cannot add outgoing particle to " << *this << endl;
	return false;
}


std::vector<complex<double> >
leptoProductionVertex::productionAmp() const
{

	const std::vector<TLorentzVector>& targetVec = target()->lzVecs();
	const std::vector<TLorentzVector>& beamVec            = beamLepton()->lzVecs();
	const std::vector<TLorentzVector>& scatteredLeptonVec = scatteredLepton()->lzVecs();
	const std::vector<TLorentzVector>& virtPhotonVec      = virtPhoton()->lzVecs();
	const std::vector<TLorentzVector>& XParticleVec       = XParticle()->lzVecs();

	std::vector<double> epsilons(targetVec.size());
	this->epsilon(epsilons);

	std::vector<double> deltas(targetVec.size());
	this->delta(epsilons, deltas);

	if(targetVec.size() != beamVec.size() || targetVec.size() != scatteredLeptonVec.size() ||
			targetVec.size() != virtPhotonVec.size() || targetVec.size() != XParticleVec.size()) {
		printErr << "size of per-event-data vectors does not match. aborting." << endl;
		throw;
	}

	std::vector<complex<double> > result(target()->numParallelEvents());
	// !! EVENT PARALLEL LOOP
	boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
	const unsigned int size = result.size();
	#pragma omp parallel for
	for(unsigned int i = 0; i < size; ++i) {

		// calculate azimuthal angle between lepton-scattering and
		// production plane in (virtual photon, target) CM system since
		TVector3             k1, k2, q, v;  // vectors of beam lepton, scattered lepton, virtual photon,
											// and X particle used to define lepton-scattering and
											// production plane
		const TLorentzVector targetLv = targetVec[i];
		if (targetLv.Vect() == TVector3(0, 0, 0)) {
			// fixed target case:
			// since going from fixed target frame into (virtual photon,
			// target) CM system involves only a boost along the virtual
			// photon direction and since the normals of both the
			// lepton-scattering and the production plane are perpendicular to
			// the virtual photon, the azimuthal angle can be calculated in
			// the lab frame as well
			k1 = beamVec[i].Vect();
			k2 = scatteredLeptonVec[i].Vect();
			q  = virtPhotonVec[i].Vect();
			v  = XParticleVec[i].Vect();
		} else {
			// general case
			// boost vectors to (virtual photon, target) CM system
			TLorentzVector       beamLeptonLv      = beamVec[i];
			TLorentzVector       scatteredLeptonLv = scatteredLeptonVec[i];
			TLorentzVector       virtPhotonLv      = virtPhotonVec[i];
			TLorentzVector       XParticleLv       = XParticleVec[i];
			const TLorentzVector photonTargetCM    = virtPhotonLv + targetLv;
			const TVector3       cmBoost           = photonTargetCM.BoostVector();
			beamLeptonLv.Boost(cmBoost);
			scatteredLeptonLv.Boost(cmBoost);
			virtPhotonLv.Boost(cmBoost);
			XParticleLv.Boost(cmBoost);
			k1 = beamLeptonLv.Vect();
			k2 = scatteredLeptonLv.Vect();  // use scattered lepton, because virtual photon is not directly measured
			q  = virtPhotonLv.Vect();
			v  = XParticleLv.Vect();
		}
		// calculate azimuthal angle in [-pi, +pi]
		const TVector3 leptonPlaneNormal     = k1.Cross(k2);
		const TVector3 productionPlaneNormal = q.Cross(v);
		const double   parProjection         = leptonPlaneNormal.Dot(productionPlaneNormal);
		const double   perpProjection        = leptonPlaneNormal.Cross(productionPlaneNormal).Mag();
		const double   phi                   = atan2(perpProjection, parProjection);

		// compute some kinematic variables
		const double epsilon = epsilons[i];
		const double delta   = deltas[i];
		//printInfo << "phi = " << phi << ", epsilon = " << maxPrecision(epsilon) << ", "
		//		  << "delta = " << maxPrecision(delta) << ", pol = " << _longPol << endl;
		// compute some intermediary terms
		const double          xi    = sqrt(epsilon * (1 + epsilon + 2 * delta));
		const double          zeta  = sqrt(xi * (1 - epsilon) / (1 + epsilon));
		const double          term  = _longPol * sqrt(1 - epsilon * epsilon);
		const complex<double> phase = exp(complex<double>(0, phi));
		//printInfo << "xi = " << xi << ", zeta = " << zeta << ", term = " << term << endl;

		// define lower triangle of virtual photon spin density matrix
		// rho_{lambda, lambda'}; lambda = -1, 0, 1 as in Schilling and
		// Wolf (common factor 1/2 is dropped):
		//   eq. 44 and 57 with 60, where alpha_2 was set to 0 (long. lepton
		//   polarization), into eq. 62
		complex<double> rho[3][3];
		// column with lambda' = -1
		rho[0][0] = 1 + term;                                     // lambda = -1
		rho[1][0] = (_longPol * zeta + xi) * phase;               // lambda =  0
		rho[2][0] = -epsilon * exp(complex<double>(0, 2 * phi));  // lambda = +1
		// column with lambda' =  0
		rho[1][1] = 2 * (epsilon + delta);                        // lambda =  0
		rho[2][1] = (_longPol * zeta - xi) * phase;               // lambda = +1
		// column with lambda' = +1
		rho[2][2] = 1 - term;                                     // lambda = +1
		// conjugated elements
		rho[0][1] = conj(rho[1][0]);
		rho[0][2] = conj(rho[2][0]);
		rho[1][2] = conj(rho[2][1]);
		//for (unsigned int j = 0; j < 3; ++j)
		//	for (unsigned int i = 0; i < 3; ++i)
		//		printInfo << "rho[" << i << "][" << j << "] = " << maxPrecisionDouble(rho[i][j]) << endl;
		// const complex<double> detRho
		// 	=   rho[0][0] * rho[1][1] * rho[2][2]
		// 	  + rho[0][1] * rho[1][2] * rho[2][0]
		// 	  + rho[1][0] * rho[2][1] * rho[0][2]
		// 	  - rho[0][2] * rho[1][1] * rho[2][0]
		// 	  - rho[1][2] * rho[0][0] * rho[2][1]
		// 	  - rho[0][1] * rho[2][2] * rho[1][0];
		//printInfo << "det[rho] = " << detRho << " vs. "
		//		  << (  (epsilon + delta) * (1 - epsilon * epsilon) * (1 - _longPol * _longPol)
		//				- (1 + epsilon) * (xi * xi + _longPol * _longPol * zeta * zeta)
		//				+ 2 * _longPol * _longPol * zeta * xi * sqrt(1 - epsilon * epsilon)) / 4
		//		  << endl;

		// perform Cholesky decomposition rho_ij = sum_r V_ir * V_jr^*, where V_ir is a lower
		// triangle matrix with real diagonal elements
		complex<double> V[3][3];
		// first column
		V[0][0] = sqrt(real(rho[0][0]));
		V[1][0] = rho[1][0] / real(V[0][0]);
		V[2][0] = rho[2][0] / real(V[0][0]);
		// second column
		V[1][1] = sqrt(real(rho[1][1]) - norm(V[1][0]));
		V[2][1] = (rho[2][1] - V[2][0] * conj(V[1][0])) / real(V[1][1]);
		// third column
		V[2][2] = sqrt(real(rho[2][2]) - norm(V[2][1]) - norm(V[2][0]));
		// zero elements
		V[0][1] = 0;
		V[0][2] = 0;
		V[1][2] = 0;
		//printInfo << "V[2][2]^2 = " << real(rho[2][2]) - norm(V[2][1]) - norm(V[2][0]) << ": "
		//		  << real(rho[2][2]) << " - " << norm(V[2][1]) << " - " << norm(V[2][0]) << endl;
		//for (unsigned int j = 0; j < 3; ++j)
		//	for (unsigned int i = 0; i < 3; ++i)
		//		printInfo << "V[" << i << "][" << j << "] = " << maxPrecisionDouble(V[i][j]) << endl;
		complex<double> rhoPrime[3][3];
		for (unsigned int j = 0; j < 3; ++j)
			for (unsigned int i = 0; i < 3; ++i) {
				rhoPrime[i][j] = 0;
				for (unsigned int r = 0; r < 3; ++r)
					rhoPrime[i][j] += V[i][r] * conj(V[j][r]);
			}
		//for (unsigned int j = 0; j < 3; ++j)
		//	for (unsigned int i = 0; i < 3; ++i)
		//		printInfo << "deltaRho[" << i << "][" << j << "] = "
		//				  << maxPrecisionDouble(rho[i][j] - rhoPrime[i][j]) << endl;


		// compute production amplitude for given photon helicity
		complex<double> prodAmp = 0;
		result[i] = prodAmp;

	}
	boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
	uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
	cout << "EPL: leptoProductionVertex::productionAmp timediff = " << timeDiff << endl;

	return result;
}


void
leptoProductionVertex::Q2(std::vector<double>& result) const
{
	const std::vector<TLorentzVector>& virtPhotonVec = virtPhoton()->lzVecs();

	if(result.size() != virtPhotonVec.size()) {
		printErr << "size of per-event-data vectors does not match. aborting." << endl;
		throw;
	}

	// !! EVENT PARALLEL LOOP
	boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
	const unsigned int size = virtPhotonVec.size();
	#pragma omp parallel for
	for(unsigned int i = 0; i < size; ++i) {
		result[i] = - virtPhotonVec[i].Mag2();
	}
	boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
	uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
	cout << "EPL: leptoProductionVertex::Q2 timediff = " << timeDiff << endl;
}


void
leptoProductionVertex::nu(std::vector<double>& result) const
{
	const std::vector<TLorentzVector>& targetVec = target()->lzVecs();
	const std::vector<TLorentzVector>& virtPhotonVec = virtPhoton()->lzVecs();
	double targetMass = target()->mass();

	if(result.size() != targetVec.size() || result.size() != virtPhotonVec.size()) {
		printErr << "size of per-event-data vectors does not match. aborting." << endl;
		throw;
	}

	// !! EVENT PARALLEL LOOP
	boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
	const unsigned int size = result.size();
	#pragma omp parallel for
	for(unsigned int i = 0; i < size; ++i) {
		result[i] = (targetVec[i] * virtPhotonVec[i]) / targetMass;
	}
	boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
	uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
	cout << "EPL: leptoProductionVertex::nu timediff = " << timeDiff << endl;
}


void
leptoProductionVertex::y(std::vector<double>& result) const
{
	//const std::vector<TLorentzVector>& beamVec = beamLepton()->lzVecs();

	const std::vector<TLorentzVector>& targetVec = target()->lzVecs();
	const std::vector<TLorentzVector>& virtPhotonVec = virtPhoton()->lzVecs();
	const std::vector<TLorentzVector>& beamVec = beamLepton()->lzVecs();

	if(result.size() != targetVec.size() || result.size() != virtPhotonVec.size() || result.size() != beamVec.size()) {
		printErr << "size of per-event-data vectors does not match. aborting." << endl;
		throw;
	}

	// !! EVENT PARALLEL LOOP
	boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
	const unsigned int size = result.size();
	#pragma omp parallel for
	for(unsigned int i = 0; i < size; ++i) {
		result[i] = (targetVec[i] * virtPhotonVec[i]) / (targetVec[i] * beamVec[i]);
	}
	boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
	uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
	cout << "EPL: leptoProductionVertex::y timediff = " << timeDiff << endl;
}


void
leptoProductionVertex::epsilon(std::vector<double>& result) const
{

	std::vector<double> Q2(result.size());
	this->Q2(Q2);

	std::vector<double> xBj(result.size());
	this->xBj(xBj);

	std::vector<double> y(result.size());
	this->y(y);

	const double targetMass2 = target()->mass2();
	const double beamLeptonMass2 = beamLepton()->mass2();

	// !! EVENT PARALLEL LOOP
	boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
	const unsigned int size = result.size();
	#pragma omp parallel for
	for(unsigned int i = 0; i < size; ++i) {

		const double xBj2  = xBj[i] * xBj[i];
		const double y2    = y[i] * y[i];

		// calculate kinematic values based on Lorentz invariants
		// see COMPASS note 2009-6
		const double gamma2         = 4 * xBj2 * targetMass2 / Q2[i];  // eq. 77d^2
		const double term1          = beamLeptonMass2 * y2 * gamma2 / Q2[i];
		const double term2          = 1 / ((1 - y[i]) * (1 - y[i]));
		const double u0             = sqrt(1 - term1 * (1 + term2) + term1 * term1 * term2);  // eq. 81
		const double E1E2           = Q2[i] * (1 - y[i]) / (y2 * gamma2);  // eq. 80
		const double Q2Min          = 2 * (E1E2 * (1 - u0) - beamLeptonMass2);  // eq. 5c and 82
		const double k1k2           = E1E2 * u0;  // eq. 82
		const double sin2ThetaHalf  = (Q2[i] - Q2Min) / (4 * k1k2);  // from eq. 5b
		const double tan2ThetaHalf  = sin2ThetaHalf / (1 - sin2ThetaHalf);
		const double term3          = 1 - Q2Min / Q2[i];
		const double oneOverEpsilon = 1 + 2 * (1 + 1 / gamma2) * tan2ThetaHalf / (term3 * term3);  // eq. 31 and 77d
		result[i] = 1 / oneOverEpsilon;

		// // calculate kinematic values based on Lorentz invariants
		// // see COMPASS note 2009-6
		// const double gamma2         = 4 * xBj2 * target()->mass2() / Q2;  // eq. 77d^2
		// const double term1          = beamLepton()->mass2() * y2 * gamma2 / Q2;
		// const double term2          = 1 / ((1 - y) * (1 - y));
		// const double u0             = sqrt(1 - term1 * (1 + term2) + term1 * term1 * term2);  // eq. 81
		// const double v0             = 1 + 2 * (  beamLepton()->mass2() / Q2
		//                                          - (1 - y) * (1 - u0) / (y2 * gamma2));  // eq. 85
		// const double v02            = v0 * v0;
		// const double term3          = (1 - y) * u0 * v0;
		// const double term4          = y2 * gamma2 / 4;
		// return (term3 - term4 * v02) / (term3 + y2 / 2 + term4 * (2 - v02));  // eq. 87

	}
	boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
	uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
	cout << "EPL: leptoProductionVertex::epsilon timediff = " << timeDiff << endl;

}


void
leptoProductionVertex::delta(std::vector<double>& result) const
{

	std::vector<double> Q2(result.size());
	this->Q2(Q2);

	std::vector<double> epsilon(result.size());
	this->epsilon(epsilon);

	const double beamMass2 = beamLepton()->mass2();

	// !! EVENT PARALLEL LOOP
	boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
	const unsigned int size = result.size();
	#pragma omp parallel for
	for(unsigned int i = 0; i < size; ++i) {
		result[i] = (2 * beamMass2 / Q2[i]) * (1 - epsilon[i]);
	}
	boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
	uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
	cout << "EPL: leptoProductionVertex::delta timediff = " << timeDiff << endl;

}


void
leptoProductionVertex::xBj(std::vector<double>& result) const
{

	std::vector<double> Q2(result.size());
	this->Q2(Q2);

	const std::vector<TLorentzVector>& targetVec = target()->lzVecs();
	const std::vector<TLorentzVector>& virtPhotonVec = virtPhoton()->lzVecs();

	if(result.size() != targetVec.size() || result.size() != virtPhotonVec.size()) {
		printErr << "size of per-event-data vectors does not match. aborting." << endl;
		throw;
	}

	// !! EVENT PARALLEL LOOP
	boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
	const unsigned int size = result.size();
	#pragma omp parallel for
	for(unsigned int i = 0; i < size; ++i) {
		result[i] = Q2[i] / (2 * (targetVec[i] * virtPhotonVec[i]));
	}
	boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
	uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
	cout << "EPL: leptoProductionVertex::xBj timediff = " << timeDiff << endl;

}


void
leptoProductionVertex::s(std::vector<double>& result) const
{
	const std::vector<TLorentzVector>& targetVec = target()->lzVecs();
	const std::vector<TLorentzVector>& virtPhotonVec = virtPhoton()->lzVecs();

	if(result.size() != targetVec.size() || result.size() != virtPhotonVec.size()) {
		printErr << "size of per-event-data vectors does not match. aborting." << endl;
		throw;
	}

	// !! EVENT PARALLEL LOOP
	boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
	const unsigned int size = result.size();
	#pragma omp parallel for
	for(unsigned int i = 0; i < size; ++i) {
		result[i] = (targetVec[i] + virtPhotonVec[i]).Mag2();
	}
	boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
	uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
	cout << "EPL: leptoProductionVertex::s timediff = " << timeDiff << endl;

}


void
leptoProductionVertex::W(std::vector<double>& result) const
{
	std::vector<double> s(result.size());
	this->s(s);

	// !! EVENT PARALLEL LOOP
	boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
	const unsigned int size = result.size();
	#pragma omp parallel for
	for(unsigned int i = 0; i < size; ++i) {
		result[i] = sqrt(s[i]);
	}
	boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
	uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
	cout << "EPL: leptoProductionVertex::W timediff = " << timeDiff << endl;
}


void
leptoProductionVertex::delta(const std::vector<double>& epsilon, std::vector<double>& result) const
{
	if(epsilon.size() != epsilon.size()) {
		printErr << "size of per-event-data vectors does not match. aborting." << endl;
		throw;
	}

	const double beamMass2 = beamLepton()->mass2();

	std::vector<double> Q2(result.size());
	this->Q2(Q2);

	// !! EVENT PARALLEL LOOP
	boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
	const unsigned int size = result.size();
	#pragma omp parallel for
	for(unsigned int i = 0; i < size; ++i) {
		result[i] = (2 * beamMass2 / Q2[i]) * (1 - epsilon[i]);
	}
	boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
	uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
	cout << "EPL: leptoProductionVertex::delta timediff = " << timeDiff << endl;
}


void
leptoProductionVertex::setXFlavorQN()
{
	//!!! check this
	particle& X    = *XParticle();
	particle& beam = *virtPhoton();
	X.setBaryonNmb  (beam.baryonNmb());
	X.setStrangeness(beam.strangeness());
	X.setCharm      (beam.charm());
	X.setBeauty     (beam.beauty());
}


bool
leptoProductionVertex::initKinematicsData(const TClonesArray& prodKinPartNames)
{
	_nmbProdKinPart = 0;

	// check production vertex data
	const string partClassName = prodKinPartNames.GetClass()->GetName();
	if (partClassName != "TObjString") {
		printWarn << "production kinematics particle names are of type '" << partClassName
		          << "' and not TObjString." << endl;
		return false;
	}
	_nmbProdKinPart = prodKinPartNames.GetEntriesFast();
	if (_nmbProdKinPart < 2) {
		printWarn << "array of production kinematics particle names has wrong size: "
		          << _nmbProdKinPart << ". need at least beam lepton (index 0) and "
		          << "scattered lepton (index 1); "
		          << "recoil (index 2) target (index 3) are optional." << endl;
		return false;
	}

	// beam lepton at index 0
	bool success = true;
	const string beamLeptonName = ((TObjString*)prodKinPartNames[0])->GetString().Data();
	if (beamLeptonName != beamLepton()->name()) {
		printWarn << "wrong particle at index 0 in production kinematics input data: "
		          << "read '" << beamLeptonName << "', "
		          << "expected beam lepton '" << beamLepton()->name() << "'" << endl;
		success = false;
	}

	// scattered lepton at index 1
	const string scatteredLeptonName = ((TObjString*)prodKinPartNames[1])->GetString().Data();
	if (scatteredLeptonName != scatteredLepton()->name()) {
		printWarn << "wrong particle at index 1 in production kinematics input data: "
		          << "read '" << scatteredLeptonName << "', "
		          << "expected scattered lepton '" << scatteredLepton()->name() << "'" << endl;
		success = false;
	}

	// recoil at index 2 (optional)
	if (_nmbProdKinPart >= 3) {
		const string recoilName = ((TObjString*)prodKinPartNames[2])->GetString().Data();
		if (recoilName != recoil()->name()) {
			printWarn << "wrong particle at index 2 in production kinematics input data: "
			          << "read '" << recoilName << "', "
			          << "expected recoil particle '" << recoil()->name() << "'" << endl;
			success = false;
		}
	}

	// target at index 3 (optional)
	if (_nmbProdKinPart >= 4) {
		const string targetName = ((TObjString*)prodKinPartNames[3])->GetString().Data();
		if (targetName != target()->name()) {
			printWarn << "wrong particle at index 3 in production kinematics input data: "
			          << "read '" << targetName << "', "
			          << "expected target particle '" << target()->name() << "'" << endl;
			success = false;
		}
	}

	return success;
}


bool
leptoProductionVertex::readKinematicsData(const vector<vector<TVector3> >& prodKinMomenta)
{
	_beamLeptonMomCache.clear();
	_scatteredLeptonMomCache.clear();
	_recoilMomCache.clear();
	_targetMomCache.clear();


	// check production vertex data
	const int nmbProdKinMom = prodKinMomenta.size();
	if (nmbProdKinMom != _nmbProdKinPart) {
		printWarn << "array of production kinematics particle momenta has wrong size: "
		          << nmbProdKinMom << " (expected " << _nmbProdKinPart << "). "
		          << "cannot read production kinematics." << endl;
		return false;
	}

	// set beam lepton
	bool success = true;
	for(unsigned int i = 0; i < prodKinMomenta[0].size(); ++i) {

		const TVector3& beamLeptonMom = prodKinMomenta[0][i];
		if (_debug)
			printDebug << "setting momentum of beam lepton '" << beamLepton()->name()
					   << "' to " << beamLeptonMom << " GeV" << endl;
		_beamLeptonMomCache.push_back(beamLeptonMom);

		// set scattered lepton
		const TVector3& scatteredLeptonMom = prodKinMomenta[1][i];
		if (_debug)
			printDebug << "setting momentum of scattered lepton '" << beamLepton()->name()
					   << "' to " << beamLeptonMom << " GeV" << endl;
		_scatteredLeptonMomCache.push_back(scatteredLeptonMom);


		// set recoil (optional)
		if (_nmbProdKinPart >= 3) {
			const TVector3& recoilMom = prodKinMomenta[2][i];
			if (_debug)
				printDebug << "setting momentum of recoil particle '" << recoil()->name()
						   << "' to " << recoilMom << " GeV" << endl;
			_recoilMomCache.push_back(recoilMom);
		}

		// set target (optional); if not defined fixed target is assumed
		if (_nmbProdKinPart >= 4) {
			const TVector3& targetMom = prodKinMomenta[3][i];
			if (_debug)
				printDebug << "setting momentum of target particle '" << target()->name()
						   << "' to " << targetMom << " GeV" << endl;
			_targetMomCache.push_back(targetMom);
		}
	}

	return success;
}


bool
leptoProductionVertex::revertMomenta()
{
	if (_debug) {
		printDebug << "resetting beam lepton momentum to " << _beamLeptonMomCache << " GeV" << endl
		           << "    resetting scattered lepton momentum to "
		           << _scatteredLeptonMomCache << " GeV" << endl
		           << "    resetting recoil momentum to " << _recoilMomCache << " GeV" << endl
		           << "    resetting target momentum to " << _targetMomCache << " GeV" << endl;
	}
	beamLepton     ()->setMomenta(_beamLeptonMomCache     );
	scatteredLepton()->setMomenta(_scatteredLeptonMomCache);
	recoil         ()->setMomenta(_recoilMomCache         );
	target         ()->setMomenta(_targetMomCache         );

	// set virtual photon

	std::vector<TLorentzVector>& beamLeptonVec = beamLepton()->mutableLzVec(); // mutable !!!
	const std::vector<TLorentzVector>& scatteredLeptonVec = scatteredLepton()->mutableLzVec(); // mutable !!!

	if(beamLeptonVec.size() != scatteredLeptonVec.size()) {
		printErr << "size of per-event-data vectors does not match. aborting." << std::endl;
		throw;
	}

	// !! EVENT PARALLEL LOOP
	boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
	const unsigned int size = beamLeptonVec.size();
	#pragma omp parallel for
	for(unsigned int i = 0; i < size; ++i) {
		beamLeptonVec[i] -= scatteredLeptonVec[i];
	}
	boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
	uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
	std::cout << "EPL: leptoProductionVertex::revertMomenta timediff = " << timeDiff << std::endl;

	return true;
}


ostream&
leptoProductionVertex::print(ostream& out) const
{
	out << name() << ": "
	    << "beam " << beamLepton()->qnSummary() << " P_L = " << _longPol << " -> "
	    << "scattered "  << scatteredLepton()->qnSummary() << " + virtual " << virtPhoton()->qnSummary()
	    << "  +  target " << target()->qnSummary() << "  --->  "
	    << XParticle()->qnSummary() << "  +  recoil " << recoil()->qnSummary();
	return out;
}


ostream&
leptoProductionVertex::dump(ostream& out) const
{
	out << name() << ": " << endl
	    << "    beam lepton ........ " << *beamLepton()      << endl
	    << "    target ............. " << *target()          << endl
	    << "    virtual photon ..... " << *virtPhoton()      << endl
	    << "    scattered lepton ... " << *scatteredLepton() << endl
	    << "    X .................. " << *XParticle()       << endl
	    << "    recoil ............. " << *recoil()          << endl;
	return out;
}


ostream&
leptoProductionVertex::printPointers(ostream& out) const
{
	out << name() << " " << this << ": "
	    << "beam lepton = "      << beamLepton()      << "; "
	    << "target particle = "  << target()          << "; "
	    << "virtual photon = "   << virtPhoton()      << "; "
	    << "scattered lepton = " << scatteredLepton() << "; "
	    << "X particle = "       << XParticle()       << "; "
	    << "recoil particle = "  << recoil()          << endl;
	return out;
}
