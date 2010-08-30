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

#include "utilities.h"
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
	  _targetMomCache         (),
	  _scatteredLeptonMomCache(),
	  _recoilMomCache         ()
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
		printInfo << "constructed " << *this << endl;
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
		_targetMomCache          = vert._targetMomCache;
		_recoilMomCache          = vert._recoilMomCache;
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
		printInfo << "cloned " << *this << "; " << this << " -> " << vertexClone << " "
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


complex<double>
leptoProductionVertex::productionAmp() const
{
	// calculate azimuthal angle between lepton-scattering and
	// production plane in (virtual photon, target) CM system since
	TVector3             k1, k2, q, v;  // vectors of beam lepton, scattered lepton, virtual photon,
	                                    // and X particle used to define lepton-scattering and
	                                    // production plane
	const TLorentzVector targetLv = target()->lzVec();
	if (targetLv.Vect() == TVector3(0, 0, 0)) {
		// fixed target case:
		// since going from fixed target frame into (virtual photon,
		// target) CM system involves only a boost along the virtual
		// photon direction and since the normals of both the
		// lepton-scattering and the production plane are perpendicular to
		// the virtual photon, the azimuthal angle can be calculated in
		// the lab frame as well
		k1 = beamLepton()->lzVec().Vect();
		k2 = scatteredLepton()->lzVec().Vect();
		q  = virtPhoton()->lzVec().Vect();
		v  = XParticle()->lzVec().Vect();
	} else {
		// general case
		// boost vectors to (virtual photon, target) CM system
		TLorentzVector       beamLeptonLv      = beamLepton()->lzVec();
		TLorentzVector       scatteredLeptonLv = scatteredLepton()->lzVec();
		TLorentzVector       virtPhotonLv      = virtPhoton()->lzVec();
		TLorentzVector       XParticleLv       = XParticle()->lzVec();
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
	const double epsilon = this->epsilon();
	const double delta   = this->delta(epsilon);
	printInfo << "phi = " << phi << ", epsilon = " << maxPrecision(epsilon) << ", "
	          << "delta = " << maxPrecision(delta) << ", pol = " << _longPol << endl;
	// compute some intermediary terms
	const double          xi    = sqrt(epsilon * (1 + epsilon + 2 * delta));
	const double          zeta  = sqrt(xi * (1 - epsilon) / (1 + epsilon));
	const double          term  = _longPol * sqrt(1 - epsilon * epsilon);
	const complex<double> phase = exp(complex<double>(0, phi));
	printInfo << "xi = " << xi << ", zeta = " << zeta << ", term = " << term << endl;
	
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
	for (unsigned int j = 0; j < 3; ++j)
		for (unsigned int i = 0; i < 3; ++i)
			printInfo << "rho[" << i << "][" << j << "] = " << maxPrecisionDouble(rho[i][j]) << endl;
	const complex<double> detRho =   rho[0][0] * rho[1][1] * rho[2][2]
		                             + rho[0][1] * rho[1][2] * rho[2][0]
		                             + rho[1][0] * rho[2][1] * rho[0][2]
		                             - rho[0][2] * rho[1][1] * rho[2][0]
		                             - rho[1][2] * rho[0][0] * rho[2][1]
		                             - rho[0][1] * rho[2][2] * rho[1][0];
	printInfo << "det[rho] = " << detRho << " vs. "
	          << (  (epsilon + delta) * (1 - epsilon * epsilon) * (1 - _longPol * _longPol)
	              - (1 + epsilon) * (xi * xi + _longPol * _longPol * zeta * zeta)
	              + 2 * _longPol * _longPol * zeta * xi * sqrt(1 - epsilon * epsilon)) / 4
	          << endl;

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
  printInfo << "V[2][2]^2 = " << real(rho[2][2]) - norm(V[2][1]) - norm(V[2][0]) << ": "
            << real(rho[2][2]) << " - " << norm(V[2][1]) << " - " << norm(V[2][0]) << endl;
  for (unsigned int j = 0; j < 3; ++j)
	  for (unsigned int i = 0; i < 3; ++i)
		  printInfo << "V[" << i << "][" << j << "] = " << maxPrecisionDouble(V[i][j]) << endl;
  complex<double> rhoPrime[3][3];
  for (unsigned int j = 0; j < 3; ++j)
	  for (unsigned int i = 0; i < 3; ++i) {
		  rhoPrime[i][j] = 0;
		  for (unsigned int r = 0; r < 3; ++r)
			  rhoPrime[i][j] += V[i][r] * conj(V[j][r]);
	  }
  for (unsigned int j = 0; j < 3; ++j)
	  for (unsigned int i = 0; i < 3; ++i)
		  printInfo << "deltaRho[" << i << "][" << j << "] = "
		            << maxPrecisionDouble(rho[i][j] - rhoPrime[i][j]) << endl;


	// compute production amplitude for given photon helicity
	complex<double> prodAmp = 0;
	return prodAmp;
}


double
leptoProductionVertex::epsilon() const
{
	const double Q2    = this->Q2();
	const double xBj   = this->xBj();
	const double xBj2  = xBj * xBj;
	const double y     = this->y();
	const double y2    = y * y;

	// calculate kinematic values based on Lorentz invariants
	// see COMPASS note 2009-6
	const double gamma2         = 4 * xBj2 * target()->mass2() / Q2;  // eq. 77d^2
	const double term1          = beamLepton()->mass2() * y2 * gamma2 / Q2;
	const double term2          = 1 / ((1 - y) * (1 - y));
	const double u0             = sqrt(1 - term1 * (1 + term2) + term1 * term1 * term2);  // eq. 81
	const double E1E2           = Q2 * (1 - y) / (y2 * gamma2);  // eq. 80
	const double Q2Min          = 2 * (E1E2 * (1 - u0) - beamLepton()->mass2());  // eq. 5c and 82
	const double k1k2           = E1E2 * u0;  // eq. 82
	const double sin2ThetaHalf  = (Q2 - Q2Min) / (4 * k1k2);  // from eq. 5b
	const double tan2ThetaHalf  = sin2ThetaHalf / (1 - sin2ThetaHalf);
	const double term3          = 1 - Q2Min / Q2;
	const double oneOverEpsilon = 1 + 2 * (1 + 1 / gamma2) * tan2ThetaHalf / (term3 * term3);  // eq. 31 and 77d
	return 1 / oneOverEpsilon;

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


bool
leptoProductionVertex::readData(const TClonesArray& prodKinParticles,
                                const TClonesArray& prodKinMomenta)
{
	_beamLeptonMomCache      = TVector3();
	_targetMomCache          = TVector3();
	_scatteredLeptonMomCache = TVector3();
	_recoilMomCache          = TVector3();
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
	if (nmbEntries < 3) {
		printWarn << "arrays of production kinematics particles and momenta have " << nmbEntries
		          << " entry/ies. need at least beam lepton (index 0), target (index 1), and "
		          << "scattered lepton (index 2); recoil (index 3) is optional." << endl;
		success = false;
	}
	if (not success)
		return false;
	// set beam lepton
	const string beamLeptonName = ((TObjString*)prodKinParticles[0])->GetString().Data();
	if (beamLeptonName != beamLepton()->name()) {
		printWarn << "cannot find entry for beam lepton '" << beamLepton()->name() << "' "
		          << "at index 0 in data." << endl;
		success = false;
	} else {
		if (_debug)
			printInfo << "setting momentum of beam lepton " << beamLeptonName
			          << " to " << *((TVector3*)prodKinMomenta[0]) << " GeV" << endl;
		beamLepton()->setMomentum(*((TVector3*)prodKinMomenta[0]));
		_beamLeptonMomCache = beamLepton()->lzVec().Vect();
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
	// set scattered lepton
	const string scatteredLeptonName = ((TObjString*)prodKinParticles[2])->GetString().Data();
	if (scatteredLeptonName != scatteredLepton()->name()) {
		printWarn << "cannot find entry for scattered lepton '" << scatteredLepton()->name() << "' "
		          << "at index 2 in data." << endl;
		success = false;
	} else {
		if (_debug)
			printInfo << "setting momentum of scattered lepton " << scatteredLeptonName
			          << " to " << *((TVector3*)prodKinMomenta[2]) << " GeV" << endl;
		scatteredLepton()->setMomentum(*((TVector3*)prodKinMomenta[2]));
		_scatteredLeptonMomCache = scatteredLepton()->lzVec().Vect();
	}
	// set recoil (optional)
	if (nmbEntries >= 4) {
		const string recoilName = ((TObjString*)prodKinParticles[3])->GetString().Data();
		if (recoilName != recoil()->name()) {
			printWarn << "cannot find entry for recoil particle '" << recoil()->name() << "' "
			          << "at index 3 in data." << endl;
			success = false;
		} else {
			if (_debug)
				printInfo << "setting momentum of recoil particle " << recoilName
				          << " to " << *((TVector3*)prodKinMomenta[3]) << " GeV" << endl;
			recoil()->setMomentum(*((TVector3*)prodKinMomenta[3]));
			_recoilMomCache = recoil()->lzVec().Vect();
		}
	}
	// set virtual photon
	virtPhoton()->setLzVec(beamLepton()->lzVec() - scatteredLepton()->lzVec());
	return success;
}


bool
leptoProductionVertex::revertMomenta()
{
	if (_debug)
		printInfo << "resetting beam lepton momentum to " << _beamLeptonMomCache << " GeV" << endl;
	beamLepton()->setMomentum(_beamLeptonMomCache);
	if (_debug)
		printInfo << "resetting target momentum to " << _targetMomCache << " GeV" << endl;
	target()->setMomentum(_targetMomCache);
	if (_debug)
		printInfo << "resetting scattered lepton momentum to " << _scatteredLeptonMomCache << " GeV"
		          << endl;
	scatteredLepton()->setMomentum(_scatteredLeptonMomCache);
	if (_debug)
		printInfo << "resetting recoil momentum to " << _recoilMomCache << " GeV" << endl;
	recoil()->setMomentum(_recoilMomCache);
	// set virtual photon
	virtPhoton()->setLzVec(beamLepton()->lzVec() - scatteredLepton()->lzVec());
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
