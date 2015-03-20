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
//      CUDA code for leptoProductionVertex.cc
//      This is only used when compiling with CUDA enabled.
//
//
// Author List:
//      Armin Gensler          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#include <cuda.h>
#include <thrust/transform.h>

#include "cudaUtils.hpp"
#include "leptoProductionVertex_cuda.h"

using namespace rpwa;

struct ThrustFunctor_leptoProductionVertex_productionAmps
{
	double _longPol;
	ThrustFunctor_leptoProductionVertex_productionAmps(double longPol):
		_longPol(longPol) {}
  
	template<typename T>
	HOST_DEVICE
	Complex operator()(T t) const
	{
		const LorentzVector& targetLv = thrust::get<0>(t);
		const LorentzVector& beamVec = thrust::get<1>(t);
		const LorentzVector& scatteredLeptonVec = thrust::get<2>(t);
		const LorentzVector& virtPhotonVec = thrust::get<3>(t);
		const LorentzVector& XParticleVec = thrust::get<4>(t);
		const double epsilon = thrust::get<5>(t);
		const double delta = thrust::get<6>(t);
	  
		//const LorentzRotation& lzRot = thrust::get<0>(t);
		//LorentzVector& lzVec = thrust::get<1>(t);
		
		// calculate azimuthal angle between lepton-scattering and
		// production plane in (virtual photon, target) CM system since
		Vector3             k1, k2, q, v;  // vectors of beam lepton, scattered lepton, virtual photon,
										   // and X particle used to define lepton-scattering and
										   // production plane
		if (targetLv.Vect() == Vector3(0, 0, 0)) {
			// fixed target case:
			// since going from fixed target frame into (virtual photon,
			// target) CM system involves only a boost along the virtual
			// photon direction and since the normals of both the
			// lepton-scattering and the production plane are perpendicular to
			// the virtual photon, the azimuthal angle can be calculated in
			// the lab frame as well
			k1 = beamVec           .Vect();
			k2 = scatteredLeptonVec.Vect();
			q  = virtPhotonVec     .Vect();
			v  = XParticleVec      .Vect();
		} else {
			// general case
			// boost vectors to (virtual photon, target) CM system
			LorentzVector       beamLeptonLv      = beamVec           ;
			LorentzVector       scatteredLeptonLv = scatteredLeptonVec;
			LorentzVector       virtPhotonLv      = virtPhotonVec     ;
			LorentzVector       XParticleLv       = XParticleVec      ;
			const LorentzVector photonTargetCM    = virtPhotonLv + targetLv;
			const Vector3       cmBoost           = photonTargetCM.BoostVector();
			beamLeptonLv.Boost     (cmBoost);
			scatteredLeptonLv.Boost(cmBoost);
			virtPhotonLv.Boost     (cmBoost);
			XParticleLv.Boost      (cmBoost);
			k1 = beamLeptonLv.Vect     ();
			k2 = scatteredLeptonLv.Vect();  // use scattered lepton, because virtual photon is not directly measured
			q  = virtPhotonLv.Vect     ();
			v  = XParticleLv.Vect      ();
		}
		// calculate azimuthal angle in [-pi, +pi]
		const Vector3 leptonPlaneNormal     = k1.Cross(k2);
		const Vector3 productionPlaneNormal = q.Cross(v);
		const double   parProjection        = leptonPlaneNormal.Dot(productionPlaneNormal);
		const double   perpProjection       = leptonPlaneNormal.Cross(productionPlaneNormal).Mag();
		const double   phi                  = atan2(perpProjection, parProjection);

		//printInfo << "phi = " << phi << ", epsilon = " << maxPrecision(epsilon) << ", "
		//		  << "delta = " << maxPrecision(delta) << ", pol = " << _longPol << endl;
		// compute some intermediary terms
		const double  xi    = sqrt(epsilon * (1 + epsilon + 2 * delta));
		const double  zeta  = sqrt(xi * (1 - epsilon) / (1 + epsilon));
		const double  term  = _longPol * sqrt(1 - epsilon * epsilon);
		const Complex phase = exp(Complex(0, phi));
		//printInfo << "xi = " << xi << ", zeta = " << zeta << ", term = " << term << endl;

		// define lower triangle of virtual photon spin density matrix
		// rho_{lambda, lambda'}; lambda = -1, 0, 1 as in Schilling and
		// Wolf (common factor 1/2 is dropped):
		//   eq. 44 and 57 with 60, where alpha_2 was set to 0 (long. lepton
		//   polarization), into eq. 62
		Complex rho[3][3];
		// column with lambda' = -1
		rho[0][0] = 1 + term;                             // lambda = -1
		rho[1][0] = (_longPol * zeta + xi) * phase;       // lambda =  0
		rho[2][0] = -epsilon * exp(Complex(0, 2 * phi));  // lambda = +1
		// column with lambda' =  0
		rho[1][1] = 2 * (epsilon + delta);                // lambda =  0
		rho[2][1] = (_longPol * zeta - xi) * phase;       // lambda = +1
		// column with lambda' = +1
		rho[2][2] = 1 - term;                             // lambda = +1
		// conjugated elements
		rho[0][1] = conj(rho[1][0]);
		rho[0][2] = conj(rho[2][0]);
		rho[1][2] = conj(rho[2][1]);
		//for (unsigned int j = 0; j < 3; ++j)
		//	for (unsigned int i = 0; i < 3; ++i)
		//		printInfo << "rho[" << i << "][" << j << "] = " << maxPrecisionDouble(rho[i][j]) << endl;
		// const Complex detRho
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
		Complex V[3][3];
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
		Complex rhoPrime[3][3];
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
		Complex prodAmp = 0;
		return prodAmp;
	}
};
		
void 
rpwa::thrust_leptoProductionVertex_productionAmps(
	const ParVector<LorentzVector>& targetVec,
	const ParVector<LorentzVector>& beamVec,
	const ParVector<LorentzVector>& scatteredLeptonVec,
	const ParVector<LorentzVector>& virtPhotonVec,
	const ParVector<LorentzVector>& XParticleVec,
	const ParVector<double>& epsilons,
	const ParVector<double>& deltas,
	ParVector<Complex>& result,
	double longPol)
{
	thrust::transform(
		thrust::make_zip_iterator(thrust::make_tuple(
			targetVec.begin(), 
			beamVec.begin(), 
			scatteredLeptonVec.begin(), 
			virtPhotonVec.begin(), 
			XParticleVec.begin(),
			epsilons.begin(), 
			deltas.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(
			targetVec.end(), 
			beamVec.end(), 
			scatteredLeptonVec.end(), 
			virtPhotonVec.end(), 
			XParticleVec.end(),
			epsilons.end(), 
			deltas.end())),
		result.begin(),
		ThrustFunctor_leptoProductionVertex_productionAmps(longPol));
}


struct ThrustFunctor_leptoProductionVertex_Q2
{
	HOST_DEVICE
	double operator()(const LorentzVector& virtPhotonVec)
	{
		return -virtPhotonVec.Mag2();
	}
};


void
rpwa::thrust_leptoProductionVertex_Q2(
	const ParVector<LorentzVector>& virtPhotonVec,
	ParVector<double>& result)
{
	thrust::transform(virtPhotonVec.begin(), virtPhotonVec.end(), result.begin(),
			ThrustFunctor_leptoProductionVertex_Q2());
}


struct ThrustFunctor_leptoProductionVertex_nu
{
	double targetMass;
	ThrustFunctor_leptoProductionVertex_nu(double targetMass):
		targetMass(targetMass){}

	HOST_DEVICE
	double operator()(const LorentzVector& targetVec, const LorentzVector& virtPhotonVec)
	{
		return (targetVec * virtPhotonVec) / targetMass;
	}
};


void
rpwa::thrust_leptoProductionVertex_nu(
	const ParVector<LorentzVector>& targetVec,
	const ParVector<LorentzVector>& virtPhotonVec,
	ParVector<double>& result,
	double targetMass)
{
	thrust::transform(targetVec.begin(), targetVec.end(), virtPhotonVec.begin(), result.begin(),
			ThrustFunctor_leptoProductionVertex_nu(targetMass));
}


struct ThrustFunctor_leptoProductionVertex_y
{
	template<typename T>
	HOST_DEVICE
	double operator()(T t)
	{
		const LorentzVector& targetVec = thrust::get<0>(t);
		const LorentzVector& virtPhotonVec = thrust::get<1>(t);
		const LorentzVector& beamVec = thrust::get<2>(t);
		return (targetVec * virtPhotonVec) / (targetVec * beamVec);
	}
};


void
rpwa::thrust_leptoProductionVertex_y(
	const ParVector<LorentzVector>& targetVec,
	const ParVector<LorentzVector>& virtPhotonVec,
	const ParVector<LorentzVector>& beamVec,
	ParVector<double>& result)
{
	thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(targetVec.begin(), virtPhotonVec.begin(), beamVec.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(targetVec.end(),   virtPhotonVec.end(),   beamVec.end())),
			result.begin(), ThrustFunctor_leptoProductionVertex_y());
}


struct ThrustFunctor_leptoProductionVertex_epsilon
{
	double targetMass2;
	double beamLeptonMass2;
	ThrustFunctor_leptoProductionVertex_epsilon(double targetMass2, double beamLeptonMass2):
		targetMass2(targetMass2), beamLeptonMass2(beamLeptonMass2) {}

	template<typename T>
	HOST_DEVICE
	double operator()(T t)
	{
		double Q2 = thrust::get<0>(t);
		double xBj = thrust::get<1>(t);
		double y = thrust::get<2>(t);

		const double xBj2  = xBj * xBj;
		const double y2    = y * y;

		// calculate kinematic values based on Lorentz invariants
		// see COMPASS note 2009-6
		const double gamma2         = 4 * xBj2 * targetMass2 / Q2;  // eq. 77d^2
		const double term1          = beamLeptonMass2 * y2 * gamma2 / Q2;
		const double term2          = 1 / ((1 - y) * (1 - y));
		const double u0             = sqrt(1 - term1 * (1 + term2) + term1 * term1 * term2);  // eq. 81
		const double E1E2           = Q2 * (1 - y) / (y2 * gamma2);  // eq. 80
		const double Q2Min          = 2 * (E1E2 * (1 - u0) - beamLeptonMass2);  // eq. 5c and 82
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
};


void
rpwa::thrust_leptoProductionVertex_epsilon(
	const ParVector<double>& Q2,
	const ParVector<double>& xBj,
	const ParVector<double>& y,
	ParVector<double>& result,
	double targetMass2,
	double beamLeptonMass2)
{
	thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(Q2.begin(), xBj.begin(), y.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(Q2.end(),   xBj.end(),   y.end())),
			result.begin(), ThrustFunctor_leptoProductionVertex_epsilon(targetMass2, beamLeptonMass2));
}


struct ThrustFunctor_leptoProductionVertex_delta
{
	double beamMass2;
	ThrustFunctor_leptoProductionVertex_delta(double beamMass2):
		beamMass2(beamMass2) {}

	HOST_DEVICE
	double operator()(double Q2, double epsilon)
	{
		return (2 * beamMass2 / Q2) * (1 - epsilon);
	}
};


void
rpwa::thrust_leptoProductionVertex_delta(
	const ParVector<double>& Q2,
	const ParVector<double>& epsilons,
	ParVector<double>& result,
	double beamMass2)
{
	thrust::transform(Q2.begin(), Q2.end(), epsilons.begin(), result.begin(),
			ThrustFunctor_leptoProductionVertex_delta(beamMass2));
}


struct ThrustFunctor_leptoProductionVertex_xBj
{
	template<typename T>
	HOST_DEVICE
	double operator()(T t)
	{
		const LorentzVector& targetVec = thrust::get<0>(t);
		const LorentzVector& virtPhotonVec = thrust::get<1>(t);
		double Q2 = thrust::get<2>(t);
		return Q2 / (2 * (targetVec * virtPhotonVec));
	}
};


void
rpwa::thrust_leptoProductionVertex_xBj(
	const ParVector<LorentzVector>& targetVec,
	const ParVector<LorentzVector>& virtPhotonVec,
	const ParVector<double>& Q2,
	ParVector<double>& result)
{
	thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(targetVec.begin(), virtPhotonVec.begin(), Q2.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(targetVec.end(),   virtPhotonVec.end(),   Q2.end())),
			result.begin(), ThrustFunctor_leptoProductionVertex_xBj());
}


struct ThrustFunctor_leptoProductionVertex_s
{
	HOST_DEVICE
	double operator()(const LorentzVector& targetVec, const LorentzVector& virtPhotonVec)
	{
		return (targetVec + virtPhotonVec).Mag2();
	}
};


void
rpwa::thrust_leptoProductionVertex_s(
	const ParVector<LorentzVector>& targetVec,
	const ParVector<LorentzVector>& virtPhotonVec,
	ParVector<double>& result)
{
	thrust::transform(targetVec.begin(), targetVec.end(), virtPhotonVec.begin(), result.begin(),
			ThrustFunctor_leptoProductionVertex_s());
}


struct ThrustFunctor_leptoProductionVertex_W
{
	HOST_DEVICE
	double operator()(double s)
	{
		return sqrt(s);
	}
};


void
rpwa::thrust_leptoProductionVertex_W(
	const ParVector<double>& s,
	ParVector<double>& result)
{
	thrust::transform(s.begin(), s.end(), result.begin(),
			ThrustFunctor_leptoProductionVertex_W());
}


void
rpwa::thrust_leptoProductionVertex_revertMomenta(
	const ParVector<LorentzVector>& scatteredLeptonVec,
	ParVector<LorentzVector>& result)
{
	thrust::transform(scatteredLeptonVec.begin(), scatteredLeptonVec.end(), result.begin(),
			thrust::negate<LorentzVector>());
}
