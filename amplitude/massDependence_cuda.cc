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
//      CUDA code for massDependence.cc
//      This is only used when compiling with CUDA enabled.
//
//
// Author List:
//      Armin Gensler          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#include <complex>
#include <cuda.h>
#include <thrust/transform.h>

#include "typedefs.h"
#include "cudaUtils.hpp"
#include "dFunction.hpp"
#include "mathUtils.hpp"
#include "massDependence_cuda.h"

using namespace rpwa;
using namespace std;


struct ThrustFunctor_massDependence_flatRangeMassDependence_amp
{
	const double parentMass;
	const double parentWidth;
	ThrustFunctor_massDependence_flatRangeMassDependence_amp(double parentMass, double parentWidth):
		parentMass(parentMass), parentWidth(parentWidth) {}

	HOST_DEVICE
	Complex operator()(const LorentzVector& parentVec)
	{
		if (fabs(parentVec.M() - parentMass) < parentWidth / 2) {
			return 1;
		} else {
			return 0;
		}
	}
};

void
rpwa::thrust_massDependence_flatRangeMassDependence_amp(
	const ParVector<LorentzVector>& parentVec,
	ParVector<Complex>& result,
	double parentMass,
	double parentWidth)
{
	thrust::transform(parentVec.begin(), parentVec.end(), result.begin(),
			ThrustFunctor_massDependence_flatRangeMassDependence_amp(parentMass, parentWidth));
}

struct ThrustFunctor_massDependence_relativisticBreitWigner_amp
{
	const double parentMass;
	const double parentWidth;
	const int L;
	ThrustFunctor_massDependence_relativisticBreitWigner_amp(double parentMass, double parentWidth, int L):
		parentMass(parentMass), parentWidth(parentWidth), L(L) {}

	template<typename T>
	HOST_DEVICE
	Complex operator()(T t)
	{
		const LorentzVector& parentVec = thrust::get<0>(t);
		const LorentzVector& daughter1Vec = thrust::get<1>(t);
		const LorentzVector& daughter2Vec = thrust::get<2>(t);

		// get Breit-Wigner parameters
		const double M      = parentVec.M();     // parent mass
		const double m1     = daughter1Vec.M();  // daughter 1 mass
		const double m2     = daughter2Vec.M();  // daughter 2 mass
		const double q      = breakupMomentum(M,  m1, m2);
		const double M0     = parentMass;        // resonance peak position
		const double q02    = breakupMomentumSquared(M0, m1, m2, true);
		// !NOTE! the following is incorrect but this is how it was done in PWA2000
		const double q0     = sqrt(fabs(q02));
		const double Gamma0 = parentWidth;       // resonance peak width

		Complex bw = breitWigner<Complex>(M, M0, Gamma0, L, q, q0);
		return bw;
	}
};

void
rpwa::thrust_massDependence_relativisticBreitWigner_amp(
	const ParVector<LorentzVector>& parentVec,
	const ParVector<LorentzVector>& daughter1Vec,
	const ParVector<LorentzVector>& daughter2Vec,
	ParVector<Complex>& result,
	double parentMass,
	double parentWidth,
	unsigned int L)
{
	thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(parentVec.begin(), daughter1Vec.begin(), daughter2Vec.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(parentVec.end(),   daughter1Vec.end(),   daughter2Vec.end())),
			result.begin(),
			ThrustFunctor_massDependence_relativisticBreitWigner_amp(parentMass, parentWidth, L));
}

struct ThrustFunctor_massDependence_constWidthBreitWigner_amp
{
	const double M0;
	const double Gamma0;
	ThrustFunctor_massDependence_constWidthBreitWigner_amp(double M0, double Gamma0):
		M0(M0), Gamma0(Gamma0) {}

	HOST_DEVICE
	Complex operator()(const LorentzVector& parentVec)
	{
		// get Breit-Wigner parameters
		const double M = parentVec.M();  // parent mass

		// A / (B - iA) = (A / (B^2 + A^2)) * (B + iA)
		const double          A  = M0 * Gamma0;
		const double          B  = M0 * M0 - M * M;
		const Complex bw = (A / (B * B + A * A)) * Complex(B, A);
		// const Complex bw = (M0 * Gamma0) / (M0 * M0 - M * M - imag * M0 * Gamma0);
		return bw;
	}
};

void
rpwa::thrust_massDependence_constWidthBreitWigner_amp(
	const ParVector<LorentzVector>& parentVec,
	ParVector<Complex>& result,
	double M0,
	double Gamma0)
{
	thrust::transform(parentVec.begin(), parentVec.end(), result.begin(),
			ThrustFunctor_massDependence_constWidthBreitWigner_amp(M0, Gamma0));
}

struct ThrustFunctor_massDependence_rhoBreitWigner_amp
{
	const double parentMass;
	const double parentWidth;
	ThrustFunctor_massDependence_rhoBreitWigner_amp(double parentMass, double parentWidth):
		parentMass(parentMass), parentWidth(parentWidth) {}

	template<typename T>
	HOST_DEVICE
	Complex operator()(T t)
	{
		const LorentzVector& parentVec = thrust::get<0>(t);
		const LorentzVector& daughter1Vec = thrust::get<1>(t);
		const LorentzVector& daughter2Vec = thrust::get<2>(t);

		// get Breit-Wigner parameters
		const double M      = parentVec.M();         // parent mass
		const double m1     = daughter1Vec.M();  // daughter 1 mass
		const double m2     = daughter2Vec.M();  // daughter 2 mass
		const double q2     = breakupMomentumSquared(M,  m1, m2);
		const double q      = sqrt(q2);
		const double M0     = parentMass;              // resonance peak position
		const double q02    = breakupMomentumSquared(M0, m1, m2);
		const double q0     = sqrt(q02);
		const double Gamma0 = parentWidth;             // resonance peak width

		const double F      = 2 * q2 / (q02 + q2);
		const double Gamma  = Gamma0 * (M0 / M) * (q / q0) * F;
		// in the original publication the width reads
		// Gamma = Gamma0 * (q / q0) * F

		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double  A  = M0 * Gamma0 * sqrt(F);
		// in the original publication A reads
		// A = sqrt(M0 * Gamma0 * (m / q0) * F)
		const double  B  = M0 * M0 - M * M;
		const double  C  = M0 * Gamma;
		const Complex bw = (A / (B * B + C * C)) * Complex(B, C);
		// return (M0 * Gamma0 * sqrt(F)) / (M0 * M0 - M * M - imag * M0 * Gamma);
		return bw;
	}
};

void
rpwa::thrust_massDependence_rhoBreitWigner_amp(
	const ParVector<LorentzVector>& parentVec,
	const ParVector<LorentzVector>& daughter1Vec,
	const ParVector<LorentzVector>& daughter2Vec,
	ParVector<Complex>& result,
	double parentMass,
	double parentWidth)
{
	thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(parentVec.begin(), daughter1Vec.begin(), daughter2Vec.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(parentVec.end(),   daughter1Vec.end(),   daughter2Vec.end())),
			result.begin(),
			ThrustFunctor_massDependence_rhoBreitWigner_amp(parentMass, parentWidth));
}

struct ThrustFunctor_massDependence_f0980BreitWigner_amp
{
	const double parentMass;
	const double parentWidth;
	ThrustFunctor_massDependence_f0980BreitWigner_amp(double parentMass, double parentWidth):
		parentMass(parentMass), parentWidth(parentWidth) {}

	template<typename T>
	HOST_DEVICE
	Complex operator()(T t)
	{
		const LorentzVector& parentVec = thrust::get<0>(t);
		const LorentzVector& daughter1Vec = thrust::get<1>(t);
		const LorentzVector& daughter2Vec = thrust::get<2>(t);

		// get Breit-Wigner parameters
		const double M      = parentVec.M();         // parent mass
		const double m1     = daughter1Vec.M();  // daughter 1 mass
		const double m2     = daughter2Vec.M();  // daughter 2 mass
		const double q      = breakupMomentum(M,  m1, m2);
		const double M0     = parentMass;              // resonance peak position
		const double q0     = breakupMomentum(M0, m1, m2);
		const double Gamma0 = parentWidth;             // resonance peak width

		const double Gamma  = Gamma0 * (q / q0);

		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double  C  = M0 * Gamma;
		const double  A  = C * M / q;
		const double  B  = M0 * M0 - M * M;
		const Complex bw = (A / (B * B + C * C)) * Complex(B, C);
		// return ((M0 * Gamma0 * M / q) / (M0 * M0 - M * M - imag * M0 * Gamma);
		return bw;
}
};

void
rpwa::thrust_massDependence_f0980BreitWigner_amp(
	const ParVector<LorentzVector>& parentVec,
	const ParVector<LorentzVector>& daughter1Vec,
	const ParVector<LorentzVector>& daughter2Vec,
	ParVector<Complex>& result,
	double parentMass,
	double parentWidth)
{
	thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(parentVec.begin(), daughter1Vec.begin(), daughter2Vec.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(parentVec.end(),   daughter1Vec.end(),   daughter2Vec.end())),
			result.begin(),
			ThrustFunctor_massDependence_f0980BreitWigner_amp(parentMass, parentWidth));
}

struct ThrustFunctor_massDependence_rhoPrimeMassDep_amp
{
	const double M01;
	const double Gamma01;
	const double M02;
	const double Gamma02;
	ThrustFunctor_massDependence_rhoPrimeMassDep_amp(double M01, double Gamma01, double M02, double Gamma02):
		M01(M01), Gamma01(Gamma01), M02(M02), Gamma02(Gamma02) {}

	HOST_DEVICE
	Complex operator()(const LorentzVector& parentVec)
	{
		// get Breit-Wigner parameters
		const double M  = parentVec.M();                 // parent mass

		// const Complex bw1 = breitWigner(M, M01, Gamma01, L, q, q0);
		// const Complex bw2 = breitWigner(M, M02, Gamma02, L, q, q0);

		// A / (B - iA) = (A / (B^2 + A^2)) * (B + iA)
		const double  A1  = M01 * Gamma01;
		const double  B1  = M01 * M01 - M * M;
		const Complex bw1 = (A1 / (B1 * B1 + A1 * A1)) * Complex(B1, A1);
		// const Complex bw = (M0 * Gamma0) / (M0 * M0 - M * M - imag * M0 * Gamma0);

		// A / (B - iA) = (A / (B^2 + A^2)) * (B + iA)
		const double  A2  = M02 * Gamma02;
		const double  B2  = M02 * M02 - M * M;
		const Complex bw2 = (A2 / (B2 * B2 + A2 * A2)) * Complex(B2, A2);
		// const Complex bw = (M0 * Gamma0) / (M0 * M0 - M * M - imag * M0 * Gamma0);

		const Complex bw = (4.0 * bw1 - 3.0 * bw2) / 7.0;
		return bw;
}
};

void
rpwa::thrust_massDependence_rhoPrimeMassDep_amp(
	const ParVector<LorentzVector>& parentVec,
	ParVector<Complex>& result,
	double M01,
	double Gamma01,
	double M02,
	double Gamma02)
{
	thrust::transform(parentVec.begin(), parentVec.end(), result.begin(),
			ThrustFunctor_massDependence_rhoPrimeMassDep_amp(M01, Gamma01, M02, Gamma02));
}
