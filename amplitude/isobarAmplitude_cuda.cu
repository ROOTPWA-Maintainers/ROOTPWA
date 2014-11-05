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
//      CUDA code for isobarAmplitude.cc
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
#include "mathUtils.hpp"
#include "isobarAmplitude_cuda.h"

using namespace rpwa;


struct ThrustFunctor_isobarAmplitude_amplitude
{
	Complex smyTermFactor;
	ThrustFunctor_isobarAmplitude_amplitude(Complex smyTermFactor):
		smyTermFactor(smyTermFactor) {}

	template<typename T>
	HOST_DEVICE
	void operator()(T t)
	{
		const Complex& permAmp = thrust::get<0>(t);
		Complex& amp = thrust::get<1>(t);
		amp += permAmp * smyTermFactor;
	}
};

void
rpwa::thrust_isobarAmplitude_amplitude(
	const ParVector<Complex>& permAmp,
	ParVector<Complex>& amp,
	Complex smyTermFactor)
{
	thrust::for_each(
			thrust::make_zip_iterator(thrust::make_tuple(permAmp.begin(), amp.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(permAmp.end(),   amp.end())),
			ThrustFunctor_isobarAmplitude_amplitude(smyTermFactor));
}

struct ThrustFunctor_isobarAmplitude_gjTransform
{
	double piHalf;
	ThrustFunctor_isobarAmplitude_gjTransform(double piHalf):
		piHalf(piHalf) {}

	HOST_DEVICE
	LorentzRotation operator()(LorentzVector beam, LorentzVector X)
	{
		const Vector3 yGjAxis = beam.Vect().Cross(X.Vect());  // y-axis of Gottfried-Jackson frame
		// rotate so that yGjAxis becomes parallel to y-axis and beam momentum ends up in (x, z)-plane
		Rotation rot1;
		rot1.RotateZ(piHalf - yGjAxis.Phi());
		rot1.RotateX(yGjAxis.Theta() - piHalf);
		beam *= rot1;
		X    *= rot1;
		// boost to X RF
		LorentzRotation boost;
		boost.Boost(-X.BoostVector());
		beam *= boost;
		// rotate about yGjAxis so that beam momentum is along z-axis
		Rotation rot2;
		rot2.RotateY(-signum(beam.X()) * beam.Theta());
		// construct total transformation
		LorentzRotation gjTransform(rot1);
		gjTransform.Transform(boost);
		gjTransform.Transform(rot2);

		return gjTransform;
	}
};

void
rpwa::thrust_isobarAmplitude_gjTransform(
	const ParVector<LorentzVector>& beamLv,
	const ParVector<LorentzVector>& XLv,
	ParVector<LorentzRotation>& result)
{
	thrust::transform(beamLv.begin(), beamLv.end(), XLv.begin(), result.begin(),
			ThrustFunctor_isobarAmplitude_gjTransform(piHalf));
}

struct ThrustFunctor_isobarAmplitude_parallelLorentzRotationInvert
{
	HOST_DEVICE
	void operator()(LorentzRotation& lzRot)
	{
		lzRot.Invert();
	}
};

void
rpwa::thrust_isobarAmplitude_parallelLorentzRotationInvert(
	ParVector<LorentzRotation>& lzRot)
{
	thrust::for_each(lzRot.begin(), lzRot.end(),
			ThrustFunctor_isobarAmplitude_parallelLorentzRotationInvert());
}

struct ThrustFunctor_isobarAmplitude_spaceInvertDecay
{
	HOST_DEVICE
	void operator()(LorentzVector& lzVec)
	{
		lzVec.SetXYZT(-lzVec.X(), -lzVec.Y(), -lzVec.Z(), lzVec.E());
	}
};

void
rpwa::thrust_isobarAmplitude_spaceInvertDecay(
	ParVector<LorentzVector>& lzVec)
{
	thrust::for_each(lzVec.begin(), lzVec.end(),
			ThrustFunctor_isobarAmplitude_spaceInvertDecay());
}

struct ThrustFunctor_isobarAmplitude_reflectDecay
{
	HOST_DEVICE
	void operator()(LorentzVector& lzVec)
	{
		lzVec.SetXYZT(lzVec.X(), -lzVec.Y(), lzVec.Z(), lzVec.E());
	}
};

void
rpwa::thrust_isobarAmplitude_reflectDecay(
	ParVector<LorentzVector>& lzVec)
{
	thrust::for_each(lzVec.begin(), lzVec.end(),
			ThrustFunctor_isobarAmplitude_reflectDecay());
}

struct ThrustFunctor_isobarAmplitude_twoBodyDecayAmplitudeSum
{
	template<typename T>
	HOST_DEVICE
	void operator()(T t)
	{
		const Complex& parentAmp = thrust::get<0>(t);
		const Complex& daughter1Amp = thrust::get<1>(t);
		const Complex& daughter2Amp = thrust::get<2>(t);
		Complex& ampSum = thrust::get<3>(t);
		ampSum += parentAmp * daughter1Amp * daughter2Amp;
	}
};

void
rpwa::thrust_isobarAmplitude_twoBodyDecayAmplitudeSum(
	const ParVector<Complex>& parentAmp,
	const ParVector<Complex>& daughter1Amp,
	const ParVector<Complex>& daughter2Amp,
	ParVector<Complex>& ampSum)
{
	thrust::for_each(
			thrust::make_zip_iterator(thrust::make_tuple(parentAmp.begin(), daughter1Amp.begin(), daughter2Amp.begin(), ampSum.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(parentAmp.end(),   daughter1Amp.end(),   daughter2Amp.end(),   ampSum.begin())),
			ThrustFunctor_isobarAmplitude_twoBodyDecayAmplitudeSum());
}
