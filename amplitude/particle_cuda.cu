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
//      CUDA code for particle.cc
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
#include "particle_cuda.h"

using namespace rpwa;

struct ThrustFunctor_particle_setMomenta
{
	double mass2;
	ThrustFunctor_particle_setMomenta(double mass2): mass2(mass2) {}

	HOST_DEVICE
	LorentzVector operator()(const Vector3& mom)
	{
		return LorentzVector(mom, rpwa::sqrt(mom.Mag2() + mass2));
	}
};

void 
rpwa::thrust_particle_setMomenta(const ParVector<Vector3>& moms, 
				 ParVector<LorentzVector>& lzVecs, 
				 double mass2)
{
	thrust::transform(moms.begin(), moms.end(), lzVecs.begin(), 
			  ThrustFunctor_particle_setMomenta(mass2));
}

struct ThrustFunctor_particle_transformRot
{
	template<typename T>
	HOST_DEVICE
	void operator()(T t) const
	{
		const LorentzRotation& lzRot = thrust::get<0>(t);
		LorentzVector& lzVec = thrust::get<1>(t);
		lzVec.Transform(lzRot);
	}
};

void 
rpwa::thrust_particle_transformRot(const ParVector<LorentzRotation>& lzTrafos,
				   ParVector<LorentzVector>& lzVecs)
{
	thrust::for_each(
		thrust::make_zip_iterator(thrust::make_tuple(lzTrafos.begin(), lzVecs.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(lzTrafos.end(),   lzVecs.end())),
		ThrustFunctor_particle_transformRot());
}

struct ThrustFunctor_particle_transformBoost
{
	template<typename T>
	HOST_DEVICE
	void operator()(T t) const
	{
		const Vector3& boost = thrust::get<0>(t);
		LorentzVector& lzVec = thrust::get<1>(t);
		lzVec.Boost(boost);
	}
};

void 
rpwa::thrust_particle_transformBoost(const ParVector<Vector3>& boosts,
				     ParVector<LorentzVector>& lzVecs)
{
	thrust::for_each(
		thrust::make_zip_iterator(thrust::make_tuple(boosts.begin(), lzVecs.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(boosts.end(),   lzVecs.end())),
		ThrustFunctor_particle_transformBoost());
}
