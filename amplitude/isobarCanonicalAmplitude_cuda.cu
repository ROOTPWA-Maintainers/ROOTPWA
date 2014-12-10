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
//      CUDA code for isobarCanonicalAmplitude.cc
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
#include "isobarCanonicalAmplitude_cuda.h"

using namespace rpwa;
using namespace std;

struct ThrustFunctor_isobarCanonicalAmplitude_transformDaughters
{
	HOST_DEVICE
	Vector3 operator()(const LorentzVector& parentVec)
	{
		return - parentVec.BoostVector();
	}
};

void
rpwa::thrust_isobarCanonicalAmplitude_transformDaughters(
	const ParVector<LorentzVector>& parentVec,
	ParVector<Vector3>& result)
{
	thrust::transform(parentVec.begin(), parentVec.end(), result.begin(),
			ThrustFunctor_isobarCanonicalAmplitude_transformDaughters());
}

struct ThrustFunctor_isobarAmplitude_twoBodyDecayAmplitude_1
{
	double LSClebsch;
	int L;
	int mL;
	const dFunctionCached<double>::cacheType* dFunctionCudaCache;
	ThrustFunctor_isobarAmplitude_twoBodyDecayAmplitude_1(double LSClebsch, int L, int mL,
			const dFunctionCached<double>::cacheType* dFunctionCudaCache):
		LSClebsch(LSClebsch), L(L), mL(mL), dFunctionCudaCache(dFunctionCudaCache) {}

	template<typename T>
	HOST_DEVICE
	void operator()(T t)
	{
		const LorentzVector& lzVec = thrust::get<0>(t);
		Complex& ampSum = thrust::get<1>(t);
		double phi = lzVec.Phi(); // use daughter1 as analyzer
		double theta = lzVec.Theta();
		ampSum += LSClebsch * sphericalHarmonic<Complex>(L, mL, theta, phi, false, dFunctionCudaCache);
	}
};

void
rpwa::thrust_isobarCanonicalAmplitude_twoBodyDecayAmplitude_1(
	const ParVector<LorentzVector>& lzVec,
	ParVector<Complex>& ampSum,
	int L,
	int mL,
	double LSClebsch)
{
	thrust::for_each(
			thrust::make_zip_iterator(thrust::make_tuple(lzVec.begin(), ampSum.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(lzVec.end(),   ampSum.end())),
			ThrustFunctor_isobarAmplitude_twoBodyDecayAmplitude_1(LSClebsch, L, mL, getDFunctionCudaCache<double>()));
}

struct ThrustFunctor_isobarCanonicalAmplitude_twoBodyDecayAmplitude_2
{
	int L;
	double norm;
	double ssClebsch;
	ThrustFunctor_isobarCanonicalAmplitude_twoBodyDecayAmplitude_2(int L, double norm, double ssClebsch):
		L(L), norm(norm), ssClebsch(ssClebsch) {}

	template<typename T>
	HOST_DEVICE
	void operator()(T t)
	{
		const LorentzVector& lzVec = thrust::get<0>(t);
		const Complex& bw = thrust::get<1>(t);
		Complex& ampProd = thrust::get<2>(t);

		// calulate barrier factor
		const double q  = lzVec.Vect().Mag();
		const double bf = barrierFactor(L, q, false);

		// calculate decay amplitude
		ampProd *= norm * ssClebsch * bf * bw;

	}
};

void
rpwa::thrust_isobarCanonicalAmplitude_twoBodyDecayAmplitude_2(
	const ParVector<LorentzVector>& lzVec,
	const ParVector<Complex>& bw,
	ParVector<Complex>& ampProd,
	int L,
	double norm,
	double ssClebsch)
{
	thrust::for_each(
			thrust::make_zip_iterator(thrust::make_tuple(lzVec.begin(), bw.begin(), ampProd.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(lzVec.end(),   bw.end(),   ampProd.end())),
			ThrustFunctor_isobarCanonicalAmplitude_twoBodyDecayAmplitude_2(L, norm, ssClebsch));
}
