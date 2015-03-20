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
//      CUDA code for isobarHelicityAmplitude.cc
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
#include "isobarHelicityAmplitude_cuda.h"

using namespace rpwa;
using namespace std;

struct ThrustFunctor_isobarHelicityAmplitude_hfTransform
{
	const double piHalf;
	ThrustFunctor_isobarHelicityAmplitude_hfTransform(double piHalf):
		piHalf(piHalf) {}

	HOST_DEVICE
	LorentzRotation operator()(const LorentzVector& daughterLv)
	{
		LorentzVector daughter = daughterLv;
		const Vector3 zAxisParent(0, 0, 1);  // take z-axis as defined in parent frame
		const Vector3 yHfAxis = zAxisParent.Cross(daughter.Vect());  // y-axis of helicity frame
		// rotate so that yHfAxis becomes parallel to y-axis and zHfAxis ends up in (x, z)-plane
		Rotation rot1;
		rot1.RotateZ(piHalf - yHfAxis.Phi());
		rot1.RotateX(yHfAxis.Theta() - piHalf);
		daughter *= rot1;
		// rotate about yHfAxis so that daughter momentum is along z-axis
		Rotation rot2;
		rot2.RotateY(-signum(daughter.X()) * daughter.Theta());
		daughter *= rot2;
		// boost to daughter RF
		rot1.Transform(rot2);
		LorentzRotation hfTransform(rot1);
		hfTransform.Boost(-daughter.BoostVector());
		return hfTransform;
	}
};

void
rpwa::thrust_isobarHelicityAmplitude_hfTransform(
	const ParVector<LorentzVector>& daughterLv,
	ParVector<LorentzRotation>& result)
{
	thrust::transform(daughterLv.begin(), daughterLv.end(), result.begin(),
			ThrustFunctor_isobarHelicityAmplitude_hfTransform(piHalf));
}

struct ThrustFunctor_isobarHelicityAmplitude_twoBodyDecayAmplitude_1
{
	const int J;
	const int Lambda;
	const int lambda;
	const int P;
	const int refl;
	const dFunctionCached<double>::cacheType* dFunctionCudaCache;
	ThrustFunctor_isobarHelicityAmplitude_twoBodyDecayAmplitude_1(int J, int Lambda, int lambda, int P, int refl,
			const dFunctionCached<double>::cacheType* dFunctionCudaCache):
		J(J), Lambda(Lambda), lambda(lambda), P(P), refl(refl), dFunctionCudaCache(dFunctionCudaCache) {}

	HOST_DEVICE
	Complex operator()(const LorentzVector& daughterLv)
	{
		double phi = daughterLv.Phi(); // use daughter1 as analyzer
		double theta = daughterLv.Theta();
		return DFunctionReflConj<Complex>(J, Lambda, lambda, P, refl, phi, theta, 0, false, dFunctionCudaCache);
	}
};

void
rpwa::thrust_isobarHelicityAmplitude_twoBodyDecayAmplitude_1(
	const ParVector<LorentzVector>& lzVec,
	ParVector<Complex>& result,
	int J,
	int Lambda,
	int lambda,
	int P,
	int refl)
{
	thrust::transform(lzVec.begin(), lzVec.end(), result.begin(),
			ThrustFunctor_isobarHelicityAmplitude_twoBodyDecayAmplitude_1(J, Lambda, lambda, P, refl,
					getDFunctionCudaCache<double>()));
}

struct ThrustFunctor_isobarHelicityAmplitude_twoBodyDecayAmplitude_2
{
	const int J;
	const int Lambda;
	const int lambda;
	const int P;
	const int refl;
	const dFunctionCached<double>::cacheType* dFunctionCudaCache;
	ThrustFunctor_isobarHelicityAmplitude_twoBodyDecayAmplitude_2(int J, int Lambda, int lambda, int P, int refl,
			const dFunctionCached<double>::cacheType* dFunctionCudaCache):
		J(J), Lambda(Lambda), lambda(lambda), P(P), refl(refl), dFunctionCudaCache(dFunctionCudaCache) {}

	HOST_DEVICE
	Complex operator()(const LorentzVector& daughterLv)
	{
		double phi = daughterLv.Phi(); // use daughter1 as analyzer
		double theta = daughterLv.Theta();
		return DFunctionConj<Complex>(J, Lambda, lambda, phi, theta, 0, false, dFunctionCudaCache);
	}
};

void
rpwa::thrust_isobarHelicityAmplitude_twoBodyDecayAmplitude_2(
	const ParVector<LorentzVector>& lzVec,
	ParVector<Complex>& result,
	int J,
	int Lambda,
	int lambda,
	int P,
	int refl)
{
	thrust::transform(lzVec.begin(), lzVec.end(), result.begin(),
			ThrustFunctor_isobarHelicityAmplitude_twoBodyDecayAmplitude_2(J, Lambda, lambda, P, refl,
					getDFunctionCudaCache<double>()));
}

struct ThrustFunctor_isobarHelicityAmplitude_twoBodyDecayAmplitude_3
{
	const int L;
	const double norm;
	const double lsClebsch;
	const double ssClebsch;
	ThrustFunctor_isobarHelicityAmplitude_twoBodyDecayAmplitude_3(int L, double norm, double lsClebsch, double ssClebsch):
		L(L), norm(norm), lsClebsch(lsClebsch), ssClebsch(ssClebsch) {}

	template<typename T>
	HOST_DEVICE
	Complex operator()(T t)
	{
		const LorentzVector& daughterLv = thrust::get<0>(t);
		const Complex& DFunc = thrust::get<1>(t);
		const Complex& bw = thrust::get<2>(t);

		// calulate barrier factor
		const double q  = daughterLv.Vect().Mag();
		const double bf = barrierFactor(L, q, false);
		return norm * DFunc * lsClebsch * ssClebsch * bf * bw;
	}
};

void
rpwa::thrust_isobarHelicityAmplitude_twoBodyDecayAmplitude_3(
	const ParVector<LorentzVector>& daughterLv,
	const ParVector<Complex>& DFunc,
	const ParVector<Complex>& bw,
	ParVector<Complex>& result,
	int L,
	double  norm,
	double lsClebsch,
	double ssClebsch)
{
	thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(daughterLv.begin(), DFunc.begin(), bw.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(daughterLv.end(),   DFunc.end(),   bw.end())),
			result.begin(),
			ThrustFunctor_isobarHelicityAmplitude_twoBodyDecayAmplitude_3(L, norm, lsClebsch, ssClebsch));
}
