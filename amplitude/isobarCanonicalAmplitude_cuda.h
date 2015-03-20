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


#ifndef ISOBARCANONICALAMPLITUDE_CUDA_H
#define ISOBARCANONICALAMPLITUDE_CUDA_H

#include "typedefs.h"


namespace rpwa {

	void 
	thrust_isobarCanonicalAmplitude_transformDaughters(
		const ParVector<LorentzVector>& parentVec,
		ParVector<Vector3>& result);

	void
	thrust_isobarCanonicalAmplitude_twoBodyDecayAmplitude_1(
		const ParVector<LorentzVector>& lzVec,
		ParVector<Complex>& ampSum,
		int L,
		int mL,
		double LSClebsch);

	void
	thrust_isobarCanonicalAmplitude_twoBodyDecayAmplitude_2(
		const ParVector<LorentzVector>& lzVec,
		const ParVector<Complex>& bw,
		ParVector<Complex>& ampProd,
		int L,
		double norm,
		double ssClebsch);

}  // namespace rpwa


#endif  // ISOBARCANONICALAMPLITUDE_CUDA_H
