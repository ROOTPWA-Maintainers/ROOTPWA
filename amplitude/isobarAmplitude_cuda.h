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


#ifndef ISOBARAMPLITUDE_CUDA_H
#define ISOBARAMPLITUDE_CUDA_H

#include "typedefs.h"


namespace rpwa {

	void 
	thrust_isobarAmplitude_amplitude(
		const ParVector<Complex>& permAmp,
		ParVector<Complex>& amp,
		Complex smyTermFactor);

	void
	thrust_isobarAmplitude_gjTransform(
		const ParVector<LorentzVector>& beamLv,
		const ParVector<LorentzVector>& XLv,
		ParVector<LorentzRotation>& result);

	void
	thrust_isobarAmplitude_parallelLorentzRotationInvert(
		ParVector<LorentzRotation>& lzRot);

	void
	thrust_isobarAmplitude_spaceInvertDecay(
		ParVector<LorentzVector>& lzVec);

	void
	thrust_isobarAmplitude_reflectDecay(
		ParVector<LorentzVector>& lzVec);

	void
	thrust_isobarAmplitude_twoBodyDecayAmplitudeSum(
		const ParVector<Complex>& parentAmp,
		const ParVector<Complex>& daughter1Amp,
		const ParVector<Complex>& daughter2Amp,
		ParVector<Complex>& ampSum);

}  // namespace rpwa


#endif  // ISOBARAMPLITUDE_CUDA_H
