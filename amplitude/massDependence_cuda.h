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


#ifndef MASSDEPENDENCE_CUDA_H
#define MASSDEPENDENCE_CUDA_H

#include "typedefs.h"


namespace rpwa {

	void
	thrust_massDependence_flatRangeMassDependence_amp(
		const ParVector<LorentzVector>& parentVec,
		ParVector<Complex>& result,
		double parentMass,
		double parentWidth);

	void
	thrust_massDependence_relativisticBreitWigner_amp(
		const ParVector<LorentzVector>& parentVec,
		const ParVector<LorentzVector>& daughter1Vec,
		const ParVector<LorentzVector>& daughter2Vec,
		ParVector<Complex>& result,
		double parentMass,
		double parentWidth,
		unsigned int L);

	void
	thrust_massDependence_constWidthBreitWigner_amp(
		const ParVector<LorentzVector>& parentVec,
		ParVector<Complex>& result,
		double M0,
		double Gamma0);

	void
	thrust_massDependence_rhoBreitWigner_amp(
		const ParVector<LorentzVector>& parentVec,
		const ParVector<LorentzVector>& daughter1Vec,
		const ParVector<LorentzVector>& daughter2Vec,
		ParVector<Complex>& result,
		double parentMass,
		double parentWidth);

	void
	thrust_massDependence_f0980BreitWigner_amp(
		const ParVector<LorentzVector>& parentVec,
		const ParVector<LorentzVector>& daughter1Vec,
		const ParVector<LorentzVector>& daughter2Vec,
		ParVector<Complex>& result,
		double parentMass,
		double parentWidth);

	void
	thrust_massDependence_piPiSWaveAuMorganPenningtonVes_amp(
		const ParVector<LorentzVector>& parentVec,
		const ParVector<Complex>& ampM,
		ParVector<Complex>& result,
		double _piChargedMass,
		double f0Mass,
		double f0Width,
		Complex coupling);

	void
	thrust_massDependence_rhoPrimeMassDep_amp(
		const ParVector<LorentzVector>& parentVec,
		ParVector<Complex>& result,
		double M01,
		double Gamma01,
		double M02,
		double Gamma02);

}  // namespace rpwa


#endif  // ISOBARCANONICALAMPLITUDE_CUDA_H
