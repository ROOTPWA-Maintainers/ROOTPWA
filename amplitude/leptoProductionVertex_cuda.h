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


#ifndef LEPTOPRODUCTIONVERTEX_CUDA_H
#define LEPTOPRODUCTIONVERTEX_CUDA_H

#include "typedefs.h"


namespace rpwa {

	void 
	thrust_leptoProductionVertex_productionAmps(
		const ParVector<LorentzVector>& targetVec,
		const ParVector<LorentzVector>& beamVec,
		const ParVector<LorentzVector>& scatteredLeptonVec,
		const ParVector<LorentzVector>& virtPhotonVec,
		const ParVector<LorentzVector>& XParticleVec,
		const ParVector<double>& epsilons,
		const ParVector<double>& deltas,
		ParVector<Complex>& result,
		double longPol);

	void
	thrust_leptoProductionVertex_Q2(
		const ParVector<LorentzVector>& virtPhotonVec,
		ParVector<double>& result);

	void
	thrust_leptoProductionVertex_nu(
		const ParVector<LorentzVector>& targetVec,
		const ParVector<LorentzVector>& virtPhotonVec,
		ParVector<double>& result,
		double targetMass);

	void
	thrust_leptoProductionVertex_y(
		const ParVector<LorentzVector>& targetVec,
		const ParVector<LorentzVector>& virtPhotonVec,
		const ParVector<LorentzVector>& beamVec,
		ParVector<double>& result);

	void
	thrust_leptoProductionVertex_epsilon(
		const ParVector<double>& Q2,
		const ParVector<double>& xBj,
		const ParVector<double>& y,
		ParVector<double>& result,
		double targetMass2,
		double beamLeptonMass2);

	void
	thrust_leptoProductionVertex_delta(
		const ParVector<double>& Q2,
		const ParVector<double>& epsilons,
		ParVector<double>& result,
		double beamMass2);

	void
	thrust_leptoProductionVertex_xBj(
		const ParVector<LorentzVector>& targetVec,
		const ParVector<LorentzVector>& virtPhotonVec,
		const ParVector<double>& Q2,
		ParVector<double>& result);

	void
	thrust_leptoProductionVertex_s(
		const ParVector<LorentzVector>& targetVec,
		const ParVector<LorentzVector>& virtPhotonVec,
		ParVector<double>& result);

	void
	thrust_leptoProductionVertex_W(
		const ParVector<double>& s,
		ParVector<double>& result);

	void
	thrust_leptoProductionVertex_revertMomenta(
		const ParVector<LorentzVector>& scatteredLeptonVec,
		ParVector<LorentzVector>& result);

}  // namespace rpwa


#endif  // LEPTOPRODUCTIONVERTEX_CUDA_H
