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


#ifndef PARTICLE_CUDA_H
#define PARTICLE_CUDA_H

#include "typedefs.h"


namespace rpwa {

	void 
	thrust_particle_setMomenta(const ParVector<Vector3>& moms, 
				   ParVector<LorentzVector>& lzVecs, 
				   double mass2);

	void 
	thrust_particle_transformRot(const ParVector<LorentzRotation>& lzTrafos,
				     ParVector<LorentzVector>& lzVecs);

	void 
	thrust_particle_transformBoost(const ParVector<Vector3>& boosts,
				       ParVector<LorentzVector>& lzVecs);


}  // namespace rpwa


#endif  // PARTICLE_CUDA_H
