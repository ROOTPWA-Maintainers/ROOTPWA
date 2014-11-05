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
//      CUDA code for isobarDecayVertex.cc
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
#include "isobarDecayVertex_cuda.h"

using namespace rpwa;

void 
rpwa::thrust_isobarDecayVertex_calcParentLzVecs(const ParVector<LorentzVector>& daughter1, 
						const ParVector<LorentzVector>& daughter2, 
						ParVector<LorentzVector>& parent)
{
	thrust::transform(daughter1.begin(), daughter1.end(), daughter2.begin(), 
			  parent.begin(), thrust::plus<LorentzVector>());
}
