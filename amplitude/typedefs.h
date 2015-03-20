///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2014
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
//      Defines the types used for vectors, rotations etc.
//
//
// Author List:
//      Armin Gensler          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef RPWA_TYPEDEFS_H
#define RPWA_TYPEDEFS_H

#include <vector>

#ifdef USE_CUDA

	#include "thrustVector.h"

	#include "cudaComplex.hpp"
	#include "cudaVector3.hpp"
	#include "cudaLorentzVector.hpp"
	#include "cudaRotation.hpp"
	#include "cudaLorentzRotation.hpp"

	#define ParVector rpwa::ThrustVector

	namespace rpwa {

		typedef double Scalar;

		typedef CudaComplex<Scalar> Complex;
		typedef CudaVector3<Scalar> Vector3;
		typedef CudaLorentzVector<Scalar> LorentzVector;
		typedef CudaRotation<Scalar> Rotation;
		typedef CudaLorentzRotation<Scalar> LorentzRotation;

		template<typename T>
		inline ParVector<T> makeParVector(const std::vector<T>& v) {
			return rpwa::ThrustVector<T>(v);
		}

		template<typename T>
		inline void copyToParVector(ParVector<T>& parVec, const std::vector<T>& v) {
			parVec.fromStdVector(v);
		}

	}

#else

	#include <complex>

	#include "TVector3.h"
	#include "TLorentzVector.h"
	#include "TRotation.h"
	#include "TLorentzRotation.h"

	#define ParVector std::vector

	namespace rpwa {

		typedef double Scalar;

		typedef std::complex<Scalar> Complex;
		typedef TVector3 Vector3;
		typedef TLorentzVector LorentzVector;
		typedef TRotation Rotation;
		typedef TLorentzRotation LorentzRotation;

		template<typename T>
		inline ParVector<T> makeParVector(const std::vector<T>& v) {
			return v;
		}

		template<typename T>
		inline void copyToParVector(ParVector<T>& parVec, const std::vector<T>& v) {
			parVec = v;
		}

	}

#endif

#endif  // RPWA_TYPEDEFS_HPP
