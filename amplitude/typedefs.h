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

namespace rpwa {

	template<typename T>
	class RpwaVec {

	public:

		typedef typename std::vector<T>::iterator iterator;
		typedef typename std::vector<T>::const_iterator const_iterator;

		RpwaVec(): vec() {}
		RpwaVec(const RpwaVec<T>& v): vec(v.vec) {}
		RpwaVec(std::size_t s, T v = T()): vec(s, v) {}
		explicit RpwaVec(const std::vector<T>& v): vec(v) {}

		const std::vector<T>& toStdVector() const { return vec; }

		T& operator[](int i) { return vec[i]; }
		const T& operator[](int i) const { return vec[i]; }

		T& front() { return vec.front(); }
		const T& front() const { return vec.front(); }

		T& back() { return vec.back(); }
		const T& back() const { return vec.back(); }

		void push_back(const T& x) { vec.push_back(x); }

		void clear() { vec.clear(); }
		void resize(std::size_t s, T v = T()) { vec.resize(s, v); }

		std::size_t size() const { return vec.size(); }

		iterator begin() { return vec.begin(); }
		iterator end() { return vec.end(); }
		const_iterator begin() const { return vec.begin(); }
		const_iterator end() const { return vec.end(); }


	private:

		std::vector<T> vec;

	};

}

#ifdef USE_CUDA



	//#include "cuda.h"
	//#include <thrust/host_vector.h>
	//#include <thrust/device_vector.h>
	#include "thrustVector.h"

	#include "cudaComplex.hpp"
	#include "cudaVector3.hpp"
	#include "cudaLorentzVector.hpp"
	#include "cudaRotation.hpp"
	#include "cudaLorentzRotation.hpp"

	//#define ParVector thrust::device_vector
	#define ParVector rpwa::ThrustVector
	//#define ParVector RpwaVec

	namespace rpwa {

		typedef double Scalar;

		typedef CudaComplex<Scalar> Complex;
		typedef CudaVector3<Scalar> Vector3;
		typedef CudaLorentzVector<Scalar> LorentzVector;
		typedef CudaRotation<Scalar> Rotation;
		typedef CudaLorentzRotation<Scalar> LorentzRotation;

		/*template<typename T>
		inline ParVector<T> toParVector(const std::vector<T>& v) {
			return RpwaVec<T>(v);
		}*/

		template<typename T>
		inline ParVector<T> makeParVector(const std::vector<T>& v) {
			return rpwa::ThrustVector<T>(v);
		}

		template<typename T>
		inline void copyToParVector(ParVector<T>& parVec, const std::vector<T>& v) {
			parVec.fromStdVector(v);
		}

		/*template<typename T>
		inline ParVector<T> toParVector(const std::vector<T>& v) {
			thrust::host_vector h(v.size());
			std::copy(v.begin(), v.end(), h.begin());
			return thrust::device_vector(h);
		}*/

	}

#else

	#include <complex>

	#include "TVector3.h"
	#include "TLorentzVector.h"
	#include "TRotation.h"
	#include "TLorentzRotation.h"

	#include "cudaComplex.hpp"
	#include "cudaVector3.hpp"
	#include "cudaLorentzVector.hpp"
	#include "cudaRotation.hpp"
	#include "cudaLorentzRotation.hpp"

	//#define ParVector std::vector
	#define ParVector RpwaVec

	namespace rpwa {

		typedef double Scalar;

		typedef CudaComplex<Scalar> Complex;
		typedef CudaVector3<Scalar> Vector3;
		typedef CudaLorentzVector<Scalar> LorentzVector;
		typedef CudaRotation<Scalar> Rotation;
		typedef CudaLorentzRotation<Scalar> LorentzRotation;

		/*
		typedef std::complex<Scalar> Complex;
		typedef TVector3 Vector3;
		typedef TLorentzVector LorentzVector;
		typedef TRotation Rotation;
		typedef TLorentzRotation LorentzRotation;
		*/

		template<typename T>
		inline ParVector<T> makeParVector(const std::vector<T>& v) {
			return RpwaVec<T>(v);
		}

		template<typename T>
		inline void copyToParVector(ParVector<T>& parVec, const std::vector<T>& v) {
			parVec = RpwaVec<T>(v);
		}

		/*
		template<typename T>
		inline ParVector<T> toParVector(const std::vector<T>& v) {
			return v;
		}
		*/

	}

#endif

#endif  // RPWA_TYPEDEFS_HPP
