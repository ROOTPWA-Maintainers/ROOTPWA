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
//      Defines the types used for vectors, rotations etc.
//
//
// Author List:
//      Armin Gensler          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef RPWA_TYPEDEFS_HPP
#define RPWA_TYPEDEFS_HPP

#include <complex>
#include <vector>

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRotation.h"
#include "TLorentzRotation.h"

#include "Complex.hpp"
#include "Vector3.hpp"
#include "LorentzVector.hpp"
#include "Rotation.hpp"
#include "LorentzRotation.hpp"

namespace rpwa {

	template<typename T>
	class RpwaVec {

	public:

		typedef typename std::vector<T>::iterator iterator;
		typedef typename std::vector<T>::const_iterator const_iterator;

		RpwaVec(): vec() {}
		RpwaVec(const RpwaVec<T>& v): vec(v.vec) {}
		RpwaVec(size_t s, T v = T()): vec(s, v) {}
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
		void resize(size_t s, T v = T()) { vec.resize(s, v); }

		size_t size() const { return vec.size(); }

		iterator begin() { return vec.begin(); }
		iterator end() { return vec.end(); }
		const_iterator begin() const { return vec.begin(); }
		const_iterator end() const { return vec.end(); }


	private:

		std::vector<T> vec;

	};

}

//#define ParVector std::vector
#define ParVector RpwaVec
//#define ParVector thrust::device_vector

namespace rpwa {

	// utility function to create single-element vectors
	template<typename T>
	inline
	ParVector<T>
	make_vector_1(const T& element) {
		ParVector<T> vec(1);
		vec[0] = element;
		return vec;
	}


	typedef double Scalar;

	typedef RpwaComplex<Scalar> Complex;

	typedef RpwaVector3<Scalar> Vector3;

	typedef RpwaLorentzVector<Scalar> LorentzVector;

	typedef RpwaRotation<Scalar> Rotation;

	typedef RpwaLorentzRotation<Scalar> LorentzRotation;

	/*
	typedef double Scalar;

	typedef std::complex<Scalar> Complex;
  
	typedef TVector3 Vector3;
	
	typedef TLorentzVector LorentzVector;
	
	typedef TRotation Rotation;
	
	typedef TLorentzRotation LorentzRotation;
	*/

}  // namespace rpwa


#endif  // RPWA_TYPEDEFS_HPP
