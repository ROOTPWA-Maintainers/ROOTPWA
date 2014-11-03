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
//      A wrapper class to allow using thrust vectors in code not compiled with nvcc.
//
//
// Author List:
//      Armin Gensler          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#ifndef RPWA_THRUST_VECTOR_H
#define RPWA_THRUST_VECTOR_H

#include <vector>

#include "cuda.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

namespace rpwa {

	template<typename T>
	class ThrustVector {

	public:

		ThrustVector();
		ThrustVector(const ThrustVector<T>& v);
		ThrustVector(std::size_t s, T v = T());
		explicit ThrustVector(const std::vector<T>& v);
		~ThrustVector();

		ThrustVector<T>& operator = (const ThrustVector<T>& other);

		thrust::host_vector<T>& toRawThrustVector();
		const thrust::host_vector<T>& toRawThrustVector() const;

		typename thrust::host_vector<T>::iterator begin();
		typename thrust::host_vector<T>::const_iterator begin() const;

		typename thrust::host_vector<T>::iterator end();
		typename thrust::host_vector<T>::const_iterator end() const;

		T& operator[](int i);
		const T& operator[](int i) const;

		void clear();
		void resize(std::size_t s, T v = T());

		std::size_t size() const;

		void fromStdVector(const std::vector<T>& vec);
		void insertInStdVector(std::vector<T>& vec, typename std::vector<T>::iterator it) const;

	private:

		// make object instead of pointer? (should work with external implementation)
		thrust::host_vector<T>* vec;

	};

}

#endif
