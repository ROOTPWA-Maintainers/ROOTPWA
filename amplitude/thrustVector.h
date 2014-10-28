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

#include "cuda.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

namespace rpwa {

	template<typename T>
	class ThrustVector {

	public:

		typedef typename thrust::host_vector<T>::iterator iterator;
		typedef typename thrust::host_vector<T>::const_iterator const_iterator;

		ThrustVector();
		ThrustVector(const ThrustVector<T>& v);
		ThrustVector(std::size_t s, T v = T());
		explicit ThrustVector(const std::vector<T>& v);
		~ThrustVector();

		std::vector<T> toStdVector() const;
		thrust::host_vector<T>& toRawThrustVector();

		T& operator[](int i);
		const T& operator[](int i) const;

		T& front();
		const T& front() const;

		T& back();
		const T& back() const;

		void push_back(const T& x);

		void clear();
		void resize(std::size_t s, T v = T());

		std::size_t size() const;

		iterator begin();
		iterator end();
		const_iterator begin() const;
		const_iterator end() const;


	private:

		thrust::host_vector<T>* vec;

	};

}

#endif
