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
//      complex template class for CUDA; interface is compatible with std::complex
//
//      inspired by cbuchner1's complex class
//      see http://forums.nvidia.com/lofiversion/index.php?t73978.html
//      and http://forums.nvidia.com/index.php?showtopic=108787
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef COMPLEXTEST_CUH
#define COMPLEXTEST_CUH


#include <iostream>

#include <cuda_runtime.h>

#include "../utilities/cudaUtils.hpp"


namespace rpwa {

	namespace cuda {


		// possible storage type
		template<typename T>
		struct
		complexStruct {
			T x;
			T y;
		};


		template<typename storageT, typename T>
		struct complexTest {

			typedef storageT storage_type;
			typedef T  value_type;
  
			storageT _complex;
  
			// accessors for real and imaginary components
			inline HOST_DEVICE T& real() { return _complex.x; };
			inline HOST_DEVICE T& imag() { return _complex.y; };

			inline HOST_DEVICE const T& real() const { return _complex.x; };
			inline HOST_DEVICE const T& imag() const { return _complex.y; };
  
			// constructors
			inline
			HOST_DEVICE
			complexTest(const T& re = (T)0,
			            const T& im = (T)0)
			{
				_complex.x = re;
				_complex.y = im;
			}

			inline
			HOST_DEVICE
			complexTest(const complexTest& a)
			{
				_complex.x = a._complex.x;
				_complex.y = a._complex.y;
			}
  
			// assignment operators
			inline
			HOST_DEVICE
			complexTest<storageT, T>&
			operator =(const complexTest& a)
			{
				_complex.x = a._complex.x;
				_complex.y = a._complex.y;
				return *this;
			};

			inline
			HOST_DEVICE
			complexTest<storageT, T>&
			operator =(const T& re)
			{
				_complex.x = re;
				_complex.y = 0;
				return *this;
			};

			inline
			HOST_DEVICE
			complexTest<storageT, T>&
			operator +=(const complexTest<storageT, T>& a)
			{
				this->_complex.x += a._complex.x;
				this->_complex.y += a._complex.y;
				return *this;
			}

			// constants: 0, 1, i
			static HOST_DEVICE complexTest<storageT, T> zero() { return complexTest<storageT, T>(0, 0); }
			static HOST_DEVICE complexTest<storageT, T> one()  { return complexTest<storageT, T>(1, 0); }
			static HOST_DEVICE complexTest<storageT, T> i()    { return complexTest<storageT, T>(0, 1); }
		};


		// contruction function
		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		complexTest<storageT, T>
		makeComplex(const T re = (T)0,
		            const T im = (T)0)
		{
			return complexTest<storageT, T>(re, im);
		}


		// unary +
		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		complexTest<storageT, T>
		operator +(const complexTest<storageT, T>& a)
		{
			return a;
		}


		// summation
		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		complexTest<storageT, T>
		operator +(const complexTest<storageT, T>& a,
		           const complexTest<storageT, T>& b)
		{
			return complexTest<storageT, T>(a._complex.x + b._complex.x, a._complex.y  + b._complex.y);
		}

		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		complexTest<storageT, T>
		operator +(const T&                    a,
		           const complexTest<storageT, T>& b)
		{
			return complexTest<storageT, T>(a + b.value.x, b.value.y);
		}

		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		complexTest<storageT, T>
		operator +(const complexTest<storageT, T>& a,
		           const T&                    b)
		{
			return complexTest<storageT, T>(a.value.x + b, a.value.y);
		}

		
		// unary -
		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		complexTest<storageT, T>
		operator -(const complexTest<storageT, T>& a)
		{
			return complexTest<storageT, T>(-a._complex.x, -a._complex.y);
		}


		// subtraction
		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		complexTest<storageT, T>
		operator -(const complexTest<storageT, T>& a,
		           const complexTest<storageT, T>& b)
		{
			return complexTest<storageT, T>(a._complex.x - b._complex.x, a._complex.y  - b._complex.y);
		}

		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		complexTest<storageT, T>
		operator -(const T&                    a,
		           const complexTest<storageT, T>& b)
		{
			return complexTest<storageT, T>(a - b.value.x, -b.value.y);
		}

		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		complexTest<storageT, T>
		operator -(const complexTest<storageT, T>& a,
		           const T&                    b)
		{
			return complexTest<storageT, T>(a.value.x - b, a.value.y);
		}


		// multiplication
		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		complexTest<storageT, T> operator *(const complexTest<storageT, T>& a,
		                                    const complexTest<storageT, T>& b)
		{
			return complexTest<storageT, T>
				(a._complex.x * b._complex.x - a._complex.y * b._complex.y,
				 a._complex.y * b._complex.x + a._complex.x * b._complex.y);
		}

		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		complexTest<storageT, T>
		operator *(const T&                    a,
		           const complexTest<storageT, T>& b)
		{
			return complexTest<storageT, T>(a * b._complex.x, a * b._complex.y);
		}

		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		complexTest<storageT, T>
		operator *(const complexTest<storageT, T>& a,
		           const T&                    b)
		{
			return complexTest<storageT, T>(a._complex.x * b, a._complex.y * b);
		}


		// division
		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		complexTest<storageT, T>
		operator /(const complexTest<storageT, T>& a,
		           const complexTest<storageT, T>& b)
		{
			T denom = (b._complex.x * b._complex.x + b._complex.y * b._complex.y );
			return complexTest<storageT, T>
				((a._complex.x * b._complex.x + a._complex.y * b._complex.y ) / denom,
				 (a._complex.y * b._complex.x - a._complex.x * b._complex.y ) / denom);
		}

		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		complexTest<storageT, T>
		operator /(const T&                    a,
		           const complexTest<storageT, T>& b)
		{
			T denom = b._complex.x * b._complex.x + b._complex.y * b._complex.y;
			return complexTest<storageT, T>(a * b._complex.x / denom, -a * b._complex.y / denom);
		}

		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		complexTest<storageT, T>
		operator /(const complexTest<storageT, T>& a,
		           const T&                    b)
		{
			return complexTest<storageT, T>(a._complex.x / b, a._complex.y / b);
		}

		// real part
		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		T
		real(const complexTest<storageT, T>& a)
		{
			return a.real();
		}

		// imaginary part
		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		T
		imag(const complexTest<storageT, T>& a)
		{
			return a.imag();
		}

		// complex norm = Re^2 + Im^2
		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		T
		norm(const complexTest<storageT, T>& a)
		{
			return a._complex.x * a._complex.x + a._complex.y * a._complex.y;
		}

		// complex absolute value
		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		T
		abs(const complexTest<storageT, T>& a)
		{
			return sqrt(a.norm());
		}


		// complex conjugate
		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		complexTest<storageT, T>
		conj(const complexTest<storageT, T>& a)
		{
			return complexTest<storageT, T>(a._complex.x, -a._complex.y);
		}


		// complex phase angle
		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		T
		arg(const complexTest<storageT, T>& a)
		{
			return atan2(a._complex.y, a._complex.x);
		}


		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		bool
		operator ==(const complexTest<storageT, T>& a,
		            const complexTest<storageT, T>& b)
		{
			return ((a._complex.x == b._complex.x) and (a._complex.y == b._complex.y));
		}


		template<typename storageT, typename T>
		inline
		HOST_DEVICE
		bool
		operator !=(const complexTest<storageT, T>& a,
		            const complexTest<storageT, T>& b)
		{
			return not(a == b);
		}


		template<typename storageT, typename T>
		HOST
		std::ostream&
		operator <<(std::ostream& out,
		            const complexTest<storageT, T>& a)
		{
			out << "(" << a.real() << ", " << a.imag() << ")";
			return out;
		}


	}  // namespace cuda

}  // namespace rpwa


#endif  // COMPLEXTEST_CUH
