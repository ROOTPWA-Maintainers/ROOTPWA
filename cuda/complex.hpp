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
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
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


#ifndef COMPLEX_HPP
#define COMPLEX_HPP


#include <iostream>

#include <cuda_runtime.h>


// boosting ALIGN from 8 to 16 doubles throughput for double scalars
// in copy and write operations, but leaves read-only bandwidth unchanged
// effective bandwitdth for floats drops to half
#ifndef ALIGN
#define ALIGN 8
#endif


// namespace rpwa {

//   namespace gpu {

// possible storage type
// without __align__ directive throughput drops to 50%
template<typename T>
struct __align__(ALIGN) complexStruct {
	T x;
	T y;
};


template<typename storageT, typename scalarT>
struct complex {

	typedef storageT storage_type;
	typedef scalarT  value_type;
  
	storageT _complex;
  
	// accessors for real and imaginary components
	inline __host__ __device__ scalarT& real() { return _complex.x; };
	inline __host__ __device__ scalarT& imag() { return _complex.y; };

	inline __host__ __device__ const scalarT& real() const { return _complex.x; };
	inline __host__ __device__ const scalarT& imag() const { return _complex.y; };
  
	// constructors
	inline __host__ __device__
	complex(const scalarT& re = (scalarT)0,
	        const scalarT& im = (scalarT)0)
	{
		_complex.x = re;
		_complex.y = im;
	}

	inline __host__ __device__
	complex(const complex& a)
	{
		_complex.x = a._complex.x;
		_complex.y = a._complex.y;
	}
  
	// assignment operators
	inline __host__ __device__
	complex<storageT, scalarT>&
	operator =(const complex& a)
	{
		_complex.x = a._complex.x;
		_complex.y = a._complex.y;
		return *this;
	};

	inline __host__ __device__
	complex<storageT, scalarT>&
	operator =(const scalarT& re)
	{
		_complex.x = re;
		_complex.y = 0;
		return *this;
	};

	inline __host__ __device__
	complex<storageT, scalarT>&
	operator +=(const complex<storageT, scalarT>& a)
	{
		this->_complex.x += a._complex.x;
		this->_complex.y += a._complex.y;
		return *this;
	}

	// constants: 0, 1, i
	static __host__ __device__ const complex<storageT, scalarT> zero() { return complex<storageT, scalarT>(0, 0); }
	static __host__ __device__ const complex<storageT, scalarT> one()  { return complex<storageT, scalarT>(1, 0); }
	static __host__ __device__ const complex<storageT, scalarT> i()    { return complex<storageT, scalarT>(0, 1); }
};


// contruction function
template<typename storageT, typename scalarT>
inline __host__ __device__
complex<storageT, scalarT>
makeComplex(const scalarT re = (scalarT)0,
            const scalarT im = (scalarT)0)
{
	return complex<storageT, scalarT>(re, im);
}


// summation
template<typename storageT, typename scalarT>
inline __host__ __device__
complex<storageT, scalarT>
operator +(const complex<storageT, scalarT>& a)
{
	return a;
}

template<typename storageT, typename scalarT>
inline __host__ __device__
complex<storageT, scalarT>
operator +(const complex<storageT, scalarT>& a,
           const complex<storageT, scalarT>& b)
{
	return complex<storageT, scalarT>(a._complex.x + b._complex.x, a._complex.y  + b._complex.y);
}

template<typename storageT, typename scalarT>
inline __host__ __device__
complex<storageT, scalarT>
operator +(const scalarT&                    a,
           const complex<storageT, scalarT>& b)
{
	return complex<storageT, scalarT>(a + b.value.x, b.value.y);
}

template<typename storageT, typename scalarT>
inline __host__ __device__
complex<storageT, scalarT>
operator +(const complex<storageT, scalarT>& a,
           const scalarT&                    b)
{
	return complex<storageT, scalarT>(a.value.x + b, a.value.y);
}


// subtraction
template<typename storageT, typename scalarT>
inline __host__ __device__
complex<storageT, scalarT>
operator -(const complex<storageT, scalarT>& a)
{
	return complex<storageT, scalarT>(-a._complex.x, -a._complex.y);
}

template<typename storageT, typename scalarT>
inline __host__ __device__
complex<storageT, scalarT>
operator -(const complex<storageT, scalarT>& a,
           const complex<storageT, scalarT>& b)
{
	return complex<storageT, scalarT>(a._complex.x - b._complex.x, a._complex.y  - b._complex.y);
}

template<typename storageT, typename scalarT>
inline __host__ __device__
complex<storageT, scalarT>
operator -(const scalarT&                    a,
           const complex<storageT, scalarT>& b)
{
	return complex<storageT, scalarT>(a - b.value.x, -b.value.y);
}

template<typename storageT, typename scalarT>
inline __host__ __device__
complex<storageT, scalarT>
operator -(const complex<storageT, scalarT>& a,
           const scalarT&                    b)
{
	return complex<storageT, scalarT>(a.value.x - b, a.value.y);
}


// multiplication
template<typename storageT, typename scalarT>
inline __host__ __device__
complex<storageT, scalarT> operator *(const complex<storageT, scalarT>& a,
                                      const complex<storageT, scalarT>& b)
{
	return complex<storageT, scalarT>(a._complex.x * b._complex.x - a._complex.y * b._complex.y,
	                                  a._complex.y * b._complex.x + a._complex.x * b._complex.y);
}

template<typename storageT, typename scalarT>
inline __host__ __device__
complex<storageT, scalarT>
operator *(const scalarT&                    a,
           const complex<storageT, scalarT>& b)
{
	return complex<storageT, scalarT>(a * b._complex.x, a * b._complex.y);
}

template<typename storageT, typename scalarT>
inline __host__ __device__
complex<storageT, scalarT>
operator *(const complex<storageT, scalarT>& a,
           const scalarT&                    b)
{
	return complex<storageT, scalarT>(a._complex.x * b, a._complex.y * b);
}


// division
template<typename storageT, typename scalarT>
inline __host__ __device__
complex<storageT, scalarT>
operator /(const complex<storageT, scalarT>& a,
           const complex<storageT, scalarT>& b)
{
	scalarT denom = (b._complex.x * b._complex.x + b._complex.y * b._complex.y );
	return complex<storageT, scalarT>((a._complex.x * b._complex.x + a._complex.y * b._complex.y ) / denom,
	                                  (a._complex.y * b._complex.x - a._complex.x * b._complex.y ) / denom);
}

template<typename storageT, typename scalarT>
inline __host__ __device__
complex<storageT, scalarT>
operator /(const scalarT&                    a,
           const complex<storageT, scalarT>& b)
{
	scalarT denom = b._complex.x * b._complex.x + b._complex.y * b._complex.y;
	return complex<storageT, scalarT>(a * b._complex.x / denom, -a * b._complex.y / denom);
}

template<typename storageT, typename scalarT>
inline __host__ __device__
complex<storageT, scalarT>
operator /(const complex<storageT, scalarT>& a,
           const scalarT&                    b)
{
	return complex<storageT, scalarT>(a._complex.x / b, a._complex.y / b);
}

// real part
template<typename storageT, typename scalarT>
inline __host__ __device__
scalarT
real(const complex<storageT, scalarT>& a)
{
	return a.real();
}

// imaginary part
template<typename storageT, typename scalarT>
inline __host__ __device__
scalarT
imag(const complex<storageT, scalarT>& a)
{
	return a.imag();
}

// complex norm = Re^2 + Im^2
template<typename storageT, typename scalarT>
inline __host__ __device__
scalarT
norm(const complex<storageT, scalarT>& a)
{
	return a._complex.x * a._complex.x + a._complex.y * a._complex.y;
}

// complex absolute value
template<typename storageT, typename scalarT>
inline __host__ __device__
scalarT
abs(const complex<storageT, scalarT>& a)
{
	return sqrt(a.norm());
}


// complex conjugate
template<typename storageT, typename scalarT>
inline __host__ __device__
complex<storageT, scalarT>
conj(const complex<storageT, scalarT>& a)
{
	return complex<storageT, scalarT>(a._complex.x, -a._complex.y);
}


// complex phase angle
template<typename storageT, typename scalarT>
inline __host__ __device__
scalarT
arg(const complex<storageT, scalarT>& a)
{
	return atan2(a._complex.y, a._complex.x);
}


template<typename storageT, typename scalarT>
inline __host__ __device__
bool
operator ==(const complex<storageT, scalarT>& a,
            const complex<storageT, scalarT>& b)
{
	return ((a._complex.x == b._complex.x) and (a._complex.y == b._complex.y));
}


template<typename storageT, typename scalarT>
inline __host__ __device__
bool
operator !=(const complex<storageT, scalarT>& a,
            const complex<storageT, scalarT>& b)
{
	return not(a == b);
}


template<typename storageT, typename scalarT>
__host__
std::ostream&
operator <<(std::ostream& out,
            const complex<storageT, scalarT>& a)
{
	out << "(" << a.real() << ", " << a.imag() << ")";
	return out;
}


//   }  // namespace gpu

// }  // namespace rpwa


#endif  // COMPLEX_HPP
