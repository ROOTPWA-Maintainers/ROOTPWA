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
//      RpwaComplex number class fully compatible to std::RpwaComplex
//      only transcedental functions are missing
//
//
// Author List:
//      Philipp Meyer          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef RPWA_COMPLEX_HPP
#define RPWA_COMPLEX_HPP

#include <cmath>
#include <sstream>

#include "cudaUtils.hpp"

namespace rpwa {

	template<typename T> 
	class
	ALIGN(2 * sizeof(double))
	RpwaComplex {

	public:

		typedef T value_type;

		T _real;
		T _imag;

		HOST_DEVICE
		RpwaComplex(const T& real = T(), const T& imag = T()): _real(real), _imag(imag) { }
			
		template<typename U>
		HOST_DEVICE
		RpwaComplex(const RpwaComplex<U>& z): _real(z._real), _imag(z._imag) { }

		template<typename U>
		HOST std::complex<U> toStdComplex() const { return std::complex<U>(_real, _imag); }

		//template<typename U>
		HOST operator std::complex<T>() const { return std::complex<T>(_real, _imag); }

		//////////////////////////////////////////////////////////////////////////
		// accessors
		HOST_DEVICE T&       real()       { return _real; }
		HOST_DEVICE const T& real() const { return _real; }
		HOST_DEVICE T&       imag()       { return _imag; }
		HOST_DEVICE const T& imag() const { return _imag; }

		HOST_DEVICE void real(const T& val) { _real = val; }
		HOST_DEVICE void imag(const T& val) { _imag = val; }

		//////////////////////////////////////////////////////////////////////////
		// assignment operator for scalars
		HOST_DEVICE RpwaComplex<T>& operator =(const T& t)
		{
			_real = t;
			_imag = T();
			return *this;
		} 
		
		HOST_DEVICE RpwaComplex<T>& operator +=(const T& t)
		{
			_real += t;
			return *this;
		}

		HOST_DEVICE RpwaComplex<T>& operator -=(const T& t)
		{
			_real -= t;
			return *this;
		}

		HOST_DEVICE RpwaComplex<T>& operator *=(const T& t)
		{
			_real *= t;
			_imag *= t;
			return *this;
		}

		HOST_DEVICE RpwaComplex<T>& operator /=(const T& t)
		{
			_real /= t;
			_imag /= t;
			return *this;
		}

		//////////////////////////////////////////////////////////////////////////
		// assignment operator for RpwaComplex numbers
		template<typename U>
		HOST_DEVICE RpwaComplex<T>& operator =(const RpwaComplex<U>& z)
		{
			_real = z.real();
			_imag = z.imag();
			return *this;
		} 
	
		template<typename U>
		HOST_DEVICE RpwaComplex<T>& operator +=(const RpwaComplex<U>& z)
		{
			_real += z.real();
			_imag += z.imag();
			return *this;
		}

		template<typename U>
		HOST_DEVICE RpwaComplex<T>& operator -=(const RpwaComplex<U>& z)
		{
			_real -= z.real();
			_imag -= z.imag();
			return *this;
		}

		template<typename U>
		HOST_DEVICE RpwaComplex<T>& operator *=(const RpwaComplex<U>& z)
		{
			const T newReal = _real * z.real() - _imag * z.imag();
			_imag = _real * z.imag() + _imag * z.real();
			_real = newReal;
			return *this;
		}

		template<typename U>
		HOST_DEVICE RpwaComplex<T>& operator /=(const RpwaComplex<U>& z)
		{
			const T newReal = _real * z.real() + _imag * z.imag();
			const T norm    = sqrt(z.real() * z.real() + z.imag() * z.imag());
			_imag = (_imag * z.real() - _real * z.imag()) / norm;
			_real = newReal / norm;
			return *this;
		}

	};

	
	//////////////////////////////////////////////////////////////////////////
	// operators with scalars
	template<typename T>
	inline
	HOST_DEVICE
	RpwaComplex<T>
	operator +(const RpwaComplex<T>& z,
		    const T&          t)
	{
		RpwaComplex<T> result = z;
		result += t;
		return result;
	}

	template<typename T>
	inline
	HOST_DEVICE
	RpwaComplex<T>
	operator +(const T&          t,
		    const RpwaComplex<T>& z)
	{
		RpwaComplex<T> result = z;
		result += t;
		return result;
	}


	template<typename T>
	inline
	HOST_DEVICE
	RpwaComplex<T>
	operator -(const RpwaComplex<T>& z,
		    const T&          t)
	{
		RpwaComplex<T> result = z;
		result -= t;
		return result;
	}

	template<typename T>
	inline
	HOST_DEVICE
	RpwaComplex<T>
	operator -(const T&          t,
		    const RpwaComplex<T>& z)
	{
		RpwaComplex<T> result(t, -z.imag());
		result -= z.real();
		return result;
	}


	template<typename T>
	inline
	HOST_DEVICE
	RpwaComplex<T>
	operator *(const RpwaComplex<T>& z,
		    const T&          t)
	{
		RpwaComplex<T> result = z;
		result *= t;
		return result;
	}

	template<typename T>
	inline
	HOST_DEVICE
	RpwaComplex<T>
	operator *(const T&          t,
		    const RpwaComplex<T>& z)
	{
		RpwaComplex<T> result = z;
		result *= t;
		return result;
	}


	template<typename T>
	inline
	HOST_DEVICE
	RpwaComplex<T>
	operator /(const RpwaComplex<T>& z,
		    const T&          t)
	{
		RpwaComplex<T> result = z;
		result /= t;
		return result;
	}

	template<typename T>
	inline
	HOST_DEVICE
	RpwaComplex<T>
	operator /(const T&          t,
		    const RpwaComplex<T>& z)
	{
		RpwaComplex<T> result = t;
		result /= z;
		return result;
	}


	template<typename T>
	inline
	HOST_DEVICE
	bool
	operator ==(const RpwaComplex<T>& z,
		    const T&          t)
	{ return (z.real() == t) and (z.imag() == T()); }

	template<typename T>
	inline
	HOST_DEVICE
	bool
	operator ==(const T&          t,
		    const RpwaComplex<T>& z)
	{ return (t == z.real()) and (T() == z.imag()); }


	template<typename T>
	inline
	HOST_DEVICE
	bool
	operator !=(const RpwaComplex<T>& z,
		    const T&          t)
	{ return (z.real() != t) or (z.imag() != T()); }

	template<typename T>
	inline
	HOST_DEVICE
	bool
	operator !=(const T&          t,
		    const RpwaComplex<T>& z)
	{ return (t != z.real()) or (T() != z.imag()); }


	//////////////////////////////////////////////////////////////////////////
	// operators with RpwaComplex numbers
	template<typename T>
	inline
	HOST_DEVICE
	RpwaComplex<T>
	operator +(const RpwaComplex<T>& x,
		    const RpwaComplex<T>& y)
	{
		RpwaComplex<T> result = x;
		result += y;
		return result;
	}

	template<typename T>
	inline
	HOST_DEVICE
	RpwaComplex<T>
	operator -(const RpwaComplex<T>& x,
		    const RpwaComplex<T>& y)
	{
		RpwaComplex<T> result = x;
		result -= y;
		return result;
	}

	template<typename T>
	inline
	HOST_DEVICE
	RpwaComplex<T>
	operator *(const RpwaComplex<T>& x,
		    const RpwaComplex<T>& y)
	{
		RpwaComplex<T> result = x;
		result *= y;
		return result;
	}

	template<typename T>
	inline
	HOST_DEVICE
	RpwaComplex<T>
	operator /(const RpwaComplex<T>& x,
		    const RpwaComplex<T>& y)
	{
		RpwaComplex<T> result = x;
		result /= y;
		return result;
	}

	template<typename T>
	inline
	HOST_DEVICE
	bool
	operator ==(const RpwaComplex<T>& x,
		    const RpwaComplex<T>& y)
	{ return (x.real() == y.real()) and (x.imag() == y.imag()); }

	template<typename T>
	inline
	HOST_DEVICE
	bool
	operator !=(const RpwaComplex<T>& x,
		    const RpwaComplex<T>& y)
	{ return (x.real() != y.real()) or (x.imag() != y.imag()); }


	//////////////////////////////////////////////////////////////////////////
	// sign operators
	template<typename T>
	inline
	HOST_DEVICE
	RpwaComplex<T>
	operator +(const RpwaComplex<T>& z)
	{ return z; }

	template<typename T>
	inline
	HOST_DEVICE
	RpwaComplex<T>
	operator -(const RpwaComplex<T>& z)
	{ return RpwaComplex<T>(-z.real(), -z.imag()); }
		

	//////////////////////////////////////////////////////////////////////////
	// accessor functions
	template<typename T> inline HOST_DEVICE T&       real(RpwaComplex<T>&       z) { return z.real(); }
	template<typename T> inline HOST_DEVICE const T& real(const RpwaComplex<T>& z) { return z.real(); }
	template<typename T> inline HOST_DEVICE T&       imag(RpwaComplex<T>&       z) { return z.imag(); }
	template<typename T> inline HOST_DEVICE const T& imag(const RpwaComplex<T>& z) { return z.imag(); }


	//////////////////////////////////////////////////////////////////////////
	// other functions
	template<typename T>
	inline
	HOST_DEVICE
	T
	abs(const RpwaComplex<T>& z)
	{
		T       real = z.real();
		T       imag = z.imag();
		const T s    = max(std::abs(real), std::abs(imag));
		if (s == T())
			return s;
		real /= s; 
		imag /= s;
		return s * std::sqrt(real * real + imag * imag);
	}

	template<typename T>
	inline
	HOST_DEVICE
	T
	arg(const RpwaComplex<T>& z)
	{	return std::atan2(z.imag(), z.real());	}

	template<typename T>
	inline
	HOST_DEVICE
	T
	norm(const RpwaComplex<T>& z)
	{
		const T real = z.real();
		const T imag = z.imag();
		return sqrt(real * real + imag * imag);
		// const T rho = cuda::abs(z);
		// return rho * rho;
	}

	template<typename T>
	inline
	HOST_DEVICE
	RpwaComplex<T>
	polar(const T& rho,
	      const T& phi)
	{ return RpwaComplex<T>(rho * std::cos(phi), rho * std::sin(phi)); }

	template<typename T>
	inline
	HOST_DEVICE
	RpwaComplex<T>
	conj(const RpwaComplex<T>& z)
	{ return RpwaComplex<T>(z.real(), -z.imag()); }

	template<typename T>
	inline
	HOST_DEVICE
	RpwaComplex<T>
	exp(const RpwaComplex<T>& z)
	{ return std::exp(z.real()) * RpwaComplex<T>(std::cos(z.imag()), std::sin(z.imag())); } /// do NOT listen to microsoft, never, please

	//////////////////////////////////////////////////////////////////////////
	// stream operators
	template<typename T, typename CharT, class Traits>
	HOST
	std::basic_istream<CharT, Traits>&
	operator >>(std::basic_istream<CharT, Traits>& in,
		    RpwaComplex<T>&                        z)
	{
		T     real, imag;
		CharT ch;
		in >> ch;
		if (ch == '(') {
			in >> real >> ch;
			if (ch == ',') {
				in >> imag >> ch;
				if (ch == ')') 
					z = RpwaComplex<T>(real, imag);
				else
					in.setstate(std::ios_base::failbit);
			}	else if (ch == ')')
				z = real;
			else
				in.setstate(std::ios_base::failbit);
		} else {
			in.putback(ch);
			in >> real;
			z = real;
		}
		return in;
	}

	template<typename T, typename CharT, class Traits>
	HOST
	std::basic_ostream<CharT, Traits>&
	operator<<(std::basic_ostream<CharT, Traits>& out,
		    const RpwaComplex<T>&                  z)
	{
		std::basic_ostringstream<CharT, Traits> s;
		s.flags(out.flags());
		s.imbue(out.getloc());
		s.precision(out.precision());
		s << '(' << z.real() << ',' << z.imag() << ')';
		return out << s.str();
	}
	

};  // namespace rpwa


#endif // RPWA_COMPLEX_HPP
