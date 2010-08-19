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
//      Complex class
//
//
// Author List:
//      Philipp Meyer          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef COMPLEX_CUH
#define COMPLEX_CUH


#include "../utilities/cudaUtils.hpp"


#define COMPLEX_ALIGN 2 * sizeof(double)


namespace rpwa {

	namespace cuda {


		template<typename T> class
		ALIGN(COMPLEX_ALIGN)
		complex {

		public:

			typedef T value_type;

			T _re;
			T _im;

			HOST_DEVICE complex(T re = 0, T im = 0) : _re(re), _im(im)        { }
			HOST_DEVICE complex(const complex<T>& z) : _re(z._re), _im(z._im) { }

			inline HOST_DEVICE T&       real()       { return _re; }
			inline HOST_DEVICE const T& real() const { return _re; }
			inline HOST_DEVICE T&       imag()       { return _im; }
			inline HOST_DEVICE const T& imag() const { return _im; }
		
			inline HOST_DEVICE complex<T>& operator =(const complex<T>& z)
				{
					_re = z.real();
					_im = z.imag();
					return *this;
				} 
		
			inline HOST_DEVICE complex<T>& operator +=(const complex<T>& z)
			{
				_re += z.real();
				_im += z.imag();
				return *this;
			}
    
			inline HOST_DEVICE friend complex<T> operator+(const complex<T>& a, const complex<T>& b)      
			{
				complex<T> result;
				result._re = a._re + b._re;
				result._im = a._im + b._im;
				return result;
			}
  
			inline HOST_DEVICE friend complex<T> operator*(const complex<T>& a, const complex<T>& b)
				{
					complex<T> result;
					result._re = (a._re * b._re) - (a._im * b._im);
					result._im = (a._re * b._im) + (a._im * b._re);
					return result;
				}

			inline HOST_DEVICE friend T abs (const complex<T>& z) { return sqrt(norm(z));                     }
			inline HOST_DEVICE friend T norm(const complex<T>& z) { return (z._re * z._re) + (z._im * z._im); }
			inline HOST_DEVICE friend T real(const complex<T>& z) { return z._re;                             }
			inline HOST_DEVICE friend T imag(const complex<T>& z) { return z._im;                             }

		};
  

	}  // namespace cuda

};  // namespace rpwa


#endif // COMPLEX_CUH
