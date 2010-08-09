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


namespace rpwa{
  template <typename T>
  class __align__(ALIGN) complex
  {
  public:
    T _re;
    T _im;
    __host__ __device__ complex(T re = 0, T im = 0) : _re(re), _im(im) {}
    /* __host__ __device__ complex <T>& operator=(const complex <T>& c) */
    /* { */
    /*   _re=c._re; _im=c._im; */
    /*   return *this; */
    /* }     */
    __host__ __device__ friend complex <T> operator+(const complex <T>& a, const complex <T>& b)      
    {
      complex <T> result;
      result._re = a._re + b._re;
      result._im = a._im + b._im;
      return result;
    }
  
    __host__ __device__ friend complex <T> operator*(const complex <T>& a, const complex <T>& b)
    {
      complex <T> result;
      result._re = (a._re * b._re) - (a._im * b._im);
      result._im = (a._re * b._im) + (a._im * b._re);
      return result;
    }
    __host__ __device__ friend T norm(const complex <T>& a) {return (sqrt((a._re * a._re) + (a._im * a._im)));}
    __host__ __device__ friend T real(const complex <T>& a) {return a._re;}
    __host__ __device__ friend T imag(const complex <T>& a) {return a._im;}
  };
  
  /*  class __align__(ALIGN) ccomplex
  {
  public:
    double2 c;
    __host__ __device__ ccomplex() : c.x(0), c.y(0) {}
    __host__ __device__ ccomplex(double re, double im) : c.x(re), c.y(im) {}
    __host__ __device__ ccomplex(double re) : c.x(re), c.y(0) {}
  
    __host__ __device__ friend ccomplex operator+(ccomplex a, ccomplex b)
    {
      return make_double2(a.x + b.x, a.y + b.y);
    }

    __host__ __device__ friend ccomplex operator*(ccomplex a, ccomplex b)
    {
      return make_double2((a.x * b.x) - (a.y * b.y),(a.x * b.y) + (a.y * b.x));
    }

    __host__ __device__ friend double norm(ccomplex a)
    {
      return (sqrt((a.x * a.x) + (a.y * a.y)));
    }

    __host__ __device__ friend double real(ccomplex a)
    {
      return (a.x);
    }
    __host__ __device__ friend double imag(ccomplex a)
    {
      return (a.y);
    }
    };*/
};
