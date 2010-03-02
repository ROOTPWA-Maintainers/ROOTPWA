#ifndef __MATRIX_H_
#define __MATRIX_H_
	

#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cassert>
#include <cstring>
#include <cmath>
#include <complex>

#include "Vec.h"


inline double conj(const double x) { return x; }

	
template<typename T> class matrix;


template<typename T> fourVec   operator *= (fourVec& v, const matrix<T>& M);
template<typename T> matrix<T> operator *  (const T& a, const matrix<T>& M);


template<typename T> class matrix {

public:

  matrix();
  matrix(const int n, const int m);
  matrix(const int n, const int m, T* p);
  matrix(const matrix& M);
  ~matrix();
  matrix& operator = (const matrix&);
   
  T&       el     (const int i, const int j);
  T&       element(const int i, const int j) { return el(i, j); }
  const T& el     (const int i, const int j) const;
  const T& element(const int i, const int j) const { return el(i, j); }

  matrix   operator *= (const double    a);
  matrix   operator += (const matrix&   M);
  matrix   operator -= (const matrix&   M);
  matrix   operator *  (const matrix&   M) const;
  matrix   operator +  (const matrix&   M) const;
  matrix   operator -  (const matrix&   M) const;
  threeVec operator *  (const threeVec& v) const;
  fourVec  operator *  (const fourVec&  v) const;
  matrix   operator *  (const T&        a) const;

  matrix conjugate() const;
  matrix transpose() const;
  matrix adjoint  () const;
  int status() const { return _status; }
  void setStatus(const int s) { _status = s; }
  T trace() const;
  T det() const;
  matrix LU() const;
  matrix inv() ;
  int nrows() const { return _nrows; }
  int ncols() const { return _ncols; }
  matrix zero();

  const matrix& print(std::ostream& os = std::cout) const;
  matrix& scan(std::istream& is = std::cin);

private:

  int _nrows;
  int _ncols;
  T*  _f;
  int _status;

  void _create(const int, const int);
  void _destroy(void);
  void _copy(const matrix&);

  matrix _LU(int*, int*) const;
  matrix _lubksb(int* indx, matrix& b);

};

	
template<typename T> class identityMatrix : public matrix<T> {

public:

  identityMatrix(const int n)
    : matrix<T>(n, n)
  {
    for (int i = 0; i < n; ++i)
      this->el(i, i) = 1;
  }
  ~identityMatrix() { }

};
	

template<typename T>
T
matrix<T>::trace() const
{
  T tr = 0;
  assert(_nrows == _ncols);
  for (int i = 0; i < _nrows; ++i)
    tr += el(i, i);
  return tr;
}


template<typename T>
T
matrix<T>::det() const
{
  T    dt   = 1.0;
  int* indx = new int[_ncols];
  int  d;
  assert(_nrows == _ncols);
  matrix r(_nrows, _ncols);
  r = _LU(&d, indx);
  for (int i = 0; i < _nrows; ++i)
    dt *= r.el(i, i);
  if (indx)
    delete[] indx;
  return d * dt;
}


template<typename T>
double
mag(T t)
{
  return abs(std::complex<double>(t));
}


template<typename T>
matrix<T>
matrix<T>::_LU(int *d, int *indx) const
{
  int i, j, k, imax = 0;
  assert(_nrows == _ncols);
  matrix<T> r(_nrows, _ncols);
#define TINY 1.e-20;
  T big, dum, sum;
  // T temp;
  T *vv = new T[_nrows];
  r = *this;
  *d = 1;
  for (i = 0; i < _nrows; ++i) {
    big = 0.0;
    for (j = 0; j < _nrows; ++j) {
      if (mag(r.el(i, j)) > mag(big))
	big = r.el(i, j);
    }
    if (mag(big) == 0.0) {
      this->print();
      std::cerr << "singular matrix " << std::endl;
      r.setStatus(0);
      return(r);
    }
    else
      r.setStatus(1);
    vv[i] = 1.0 / big;
  }
  for (j = 0; j < _nrows; ++j) {
    for (i = 0; i < j; ++i) {
      sum = r.el(i, j);
      for (k = 0; k < i; ++k)
	sum -= r.el(i, k) * r.el(k, j);
      r.el(i, j) = sum;
    }
    big = 0.0;
    for (i = j; i < _nrows; ++i) {
      sum = r.el(i, j);
      for (k = 0; k < j; ++k) {
	sum -= r.el(i, k) * r.el(k, j);
      }
      r.el(i, j) = sum;
      if (mag(vv[i] * sum) >= mag(big)) {
	big = mag(vv[i] * sum);
	imax = i;
      }
    }
    if (j != imax) {      /* do we need to interchange rows? */
      for (k = 0; k < _nrows; ++k) {    /* Yes, make it so */
	dum = r.el(imax, k);
	r.el(imax, k) = r.el(j, k);
	r.el(j, k) = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (r.el(j, j) == 0.0)
      r.el(j, j) = TINY;
    if (j != _nrows) {
      dum = 1.0 / r.el(j, j);
      for (i = j + 1; i < _nrows; ++i)
	r.el(i, j) *= dum;
    }
  }
  if (vv)
    delete[] vv;
  r.setStatus(1);
  return(r);
}


template<typename T>
matrix<T>
matrix<T>::LU() const
{
  int d;
  int *indx = new int[_nrows];
  matrix<T> r(_nrows, _ncols);
  r = _LU(&d, indx);
  if (indx)
    delete[] indx;
  return(r);
}


template<typename T>
matrix<T>
matrix<T>::_lubksb(int*    indx,
		   matrix& b) 
{
  int i, ii = -1, ip, j;
  T sum;
  matrix r(b.nrows(), 1);
  r = b;
  for (i = 0; i < _nrows; ++i) {
    ip = indx[i];
    sum = r.el(ip, 0);
    r.el(ip, 0) = r.el(i, 0);
    if (ii > -1) {
      for (j = ii; j < _nrows; ++j)
	sum -= element(i, j) * r.element(j, 0);
    }
    else if (mag(sum)) 
      ii = i;
    r.el(i, 0) = sum;
  }
  for (i = _nrows-1; i > -1; --i) {
    sum = r.el(i, 0);
    for (j = i + 1; j < _nrows; ++j)
      sum -= element(i, j) * r.element(j, 0);
    r.el(i, 0) = sum / element(i, i);
  }
  return(r);
}
		

template<typename T>
matrix<T>
matrix<T>::inv()
{
  int d;
  int i, j;
  int *indx = new int[_ncols];
  matrix r (_nrows, _ncols);
  matrix lu(_nrows, _ncols);
  lu = _LU(&d, indx);
  if (lu.status()) {
    for (j = 0; j < _nrows; ++j) {
      matrix col(_nrows, 1);
      matrix x(_nrows, 1);
      col.el(j, 0) = 1.0;
      x = lu._lubksb(indx, col);
      for (i = 0; i < _nrows; ++i)
	r.el(i, j) = x.el(i, 0);
    }
    r.setStatus(1);
    if (indx)
      delete[] indx;
  } else {
    r= lu;
    throw "singular matrix";
  }
  return(r);
}
		

template<typename T>
matrix<T>
matrix<T>::conjugate() const
{
  matrix<T> r(_nrows, _ncols);
  for (int i = 0; i < _nrows; ++i)
    for (int j = 0; j < _ncols; ++j)
      r.el(i, j) = conj(el(i, j));
  return(r);
}


template<typename T>
matrix<T>
matrix<T>::transpose() const
{
  matrix<T> r(_ncols, _nrows);
  for (int i = 0; i < _nrows; ++i)
    for (int j = 0; j < _ncols; ++j)
      r.el(j, i) = el(i, j);
  return(r);
}


template<typename T>
matrix<T>
matrix<T>::adjoint() const
{
  matrix<T> r(_ncols, _nrows);
  r = transpose().conjugate();
  return(r);
}


template<typename T>
matrix<T>
matrix<T>::operator + (const matrix& M) const
{
  assert((_ncols == M._ncols) && (_nrows == M._nrows));
  matrix<T> r(_nrows, _ncols);
  for (int i = 0; i < _nrows; ++i)
    for (int j = 0; j < M._ncols; ++j)
      r.el(i, j) = el(i, j) + M.el(i, j);
  return(r);
}


template<typename T>
matrix<T>
matrix<T>::operator += (const matrix& M)
{
  assert((_ncols == M._ncols) && (_nrows == M._nrows));
  for (int i = 0; i < _nrows; ++i)
    for (int j = 0; j < M._ncols; ++j)
      el(i, j) += M.el(i, j);
  return *this;
}

template<typename T>
matrix<T>
matrix<T>::operator *= (const double k)
{
  for (int i = 0; i < _nrows; ++i)
    for (int j = 0; j < _ncols; ++j)
      el(i, j) *= k;
  return *this;
}


template<typename T>
matrix<T>
matrix<T>::operator - (const matrix& M) const
{
  assert((_ncols == M._ncols) && (_nrows == M._nrows));
  matrix<T> r(_nrows, _ncols);
  for (int i = 0; i < _nrows; ++i)
    for (int j = 0; j < M._ncols; ++j)
      r.el(i, j) = el(i, j) - M.el(i, j);
  return(r);
}


template<typename T>
matrix<T>
matrix<T>::operator -= (const matrix& M)
{
  assert((_ncols == M._ncols) && (_nrows == M._nrows));
  for (int i = 0; i < _nrows; ++i)
    for (int j = 0; j < M._ncols; ++j)
      el(i, j) -= M.el(i, j);
  return *this;
}


template<typename T>
matrix<T>
matrix<T>::operator * (const matrix& M) const
{
  assert(_ncols == M._nrows);
  matrix<T> r(_nrows, M._ncols);
  for (int i = 0; i < _nrows; ++i)
    for (int j = 0; j < M._ncols; ++j)
      for (int k = 0; k < _ncols; ++k)
	r.el(i, j) += el(i, k) * M.el(k, j);
  return(r);
}


template<typename T>
matrix<T>
operator * (const T&         a,
	    const matrix<T>& M)
{
  return(M * a);
}


template<typename T>
matrix<T>
matrix<T>::operator * (const T& a) const
{
  matrix<T> r(nrows(), ncols());
  for (int i = 0; i < nrows(); ++i)
    for (int j = 0; j < ncols(); ++j)
      r.el(i, j) = a * element(i, j);
  return(r);
}


template<typename T>
T&
matrix<T>::el(const int i,
	      const int j)
{
  assert((i < _nrows) && (j < _ncols));
  return _f[j + i * _ncols];
}


template<typename T>
const T&
matrix<T>::el(const int i,
	      const int j) const
{
  assert((i < _nrows) && (j < _ncols));
  return _f[j + i * _ncols];
}


template<typename T>
void
matrix<T>::_create(const int n,
		   const int m)
{	
  _nrows = n;
  _ncols = m;
  _f = NULL;
  if (n * m != 0) {
    _f = new T[n * m];
    memset(_f, 0, n * m * sizeof(T));
  }
}


template<typename T>
void
matrix<T>::_copy(const matrix<T>& src)
{
  _nrows = src._nrows;
  _ncols = src._ncols;
  _status = src.status();
  if (_f)
    delete[] _f;
  _f = new T[src._nrows * src._ncols];
  memcpy(_f, src._f, src._nrows * src._ncols * sizeof(T));
}


template<typename T>
void
matrix<T>::_destroy(void)
{
  if (_f)
    delete[] _f;
  _f = NULL;
  _nrows = 0;
  _ncols = 0;
}


template<typename T>
matrix<T>::matrix()
{
  _create(0, 0);
}


template<typename T>
matrix<T>::matrix(const int n,
		  const int m)
{
  _create(n, m);
}


template<typename T>
matrix<T>::matrix(const int n,
		  const int m,
		  T*        p)
{
  _create(n, m);
  _f = p;
}


template<typename T>
matrix<T>::matrix(const matrix<T>& M)
{
  _create(M._nrows, M._ncols);
  _copy(M);
}


template<typename T>
matrix<T>::~matrix()
{
  _destroy();
}


template<typename T>
matrix<T>&
matrix<T>::operator = (const matrix<T>& M)
{
  _copy(M);
  return *this;
}


template<typename T>
const matrix<T>&
matrix<T>::print(std::ostream& os) const
{
  const unsigned int nmbDigits = std::numeric_limits<double>::digits10 + 1;
  std::ostringstream s;
  s.precision(nmbDigits);
  s.setf(std::ios_base::scientific, std::ios_base::floatfield);
  s << _nrows << " " << _ncols << std::endl;
  for (int i = 0; i < _nrows; ++i) {
    for (int j = 0; j < _ncols; ++j)
      s << _f[j + i * _ncols] << "\t";
    s << std::endl;
  }
  os << s.str();
  return *this;
}


template<typename T>
matrix<T>&
matrix<T>::scan(std::istream& is)
{
  int i = 0, j = 0;
  is >> i;
  is >> j;
  _create(i, j);
  for (int i = 0; i < _nrows; ++i)
    for (int j = 0; j < _ncols; ++j) {
      is >> _f[j + i * _ncols];
      if (is.eof())
	throw "trunc matrix input";
    }
  return *this;
}


template<typename T>
std::ostream&
operator << (std::ostream&    os,
	     const matrix<T>& m)
{
  m.print(os);
  return os;
}


template<typename T>
std::istream&
operator >> (std::istream& is,
	     matrix<T>&    m)
{
  m.scan(is);
  return is;
}


template<typename T>
matrix<T>
matrix<T>::zero()
{
  memset(_f, 0, _nrows * _ncols * sizeof(T));
  return *this;
}


template<typename T>
threeVec
matrix<T>::operator * (const threeVec& V) const
{
  threeVec R;
  assert((_nrows == 3) && (_ncols == 3));
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      R.el(i) += el(i, j) * V.el(j); 
  return(R);
}


template<typename T>
fourVec
matrix<T>::operator * (const fourVec& v) const
{
  fourVec r;
  assert((_nrows == 4) && (_ncols == 4));
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      r.el(i) += el(i, j) * v.el(j); 
  return(r);
}


template<typename T>
fourVec
operator *= (fourVec& v, const matrix<T>& M)
{
  v = M * v;
  return v;
}


#endif
