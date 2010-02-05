#ifndef __MATRIX_H_
#define __MATRIX_H_
	

#include <iostream>
#include <cassert>
#include <cstring>
#include <cmath>
#include <complex>
#include <Vec.h>


double conj(double x);

	
template <class Type> class matrix;

template <class Type> fourVec operator *= (fourVec&, const matrix<Type>&);

template <class Type> matrix<Type> operator * (Type, matrix<Type>&);

template <class Type> class matrix {

private:

  int   _nrows;
  int   _ncols;
  Type* _f;
  int   _status;

  void _create(int, int);
  void _destroy(void);
  void _copy(const matrix&);

  matrix _LU(int*, int*) const;
  matrix _lubksb(int* indx, matrix& b);

public:

  matrix();
  matrix(int n, int m);
  matrix(int, int, Type*);
  matrix(const matrix&);
  ~matrix();
  matrix& operator = (const matrix&);
   
  Type& el(int, int);
  Type& element(int i,int j) { return this->el(i, j); }
  const Type& el(int, int) const;
  const Type& element(int i, int j) const { return this->el(i, j); }
  matrix operator * (const matrix&) const;
  matrix operator *= (const double);
  matrix operator + (const matrix&) const;
  matrix operator += (const matrix&);
  matrix operator - (const matrix&) const;
  matrix operator -= (const matrix&);
  matrix conjugate() const;
  matrix transpose() const;
  matrix adjoint() const;
  int status() const { return _status; }
  void setStatus(int s) {this->_status = s; }
  Type trace() const;
  Type det() const;
  matrix LU() const;
  matrix inv() ;
  int nrows() const { return this->_nrows; }
  int ncols() const { return this->_ncols; }
  matrix zero();
  threeVec operator * (const threeVec&) const;
  fourVec operator * (const fourVec&) const;
  matrix operator * (Type) const;

  //friend fourVec operator *=<> (fourVec&, const matrix&);
  // friend matrix operator *<> (Type, matrix&);

  const matrix& print(std::ostream& os = std::cout) const;
  matrix& scan(std::istream& is = std::cin);

};

	
template <class Type> class identityMatrix : public matrix<Type> {

public:

  identityMatrix(int n)
    : matrix<Type>(n,n)
  {
    for (int i = 0; i < n; ++i)
      this->el(i, i) = 1;
  }
  ~identityMatrix() { }

};
	

template <class Type>
Type
matrix<Type>::trace() const
{
  Type tr = 0;
  assert(this->_nrows == this->_ncols);
  for (int i = 0; i < this->_nrows; ++i)
    tr += (const_cast<matrix<Type>*>(this))->el(i, i);
  return(tr);
}


template <class Type>
Type
matrix<Type>::det() const
{
  Type dt = 1.0;
  int *indx = new int[this->_ncols];
  int d;
  assert(this->_nrows == this->_ncols);
  matrix r(this->_nrows, this->_ncols);
  r = this->_LU(&d, indx);
  for (int i = 0; i < this->_nrows; ++i)
    dt *= r.el(i, i);
  if (indx)
    delete[] indx;
  return(d * dt);
}


template <class Type>
double
mag(Type t)
{
  return abs(std::complex<double>(t));
}


template <class Type>
matrix<Type>
matrix<Type>::_LU(int *d, int *indx) const
{
  int i, j, k, imax = 0;
  assert(this->_nrows == this->_ncols);
  matrix<Type> r(this->_nrows, this->_ncols);
#define TINY 1.e-20;
  Type big, dum, sum;
  // Type temp;
  Type *vv = new Type[this->_nrows];
  r = *this;
  *d = 1;
  for (i = 0; i < this->_nrows; ++i) {
    big = 0.0;
    for (j = 0; j < this->_nrows; ++j) {
      //			if ((temp = fabs(r.el(i,j))) > big) 
      //				big = temp;
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
  for (j = 0; j < this->_nrows; ++j) {
    for (i = 0; i < j; ++i) {
      sum = r.el(i, j);
      for (k = 0; k < i; ++k)
	sum -= r.el(i, k) * r.el(k, j);
      r.el(i, j) = sum;
    }
    big = 0.0;
    for (i = j; i < this->_nrows; ++i) {
      sum = r.el(i, j);
      for (k = 0; k < j; ++k) {
	sum -= r.el(i, k) * r.el(k, j);
      }
      r.el(i, j) = sum;
      //			if (( dum=vv[i] * fabs(sum)) >= big) {
      //				big = dum;
      //				imax = i;
      //			}
      if (mag(vv[i] * sum) >= mag(big)) {
	big = mag(vv[i] * sum);
	imax = i;
      }
    }
    if (j != imax) {      /* do we need to interchange rows? */
      for (k = 0; k < this->_nrows; ++k) {    /* Yes, make it so */
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
    if (j != this->_nrows) {
      dum = 1.0 / r.el(j, j);
      for (i = j + 1; i < this->_nrows; ++i)
	r.el(i, j) *= dum;
    }
  }
  if (vv)
    delete[] vv;
  r.setStatus(1);
  return(r);
}


template <class Type>
matrix<Type>
matrix<Type>::LU() const
{
  int d;
  int *indx = new int[this->_nrows];
  matrix<Type> r(this->_nrows, this->_ncols);
  r = this->_LU(&d, indx);
  if (indx)
    delete[] indx;
  return(r);
}


template <class Type>
matrix<Type>
matrix<Type>::_lubksb(int *indx, matrix& b) 
{
  int i, ii = -1, ip, j;
  Type sum;
  matrix r(b.nrows(), 1);
  r = b;
  for (i = 0; i < this->_nrows; ++i) {
    ip = indx[i];
    sum = r.el(ip, 0);
    r.el(ip, 0) = r.el(i, 0);
    if (ii > -1) {
      for (j = ii; j < this->_nrows; ++j)
	sum -= this->element(i, j) * r.element(j, 0);
    }
    else if (mag(sum)) 
      ii = i;
    r.el(i, 0) = sum;
  }
  for (i = this->_nrows-1; i > -1; --i) {
    sum = r.el(i, 0);
    for (j = i + 1; j < this->_nrows; ++j)
      sum -= this->element(i, j) * r.element(j, 0);
    r.el(i, 0) = sum / this->element(i, i);
  }
  return(r);
}
		

template <class Type>
matrix<Type>
matrix<Type>::inv()
{
  int d;
  int i, j;
  int *indx = new int[this->_ncols];
  matrix r (this->_nrows, this->_ncols);
  matrix lu(this->_nrows, this->_ncols);
  lu = this->_LU(&d, indx);
  if (lu.status()) {
    for (j = 0; j < this->_nrows; ++j) {
      matrix col(this->_nrows, 1);
      matrix x(this->_nrows, 1);
      col.el(j, 0) = 1.0;
      x = lu._lubksb(indx, col);
      for (i = 0; i < this->_nrows; ++i)
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
		

template <class Type>
matrix<Type>
matrix<Type>::conjugate() const
{
  matrix<Type> r(this->_nrows, this->_ncols);
  for (int i = 0; i < this->_nrows; ++i)
    for (int j = 0; j < this->_ncols; ++j)
      r.el(i, j) = conj((const_cast<matrix<Type>*>(this))->el(i, j));
  return(r);
}


template <class Type>
matrix<Type>
matrix<Type>::transpose() const
{
  matrix<Type> r(this->_ncols, this->_nrows);
  for (int i = 0; i < this->_nrows; ++i)
    for (int j = 0; j < this->_ncols; ++j)
      r.el(j, i) = (const_cast<matrix<Type>*>(this))->el(i, j);
  return(r);
}


template <class Type>
matrix<Type>
matrix<Type>::adjoint() const
{
  matrix<Type> r(this->_ncols, this->_nrows);
  r = (this->transpose()).conjugate();
  return(r);
}


template <class Type>
matrix<Type>
matrix<Type>::operator + (const matrix& M) const
{
  assert((this->_ncols == M._ncols) && (this->_nrows == M._nrows));
  matrix<Type> r(this->_nrows, this->_ncols);
  for (int i = 0; i < this->_nrows; ++i)
    for (int j = 0; j < M._ncols; ++j)
      r.el(i, j) =(const_cast<matrix<Type>*>(this))->el(i, j)
	          + (const_cast<matrix<Type>*>(&M))->el(i, j);
  return(r);
}


template <class Type>
matrix<Type>
matrix<Type>::operator += (const matrix& M)
{
  assert((this->_ncols == M._ncols) && (this->_nrows == M._nrows));
  for (int i = 0; i < this->_nrows; ++i)
    for (int j = 0; j < M._ncols; ++j)
      this->el(i, j) += (const_cast<matrix<Type>*>(&M))->el(i, j);
  return(*this);
}

template <class Type>
matrix<Type>
matrix<Type>::operator *= (const double k)
{
  for (int i = 0; i < this->_nrows; ++i)
    for (int j = 0; j < this->_ncols; ++j)
      this->el(i, j) *= k;
  return(*this);
}


template <class Type>
matrix<Type>
matrix<Type>::operator - (const matrix& M) const
{
  assert((this->_ncols == M._ncols) && (this->_nrows == M._nrows));
  matrix<Type> r(this->_nrows, this->_ncols);
  for (int i = 0; i < this->_nrows; ++i)
    for (int j = 0; j < M._ncols; ++j)
      r.el(i, j) = (const_cast<matrix<Type>*>(this))->el(i, j)
	           - (const_cast<matrix<Type>*>(&M))->el(i, j);
  return(r);
}


template <class Type>
matrix<Type>
matrix<Type>::operator -= (const matrix& M)
{
  assert((this->_ncols == M._ncols) && (this->_nrows == M._nrows));
  for (int i = 0; i < this->_nrows; ++i)
    for (int j = 0; j < M._ncols; ++j)
      this->el(i, j) -= (const_cast<matrix<Type>*>(&M))->el(i, j);
  return(*this);
}


template <class Type>
matrix<Type>
matrix<Type>::operator * (const matrix& M) const
{
  assert(this->_ncols == M._nrows);
  matrix<Type> r(this->_nrows, M._ncols);
  for (int i = 0; i < this->_nrows; ++i)
    for (int j = 0; j < M._ncols; ++j)
      for (int k = 0; k < this->_ncols; ++k)
	r.el(i, j) += (const_cast<matrix<Type>*>(this))->el(i, k)
	              * (const_cast<matrix<Type>*>(&M))->el(k, j);
  return(r);
}


template <class Type>
matrix<Type>
operator * (Type a, matrix<Type>& M)
{
  return(M * a);
}


template <class Type>
matrix<Type>
matrix<Type>::operator * (Type a) const
{
  matrix<Type> r(this->nrows(), this->ncols());
  for (int i = 0; i < this->nrows(); ++i)
    for (int j = 0; j < this->ncols(); ++j)
      r.el(i, j) = a * this->element(i, j);
  return(r);
}


template <class Type>
Type&
matrix<Type>::el(int i, int j)
{
  assert((i < this->_nrows) && (j < this->_ncols));
  return _f[j + i * _ncols];
}


template <class Type>
const Type&
matrix<Type>::el(int i, int j) const
{
  assert((i < this->_nrows) && (j < this->_ncols));
  return _f[j + i * _ncols];
}


template <class Type>
void
matrix<Type>::_create(int n, int m)
{	
  this->_nrows = n;
  this->_ncols = m;
  this->_f = NULL;
  if (n * m != 0) {
    this->_f = new Type[n * m];
    memset(this->_f, 0, n * m * sizeof(Type));
  }
}


template <class Type>
void
matrix<Type>::_copy(const matrix<Type>& src)
{
  this->_nrows = src._nrows;
  this->_ncols = src._ncols;
  this->_status = src.status();
  if (this->_f)
    delete[] this->_f;
  this->_f = new Type[src._nrows * src._ncols];
  memcpy(this->_f, src._f, src._nrows * src._ncols * sizeof(Type));
}


template <class Type>
void
matrix<Type>::_destroy(void)
{
  if (this->_f)
    delete[] this->_f;
  this->_f = NULL;
  this->_nrows = 0;
  this->_ncols = 0;
}


template <class Type>
matrix<Type>::matrix()
{
  this->_create(0, 0);
}


template <class Type>
matrix<Type>::matrix(int i, int j)
{
  this->_create(i, j);
}


template <class Type>
matrix<Type>::matrix(int i, int j, Type* p)
{
  this->_create(i, j);
  this->_f = p;
}


template <class Type>
matrix<Type>::matrix(const matrix<Type>& M)
{
  this->_create(M._nrows, M._ncols);
  this->_copy(M);
}


template <class Type>
matrix<Type>::~matrix()
{
  this->_destroy();
}


template <class Type>
matrix<Type>&
matrix<Type>::operator = (const matrix<Type>& M)
{
  this->_copy(M);
  return(*this);
}


template <class Type>
const matrix<Type>&
matrix<Type>::print(std::ostream& os) const
{
  os << _nrows << " " << _ncols << std::endl;
  for (int i = 0; i < _nrows; ++i) {
    for (int j = 0; j < _ncols; ++j)
      os << _f[j + i * _ncols] << "\t";
    os << std::endl;
  }
  return(*this);
}


template <class Type>
matrix<Type>&
matrix<Type>::scan(std::istream& is)
{
  int i = 0, j = 0;
  is >> i;
  is >> j;
  this->_create(i, j);
  for (int i = 0; i < _nrows; ++i)
    for (int j = 0; j < _ncols; ++j) {
      is >> this->_f[j + i * _ncols];
      if (is.eof())
	throw "trunc matrix input";
    }
  return(*this);
}


template <class Type>
std::ostream&
operator << (std::ostream& os, const matrix<Type>& m)
{
  m.print(os);
  return os;
}


template <class Type>
std::istream&
operator >> (std::istream& is, matrix<Type>& m)
{
  m.scan(is);
  return is;
}


template <class Type>
matrix<Type>
matrix<Type>::zero()
{
  memset(this->_f, 0, this->_nrows * this->_ncols * sizeof(Type));
  return(*this);
}


template <class Type>
threeVec
matrix<Type>::operator * (const threeVec& V) const
{
  threeVec R;
  assert((this->_nrows == 3) && (this->_ncols == 3));
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      R.el(i) += (const_cast<matrix<Type>*>(this))->el(i, j)
	         * (const_cast<threeVec*>(&V))->el(j); 
  return(R);
}


template <class Type>
fourVec
matrix<Type>::operator * (const fourVec& v) const
{
  fourVec r;
  assert((this->_nrows == 4) && (this->_ncols == 4));
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      r.el(i) += (const_cast<matrix<Type>*>(this))->el(i, j)
	         * (const_cast<fourVec*>(&v))->el(j); 
  return(r);
}


template <class Type>
fourVec
operator *= (fourVec& v, const matrix<Type>& M)
{
  v = M * v;
  return v;
}


#endif
