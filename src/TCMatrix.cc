///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
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
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////


#include "TCMatrix.h"


using namespace std;


ClassImp(TCMatrix)


TCMatrix::TCMatrix(const int i,
		   const int j)
: _re(i, j),
  _im(i, j)
{ }


void
TCMatrix::set(const int              i,
	      const int              j,
	      const complex<double>& c)
{
  _re[i][j] = c.real();
  _im[i][j] = c.imag();
}
void
TCMatrix::set(const int              i,
	      const int              j,
	      const TComplex& c)
{
  _re[i][j] = c.Re();
  _im[i][j] = c.Im();
}

TCMatrix
TCMatrix::t() const {
  TMatrixD ReT(TMatrixD::kTransposed,_re);
  TMatrixD ImT(TMatrixD::kTransposed,_im);
  return TCMatrix(ReT,ImT);
}

TCMatrix
TCMatrix::dagger() const {
  TMatrixD ReT(TMatrixD::kTransposed,_re);
  TMatrixD ImT(_im.GetNcols(),_im.GetNrows());
  for(int i=0;i<_im.GetNrows();++i){
    for(int j=0;j<_im.GetNcols();++j){
      ImT[j][i]= -(_im[i][j]);
    }
  }
  return TCMatrix(ReT,ImT);
}


TCMatrix operator *(const TCMatrix& c1, const TCMatrix& c2){
  TCMatrix result(c1.nrows(),c2.ncols());
  for(int i=0;i<result.nrows();++i){
    for(int j=0;j<result.ncols();++j){
      TComplex elem(0,0);
      for(int k=0;k<c1.ncols();++k){
	elem+=(c1(i,k))*(c2(k,j));
      }
      result.set(i,j,elem);
    }
  }
  return result;
}

TCMatrix operator -(const TCMatrix& c1, const TCMatrix& c2){
  TMatrixD Re(c1._re-c2._re);
  TMatrixD Im(c1._im-c2._im);
  return TCMatrix(Re,Im);
}

TCMatrix operator +(const TCMatrix& c1, const TCMatrix& c2){
  TMatrixD Re(c1._re+c2._re);
  TMatrixD Im(c1._im+c2._im);
  return TCMatrix(Re,Im);
}
