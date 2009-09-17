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
// a matrix with complex elements
#ifndef TCMATRIX_HH
#define TCMATRIX_HH


#include "TObject.h"
#include <complex>
#include "TComplex.h"
#include "TMatrixD.h"

using std::complex;

class TCMatrix : public TObject {
public:
  TCMatrix(){}
  TCMatrix(int i, int j);
  ~TCMatrix(){};

  void ResizeTo(int i, int j) { _re.ResizeTo(i,j); _im.ResizeTo(i,j);}
  void set(int i, int j, complex<double> c);
  TComplex get(int i, int j);
  TComplex operator() (int i, int j);
  int nrows() const {return _re.GetNrows();}
  int ncols() const {return _re.GetNcols();}
  void Print(){_re.Print(),_im.Print();}

private:
  TMatrixD _re;
  TMatrixD _im;


public:
  ClassDef(TCMatrix,2);


};

#endif
