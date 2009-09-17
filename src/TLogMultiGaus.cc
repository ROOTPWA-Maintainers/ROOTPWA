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
//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class TLogMultiGaus
//      see TLogMultiGaus.hh for details
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

// This Class' Header ------------------
#include "TLogMultiGaus.h"

// C/C++ Headers ----------------------
#include <fstream>
#include <iostream>
#include <sstream>
using std::cout;
using std::endl;
using std::ifstream;
#include <assert.h>

// Collaborating Class Headers --------
#include "TString.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TSystem.h"
#include "TStopwatch.h"

// Class Member definitions -----------

TLogMultiGaus::TLogMultiGaus() : _dim(0)
{}

TLogMultiGaus::~TLogMultiGaus()
{}




// Likelyhood and first derivatives in one go: *****************************

void 
TLogMultiGaus::FdF(const double* x, double& f, double* df) const {
  TVectorD v(_dim);
  for(unsigned int i=0;i<_dim;++i){
    v[i]=x[i];
  }

  TVectorD vp=_precs*v;
  double res=v*vp;

  for(unsigned int i=0;i<_dim;++i){
    df[i]=vp[i];
  }

  res*=0.5;
  f=res;
  return;
}


//************************************************************************

double 
TLogMultiGaus::DoEval(const double* x) const
{

  // call FdF
  double df[_dim];
  double L;
  FdF(x, L,df);
  return L;

}



unsigned int 
TLogMultiGaus::NDim() const {
  return _dim;
}


double 
TLogMultiGaus::DoDerivative(const double*, unsigned int) const {
  return 0;
}
