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
//      Likelihood function Object to use with ROOT minimizers
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef TLOGMULTIGAUS_HH
#define TLOGMULTIGAUS_HH

// Base Class Headers ----------------


// Collaborating Class Headers -------
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <complex>
using std::complex;


#include "TString.h"
#include "TMatrixD.h"
#include "Math/IFunction.h"
using namespace ROOT;


// Collaborating Class Declarations --


class TLogMultiGaus : public Math::IGradientFunctionMultiDim {
public:

  // Constructors/Destructors ---------
  TLogMultiGaus();
  ~TLogMultiGaus();
  TLogMultiGaus* Clone() const {return new TLogMultiGaus(*this);}
  // Accessors -----------------------
  
  unsigned int NDim() const;


  // Modifiers -----------------------
  // Precisions = 1/sigma^2
  void Set(const TMatrixT<double>& prec){_dim=prec.GetNrows();_precs.ResizeTo(_dim,_dim);_precs=prec;}

  // Operations ----------------------
  // Load Amplitudes into memory


  virtual double DoEval(const double*) const;
  virtual double DoDerivative(const double*, unsigned int) const;

  // Calculate Function f and Derivatives df in one go
  virtual void FdF(const double* x, double& f, double* df) const;

private:
  unsigned int _dim;
  TMatrixT<double> _precs;
 
};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
