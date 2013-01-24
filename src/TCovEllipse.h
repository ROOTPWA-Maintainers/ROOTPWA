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
//
// Description:
//      A covariance ellipse
//
//
// Environment:
//      Software developed for the PANDA Detector at FAIR.
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef TCOVELLIPSE_HH
#define TCOVELLIPSE_HH

// Base Class Headers ----------------
#include "TEllipse.h"

// Collaborating Class Headers -------
#include "TMatrixT.h"

// Collaborating Class Declarations --



class TCovEllipse : public TEllipse{
public:

  // Constructors/Destructors ---------
  TCovEllipse();
  TCovEllipse(const TMatrixT<double>& cov, 
	      double x0=0, double y0=0,
	      int i=0, int j=1);
  virtual ~TCovEllipse(){}

  // Operators
  

  // Accessors -----------------------
  TMatrixT<double> getCov() const {return _cov;}

  // Modifiers -----------------------
  void setCov(const TMatrixT<double>& cov);
  void setMean(double x, double y);
  void select(int i, int j){_i=i;_j=j;recalc();}

  // Operations ----------------------

private:

  // Private Data Members ------------
  TMatrixT<double> _cov;
  double _sig;
  double _x0;
  double _y0;

  int _i,_j; // row and column index to be plotted

  // Private Methods -----------------
  void recalc();

public: 
  ClassDef(TCovEllipse,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
