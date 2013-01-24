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
//      Implementation of class TCovEllipse
//      see TCovEllipse.hh for details
//
// Environment:
//      Software developed for the PANDA Detector at FAIR.
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

// Panda Headers ----------------------

// This Class' Header ------------------
#include "TCovEllipse.h"

// C/C++ Headers ----------------------
#include <iostream>

// Collaborating Class Headers --------
#include "TMath.h"

// Class Member definitions -----------
ClassImp(TCovEllipse)

TCovEllipse::TCovEllipse()
  : TEllipse(), _sig(1)
{
  SetFillStyle(0);
}

TCovEllipse::TCovEllipse(const TMatrixT<double>& cov,
			 double x0, double y0,
			 int i, int j)
  : TEllipse(),_cov(cov), _x0(x0), _y0(y0), _i(i), _j(j)
{
  recalc(); 
  SetFillStyle(0);
}

void
TCovEllipse::setCov(const TMatrixT<double>& cov)
{
  _cov=cov;
  recalc();
}

void
TCovEllipse::setMean(double x, double y)
{
  SetX1(x);
  SetY1(y);
}

void
TCovEllipse::recalc()
{
  double sig1=sqrt(_cov[_i][_i]);
  double sig2=sqrt(_cov[_j][_j]);
  double rho=_cov[_i][_j]/(sig1*sig2);
  double theta=0;
  if(_cov[_i][_i]!=_cov[_j][_j]){
    double term=2.*rho*sig1*sig2/(_cov[_i][_i]-_cov[_j][_j]);
    theta=TMath::ATan(term)*0.5;
  }
  //std::cout<<theta<<std::endl;
  double num=_cov[_i][_i]*_cov[_j][_j]*(1-rho*rho);
  double mixed=2*rho*sig1*sig2*sin(theta)*cos(theta);
  double sin2=sin(theta)*sin(theta);
  double cos2=cos(theta)*cos(theta);
  double p12_=_cov[_j][_j]*cos2-mixed+_cov[_i][_i]*sin2;
  double p22_=_cov[_j][_j]*sin2+mixed+_cov[_i][_i]*cos2; 
  double p1=sqrt(num/p12_);
  double p2=sqrt(num/p22_);
  SetTheta(theta*57.2958);
  SetR1(p1);
  SetR2(p2);
  SetX1(_x0);
  SetY1(_y0);
}
