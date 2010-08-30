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

/** @brief numerical evaluation of cauchy integrals 
    following D.B.Hunter, Numer. Math. 19,419-424 (1972)
 */


#ifndef CAUCHYINTEGRAL_HH
#define CAUCHYINTEGRAL_HH

#include "TF1.h"
#include <vector>


// Base Class Headers ----------------


// Collaborating Class Declarations --

struct realPole {
  realPole(double xv,double r):x(xv),residual(r){}
  double x;
  double residual;
};


class cauchyIntegral {
public:

  cauchyIntegral(TF1* func, const std::vector<realPole>& poles,
		 double rangelow, double rangeup);

  /// calculates the principal value integral
  double eval_Hunter(unsigned int degree=4); /// degree>0

  void setRange(double low, double up){_rlow=low;_rup=up;}

private:
  std::vector<realPole> _poles; /// Poles of the function
  double _rlow;               /// lower boundary of integration range
  double _rup;                /// upper boundary of integration range
  double _range;              /// _rlow+_rup
  double _diff;               /// __rup-rlow
  TF1* _func;

  double H(unsigned int r,     /// index of zero of legendre polynomial
	   unsigned int degree);

  std::vector<std::vector<double> > _roots; ///  roots of the first 5 legendre polys 

  double trafo(double t);
  double ofart(double x);

};


#endif
