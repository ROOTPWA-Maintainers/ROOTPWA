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
//test  program for rootpwa

#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <map>
#include <complex>

#include "TPWAAmp.h"
#include "TPWAPhaseConstraint.h"

using namespace std;

int main(int argc, char** argv){

  TPWAAmp a("1-1++1+blub",1);
  if(a.reflectivity()!=1)return 1;
  double par[2];
  par[0]=2;par[1]=4;
  cout << a.name() << a.setPar(par) << endl;
  cout << "Npar=" << a.npar() << endl;
  cout << "Par0=" << a.par(0) << "   Par1=" << a.par(1) << endl;

  TPWAAmp b("1-1++1-bla",1);
  if(b.reflectivity()!=-1)return 1;
  b.setConstraint(new TPWARealConstraint());
  
  par[0]=5;par[1]=6;
  cout << b.name() << b.setPar(par) << endl;
  cout << "Npar=" << b.npar() << endl;
  if(b.npar()!=1)return 1;
  cout << "Par0=" << b.par(0) << "   Par1=" << b.par(1) << endl;
  if(b.amp().imag()!=0)return 2;


  TPWAPhaseConstraint* c1= new TPWAPhaseConstraint(0.3,&b);
  a.setConstraint(c1);
  cout << a.name() << a.amp() << endl;
  
  if(abs(a.amp())!=a.par(0))return 3;
  return 0;
}

