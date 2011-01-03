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
//#include <fitlog.h>
//#include <integral.h>
#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <map>
#include <complex>
#include "TBranch.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TComplex.h"
#include "TRandom.h"
#include "TPWALikelihood.h"
#include "TFitBin.h"

using namespace std;


int lineno = 1; // global variables needed for lex (not understood)
char *progname;

int main(int argc, char** argv){
  
 	// TPWALikelihood<complex<double> > L;
 	TPWALikelihood<double> L;
  L.useNormalizedAmps();
  //if(quiet)L.SetQuiet();
  int rank=2; // TODO: make this an option
  L.init(rank, argv[1], "norm.int", "norm.int");
  cout<<L.NDim()<<endl;

  // 12 parameters + flat
  double x[13]={0.52707,0.21068,-0.604365,0.17596,-0.216668,-0.0990815,-0.348459,0.208961,0,0,0,0,0};
  //double x[13]; for(int i=0;i<13;++i)x[i]=0.001;
  //string a[13]={"a","b","c","d","e","f","g","h","i","j","k","l","flat"};
  std::cout<<L.DoEval(x)<<std::endl;


 // check derivatives:
  // first numerical:
  double L1=L.DoEval(x);
  double h=1E-8;
  double dxNum[13];
  double dxAna[13];
  for(unsigned int i=0; i<L.NDim();++i){
    x[i]+=h;
    double L2=L.DoEval(x);
    dxNum[i]=(L2-L1)/h;
    x[i]-=h;
  }
  double F;
  L.FdF(x,F,dxAna);
  for(unsigned int i=0; i<L.NDim();++i){
    cout<< "dL/d"<<i<<"(num)="<<dxNum[i]<<endl;
    cout<< "dL/d"<<i<<"(ana)="<<dxAna[i]<<endl;
  }
  
    

  return 0;
}

// dummy function needed since we link to but do not use minuit.o
int mnparm(int, string, double, double, double, double) {
    cerr << "this is impossible" << endl;
    throw "aFit";
    return 0;
}
