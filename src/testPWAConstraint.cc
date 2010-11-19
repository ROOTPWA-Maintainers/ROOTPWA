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
#include "TPWALikelihoodC.h"

using namespace std;


int lineno = 1; // global variables needed for lex (not understood)
char *progname;

int main(int argc, char** argv){
  
  TPWALikelihoodC LC;
  TPWALikelihoodC L;
  unsigned int rank=1; 
  // There is one phase constraint in phase.constr!
  LC.Init("wavelist2",rank,"norm.int","norm.int",20000,"phase.constr");
  L.Init("wavelist2",rank,"norm.int","norm.int",20000);
  LC.LoadAmplitudes();
  L.LoadAmplitudes();
  // There is one phase constraint in phase.constr!
  cout << "LC.NDim()=" << LC.NDim() << "   L.NDim()=" << L.NDim() << endl;
  if(LC.NDim()!=L.NDim()-1)return 1;

  // 12 parameters + flat
  // wavelist2 has 4 waves:
  // 1-0-+0+rho770_11_pi-.amp
  // 1-0-+0+sigma_00_pi-.amp
  // 1-1++1-rho770_01_pi-.amp
  // 1-1++0+sigma_10_pi-.amp

  //             Re0      Re1      Im1     Re2      Im2        Re3       Im3
  double xL[13]={0.52707,0.21068,-0.604365,0.17596,-0.216668,-0.0990815,-0.348459,0.208961,0,0,0,0,0};
  // flat

  // now we contrain second and third wave:
  // Phase {
  // 1-0-+0+sigma_00_pi-.amp
  // 1-1++1-rho770_01_pi-.amp
  // 0.3
  // }

  //              Re0      Mag1      Re2      Im2        Re3       Im3    Flat
  double x[13]={0.52707,0.21068,0.17596,-0.216668,-0.0990815,-0.348459,0.208961,0,0,0,0,0};
  //double x[13]; for(int i=0;i<13;++i)x[i]=0.001;
  //string a[13]={"a","b","c","d","e","f","g","h","i","j","k","l","flat"};
  double LLC=LC(x);
  double LL=L(xL);
  std::cout<<"L(xL)="<<LL<<std::endl;
  std::cout<<"LC(x)="<<LLC<<std::endl;
  
  LC.Print();

 // check derivatives:
  // first numerical:
  double LC1=LC(x);
  double h=1E-8;
  double dxNumC[13];
  double dxAnaC[13];
  bool problem=false;
  for(unsigned int i=0; i<LC.NDim();++i){
    x[i]+=h;
    double LC2=LC(x);
    dxNumC[i]=(LC2-LC1)/h;
    x[i]-=h;
  }
  // Then analytic
   double FC;
  LC.FdF(x,FC,dxAnaC);
  for(unsigned int i=0; i<LC.NDim();++i){
    if(2*fabs(dxNumC[i]-dxAnaC[i])/(fabs(dxNumC[i])+fabs(dxAnaC[i]))>0.0001){
      problem=true;  
      cout << "ERR>>>" << endl;
    }
    cout<< "dLC/d"<<i<<"(num)="<<dxNumC[i]<<endl;
    cout<< "dLC/d"<<i<<"(ana)="<<dxAnaC[i]<<endl;
    
  }
  if(problem)return 12;
    

  // check if the contraint works
  vector<complex<double> > V;
  vector<pair<int,int> > indices;
  vector<string> names;
  TMatrixD cov(LC.NDim(),LC.NDim());
  for(unsigned int i=0;i<LC.NDim();++i){
    cov[i][i]=1.2;
     if(i<LC.NDim()-1){
       cov[i][i+1]=0.5;
       cov[i+1][i]=0.5;
     }
  }
  cout << "Assume Covariance Matrix:" << endl;
  cov.Print();
  TMatrixD cova(2,2);
  LC.buildCAmps(x,V,indices,names,cov,cova,true);
  cout << "Returned Covariance Matrix from TPWALikelihoodC" << endl;
  cova.Print();

  // compare the two constraint amps:
  double phasediff=std::arg(V[1])-std::arg(V[2]);
  cout << "PhaseDiff= " << phasediff << endl;
  if(phasediff!=0.3)return 2;

  return 0;
 
}

// dummy function needed since we link to but do not use minuit.o
int mnparm(int, string, double, double, double, double) {
    cerr << "this is impossible" << endl;
    throw "aFit";
    return 0;
}
