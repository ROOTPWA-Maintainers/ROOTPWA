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
#include <time.h>
#include "utilities.h"
#include "TBranch.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TComplex.h"
#include "TRandom3.h"
#include "TFitBin.h"
#include "TPWALikelihood.h"
#include "TPWALikelihoodC.h"

using namespace std;


int lineno = 1; // global variables needed for lex (not understood)
char *progname;

int main(int argc, char** argv){
  
  time_t seed=1254410383;
  //uint seed=1254410383;
  time(&seed);
  gRandom->SetSeed(seed);

  cout << "Seed=" << (uint)seed << endl;

  TPWALikelihood<double> L;
  TPWALikelihoodC LC;

  L.UseNormalizedAmps(true);
  LC.UseNormalizedAmps(true);


  //if(quiet)L.SetQuiet();
  L.SetWavelist(argv[1]);
  int rank=2; // TODO: make this an option
  L.SetRank(rank);
  LC.Init(argv[1],rank,"norm.int","norm.int",20000);

  cout<<"L.NDIM()="<<L.NDim()<<"   LC.NDim()="<<LC.NDim()<<endl;
  if(L.NDim()!=LC.NDim())return 1;

  

  LC.LoadAmplitudes();

  double par[LC.NDim()];
  double par2[LC.NDim()];
  for(unsigned int i=0;i<LC.NDim();++i){
    par[i]=(double)i;
    par2[i]=0;
  }
  LC.partoamp(par);
  LC.Print();
  LC.amptopar(par2);
  for(unsigned int i=0;i<LC.NDim();++i){
    if(par[i]!=par2[i])return 5;
  }


  //L.SetMaxSampDL(1000);
  L.LoadIntegrals("norm.int","norm.int");
  L.LoadAmplitudes();
 


  // 12 parameters + flat
  //double x[13]={0.52707,0.21068,-0.604365,0.17596,-0.216668,-0.0990815,-0.348459,0.208961,0.02,0.03,0.44,0,0};
  double x[L.NDim()]; for(unsigned int i=0;i<L.NDim();++i)x[i]=gRandom->Uniform(-5,5);
  //string a[13]={"a","b","c","d","e","f","g","h","i","j","k","l","flat"};
  double LL=L(x);
  double LLC=LC(x);
    
  std::cout<<"L(x)="<< maxPrecision(LL)<<std::endl;
  std::cout<<"LC(x)="<< maxPrecision(LLC)<<std::endl;
  for(unsigned int i=0;i<L.NDim();++i) cout << x[i] << endl;
  if(fabs(LL-LLC)>1E-8) return 10;

  

 // check derivatives:
  // first numerical:
  double L1=L(x);
  double LC1=LC(x);
  double h=1E-4;
  double dxNum[L.NDim()];
  double dxAna[L.NDim()];
  double dxNumC[L.NDim()];
  double dxAnaC[L.NDim()];
  bool problem=false;
  std::cout << " *** Check Numerical Gradient Old/New Likelihood" << std::endl;
  for(unsigned int i=0; i<L.NDim();++i){
    x[i]+=h;
    double L2=L(x);
    double LC2=LC(x);
    dxNum[i]=(L2-L1)/h;
    dxNumC[i]=(LC2-LC1)/h;
    if(2*fabs(dxNum[i]-dxNumC[i])/(fabs(dxNum[i])+fabs(dxNumC[i]))>0.01){
      cout << "Numerical dL/d" << i << "=" << maxPrecision(dxNum[i])
	   << " NOT EQUAL TO (Threshold 1%) "
           << "dLC/d" << i << "=" << maxPrecision(dxNumC[i]) << endl;
      //problem=true;
    }
    x[i]-=h;
  }
  if(problem)return 11;
     // Then analytic
  std::cout << " *** Check Numerical vs. Analytical Gradient" << std::endl;
  double F;
  double FC;
  L.FdF(x,F,dxAna);
  LC.FdF(x,FC,dxAnaC);
  for(unsigned int i=0; i<L.NDim();++i){
    if(2*fabs(dxNum[i]-dxAna[i])/(fabs(dxNum[i])+fabs(dxAna[i]))>0.0001){
      //problem=true;
      cout << "ERR>>>" << endl;
    }
    cout<< "dL/d"<<i<<"(num)="<<maxPrecision(dxNum[i])<<endl;
    cout<< "dL/d"<<i<<"(ana)="<<maxPrecision(dxAna[i])<<endl;

    if(2*fabs(dxNumC[i]-dxAnaC[i])/(fabs(dxNumC[i])+fabs(dxAnaC[i]))>0.0001){
      //problem=true;  
      cout << "ERR>>>" << endl;
    }
    cout<< "dLC/d"<<i<<"(num)="<<maxPrecision(dxNumC[i])<<endl;
    cout<< "dLC/d"<<i<<"(ana)="<<maxPrecision(dxAnaC[i])<<endl;
    if(2*fabs(dxAna[i]-dxAnaC[i])/(fabs(dxAna[i])+fabs(dxAnaC[i]))>0.0001){
      cout << "Analytical dL/d" << i << "=" << maxPrecision(dxNum[i])
	   << " NOT EQUAL TO (Threshold 0.01%) "
           << "dLC/d" << i << "=" << maxPrecision(dxNumC[i]) << endl;
      problem=true;
    }
  }
  if(problem)return 12;
    
  // value extraction
  vector<complex<double> > V;
  vector<pair<int,int> > indices;
  vector<string> names;
  TMatrixD cov(L.NDim(),L.NDim());
  for(unsigned int i=0;i<L.NDim();++i){
    cov[i][i]=1.2;
    if(i<L.NDim()-1){
      cov[i][i+1]=0.5;
      cov[i+1][i]=0.5;
    }
  }
  cout << " *** Check Covariance Matrix: \n Assume Covariance Matrix:" << endl;
  cov.Print();
  vector<complex<double> > VC;
  vector<pair<int,int> > indicesC;
  vector<string> namesC;
  TMatrixD covC(LC.NDim(),LC.NDim());
  TMatrixD covaC(2,2);
  L.buildCAmps(x,V,indices,names,true);
  LC.buildCAmps(x,VC,indicesC,namesC,cov,covaC,true);
  cout << "#Names: "<< names.size() << endl;
  cout << "#NamesC: "<< namesC.size() << endl;
  for(unsigned int in=0;in<names.size();++in){
    cout << "L : " <<names[in] << endl;
    cout << "LC: " << namesC[in] << endl;
  }
  cout << "Returned Covariance Matrix from TPWALikelihoodC" << endl;
  covaC.Print();
  cout.flush();
  // test some random entries in the matrices
  for(unsigned int i=0;i<100;++i){
    int a=gRandom->Integer(V.size());
    int b=gRandom->Integer(V.size());
    cout<<"a="<<a<<"   b="<<b<<endl;
    double errRE=0;
    double errIM=0;
    if(indices[a].first>-1 && indices[b].first>-1 )
      errRE=cov[indices[a].first][indices[b].first];
    if(indices[a].second>-1 && indices[b].second>-1)
      errIM=cov[indices[a].second][indices[b].second];

    cout << "L: Compare covariances\n"
	 <<"cov("<< names[a] << "," << names[b] << ")=" 
	 <<" ("<<errRE<<","<<errIM<<")"<<endl
	 <<"with"<<endl;
    
    
    // find corresponding waves;
    unsigned int aC=0;
    unsigned int bC=0;
    int bothfound=0;
    for(unsigned int iC=0;iC<namesC.size();++iC){
      if(names[a]==namesC[iC]){aC=iC;++bothfound;}
      if(names[b]==namesC[iC]){bC=iC;++bothfound;}
    }
    if(bothfound!=2){
      cout << "Did not find corresponding constraint amplitude pair "<< endl;
      return 14;
    }

    pair<int, int> ia,ib;
    ia=indicesC[aC];
    ib=indicesC[bC];
    double errCIM=covaC[ia.second][ib.second];
    double errCRE=covaC[ia.first][ib.first];
    
    cout <<"covC("<< namesC[aC] << "," << namesC[bC] << ")="
	 <<" ("<<errCRE<<","<<errCIM<<")"<<endl;
    
    if(errRE!=errCRE || errIM!=errCIM)return 14;
  }


  return 0;
}

// dummy function needed since we link to but do not use minuit.o
int mnparm(int, string, double, double, double, double) {
    cerr << "this is impossible" << endl;
    throw "aFit";
    return 0;
}
