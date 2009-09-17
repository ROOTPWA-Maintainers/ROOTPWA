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
#include "TTree.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include <iostream>
using namespace std;



void checkMCMC(TTree* tree, unsigned int ndim, unsigned int burnin=0, unsigned int step=100){
  double X[ndim];
  double Grad[ndim];
  double Xmean[ndim];
  double Var2[ndim]; // variances
  double Var3[ndim]; // sum_k (x_k - x_mean)^3*dL/dx_k
  double C[ndim];    // Convergence Ratio
  vector<TGraph*> graphs(ndim);
  unsigned int nentries=tree->GetEntriesFast();
  TMultiGraph* mgraph=new TMultiGraph();
  
  for(unsigned int ip=0;ip<ndim;++ip){
	Xmean[ip]=0;
	Grad[ip]=0;
	Var2[ip]=0;
	Var3[ip]=0;
	graphs[ip]=new TGraph((nentries-burnin)/step);
	mgraph->Add(graphs[ip]);
       }// end loop over parameters

  tree->SetBranchAddress("dL",&Grad);
  tree->SetBranchAddress("X",&X);

  
 

  for(unsigned int i=burnin; i<nentries; i+=step){
    // calculate mean
    for(unsigned int j=burnin; j<=i; ++j){
      tree->GetEntry(j);
      for(unsigned int ip=0;ip<ndim;++ip){
	Xmean[ip]+=X[ip];
      }// end loop over parameters
    }
    for(unsigned int ip=0;ip<ndim;++ip){
      Xmean[ip]/=(double)(i-burnin+1);
      //cout << Xmean[ip] << endl;
    }// end loop over parameters
    for(unsigned int j=burnin; j<=i; ++j){
      tree->GetEntry(j);
      // caluclate variance terms
      for(unsigned int ip=0;ip<ndim;++ip){
	double diff=X[ip]-Xmean[ip];
	///cout << "diff=" <<diff<< endl;
	Var2[ip]+=diff*diff;
	Var3[ip]+=diff*diff*diff*Grad[ip];
      }// end loop over parameters
    }
    
    double Rmin=10000000;
    double Rmax=-10000000;
    double Rmean=0;
    unsigned int Imin, Imax;


    if(i>burnin){
      for(unsigned int ip=0;ip<ndim;++ip){
	if(Var2[ip]==0) continue;
	double R=Var3[ip]/(3*Var2[ip]);
	graphs[ip]->SetPoint(i-burnin,i-burnin,R);
	if(R<Rmin){Rmin=R;Imin=ip;}
	if(R>Rmax){Rmax=R;Imax=ip;}
	Rmean+=R;
	
	Xmean[ip]=0;
	Grad[ip]=0;
	Var2[ip]=0;
	Var3[ip]=0;
      }// end loop over parameters
      Rmean/=(double)ndim;
      cout << "Rmin("<<Imin<<")="<<Rmin
	   <<" ... Rmean="<<Rmean
	   <<" ... Rmax("<<Imax<<")="<<Rmax<<endl;
    }

  }

  mgraph->Draw("AP");

}
