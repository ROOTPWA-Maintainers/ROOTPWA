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



void mcmcglidingmean(TTree* tree, unsigned int ndim, unsigned int burnin=0, unsigned int win=10){
  double X[ndim];
  double Grad[ndim];
  double Xsum[ndim];
  double Var2[ndim]; // variances
  double Var3[ndim]; // sum_k (x_k - x_mean)^3*dL/dx_k
  double C[ndim];    // Convergence Ratio
  vector<TGraph*> graphs(ndim);
  unsigned int nentries=tree->GetEntriesFast();
  TMultiGraph* mgraph=new TMultiGraph();

  for(unsigned int ip=0;ip<ndim;++ip){
	Xsum[ip]=0;
	Grad[ip]=0;
	Var2[ip]=0;
	Var3[ip]=0;
	graphs[ip]=new TGraph(nentries-win);
	mgraph->Add(graphs[ip]);
       }// end loop over parameters

  tree->SetBranchAddress("dL",&Grad);
    tree->SetBranchAddress("X",&X);

  // calculate first mean value
  for(unsigned int i=0; i<win; i+=1){
    tree->GetEntry(i);
    for(unsigned int ip=0;ip<ndim;++ip){
	Xsum[ip]+=X[ip];
    }// end loop over parameters
  }


  // start sliding:
  for(unsigned int i=win; i<nentries; i+=1){
    // calculate sliding mean
    tree->GetEntry(i-win);
    // remove first point
    for(unsigned int ip=0;ip<ndim;++ip){
      Xsum[ip]-=X[ip];
    }// end loop over parameters
    tree->GetEntry(i);
    // add new point
    for(unsigned int ip=0;ip<ndim;++ip){
      Xsum[ip]+=X[ip];
      // Fill Graph
      graphs[ip]->SetPoint(i-win,i,Xsum[ip]/(double)win);
    }// end loop over parameters

  }
  mgraph->Draw("AP");
}
