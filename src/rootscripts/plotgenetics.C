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

#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "TTree.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TPostScript.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"
#include "TPad.h"
#include "TGraph.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"

#include "TPrincipal.h"


#include "wset.h"



void plotgenetics(TString resultfile, TString wavepool){

  // load list of all possible waves in this run
  
  map<TString,unsigned int> waveindex;

  ifstream file(wavepool.Data());
  string line;
  // note that the index is starting at 1 !!!
  unsigned int index=1;
  while(file.good()){
    getline(file,line);
    unsigned int pos=line.find(" ");
    string name=line.substr(0,pos);
    double thres;
    if(pos<line.length()){
      thres=atof(line.substr(pos,line.length()).c_str());
    } else thres=0;
    if(line.length()>1){
      waveindex[name.c_str()]=index;
      ++index;
    }
  }
 
  unsigned int d=waveindex.size();
  cout << d << " waves in wavepool." << endl;

  // read results
  vector<double> fitness;


  // setup matrix for mean value
  TMatrixD mean(d,1);
  // list of chromosomes:
  vector<TMatrixD> csomes;

  // read in the wavelist of the individual fits
  // and analyse the waveset of this fit
  ifstream results(resultfile.Data());
  while(results.good()){
    getline(results,line);
    unsigned int pos=line.find(" ");
    string name=line.substr(0,pos);
    double fit=0;
    if(pos<line.length()){
      fit=atof(line.substr(pos,line.length()).c_str());
    } 
    if(line.length()>1){
      fitness.push_back(fit);
      set<wsetentry> myset;
      readWavelist(myset,name);
      TMatrixD ch(d,1);
      // build chromosome
      unsigned int iw=0;
      set<wsetentry>::iterator it=myset.begin();
      while(it!=myset.end()){
	unsigned int index=waveindex[it->name];
	ch(index,0)=1;
	++it;++iw;
      }
      csomes.push_back(ch);
      ch.Print();
      mean+=ch;
      
    }
  }// end loop over results files
  

    


  //pca.Test();

  
  // renormalize mean
  mean*=1./(double)csomes.size();



  

  // calculate covariances:
  TMatrixD cov(d,d);
  for(unsigned int i=0;i<csomes.size();++i){
    TMatrixD u=csomes[i]-mean;
    TMatrixD c(u,TMatrixD::kMultTranspose,u);
    cov+=c;
  }

   cov*=1./(double)csomes.size();

   //cov.Print();
   
   // calculate eigenvectors

   TMatrixDSym sym; sym.Use(cov.GetNrows(),cov.GetMatrixArray());
   TMatrixDSymEigen eigen(sym);
   TMatrixD EigenVectors = eigen.GetEigenVectors();
   TVectorD EigenValues  = eigen.GetEigenValues();

   EigenValues.Print();

   // TMatrixD uc=cov.EigenVectors(lamda);
  
   // projections onto first and second component:
   TVectorD u1=TMatrixDColumn(EigenVectors,0);
   TVectorD u2=TMatrixDColumn(EigenVectors,1);

   TCanvas* c=new TCanvas("c","PCA",10,10,1200,800);
   c->Divide(2,1);
   c->cd(1);
   TH2D* hpca=new TH2D("hpca","PCA",100,-5,5,100,-5,5);
   hpca->Draw();
   
   unsigned int n=csomes.size();
   TMatrixD* dist=new TMatrixD(n,n);

   for(unsigned int i=0;i<n;++i){
     TVectorD x=TMatrixDColumn(csomes[i],0);
     for(unsigned int j=0;j<n;++j){
       TVectorD y=TMatrixDColumn(csomes[j],0);
       TVectorD dis=x-y;
       (*dist)(i,j)=dis.Norm1();
     }
     double x1=x*u1;
     double x2=x*u2;
     cout << x1 << "    "  << x2 << endl;
     TString mark;
     mark+=i;
     int col=4-i;
     if(col<-10)col=-10;
     col+=kBlue;
     TLatex* maker=new TLatex(x1,x2,mark);
     maker->SetTextColor(col);
     maker->SetTextSize(0.02);
     maker->Draw();
   }

   c->cd(2);
   dist->DrawClone("COL");
   dist->DrawClone("TEXT SAME");
   // hpca->Draw();

}
