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
#include "TFile.h"
#include "TTree.h"
#include "TFileCollection.h"
#include "TChain.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TFitBin.h"
#include <string>
#include <vector>
#include <iostream>
using namespace std;

TChain*
loadFit(TString select){
  cout << "Loading " << select << endl;
  TChain* pwa=new TChain(select,select);
  TString kap=select;
  kap+="/pwa";
  pwa->Add(kap);

  pwa->Scan("nwaves()");


  return pwa;


}



void
compareFits(){

  vector<TString> files;
files.push_back("/afs/e18/compass/analysis/sneubert/Q3PiData/FITS/testseries/1020.1220.test6.grad.result.root");
files.push_back("/afs/e18/compass/analysis/sneubert/Q3PiData/FITS/testseries/1020.1220.test7.grad.result.root");
  files.push_back("/afs/e18/compass/analysis/sneubert/Q3PiData/FITS/testseries/1020.1220.test8.grad.result.root");
  files.push_back("/afs/e18/compass/analysis/sneubert/Q3PiData/FITS/testseries/1020.1220.test9.grad.result.root");
  files.push_back("/afs/e18/compass/analysis/sneubert/Q3PiData/FITS/testseries/1020.1220.test10.grad.result.root");
files.push_back("/afs/e18/compass/analysis/sneubert/Q3PiData/FITS/testseries/1020.1220.test11.grad.result.root");
//files.push_back("/afs/e18/compass/analysis/sneubert/Q3PiData/FITS/fit5/rank2/1020.1220.result.root");

  int b=1;


  // load chains
  vector<TChain* > chains;
  for(int k=0;k<files.size();++k){
    chains.push_back(loadFit(files[k]));
  }

  // make error plots
  
  unsigned int nc=chains.size();

  TFitBin* bin=new TFitBin();
  chains[nc-1]->SetBranchAddress("fitbin",&bin);
  chains[nc-1]->GetEntry(b);

  // now build graphs;
  int ng=bin->nwaves();
  cout << "Installing "<<ng<<" error graphs:"<<endl;
  TMultiGraph* mg=new TMultiGraph();
  vector<TGraph*> graphs(ng);
  vector<TString> names(ng);
  for(int ig=0;ig<ng;++ig){
    names[ig]=bin->waveDesignator(ig);
    cout<<names[ig]<<endl;
    TGraph* gerr=new TGraph(nc);
    graphs[ig]=gerr;
    gerr->SetMarkerStyle(22);
    if(names[ig].Contains("1++0+rho770_01"))gerr->SetLineColor(kRed);
    //gerr->SetTitel(names[ng]);
    gerr->SetName(names[ig]);
    mg->Add(gerr);
  }
  

  // loop through fits
  for(unsigned int i=0; i<nc; ++i){
    TChain* mychain=chains[i];
    mychain->SetBranchAddress("fitbin",&bin);
    mychain->GetBranch("fitbin")->SetAutoDelete();
    mychain->GetEntry(b);
    for(unsigned int ig=0;ig<ng;++ig){
      // get error of wave in first bin
      double intens=bin->intens(names[ig]);
      if(intens>0)graphs[ig]->SetPoint(i,bin->nwaves(),bin->err(names[ig])/intens);
      else graphs[ig]->SetPoint(i,bin->nwaves(),0);
   }
  }
  
  mg->Draw("APL");
  cout << "there should be "<<nc<<" wavesets with these numbers of waves:" << endl;
  for(unsigned int i=0; i<nc; ++i){
    TChain* mychain=chains[i];
    mychain->GetEntry(0);
    cout << bin->nwaves() << endl;
  }
  
}
