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
#include <vector>
#include "../TFitBin.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TPostScript.h"
#include "TMath.h"
#include "TPad.h"
#include "TGraph.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TROOT.h"
#include "TH1.h"

void plotwaves(TTree* pwa, bool dops=false, TString Outpath="./", bool mcmc=false){

TFitBin* bin=new TFitBin();

pwa->SetBranchAddress("fitbin",&bin);


// determin relative intensities


  pwa->SetBranchAddress("fitbin",&bin);
  
  pwa->GetEntry(0);
  int nw=bin->nwaves();
  std::vector<double> intens(nw,0);

  int n=pwa->GetEntries();
  int usedbins=0;
  double I=0;  // total intensity;

  for(int i=0;i<n;++i){
    pwa->GetEntry(i);
    double Ib=bin->intens();
    I+=Ib;
    if(Ib>500){
      for(int j=0;j<nw;++j){
	// calculate relative intensity
	intens[j]+=bin->intens(j);
      } // end loop over waves
      ++usedbins;
    }
  }// end loop over mnass bins
 

  // sort the wavenames
  std::map<TString,double> waves;
  std::map<double, TString> waves2;
  for(int j=0;j<nw;++j){/// loop over waves
    // normalize relatve intensities to total intensity!
    waves[bin->waveDesignator(j)]=intens[j]/I;
    waves2[intens[j]/I]=bin->waveDesignator(j);
  }// end loop over waves
  std::map<TString,int> wstrength;
  std::map<double, TString>::iterator sit=waves2.begin();
  int count=1;
  while(sit!=waves2.end()){
    wstrength[sit->second]=count++;
    ++sit;
  }




pwa->GetEntry(0);
nw=bin->nwaves();

cout << "Number of waves: "<<nw<<endl;

// int n1=5;//(int)TMath::Ceil(TMath::Sqrt(nw));
//int n2=4;//(int)TMath::Ceil((double)nw/(double)n1);
// int n=n1*n2;
// int nc=(int)TMath::Ceil((double)nw/(double)n); 
//cout << "n1="<<n1<<"  n2="<<n2<<endl;

// cout << can.size() << endl;

// sort plots according to wavename
std::map<string,int> wavemap;
std::map<TString,TCanvas*> jpcmap;
 std::map<TString,int> jpcnmap; 
 std::map<TString,int> jpccounter; 

for(int i=0;i<nw;++i){
  TString title=bin->waveDesignator(i);
  wavemap[title.Data()]=i;
  cout << title <<endl;
  TString jpc=title(2,3);
  cout << "JPC==" << jpc << endl;
  jpcnmap[jpc]+=1;
}

 int npercmin=4;
 int n1=(int)TMath::Ceil(TMath::Sqrt(npercmin));
 int n2=(int)TMath::Ceil((double)npercmin/(double)n1);
 int nperc=n1*n2;

std::map<TString,int>::iterator nit=jpcnmap.begin();
while(nit!=jpcnmap.end()){
  jpccounter[nit->first]=0;
  int n=nit->second;
  // number of canvases for this jpc
  int nc=(int)TMath::Ceil((double)n/(double)nperc);
  for(int ic=0; ic<nc;++ic){
    TString name(nit->first);
    name+=ic;
    name.ReplaceAll(" ","");
    if(jpcmap[name]==NULL){
      jpcmap[name]=new TCanvas(name,name,10,10,1000,800);
      jpcmap[name]->Divide(n1,n2);
    }
  }
  ++nit;
}

TString filename=Outpath;filename+="waves.root";
TFile* file=TFile::Open(filename,"RECREATE");
 gROOT->cd();

std::map<double,TString>::iterator it=waves2.end();
int i=0;
 TCanvas* myc;
while(it!=waves2.begin()){
  --it;
  // select canvas
  TString title=it->second;
  TString jpc=title(2,3);
  TString cname(jpc);
  int njpc=jpcnmap[jpc];
  cout << njpc << " waves for " << jpc << endl;
  cname+=(int)TMath::Floor((double)jpccounter[jpc]/(double)nperc);
  cname.ReplaceAll(" ","");
  cout << " Selecting Canvas "<< cname << endl;
  myc=jpcmap[cname];
  int ipad=(jpccounter[jpc]%nperc)+1;
  myc->cd(ipad);
  jpccounter[jpc]+=1;
  int id=wavemap[title.Data()];
  TString com;
  if(!mcmc)com=".x plotwave.C(\"";
  else com=".x plotmcmc.C(\"";
  com+=title;
  com+="\",\"\",\"";
  com+=title;
  com+="\");";
  gROOT->ProcessLine(com);
  TGraph* g=(TGraph*)gPad->FindObject(title);
  double m=g->GetMaximum();
  g->SetName(title);
  file->cd();
  g->Write();
  gROOT->cd();
  // draw additional info
  TString info("I=");
  double I=waves[title]*1000;
  I=TMath::Floor(I+0.5);
  I/=10;
  info+=I;
  info+="%(";
  info+=(wstrength.size()-wstrength[title]+1);
  info+=")";
  info.ReplaceAll(" ","");
  info.ReplaceAll("%(","% (");
  TLatex* text=new TLatex(1.0,0.9*m,info);
  text->Draw();


  ++i;
}

 TList* Hlist=gDirectory->GetList();
 Hlist->Remove(pwa);

 file->cd();

    int nobj=Hlist->GetEntries();
    std::cout<<"Found "<<nobj<<" Objects in HList"<<std::endl;
    Hlist->Print();
    Hlist->Write();


  
   
 file->Close();
 gROOT->cd();



// now plot this to ps file
 if(dops){
   std::map<TString,TCanvas*>::iterator cai=jpcmap.begin();
   Int_t type = 112;
   // create a postscript file and set the paper size
   TString com=Outpath;
   com+="waves.ps";
   TPostScript ps(com,type);
   //ps.Range(16,24);  //set x,y of printed page
   while(cai!=jpcmap.end()){
     ps.NewPage();
     cai->second->Update();
     delete cai->second;
     ++cai;
   }
   ps.Close();
   TString com2="gv ";com2+=com;
   gSystem->Exec(com2);
   
 }



}
