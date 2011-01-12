#include "TFile.h"
#include "TList.h"
#include "TKey.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TAxis.h"

#include <iostream>
using namespace std;

void plotMassDepFitResult(TString infilename){
  TFile* infile=TFile::Open(infilename);

  TList* keylist=infile->GetListOfKeys();
  unsigned int num=keylist->GetEntries();

  // loop over keys and count waves
  vector<TString> wavenames;

  for(unsigned int i=0;i<num;++i){
    TString keyname(((TKey*)keylist->At(i))->GetName());
    if(keyname.Contains("dPhi") || keyname.Contains("Re") || keyname.Contains("Im"))continue;
    wavenames.push_back(keyname);
  }


  unsigned int nwaves=wavenames.size();
  std::cout << nwaves << " waves used in fit" << endl;

  TCanvas* c=new TCanvas("c","Spin Density Matrix",10,10,1000,1000);
  c->Divide(nwaves,nwaves,0,0);
  
  // do plotting
  for(unsigned int ip=0;ip<nwaves;++ip){
    for(unsigned int jp=ip;jp<nwaves;++jp){
      c->cd(jp+ip*nwaves+1);
      if(ip==jp){
	TGraph* g=(TGraph*)infile->Get(wavenames[ip]);
	g->Draw("AP");
      }
      else{
	TString key="dPhi_"+wavenames[ip]+"---"+wavenames[jp];
	TMultiGraph* g=(TMultiGraph*)infile->Get(key);
	if(g!=NULL){
	  g->Draw("AN");
	  TAxis* a=g->GetYaxis();
	  if(a!=NULL)a->SetRangeUser(-240,240);
	  g->Draw("A");
	}

      }
    } // end inner loop
  } // end plotting loop

}
