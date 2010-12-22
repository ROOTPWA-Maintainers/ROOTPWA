#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TList.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TKey.h"
#include "TStyle.h"
#include "TH2D.h"
#include "TMultiGraph.h"

#include <vector>
#include <iostream>
#include <cmath>

using namespace std;


void roottops(const TString& infilename, const TString& normfilename=""){

  gStyle->SetFrameFillStyle(4000);
  
  TString name=infilename;
  name.ReplaceAll(".root","");

  const int    nmbPadsPerCanvMin = 6;            // minimum number of pads each canvas is subdivided into

  TFile* infile=TFile::Open(infilename,"READ");

  TFile* normfile=0;
  if(normfilename.Length()>1)normfile=TFile::Open(normfilename,"READ");

  TList* keylist=infile->GetListOfKeys();
  unsigned int num=keylist->GetEntries();

  cout << num << " entries in keylist." << endl;

  const int nmbPadsHor     = (int)floor(sqrt(nmbPadsPerCanvMin));
  const int nmbPadsVert    = (int)ceil((double)nmbPadsPerCanvMin / (double)nmbPadsHor);
  const int nmbPadsPerCanv = nmbPadsHor * nmbPadsVert;
  

  vector<TCanvas*> canvases;
  TCanvas* current=NULL;
  unsigned int padcounter=0;
  for(unsigned int i=0;i<num;++i){ // loop over keys

    if(padcounter%(nmbPadsPerCanv+1) ==0){ // create new canvas
      TString cname=name;cname+=canvases.size();
      current=new TCanvas(cname,cname, 10, 10, 800, 1000);
      current->Divide(nmbPadsHor,nmbPadsVert);
      canvases.push_back(current);
      padcounter=1;
    }
    current->cd(padcounter);
    
    // TString type("TMultiGraph");
    TString type("TH2D");
    

    if(TString(((TKey*)keylist->At(i))->GetClassName())==type){
      cout << "Found " << type << endl; 
      TString keyname(((TKey*)keylist->At(i))->GetName());
      if(keyname.Contains("PHI"))continue;
      if(type=="TH2D"){
	((TKey*)keylist->At(i))->ReadObj()->Draw("COLZ");
	if(keyname.Contains("rho1")){
	  ((TH2D*)(((TKey*)keylist->At(i))->ReadObj()))->GetYaxis()->SetRangeUser(-100,15000);}
      }
      else ((TKey*)keylist->At(i))->ReadObj()->Draw("AP");
      // plot normgraphs if available
      if(normfile!=0){
	TString name=((TKey*)keylist->At(i))->ReadObj()->GetName();
	if(name.Contains("amp")){
	  TH2D* h=(TH2D*)((TKey*)keylist->At(i))->ReadObj();
	  double xmin=h->GetXaxis()->GetXmin()*1000;
	  double xmax=h->GetXaxis()->GetXmax()*1000;
	  cout << xmin << "   " << xmax << endl;
	  TString wavename=name(1,name.Length());
	  TGraph* g=(TGraph*)normfile->Get(wavename);
	  
	  TPad* p2=new TPad("pad","PS",gPad->GetX1(),gPad->GetY1(),gPad->GetX2(),gPad->GetY2());
	  //p2->RangeAxis(xmin,0,xmax,1);
	  p2->SetFillStyle(4000);
	  p2->Draw();
	  p2->cd();
	  g->Draw("APC");
	  g->GetXaxis()->SetRangeUser(xmin,xmax);
	  g->GetXaxis()->SetTitle("Mass (GeV/c^{2}");
	}
      }
      ++padcounter;
    }

  }// end loop over keys
  int nc1=canvases.size();
  // loop again to get 3-plots for phases
  for(unsigned int i=0;i<num;++i){ // loop over keys
    TString keyname(((TKey*)keylist->At(i))->GetName());
    if(keyname.Contains("PHI")){
      
      TMultiGraph* g=(TMultiGraph*)((TKey*)keylist->At(i))->ReadObj();
      if(g==NULL)continue;
      cout << "found Phase Graph!" << TString(((TKey*)keylist->At(i))->GetName()) << endl;
      int divider=keyname.Index("---");
      TString key1=keyname(3,divider-3);
      TString key2=keyname(divider+6,keyname.Length());
      cout << key1 << endl;
      cout << key2 << endl;
      TMultiGraph* g1=(TMultiGraph*) infile->Get(key1);
      TMultiGraph* g2=(TMultiGraph*) infile->Get(key2);
      if(g1==NULL || g2==NULL) continue;
      current=new TCanvas("C"+keyname,"C"+keyname, 10, 10, 800, 1000);
      current->Divide(1,3);
      canvases.push_back(current);
      current->cd(1);
      gPad->SetGridy();
      // plot phasegraph
      
      if(g->GetListOfGraphs()->GetSize()==0){}
      else {

	g->Draw("AN");
	//g->GetYaxis()->Set(8,-720,720);
	//g->GetYaxis()->SetRangeUser(-720,720);
	g->GetYaxis()->SetRangeUser(-200,200);
	
	g->Draw("A");
	}
      // get waves
    
      current->cd(2);if(g1!=NULL)g1->Draw("APC");
      current->cd(3);if(g2!=NULL)g2->Draw("APC");
    } // endif found PHI hist
  }


  TString psFileName = name; psFileName += ".ps";
  TCanvas      dummyCanv("dummy", "dummy");
  TString option="Portrait";
  dummyCanv.Print((psFileName + "["),option);
  for(unsigned int ic=0;ic<canvases.size();++ic){
    if(ic>=nc1)option="Portrait";
    canvases[ic]->Print(psFileName,option);
    delete canvases[ic];
  }
 dummyCanv.Print((psFileName + "]"),option);
 gSystem->Exec(("gv " + psFileName));
}





//     const string psFileName = outPath + "waveIntensities.ps";
//     TCanvas      dummyCanv("dummy", "dummy");
//     dummyCanv.Print((psFileName + "[").c_str());
//     for (map<string, TCanvas*>::iterator i = canvases.begin(); i != canvases.end(); ++i) {
//       i->second->Print(psFileName.c_str());
//       delete i->second;
//       i->second = 0;
//     }
//     dummyCanv.Print((psFileName + "]").c_str());
//     gSystem->Exec(("gv " + psFileName).c_str());
//   }

//   return wavePads;




// }
