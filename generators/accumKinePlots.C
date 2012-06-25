#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TLorentzRotation.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TMath.h"
#include "TPad.h"
#include "TColor.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TString.h"
#include <iostream>
#include <vector>

using namespace std;

#define MCOLOR kOrange-3

TGraphErrors* buildGraph(vector<TH1D*> histo, unsigned int n){
  unsigned int nbins=histo[0]->GetNbinsX();

  TGraphErrors* g=new TGraphErrors(nbins);
  for(unsigned int ib=0;ib<nbins;++ib){
    double x=histo[0]->GetBinCenter(ib+1);
    double val=histo[0]->GetBinContent(ib+1);
    double dx=histo[0]->GetBinWidth(ib+1)*0.5;
    g->SetPoint(ib,x,val);
    // calculate error by looping over all histos
    double sigma=0;
    for(unsigned int ih=0;ih<n;++ih){
      double xi=val-histo[ih]->GetBinContent(ib+1);
      sigma+=xi*xi;
    }
    sigma/=n;
    double sigmastat=histo[0]->GetBinError(ib+1); // add error from event sampling
    sigmastat*=sigmastat;
    g->SetPointError(ib,dx,sqrt(sigma+sigmastat));
  }
  
  return g;

}


// routine to make a nice, accumulated plot
TCanvas* plotAccu(TString histoname, TString xtitle, vector<TString> bins, unsigned int nsamples,TFile* infile, TString plotdir){
  
  TString mcname=histoname+"MC";
  TString dataname=histoname+"Data";
  TString apsname=histoname+"APS";
  TString range=bins[0](2,4)+"."+bins.back()(7,4);


  // vector of samples
  vector<TH1D*> mcplots(nsamples);
  TH1D* dataplot;
  TH1D* apsplot;

  // loop over bins
  for(unsigned int ib=0;ib<bins.size();++ib){
    // data object
    TString name=dataname+bins[ib];
    //cout << "Searching for " << name << endl;
    TH1D* datahisto=(TH1D*)infile->Get(name);
    name=apsname+bins[ib];
    //cout << "Searching for " << name << endl;
    TH1D* apshisto=(TH1D*)infile->Get(name);
    
    if(datahisto==NULL || apshisto==NULL){
      

      return NULL;
    }
    if(ib==0){ // create clone, rename and acccumulate into this
      TString accname=dataname+range;
      dataplot=(TH1D*)datahisto->Clone(accname);
      TString accuapsname=apsname+range;
      apsplot=(TH1D*)apshisto->Clone(accuapsname);
    }
    else {
      dataplot->Add(datahisto);
      apsplot->Add(apshisto);
    }

    // loop over mc samples
    for(unsigned int is=0;is<nsamples;++is){
      // data object
      TString name=mcname+bins[ib]+"_";name+=is;
      //cout << "Searching for " << name << endl;
      TH1D* mchisto=(TH1D*)infile->Get(name);
      if(mchisto==NULL)return 0;
      if(ib==0){ // create clone, rename and acccumulate into this
	TString accname=mcname+range+"_";accname+=is;
	mcplots[is]=(TH1D*)mchisto->Clone(accname);
      }
      else {
	mcplots[is]->Add(mchisto);
      }
      mcplots[is]->SetLineColor(kRed);
    } // end loop over samples
  }// end loop over bins

  TGraphErrors* g=buildGraph(mcplots,nsamples);
  g->SetLineColor(MCOLOR);
  g->SetFillColor(MCOLOR);
  g->SetMarkerColor(MCOLOR);

  TCanvas* c=new TCanvas(histoname+range,histoname+range,10,10,700,700);
  dataplot->Draw("E");
  g->Draw("same p2");
  //for(unsigned int is=0;is<nsamples;++is)mcplots[is]->Draw("same");
  //g->Draw("p2");
  dataplot->Draw("same E");
  dataplot->GetYaxis()->SetRangeUser(0,dataplot->GetMaximum()*1.5);
  dataplot->GetXaxis()->SetTitle(xtitle);

  double xcenter=0.5;
  double ycenter=0.5;
   // compass 2004
  double xc=xcenter+0.05;
  //if(right)xc=xcenter+0.1;
  double yc=ycenter+0.35;
  TLatex* com04=new TLatex(xc,yc,"COMPASS 2004");
  com04->SetNDC();
  com04->SetTextSize(0.05);
  com04->Draw();
  
  // 5 pi on pb
  xc=xcenter+0.05;
  //xc=xcenter+0.1;
  yc=ycenter+0.31;
  TLatex* react=new TLatex(xc,yc,"#pi^{-} Pb #rightarrow #pi^{-}#pi^{+}#pi^{-}#pi^{+}#pi^{-} Pb");
  react->SetNDC();
  react->SetTextSize(0.039);
  react->Draw();
  yc=ycenter+0.27;
  TLatex* dat=new TLatex(xc,yc,"Data vs");
  dat->SetNDC();
  dat->SetTextSize(0.039);
  dat->Draw();
  //yc=ycenter+0.23;
  //TLatex* mc=new TLatex(xc+0.115,yc,"weighted MC");
  //mc->SetNDC();
  //mc->SetTextSize(0.039);
  //mc->SetTextColor(MCOLOR);
  //mc->Draw();
  TPaveText* mc= new TPaveText(xc+0.114,yc-0.009,xc+0.325,yc+0.031,"NDC");
  //mc->SetX1NDC(xc+0.112);mc->SetX2NDC(xc+0.48);mc->SetY1NDC(yc);mc->SetY2NDC(yc+0.04);
  mc->SetBorderSize(0);
  mc->SetTextSize(0.039);
  //mc->SetTextAlign(1);
  mc->SetFillColor(MCOLOR);
  mc->AddText("weighted MC");
  mc->Draw();

  yc=ycenter+0.23;
  TString rangelable("m_{5#pi} #in [");rangelable+=range;rangelable+="] MeV/c^{2}";rangelable.ReplaceAll(".",",");
  TLatex* bin=new TLatex(xc,yc,rangelable);
  bin->SetNDC();
  bin->SetTextSize(0.032);
  bin->Draw();

  range.ReplaceAll(".","_");
  c->SaveAs(plotdir+histoname+range+".eps");

  TCanvas* c2=new TCanvas(apsname+range,apsname+range,10,10,700,700);
  apsplot->Draw();
  apsplot->SetFillColor(MCOLOR);
  apsplot->GetYaxis()->SetRangeUser(0,apsplot->GetMaximum()*1.5);
  apsplot->GetXaxis()->SetTitle(xtitle);
  com04->Draw();
  react->Draw();
  yc=ycenter+0.27;
  TLatex* dat2=new TLatex(xc,yc,"accepted phase-space MC");
  dat2->SetNDC();
  dat2->SetTextSize(0.034);
  dat2->Draw();
  yc=ycenter+0.23;
  bin->Draw();
  c2->SaveAs(plotdir+apsname+range+".eps");

  return c;
}









void accumKinePlots(TString plotsfile, TString outdir){

  TFile* infile=TFile::Open(plotsfile.Data(),"READ");
  Int_t font=132;
  gStyle->SetTextFont(font);
  gStyle->SetLabelFont(font,"xy");
  gStyle->SetTitleFont(font,"xy");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetStripDecimals(1);
  TGaxis::SetMaxDigits(4);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameBorderMode(0);
  
  gROOT->ForceStyle();

  // setup binning
  unsigned int nbins=7;
  vector<vector<TString> > bins(nbins);
  int dm=60; int mstart=1360;
  unsigned int nacc=4;
  int mass=mstart;
  for(unsigned int ib=0;ib<nbins;++ib){
    cout << "Bin " << ib << " : " << endl;
    for(unsigned int ia=0;ia<nacc;++ia){
      TString masslabel("_m");
      masslabel+=mass;masslabel+=".";masslabel+=(mass+dm);
      bins[ib].push_back(masslabel);
      cout << masslabel << endl;
      mass+=dm;
    }
  }

  for(unsigned int bin=0;bin<nbins;++bin){
    plotAccu("hMIsobar","invariant mass of #pi^{-}#pi^{+}#pi^{-}#pi^{+} system (GeV/c^{2})", bins[bin],100,infile,outdir);
    plotAccu("hMIsobar2","invariant mass of #pi^{-}#pi^{+}#pi^{-} system (GeV/c^{2})", bins[bin],100,infile,outdir);
    plotAccu("hMIsobar3","invariant mass of #pi^{-}#pi^{+} system (GeV/c^{2})", bins[bin],100,infile,outdir);
    
    plotAccu("hGJ","cos #theta_{GJ}(#pi^{-}#pi^{+}#pi^{-}#pi^{+})", bins[bin],100,infile,outdir);
  plotAccu("hTY","#phi_{TY}(#pi^{-}#pi^{+}#pi^{-}#pi^{+})", bins[bin],100,infile,outdir);

    plotAccu("hGJ2","cos #theta_{GJ}(#pi^{+}#pi^{-}#pi^{-})" ,bins[bin],100,infile,outdir);
    plotAccu("hHe22Th","cos #theta_{Hel}(#pi^{-}#pi^{+})(#pi^{-}#pi^{+})" ,bins[bin],100,infile,outdir);
    plotAccu("hHe21Th","cos #theta_{Hel}(#pi^{-}#pi^{+})(#pi^{-})" ,bins[bin],100,infile,outdir);
    plotAccu("hHe31Th","cos #theta_{Hel}(#pi^{-}#pi^{+}#pi^{-})(#pi^{+})" ,bins[bin],100,infile,outdir);
    plotAccu("hHe22Phi","#phi_{Hel}(#pi^{-}#pi^{+})(#pi^{-}#pi^{+})" ,bins[bin],100,infile,outdir);
    plotAccu("hHe21Phi","#phi_{Hel}(#pi^{-}#pi^{+})(#pi^{-})" ,bins[bin],100,infile,outdir);
    plotAccu("hHe31Phi","#phi_{Hel}(#pi^{-}#pi^{+}#pi^{-})(#pi^{+})" ,bins[bin],100,infile,outdir);
    //break;
  }  

  

}
