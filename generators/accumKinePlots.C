#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TLatex.h"

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
TCanvas* plotAccu(TString histoname, TString xtitle, vector<TString> bins, unsigned int nsamples,TFile* infile){
  
  TString mcname=histoname+"MC";
  TString dataname=histoname+"Data";
  TString range=bins[0](2,4)+"."+bins.back()(7,4);


  // vector of samples
  vector<TH1D*> mcplots(nsamples);
  TH1D* dataplot;
  
  // loop over bins
  for(unsigned int ib=0;ib<bins.size();++ib){
    // data object
    TString name=dataname+bins[ib];
    //cout << "Searching for " << name << endl;
    TH1D* datahisto=(TH1D*)infile->Get(name);
    if(datahisto==NULL)return NULL;
    if(ib==0){ // create clone, rename and acccumulate into this
      TString accname=dataname+range;
      dataplot=(TH1D*)datahisto->Clone(accname);
    }
    else {
      dataplot->Add(datahisto);
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
  g->SetLineColor(kOrange);
  g->SetFillColor(kOrange);
  g->SetMarkerColor(kOrange);

  TCanvas* c=new TCanvas(histoname+range,histoname+range,10,10,800,800);
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
  double xc=xcenter+0.1;
  //if(right)xc=xcenter+0.1;
  double yc=ycenter+0.35;
  TLatex* com04=new TLatex(xc,yc,"COMPASS 2004");
  com04->SetNDC();
  com04->SetTextSize(0.05);
  com04->Draw();
  
  // 5 pi on pb
  xc=xcenter+0.1;
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
  TLatex* mc=new TLatex(xc+0.095,yc,"weighted MC");
  mc->SetNDC();
  mc->SetTextSize(0.039);
  mc->SetTextColor(kOrange);
  mc->Draw();
  yc=ycenter+0.23;
  TString rangelable("m_{5#pi} #in [");rangelable+=range;rangelable+="] MeV/c^{2}";rangelable.ReplaceAll(".",",");
  TLatex* bin=new TLatex(xc,yc,rangelable);
  bin->SetNDC();
  bin->SetTextSize(0.032);
  bin->Draw();




  return c;
}









void accumKinePlots(TString plotsfile, unsigned int bin){

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

  plotAccu("hMIsobar","invariant mass of #pi^{-}#pi^{+}#pi^{-}#pi^{+} system (GeV/c^{2})", bins[bin],100,infile);
  plotAccu("hMIsobar2","invariant mass of #pi^{-}#pi^{+}#pi^{-} system (GeV/c^{2})", bins[bin],100,infile);
plotAccu("hMIsobar3","invariant mass of #pi^{-}#pi^{+} system (GeV/c^{2})", bins[bin],100,infile);

  plotAccu("hGJ","cos #theta_{GJ}^{14}", bins[bin],100,infile);
  plotAccu("hGJ2","cos #theta_{GJ}^{23}" ,bins[bin],100,infile);
  plotAccu("hHe22Th","cos #theta_{Hel}^{22}" ,bins[bin],100,infile);
  plotAccu("hHe21Th","cos #theta_{Hel}^{21}" ,bins[bin],100,infile);
  plotAccu("hHe31Th","cos #theta_{Hel}^{31}" ,bins[bin],100,infile);
  plotAccu("hHe22Phi","#phi_{Hel}^{22}" ,bins[bin],100,infile);
  plotAccu("hHe21Phi","#phi_{Hel}^{21}" ,bins[bin],100,infile);
  plotAccu("hHe31Phi","#phi_{Hel}^{31}" ,bins[bin],100,infile);
  

  

}
