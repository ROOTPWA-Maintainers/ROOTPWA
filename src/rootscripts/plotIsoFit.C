#include "TString.h"
#include "TCanvas.h"
#include <iostream>
#include <cstdlib>
#include <complex>
#include <vector>
#include <string>

#include "TMatrixT.h"

#include "TTree.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TAxis.h"

#include "TGraphErrors.h"
#include "TFile.h"
#include "fitResult.h"

using namespace std;
using namespace rpwa;

void plotIsoFit(TString infilename, unsigned int bin=1) {
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

   // load fitResult tree
  TFile* infile=TFile::Open(infilename,"READ");
  if(infile==NULL){
    cerr << "File " << infilename << " not found!"<< endl;
    return;
  }
  TTree* pwa=(TTree*)infile->FindObjectAny("pwa");
  if(pwa==NULL){
    cerr << "Tree not found!"<< endl;
    return;  
  }

  fitResult* res=0;//new fitResult();
  string branchname("fitResult_v2");
  if(pwa->FindBranch(branchname.c_str())==NULL){
    cerr << "Invalid branch ."<<branchname<<endl;
    throw;
   }
  

  unsigned int nbins=13;
  TGraphErrors* gArgand=new TGraphErrors(0);//nbins);
  TGraphErrors* gIntens=new TGraphErrors(0);//nbins);
  TGraphErrors* gPhase=new TGraphErrors(0);//nbins);
 
  TGraph* gPS=new TGraph(nbins);

  vector<TString> labels(nbins);
  vector<double> xlabels(nbins);
  vector<double> ylabels(nbins);


  pwa->SetBranchAddress(branchname.c_str(),&res);
  cout << "Entries: " <<pwa->GetEntries() << endl;
  pwa->GetEntry(bin);

  cout << "5pi Mass=" << res->massBinCenter() << endl;

 std::vector<std::string> names=res->prodAmpNames();
 for(unsigned int i=0;i<names.size();++i)cout << names[i]<<endl;

  // now loop through know mass binning (hardcoded here)
  // to collect isobar amplitude info
 unsigned int waa=res->waveIndex("1-2-+0+pi-_02_f21270=pi-+_1_a11269=pi+-_0_rho770.amp");
 
 //unsigned int waa=res->waveIndex("1-1++0+pi-_11_f11285=pi-+_11_a11269=pi+-_0_rho770.amp");
  typedef double tnum;
  tnum m0=700;
  tnum dm=100;
  unsigned int i=0;
  for(unsigned int c=0;c<nbins;++c){
    tnum ml=m0+dm*c;
    tnum mu=ml+dm;
    
    if(ml<1050)continue;
    
    //TString Wave("V0_1-1++0+pi-_11_f1");
    //TString Wave("V0_1-2-+0+pi-_02_f2");
    //TString Wave("V0_1-1++0+pi-_01_rho");
    TString Wave("V0_1-1++0+pi-_01_eta1");
    //TString Wave("V0_1-0-+0+pi-_00_f0");
    
    if(ml<1000)Wave.Append("0");
    //itoa(ml,buf,4);
    Wave+=(ml);
    Wave.Append("t");
    if(mu<1000)Wave.Append("0");
    //itoa(mu,buf,4);
    Wave+=(mu);
    Wave.Append("=pi-+_01_a11269=pi+-_01_rho770.amp");
    //Wave.Append("=rho770_01_sigma.amp");
    //Wave.Append("=rho770_00_rho770.amp");
    TString wav=Wave; wav.ReplaceAll("V0_","");


    //cout << Wave << endl;

    // get amplitude 
    unsigned int wi1=res->prodAmpIndex(Wave.Data());
    unsigned int wa1=res->waveIndex(wav.Data());
    complex<double> amp=res->prodAmp(wi1);
    double integral=res->phaseSpaceIntegral(wa1); // returns sqrt!
    gPS->SetPoint(i,ml+0.5*dm,integral*integral);
    if(amp.real()==0 && amp.imag()==0)continue;
    TMatrixT<double> ampCov=res->prodAmpCov(wi1);
    amp*=1./integral;
    ampCov*=1./(integral*integral);
    
    //if(sqrt(ampCov[0][0])>amp.real() || sqrt(ampCov[1][1])>amp.imag())continue;

    gArgand->SetPoint(i,amp.real(),amp.imag());
    gArgand->SetPointError(i,sqrt(ampCov[0][0]),sqrt(ampCov[1][1]));
    
    // build labels;
    xlabels[i]=amp.real();
    ylabels[i]=amp.imag();
    labels[i]+=ml+0.5*dm;

    double intens=norm(amp); //  res->intensity(wa1);
    double intensE=res->intensityErr(wa1)/integral;
    gIntens->SetPoint(i,ml+0.5*dm,intens);
    gIntens->SetPointError(i,0.5*dm,intensE);

    double ph=180/TMath::Pi()*arg(amp);//res->phase(wi1,waa);
    if(ph<0)ph+=360;
    if(ph>300)ph-=360;
    // error propagation:
    double tn=amp.real()/amp.imag();
    double dacot=-1./(1+tn*tn);
    double dacotRe=dacot/amp.imag();
    double dacotIm=dacot*amp.real();
    
    //ampCov.Print();

    //double dphi=sqrt(dacotRe*dacotRe*ampCov[0][0]+dacotIm*dacotIm*ampCov[1][1]+dacotRe*dacotIm*ampCov[0][1]+dacotRe*dacotIm*ampCov[1][0]   );


    
    double pherr=res->phaseErr(wi1,waa);
     // check if we should make a transformation by 2pi
     // this is needed because of cyclical variable phi
     if(i>10){
       double mpre;
       double phpre;
       gPhase->GetPoint(i-1,mpre,phpre);
       double diff1=fabs(ph-phpre);
       double diff2=fabs(ph+360-phpre);
       double diff3=fabs(ph-360-phpre);
       if(diff2<diff1 && diff2<diff3)ph+=360;
       else if(diff3<diff1 && diff3<diff2)ph-=360;
     }
     
     
     
     gPhase->SetPoint(i,ml+0.5*dm,ph);
     gPhase->SetPointError(i,dm*0.5,pherr);
     
     ++i;
  }// end loop
  
  TCanvas* c= new TCanvas("c","c",10,10,1000,600);
  c->Divide(3,1);
  c->cd(1);
  gIntens->Draw("APL");
  gIntens->GetXaxis()->SetRangeUser(800,2000);
  gIntens->GetXaxis()->SetTitle("4#pi mass (GeV/c^{2})");
  gIntens->GetYaxis()->SetTitle("Intensity");
  gIntens->GetYaxis()->SetTitleOffset(1.2);
  //gPS->Draw("same");
  c->cd(2);
  gPhase->Draw("APL");
  gPhase->GetXaxis()->SetRangeUser(800,2000);
  gPhase->GetXaxis()->SetTitle("4#pi mass (GeV/c^{2})");
  gPhase->GetYaxis()->SetTitle("Phase");
  gPhase->GetYaxis()->SetTitleOffset(1.2);
  c->cd(3);
  gArgand->Draw("APL");
   gArgand->GetXaxis()->SetTitle("Re");
  gArgand->GetYaxis()->SetTitle("Im");
 gArgand->GetYaxis()->SetTitleOffset(1.2);
  // plot labels

  for(unsigned int j=0;j<nbins;++j){
    cout << xlabels[j] << "   " 
	 << ylabels[j] << "   " << labels[j].Data() << endl;
    TLatex* lab=new TLatex(xlabels[j],ylabels[j],labels[j].Data());
    //lab->SetNDC();
    lab->Draw();
  }

  infile->Close();

}
