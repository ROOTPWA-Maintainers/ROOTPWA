#include "TFile.h"
#include "TList.h"
#include "TKey.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TMultiGraph.h"
#include "TAxis.h"

#include <iostream>
using namespace std;

void plotMassDepFitResult(TString infilename, TString plotdir="plots/"){
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
 TCanvas* cRe=new TCanvas("cRe","Spin Density Matrix - RealPart",10,10,1000,1000);
  cRe->Divide(nwaves,nwaves,0,0);
   TCanvas* cIm=new TCanvas("cIm","Spin Density Matrix - ImagPart",10,10,1000,1000);
  cIm->Divide(nwaves,nwaves,0,0);
  
  // do plotting
  for(unsigned int ip=0;ip<nwaves;++ip){
    for(unsigned int jp=ip;jp<nwaves;++jp){
      c->cd(jp+ip*nwaves+1);
      if(ip==jp){
	TMultiGraph* g=(TMultiGraph*)infile->Get(wavenames[ip]);
	g->Draw("APC");
	// rescale
	TGraphErrors* datag=(TGraphErrors*)g->GetListOfGraphs()->At(0);
	double* y=datag->GetY();
	double max=-1E6;
	double min=1E6;
	for(unsigned int i=0;i<datag->GetN();++i){
	  if(max<y[i])max=y[i];
	  if(min>y[i])min=y[i];
	}
	cout << min << "     " << max << endl;
	g->GetYaxis()->SetRangeUser(0 < min ? -0.8*min : 1.2*min,1.2*max);
	g->GetYaxis()->SetTitle("blubler");
	//g->GetHistogram()->Draw();//gPad->Update();
	//g->Draw("APC");
	cRe->cd(jp+ip*nwaves+1);
	g->Draw("APC");

	TCanvas* c2=new TCanvas();
	g->Draw("APC");
	g->GetXaxis()->SetTitle("mass (MeV/c^{2})");
	g->GetYaxis()->SetTitle("Intensity");
	g->GetYaxis()->SetTitleOffset(1.2);
	g->SaveAs(TString(plotdir+wavenames[ip]+".eps"));
	delete c2;
      }
      else {
	TString key="dPhi_"+wavenames[ip]+"---"+wavenames[jp];
	TMultiGraph* g=(TMultiGraph*)infile->Get(key);
	if(g!=NULL){
	  g->Draw("AN");
	  double max=-1E6;
	  double min=1E6;
	  TGraphErrors* fitg=(TGraphErrors*)g->GetListOfGraphs()->At(1);
	  double* y=fitg->GetY();
	  for(unsigned int i=0;i<fitg->GetN();++i){
	    if(max<y[i])max=y[i];
	    if(min>y[i])min=y[i];
	  }
	  TAxis* a=g->GetYaxis();
	  if(a!=NULL)a->SetRangeUser(0.5*(max+min)-200,0.5*(max+min)+200);
	  g->Draw("A");
	  TCanvas* c2=new TCanvas();
	  g->Draw("A");
	  g->GetXaxis()->SetTitle("mass (MeV/c^{2})");
	  g->GetYaxis()->SetTitle("#Delta #phi");
	  //g->GetYaxis()->SetTitleOffset(1.2);
	  c2->SaveAs(TString(plotdir+key+".eps"));
	  delete c2;
	  key.ReplaceAll("dPhi","Re");
	  TMultiGraph* g2=(TMultiGraph*)infile->Get(key);
	  cRe->cd(jp+ip*nwaves+1);
	  g2->Draw("A");
	  g2->GetXaxis()->SetTitle("mass (MeV/c^{2})");
	  g2->GetYaxis()->SetTitle("Re(#rho_{ij})");
	  c2=new TCanvas();
	  g2->Draw("A");
	  g2->GetXaxis()->SetTitle("mass (MeV/c^{2})");
	  g2->GetYaxis()->SetTitle("Re(#rho_{ij})");
	  //g2->GetYaxis()->SetTitleOffset(1.2);
	  c2->SaveAs(TString(plotdir+key+".eps"));
	  delete c2;
	  key.ReplaceAll("Re","Im");
	  TMultiGraph* g3=(TMultiGraph*)infile->Get(key);
	  cIm->cd(jp+ip*nwaves+1);
	  g3->Draw("A");
	  g3->GetXaxis()->SetTitle("mass (MeV/c^{2})");
	  g3->GetYaxis()->SetTitle("Im(#rho_{ij})");
	  c2=new TCanvas();
	  g3->Draw("A");
	  g3->GetXaxis()->SetTitle("mass (MeV/c^{2})");
	  g3->GetYaxis()->SetTitle("Im(#rho_{ij})");
	  //g2->GetYaxis()->SetTitleOffset(1.2);
	  c2->SaveAs(TString(plotdir+key+".eps"));
	  delete c2;
	  
	}

      }
    } // end inner loop
  } // end plotting loop
  c->SaveAs(TString(plotdir+"/spindensitymatrix.eps"));
  cRe->SaveAs(TString(plotdir+"/spindensitymatrixRe.eps"));
  cIm->SaveAs(TString(plotdir+"/spindensitymatrixIm.eps"));
}
