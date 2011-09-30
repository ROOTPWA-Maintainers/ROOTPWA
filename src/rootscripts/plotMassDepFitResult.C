#include "TFile.h"
#include "TList.h"
#include "TStyle.h"
#include "TKey.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TFrame.h"
#include "TObjArray.h"
#include "TROOT.h"
//#include "TQSender.h"

#include <iostream>
#include <sstream>
#include <map>
using namespace std;

TString fitname;


TString parseTitle(TString l, unsigned int level=10){
 // setup isobar dictionary key->tex
  map<TString, TString> isobars;
  isobars["pi+"     ] = "\\pi^{+}";
  isobars["pi-"     ] = "\\pi^{-}";
  isobars["pi+-"    ] = "\\pi^{\\pm}";
  isobars["pi-+"    ] = "\\pi^{\\mp}";
  isobars["sigma"   ] = "\\sigma";
  isobars["rho770"  ] = "\\rho(770)";
  isobars["a11269"  ] = "a_{1}(1269)";
  isobars["a21320"  ] = "a_{2}(1320)";
  isobars["rho1450" ] = "\\rho(1450)";
  isobars["rho1700" ] = "\\rho(1700)";
  isobars["pi1300"  ] = "\\pi(1300)";
  isobars["pi1800"  ] = "\\pi(1800)";
  isobars["pi21670" ] = "\\pi_2(1670)";
  isobars["f01370"  ] = "f_{0}(1370)";
  isobars["f01500"  ] = "f_{0}(1500)";
  isobars["f01700"  ] = "f_{0}(1700)";
  isobars["f11285"  ] = "f_{1}(1285)";
  isobars["f11420"  ] = "f_{1}(1420)";
  isobars["b11235"  ] = "b_{1}(1235)";
  isobars["b11800"  ] = "b_{1}(1800)";
  isobars["b11500"  ] = "b_{1}(1500)";
  isobars["f21270"  ] = "f_{2}(1270)";
  isobars["f21950"  ] = "f_{2}(1950)";
  isobars["f21565"  ] = "f_{2}(1565)";
  isobars["f21270"  ] = "f_{2}(1270)";
  isobars["f22010"  ] = "f_{2}(2010)";
  isobars["f11420"  ] = "f_{1}(1420)";
  isobars["eta1440" ] = "\\eta(1420)";
  isobars["eta21645"] = "\\eta_{2}(1645)";
  isobars["rho31690"] = "\\rho_{3(1690)";

   // remove file extension
    l.Remove(l.Length() - 4);
    // extract X quantum numbers
    const TString head = l(0, 7);
    const TString I    = head(0, 1); 
    const TString G    = head(1, 1);
    const TString J    = head(2, 1);
    const TString P    = head(3, 1);
    const TString C    = head(4, 1);
    const TString M    = head(5, 1);
    const TString refl = head(6, 1);
    l.Remove(0, 7);
    // print X quantum numbers
    
    stringstream res;
    res << I << "^{" << G <<"}(" << J << "^{" << P << C << "}" << M << "^{" << refl << "})";

     // tokenize input
    TObjArray* tokens = l.Tokenize("_=");
    int        mode   = 0;
    unsigned int counter=0;
    for (int i = 0; i < tokens->GetEntries(); ++i) {
      const TString token = ((TObjString*)tokens->At(i))->GetString();
      cerr << "    " << mode << ": '" << token << "'" << endl;
      if (mode == 0) {  // isobar mode
	if (isobars.find(token) != isobars.end())
	  res << isobars[token];
	else
	  res << token;
	res << " ";
	// check which mode to switch to depending whether we get _ or =
	l.Remove(0, token.Length());
	if (l(0, 1) == "_")
	  mode = 1;
	else
	  mode = 2;
      } else if (mode == 1) {  // ls mode
	if (token.Length() == 1)  // only l
	  res << "[" << token << "] ";
	else
	  res << "\\left[#splitline{" << token(0, 1) << "}{"
	       << token(1, 1) << "}\\right]";
	l.Remove(0, token.Length());
	mode = 0;
      } else if (mode == 2) {
	if(level<=counter) break;
	++counter;
	res << "\\rightarrow ";
	if (isobars.find(token) != isobars.end())
	  res << isobars[token];
	else
	  res << token;
	res << " ";
	l.Remove(0, token.Length());
	if (l(0, 1) == "_" )
	  mode = 1;
	else
	  mode = 2;
      }
      l.Remove(0, 1); // remove delimiter
    }
    res;    

    tokens->Delete();
    delete tokens;
    TString result(res.str().c_str());
    return result;
}

void plotNice(TVirtualPad* pad, TString plotDir=""){
   TCanvas* cpopup=(TCanvas*)gROOT->FindObjectAny("cpopup");
    cpopup->Clear("");
    cpopup->cd();
    TVirtualPad* clone=(TVirtualPad*)pad->Clone();
    clone->Draw();
    clone->SetPad("pu","PopUp",0,0,1,1,0);
    cpopup->cd(1);
    cpopup->Update();
    clone->cd();
    
   
    double xcenter=0.5;
    double ycenter=0.5;
    
    double x=xcenter-0.25;
    double y=xcenter-0.2;
    
    // preliminary
    TLatex* prelim=new TLatex(x,y,"preliminary");
    prelim->SetNDC();
    prelim->SetTextColor(kGray);
    prelim->SetTextSize(0.1);
    prelim->SetTextAngle(20);
    prelim->Draw();
    
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
    
    // add waves 
    cout << "####### Title:" << endl;
    TMultiGraph* gr=(TMultiGraph*)clone->GetListOfPrimitives()->At(1);

    // gr->GetListOfGraphs()->RemoveAt(2);
    gr->GetXaxis()->SetTitle("mass (MeV/c^{2})");


    TString title=gr->GetName();
    cout << title << endl;
    // check if this is intensity or off diagonal
    TString wave;
    if(!title.Contains("---")){
      cout << "Processing Intensity:" << endl;
      wave=parseTitle(title,0);
      cout << wave << endl;
      yc=ycenter+0.24;
    }
    else {
      // Split
      unsigned int i = title.Index("---");
      TString wave1=title(3,i-3);
      TString wave2=title(i+3,title.Length());
      
      cout << parseTitle(wave1,0) << endl;
      cout << parseTitle(wave2,0) << endl;
      
      

      wave="#splitline{Interference - ";
      if(title.Contains("Re"))wave+= "real part";
      else if(title.Contains("Im"))wave+="imaginary part"; 
      else wave+="phase difference";
      wave+="}{#splitline{";
      wave+=parseTitle(wave1,0);
      wave+="}{";
      wave+=parseTitle(wave2,0);
      wave+="}}";
	   // check if we put text up or down
     double max=gr->GetYaxis()->GetXmax();
     double min=gr->GetYaxis()->GetXmin();
     // get last data point
     TGraph* g=(TGraph*)gr->GetListOfGraphs()->At(0);
     double yg,xg;
     g->GetPoint((int)(g->GetN()*0.6),xg,yg);
     if(fabs(max-yg)>fabs(min-yg)){
	yc=ycenter+0.24;
     }
     else{
       yc=ycenter-0.25;
     }


    }

     xc=xcenter+0.05;
    //xc=xcenter+0.1;
  
    TLatex* waveL=new TLatex(xc,yc,wave);
    waveL->SetNDC();
    waveL->SetTextSize(0.03);
    waveL->Draw();


     // fitname
    xc=xcenter-0.4;
    //if(right)xc=xcenter+0.1;
     yc=ycenter+0.42;
    TLatex* fitNameL=new TLatex(xc,yc,fitname);
    fitNameL->SetNDC();
    fitNameL->SetTextSize(0.025);
    fitNameL->Draw();


    cpopup->UseCurrentStyle();
    cpopup->Update();

    if(plotDir.Length()>1){
      TString filename=plotDir+title+".eps";
      filename.ReplaceAll(".amp","");
      filename.ReplaceAll("+","p");
      filename.ReplaceAll("---","___");
      filename.ReplaceAll("-","m");
      filename.ReplaceAll("=","_to_");
      cpopup->SaveAs(filename);
    }

    //cpopup->Flush();
}// end plotnice


void exec3event(Int_t event, Int_t x, Int_t y, TObject *selected)
{
  TCanvas *c = (TCanvas *) gTQSender;
  //if(selected->IsA()->GetName())
  //printf("Canvas %s: event=%d, x=%d, y=%d, selected=%s\n", c->GetName(),event, x, y, selected->IsA()->GetName());
  
  if(event==1){ // clicked

    // figure out which pad we clicked onto
    TVirtualPad* pad=c->GetClickSelectedPad();
    
    plotNice(pad);
  }
}

//------------------------------------------------------
//------------------------------------------------------
//------------------------------------------------------
//------------------------------------------------------
//------------------------------------------------------


void plotMassDepFitResult(TString infilename, TString plotdir="plots/", TString fittitle="", double mmin=0, double mmax=0 ){
  if(fittitle.Length()<=1)fitname=infilename;
  else fitname=fittitle;
  TFile* infile=TFile::Open(infilename);
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

  TCanvas* cS=new TCanvas("cS","Spin Density Matrix",10,10,1000,1000);
  cS->Divide(nwaves,nwaves,0,0);
 TCanvas* cRe=new TCanvas("cRe","Spin Density Matrix - RealImagPart",10,10,1000,1000);
  cRe->Divide(nwaves,nwaves,0,0);
  cRe->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,"exec3event(Int_t,Int_t,Int_t,TObject*)");
  cS->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,"exec3event(Int_t,Int_t,Int_t,TObject*)");

  TCanvas* cpopup=new TCanvas("cpopup","PopUp",50,50,700,700);


  // TCanvas* cIm=new TCanvas("cIm","Spin Density Matrix - ImagPart",10,10,1000,1000);
  //cIm->Divide(nwaves,nwaves,0,0);
  
  // cRe->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
  //		      "exec3event(Int_t,Int_t,Int_t,TObject*)");
  //c->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
  //		      "exec3event(Int_t,Int_t,Int_t,TObject*)");

  // do plotting
  for(unsigned int ip=0;ip<nwaves;++ip){
    for(unsigned int jp=ip;jp<nwaves;++jp){
      cS->cd(jp+ip*nwaves+1);
      if(ip==jp){
	TMultiGraph* g=(TMultiGraph*)infile->Get(wavenames[ip]);
	//g->GetListOfGraphs()->RemoveAt(3); // remove black line
	//g->GetListOfGraphs()->RemoveAt(2); // remove fit
	g->Draw("APC");

	// rescale
	TGraphErrors* datag=(TGraphErrors*)g->GetListOfGraphs()->At(1);
	double* y=datag->GetY();
	double max=-1E6;
	double min=1E6;
	for(unsigned int i=0;i<datag->GetN();++i){
	  if(max<y[i])max=y[i];
	  if(min>y[i])min=y[i];
	}
	cout << min << "     " << max << endl;
	g->GetYaxis()->SetRangeUser(0 < min ? -0.8*min : 1.2*min,1.2*max);
	g->GetYaxis()->SetTitle("intensity");
	g->GetYaxis()->SetTitleOffset(1.2);
	
	if(mmin!=0 || mmax!=0){
	  g->GetXaxis()->SetRangeUser(mmin,mmax);
        }
	//g->GetHistogram()->Draw();//gPad->Update();
	//g->Draw("APC");
	cRe->cd(jp+ip*nwaves+1);
	g->Draw("APC");


	/*
	TCanvas* c2=new TCanvas();
	g->Draw("APC");
	g->GetXaxis()->SetTitle("mass (MeV/c^{2})");
	g->GetYaxis()->SetTitle("Intensity");
	g->GetYaxis()->SetTitleOffset(1.2);
	g->SaveAs(TString(plotdir+wavenames[ip]+".eps"));
	delete c2;
	*/
	plotNice(cRe->GetPad(jp+ip*nwaves+1),plotdir);
      }
      else {
	TString key="dPhi_"+wavenames[ip]+"---"+wavenames[jp];
	TMultiGraph* g=(TMultiGraph*)infile->Get(key);
	if(g!=NULL){
       	  g->Draw("AN");
	  if(mmin!=0 || mmax!=0){
	    g->GetXaxis()->SetRangeUser(mmin,mmax);
	  }
	  double max=-1E6;
	  double min=1E6;
	  TGraphErrors* fitg=(TGraphErrors*)g->GetListOfGraphs()->At(2);
	  double* y=fitg->GetY();
	  for(unsigned int i=0;i<fitg->GetN();++i){
	    if(max<y[i])max=y[i];
	    if(min>y[i])min=y[i];
	  }
	  TAxis* a=g->GetYaxis();
	  if(a!=NULL)a->SetRangeUser(0.5*(max+min)-220,0.5*(max+min)+220);
	  a->SetTitle("#Delta#Phi");
	  a->SetTitleOffset(1.2);
	  g->Draw("A");
	  plotNice(cS->GetPad(jp+ip*nwaves+1),plotdir);
	  //TCanvas* c2=new TCanvas();
	  //g->Draw("A");
	  //g->GetXaxis()->SetTitle("mass (MeV/c^{2})");
	  //g->GetYaxis()->SetTitle("#Delta #phi");
	  //g->GetYaxis()->SetTitleOffset(1.2);
	  //c2->SaveAs(TString(plotdir+key+".eps"));
	  //delete c2;
	  key.ReplaceAll("dPhi","Re");
	  TMultiGraph* g2=(TMultiGraph*)infile->Get(key);
	  //g2->GetListOfGraphs()->RemoveAt(2);
	  TVirtualPad* pa= cRe->cd(jp+ip*nwaves+1);
	  pa->SetFillColor(kYellow-9);
	  g2->Draw("A");
	  if(mmin!=0 || mmax!=0){
	    g2->GetXaxis()->SetRangeUser(mmin,mmax);
	  }
	  g2->GetXaxis()->SetTitle("mass (MeV/c^{2})");
	  g2->GetYaxis()->SetTitle("real part");
	  g2->GetYaxis()->SetTitleOffset(1.2);
	  plotNice(cRe->GetPad(jp+ip*nwaves+1),plotdir);
	  //c2=new TCanvas();
	  //g2->Draw("A");
	  //g2->GetXaxis()->SetTitle("mass (MeV/c^{2})");
	  //g2->GetYaxis()->SetTitle("Re(#rho_{ij})");
	  //g2->GetYaxis()->SetTitleOffset(1.2);
	  //c2->SaveAs(TString(plotdir+key+".eps"));
	  //delete c2;
	  key.ReplaceAll("Re","Im");
	  TMultiGraph* g3=(TMultiGraph*)infile->Get(key);
	  //g3->GetListOfGraphs()->RemoveAt(2);
	  //cIm->cd(jp+ip*nwaves+1);
	  pa=cRe->cd(ip+jp*nwaves+1);
	  pa->SetFillColor(kSpring+6);
	  g3->Draw("A");
	 if(mmin!=0 || mmax!=0){
	    g3->GetXaxis()->SetRangeUser(mmin,mmax);
	  }
	  g3->GetXaxis()->SetTitle("mass (MeV/c^{2})");
	  g3->GetYaxis()->SetTitle("imaginary part");
	  g3->GetYaxis()->SetTitleOffset(1.2);
	  plotNice(cRe->GetPad(ip+jp*nwaves+1),plotdir);
	  //c2=new TCanvas();
	  //g3->Draw("A");
	  //g3->GetXaxis()->SetTitle("mass (MeV/c^{2})");
	  // g3->GetYaxis()->SetTitle("Im(#rho_{ij})");
	  //g2->GetYaxis()->SetTitleOffset(1.2);
	  //c2->SaveAs(TString(plotdir+key+".eps"));
	  //delete c2;
	  
	}// end if g!=NULL
      } // end else
    } // end inner loop
  } // end plotting loop

  cS->SaveAs(TString(plotdir+"/spindensitymatrix.eps"));
  cRe->SaveAs(TString(plotdir+"/spindensitymatrixRe.eps"));
  //cIm->SaveAs(TString(plotdir+"/spindensitymatrixIm.eps"));
}



