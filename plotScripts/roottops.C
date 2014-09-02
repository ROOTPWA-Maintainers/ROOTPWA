#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TList.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TROOT.h"
#include "TKey.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TH2D.h"
#include "TMultiGraph.h"

#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>

using namespace std;





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
  isobars["b21800"  ] = "b_{2}(1800)";
  isobars["b11800"  ] = "b_{1}(1800)";
  isobars["b11500"  ] = "b_{1}(1500)";
  isobars["f21270"  ] = "f_{2}(1270)";
  isobars["f21950"  ] = "f_{2}(1950)";
  isobars["f21565"  ] = "f_{2}(1565)";
  isobars["f21270"  ] = "f_{2}(1270)";
  isobars["f22010"  ] = "f_{2}(2010)";
  isobars["f11420"  ] = "f_{1}(1420)";
  isobars["eta1440" ] = "\\eta(1420)";
  isobars["eta11600" ] = "\\eta_{1}(1600)";
  isobars["eta21645"] = "\\eta_{2}(1645)";
  isobars["rho31690"] = "\\rho_{3}(1690)";
  isobars["rho1600"] = "\\rho(1600)";

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
    //res;

    tokens->Delete();
    delete tokens;
    TString result(res.str().c_str());
    return result;
}




void annotatePlot(TVirtualPad* pad){

  pad->cd();

    double xcenter=0.5;
    double ycenter=0.5;

    double x=xcenter-0.2;
    double y=xcenter-0.2;

    // preliminary
    TLatex* prelim=new TLatex(x,y,"preliminary");
    prelim->SetNDC();
    prelim->SetTextColor(kGray);
    prelim->SetTextSize(0.1);
    prelim->SetTextAngle(20);
    prelim->Draw();

    // compass 2004
    double xc=xcenter+0.08;
    //if(right)xc=xcenter+0.1;
    double yc=ycenter+0.35;
    TLatex* com04=new TLatex(xc,yc,"COMPASS 2004");
    com04->SetNDC();
    com04->SetTextSize(0.05);
    com04->Draw();

    // 5 pi on pb
    xc=xcenter+0.08;
    //xc=xcenter+0.1;
    yc=ycenter+0.31;
    TLatex* react=new TLatex(xc,yc,"#pi^{-} Pb #rightarrow #pi^{-}#pi^{+}#pi^{-}#pi^{+}#pi^{-} Pb");
    react->SetNDC();
    react->SetTextSize(0.039);
    react->Draw();

    // add waves
    cout << "####### Title:" << endl;
    TMultiGraph* gr=dynamic_cast<TMultiGraph*>(pad->GetListOfPrimitives()->At(0));
    TH2D* h2=dynamic_cast<TH2D*>(pad->GetListOfPrimitives()->At(0));
    TString title;
    if(gr!=NULL){
      title=gr->GetName();
    }
    else {
       title=h2->GetName();
       title.Remove(0,1);
    }
    cout << title << endl;
    // check if this is intensity or off diagonal
    TString wave;
    if(title.Contains("amp")){
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
	wave+= title.Contains("Re") ? "real part" : "imaginary part";
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
    } // end if title.Contains amp
    else {
      wave=title;
      yc=ycenter+0.24;
    }

     xc=xcenter+0.08;
    //xc=xcenter+0.1;

    TLatex* waveL=new TLatex(xc,yc,wave);
    waveL->SetNDC();
    waveL->SetTextSize(0.03);
    waveL->Draw();


     // fitname
    // xc=xcenter-0.4;
    //if(right)xc=xcenter+0.1;
    // yc=ycenter+0.42;
    // TLatex* fitNameL=new TLatex(xc,yc,fitname);
    //fitNameL->SetNDC();
    //fitNameL->SetTextSize(0.025);
    //fitNameL->Draw();
    double xmax, xmin;
    if(gr!=NULL){
      xmax=gr->GetXaxis()->GetXmax();
      xmin=gr->GetXaxis()->GetXmin();
    }
    else {
       xmax=h2->GetXaxis()->GetXmax();
       xmin=h2->GetXaxis()->GetXmin();
    }


    TLine* zeroline=new TLine(xmin,0,xmax,0);
    zeroline->SetLineStyle(7);
    zeroline->Draw();

    if(gr!=NULL){
      gr->GetXaxis()->SetTitle("mass (GeV/c^{2})");
      gr->Draw();
    }
    if(h2!=NULL){
      h2->GetXaxis()->SetTitle("mass (GeV/c^{2})");
      //h2->Draw("COL");
    }
    pad->UseCurrentStyle();
    pad->Update();
    //cpopup->Flush();
}









void roottops(const TString& infilename, const TString& plotdir="", const TString& normfilename=""){

 gStyle->SetFrameFillStyle(4000);

  bool savePlots=plotdir.Length()>1;


  TString name=infilename;
  name.ReplaceAll(".root","");

  const int    nmbPadsPerCanvMin = 6;            // minimum number of pads each canvas is subdivided into

  TFile* infile=TFile::Open(infilename,"READ");

 Int_t font=132;
  gStyle->SetTextFont(font);
  gStyle->SetLabelFont(font,"xyz");
  gStyle->SetTitleFont(font,"xyz");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetStripDecimals(1);
    TGaxis::SetMaxDigits(4);

    gROOT->ForceStyle();




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
	((TKey*)keylist->At(i))->ReadObj()->Draw("COL");
	if(keyname.Contains("rho1")){
	  ((TH2D*)(((TKey*)keylist->At(i))->ReadObj()))->GetYaxis()->SetRangeUser(-100,15000);}
      }
      else ((TKey*)keylist->At(i))->ReadObj()->Draw("AP");
       annotatePlot(gPad);
       if(savePlots){
	 TPad* mypad=(TPad*)gPad;
	 mypad->Dump();
	 TString filename=plotdir+"plot_"+keyname;
	 filename.ReplaceAll(".amp",".eps");
	 filename.ReplaceAll("+","p");
	 filename.ReplaceAll("-","m");
	 filename.ReplaceAll("=","_to_");
	 TCanvas* cplot=new TCanvas("cplot","cplot",10,10,700,700);
	 cplot->cd();
	 TVirtualPad* clone=(TVirtualPad*)mypad->Clone();
	 clone->Draw();
	 clone->SetPad("pu","PopUp",0,0,1,1,0);
	 cplot->cd(1);
	 cplot->Update();
	 cplot->SaveAs(filename);
	 cplot;
       }
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
      continue;
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
