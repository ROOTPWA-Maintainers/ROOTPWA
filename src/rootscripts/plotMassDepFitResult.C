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
#include <fstream>
#include <sstream>
#include <map>
using namespace std;

#define RECOLOR kYellow-9
#define IMCOLOR kSpring-9


TString fitname;

TString parseTitle(TString l, unsigned int level=10){
 // setup isobar dictionary key->tex
  map<TString, TString> isobars;
  isobars["pi+"     ] = "\\pi^{+}";
  isobars["pi-"     ] = "\\pi^{-}";
  isobars["pi+-"    ] = "\\pi^{\\pm}";
  isobars["pi-+"    ] = "\\pi^{\\pm}";
  isobars["sigma"   ] = "\\sigma";
  isobars["rho770"  ] = "\\rho(770)";
  isobars["a11269"  ] = "a_{1}(1269)";
  isobars["a21320"  ] = "a_{2}(1320)";
  isobars["rho1450" ] = "\\rho(1450)";
  isobars["rho1600" ] = "\\rho(1600)";
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
 isobars["eta11600" ] = "\\eta_{1}(1600)";
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
   if(cpopup==NULL){
     cerr << "Popupo Window not found! RESPAWNING" << endl;
     cpopup=new TCanvas("cpopup","PopUp",50,50,700,700);
   }

    cpopup->Clear("");
    cpopup->cd();
    TVirtualPad* clone=(TVirtualPad*)pad->Clone();
    if(clone==NULL){
        cerr << "Cloning pad failed!" << endl;
     return;
    }


    clone->SetFillStyle(1);
    clone->Draw();
    clone->SetPad("pu","PopUp",0,0,1,1,0);
    cpopup->cd(1);
    cpopup->Update();

    TMultiGraph* gr=(TMultiGraph*)clone->GetListOfPrimitives()->At(1);
    if(gr==NULL){
        cerr << "Cloning graphs failed!" << endl;
     return;
    }




    TString title=gr->GetName();
    TH1* hx=NULL;
    if(!title.Contains("dPhi") && !title.Contains("Re") && !title.Contains("Im"))hx=(TH1*)clone->GetListOfPrimitives()->At(2);



    TGraphErrors* gdata=(TGraphErrors*)gr->GetListOfGraphs()->At(1);
    if(gdata==NULL){
        cerr << "Cloning data failed!" << endl;
     return;
    }
    gdata->SetLineWidth(2);


    double xmax=gr->GetXaxis()->GetBinCenter(gr->GetXaxis()->GetLast());
    double xmin=gr->GetXaxis()->GetBinCenter(gr->GetXaxis()->GetFirst());
    double max=-1E6;
    double min=1E6;
    TGraphErrors* fitg=(TGraphErrors*)gr->GetListOfGraphs()->At(2);

    if(fitg!=NULL){
      fitg->SetLineWidth(2);
      double* yp=fitg->GetY();
      for(unsigned int i=0;i<fitg->GetN();++i){
	if(max<yp[i])max=yp[i];
	if(min>yp[i])min=yp[i];
      }
    }
    else {
      max=0;
      min=0;
    }
    // this works nly for phase plots
    double ymin=0.5*(max+min)-220;
    double ymax=0.5*(max+min)+220;

    if(title.Contains("dPhi")){
      cerr << "Ymin: " << ymin << "   Ymax: " << ymax << endl;
      // for phase plots manually cut the systematic errors!
      TGraphErrors* gsys=(TGraphErrors*)gr->GetListOfGraphs()->At(0);
      //gsys->GetYaxis()->SetRangeUser(ymin,ymax);
      unsigned int n=gsys->GetN();
      for(unsigned int i=0;i<n;++i){
	double y, x, ey, ex;
	gsys->GetPoint(i,x,y);
	ey=gsys->GetErrorY(i);
	ex=gsys->GetErrorX(i);

	// check if error band is in limits
	if(y+ey>ymax){
	  cerr << "correcting upper" << endl;
	  double dy=0.5*(y+ey-ymax);
	  y-=dy;
	  ey-=dy;
	}
	if(y-ey<ymin){
	  double dy=0.5*(ymin-y+ey);
	  y+=dy;
	  ey-=dy;
	}
	gsys->SetPoint(i,x,y);
	gsys->SetPointError(i,ex,ey);
      }
    }
    if(hx!=NULL)hx->SetLineColor(kMagenta);
    clone->Draw();
    clone->cd();
    clone->GetFrame()->Draw();
    clone->SetBorderMode(0);




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
    //double xc=xcenter+0.05;
    double xc=0.105;//xcenter+0.05;

    //if(right)xc=xcenter+0.1;
    //double yc=ycenter+0.35;
    double yc=ycenter+0.45;
    TLatex* com04=new TLatex(xc,yc,"COMPASS 2004");
    com04->SetNDC();
    com04->SetTextSize(0.05);
    com04->Draw();

    // 5 pi on pb
    yc=yc-0.04;
    TLatex* react=new TLatex(xc,yc,"#pi^{-} Pb #rightarrow #pi^{-}#pi^{+}#pi^{-}#pi^{+}#pi^{-} Pb");
    react->SetNDC();
    react->SetTextSize(0.039);
    react->Draw();

    // add waves
    //cout << "####### Title:" << endl;



    if(ymax > 0 && ymin < 0){
      TLine* zeroline=new TLine(xmin,0,xmax,0);
      zeroline->SetLineStyle(7);
      zeroline->Draw();
    }

    // gr->GetListOfGraphs()->RemoveAt(2);
    gr->GetXaxis()->SetTitle("mass (GeV/c^{2})");
    gr->GetXaxis()->Draw();
    gr->GetYaxis()->Draw();


    //cout << title << endl;
    // check if this is intensity or off diagonal
    TString wave;
    bool isIntens=false;
    if(!title.Contains("---")){
      cout << "Processing Intensity:" << endl;
      wave=parseTitle(title,1);
      //cout << wave << endl;
      yc=ycenter+0.44;
      isIntens=true;
    }
    else {
      // Split
      //cout << title << endl;

      unsigned int offset=3;
      wave="#splitline{Interference - ";
      if(title.Contains("Re"))wave+= "real part";
      else if(title.Contains("Im"))wave+="imaginary part";
      else {
	wave+="phase difference";
	offset=5;
      }

      unsigned int i = title.Index("---");
      TString wave1=title(offset,i-offset);
      TString wave2=title(i+3,title.Length());

      //cout << parseTitle(wave1,1) << endl;
      //cout << parseTitle(wave2,1) << endl;



      wave="#splitline{";
      wave+=parseTitle(wave1,1);
      wave+="}{";
      wave+=parseTitle(wave2,1);
      wave+="}";
	   // check if we put text up or down
     double max=gr->GetYaxis()->GetXmax();
     double min=gr->GetYaxis()->GetXmin();
     // get last data point
     TGraph* g=(TGraph*)gr->GetListOfGraphs()->At(1);
     double yg,xg;
     g->GetPoint((int)(g->GetN()*0.6),xg,yg);
     if(fabs(max-yg)>fabs(min-yg)){
	yc=ycenter+0.44;
     }
     else{
       yc=ycenter+0.44;
     }


    }

    xc=xcenter-0.05;
    //xc=xcenter+0.1;

    TLatex* waveL=new TLatex(xc,yc,wave);
    waveL->SetNDC();
    if(isIntens)waveL->SetTextSize(0.03);
    else waveL->SetTextSize(0.02);
    waveL->Draw();


     // fitname
    xc=xcenter-0.03;
    //if(right)xc=xcenter+0.1;
     yc=ycenter+0.35;
    TLatex* fitNameL=new TLatex(xc,yc,fitname);
    fitNameL->SetNDC();
    fitNameL->SetTextSize(0.05);
    fitNameL->Draw();

    cpopup->UseCurrentStyle();
    if(hx!=NULL)hx->SetLineColor(kMagenta);
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
	cerr << "###### end of plotnice" << endl;
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

// plotlevel:
// 0 = data + fit + component
// 1 = data + fit
// 2 = data only

void plotMassDepFitResult(TString infilename, TString plotdir="plots/", TString fittitle="", double mmin=0, double mmax=0, unsigned int plotLevel=0 , TString xcheckfile="", bool onlyDiag=false){
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
  gStyle->SetFrameBorderMode(0);

    gROOT->ForceStyle();

    TFile* xcheck=NULL;
    ifstream xcheckmapfile;
    map<TString,TString> xcheckmap;
    if(xcheckfile.Length()>2){
      cerr << "Found XCheck File " << xcheckfile << endl;
      xcheck=TFile::Open(xcheckfile);
      xcheckfile.ReplaceAll(".root",".map");
      xcheckmapfile.open(xcheckfile.Data());
      TString wave, key;
      while(xcheckmapfile.good()){
	xcheckmapfile >> wave >> key;
	xcheckmap[wave]=key;
	cerr << wave << " ---> " << key << endl;
      }
    }




    TList* keylist=infile->GetListOfKeys();
    unsigned int num=keylist->GetEntries();

  // loop over keys and count waves
  vector<TString> wavenames;

  for(unsigned int i=0;i<num;++i){
    TString keyname(((TKey*)keylist->At(i))->GetName());
    if(keyname.Contains("dPhi") || keyname.Contains("Re") || keyname.Contains("Im") || keyname.Contains("fPS"))continue;
    wavenames.push_back(keyname);
  }



  unsigned int nwaves=wavenames.size();
  std::cout << nwaves << " waves used in fit" << endl;

  map<TString,TString> xcheckmapReIm;
  for(unsigned int i=0;i<nwaves;++i){
TString map1=xcheckmap[wavenames[i]];
 unsigned int index1=atoi(map1(1,2).Data());
    for(unsigned int j=i+1;j<nwaves;++j){
      TString keyRe="Re_"+wavenames[i]+"---"+wavenames[j];
      TString keyIm="Im_"+wavenames[i]+"---"+wavenames[j];
      TString map2=xcheckmap[wavenames[j]];
      // extract index number
      unsigned int index2=atoi(map2(1,2).Data());
      TString ReMapper="h";ReMapper+= 10000+100*index1+index2;
      TString ImMapper="h";ImMapper+= 20000+100*index1+index2;
      xcheckmapReIm[keyRe]=ReMapper;
      xcheckmapReIm[keyIm]=ImMapper;
    }
  } // end loop to build xcheckmap



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
      cerr << " ##### " << ip << "-" << jp << " ###### " << endl;
      if(ip==jp){
	TMultiGraph* g=(TMultiGraph*)infile->Get(wavenames[ip]);

	// remove components and phase space graphs
	// plotlevel:
	// 0 = data + fit + component
	// 1 = data + fit
	// 2 = data only
	if(plotLevel>0){
	  for(unsigned int i=g->GetListOfGraphs()->GetSize()-1;i>2;--i){
	    g->GetListOfGraphs()->RemoveAt(i);
	  }
	}
	// remove only ps
	else g->GetListOfGraphs()->RemoveAt(3);



	if(plotLevel>1)g->GetListOfGraphs()->RemoveAt(2); // remove fit
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
	//g->GetYaxis()->SetRangeUser(0 < min ? -0.8*min : 1.2*min,1.2*max);
	g->GetYaxis()->SetTitle("intensity");
	g->GetYaxis()->SetTitleOffset(1.3);

	if(mmin!=0 || mmax!=0){
	  g->GetXaxis()->SetRangeUser(mmin,mmax);
        }
	//g->GetHistogram()->Draw();//gPad->Update();
	//g->Draw("APC");
	cRe->cd(jp+ip*nwaves+1);
	g->Draw("APC");

	if(xcheck!=NULL){
	  // get key
	  TString xcheckkey=xcheckmap[wavenames[ip]];
	  cerr << "Adding xcheckplot " << xcheckkey << endl;
	  if(xcheckkey!=""){
	    TH1* xh=(TH1*)xcheck->Get(xcheckkey);
	    if(xh!=NULL){
	      xh->SetLineColor(kMagenta);
	      xh->Draw("same");
	    }
	    else cerr << "Did not find xcheckPlot!" << endl;
	  }
	}
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
	if(onlyDiag)continue;
	TString key="dPhi_"+wavenames[ip]+"---"+wavenames[jp];
	cerr << key << endl;
	TMultiGraph* g=(TMultiGraph*)infile->Get(key);
	if(g!=NULL){
	  TString title=g->GetName();
	  unsigned int i = title.Index("---");
	  TString wave1=title(5,i-5);
	  TString wave2=title(i+3,title.Length());
	  // check for same reflectivity
	  if(wave1(6)!=wave2(6)) continue;

	  if(plotLevel>1)g->GetListOfGraphs()->RemoveAt(2); // remove fit
	    g->Draw("AN");
	  if(mmin!=0 || mmax!=0){
	    g->GetXaxis()->SetRangeUser(mmin,mmax);
	  }
	  double max=-1E6;
	  double min=1E6;
	  TGraphErrors* fitg=(TGraphErrors*)g->GetListOfGraphs()->At(2);
	  if(fitg!=NULL){
	    double* y=fitg->GetY();
	    for(unsigned int i=0;i<fitg->GetN();++i){
	      if(max<y[i])max=y[i];
	      if(min>y[i])min=y[i];
	    }
	  }
	  TAxis* a=g->GetYaxis();
	  if(a!=NULL)a->SetRangeUser(0.5*(max+min)-220,0.5*(max+min)+220);
	  a->SetTitle("#Delta#phi");
	  a->SetTitleOffset(1.3);
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
	  if(plotLevel>1)g2->GetListOfGraphs()->RemoveAt(2);
	  TVirtualPad* pa= cRe->cd(jp+ip*nwaves+1);
	  pa->SetFillColor(RECOLOR);
	  g2->Draw("A");
	  if(mmin!=0 || mmax!=0){
	    g2->GetXaxis()->SetRangeUser(mmin,mmax);
	  }
	  g2->GetXaxis()->SetTitle("mass (MeV/c^{2})");
	  g2->GetYaxis()->SetTitle("real part");
	  g2->GetYaxis()->SetTitleOffset(1.3);

	if(xcheck!=NULL){
	  // get key
	  TString xcheckkey=xcheckmapReIm[key];
	  cerr << "Adding xcheckplot " << xcheckkey << endl;
	  if(xcheckkey!=""){
	    TH1* xh=(TH1*)xcheck->Get(xcheckkey);
	    if(xh!=NULL){
	      xh->SetLineColor(kMagenta);
	      xh->Draw("same");
	    }
	    else cerr << "Did not find xcheckPlot!" << endl;
	  }
	}


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
	  if(plotLevel>1)g3->GetListOfGraphs()->RemoveAt(2);
	  //cIm->cd(jp+ip*nwaves+1);
	  pa=cRe->cd(ip+jp*nwaves+1);
	  pa->SetFillColor(IMCOLOR);
	  g3->Draw("A");
	 if(mmin!=0 || mmax!=0){
	    g3->GetXaxis()->SetRangeUser(mmin,mmax);
	  }
	  g3->GetXaxis()->SetTitle("mass (MeV/c^{2})");
	  g3->GetYaxis()->SetTitle("imaginary part");
	  g3->GetYaxis()->SetTitleOffset(1.3);

if(xcheck!=NULL){
	  // get key
	  TString xcheckkey=xcheckmapReIm[key];
	  cerr << "Adding xcheckplot " << xcheckkey << endl;
	  if(xcheckkey!=""){
	    TH1* xh=(TH1*)xcheck->Get(xcheckkey);

	    if(xh!=NULL){
	      xh->Scale(-1);
	      xh->SetLineColor(kMagenta);
	      xh->Draw("same");
	    }
	    else cerr << "Did not find xcheckPlot!" << endl;
	  }
	}


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
