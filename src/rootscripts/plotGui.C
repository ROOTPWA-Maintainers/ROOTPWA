//
// Author: Ilka Antcheva   1/12/2006

// This macro gives an example of how to create a list box
// and how to set and use its multiple selection feature.
// To run it do either:
// .x listBox.C
// .x listBox.C++

#include <TApplication.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TGListBox.h>
#include <TList.h>
#include <TFile.h>
#include <TTree.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

#include "fitResult.h"
#include <string>
#include <vector>
#include <iostream>
#include <limits>
//#include <stringstream>

using namespace rpwa;
using namespace std;

class MyMainFrame : public TGMainFrame {

private:
   TGListBox           *fListBox;
  TGListBox           *fListBox2;
   TGCheckButton       *fCheckMulti;
   TList               *fSelected;   
  TTree               *fTree;
  std::vector<std::string> fwavenames;

public:
   MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h, TTree* tree);
   virtual ~MyMainFrame();
   void DoExit();
   void DoSelect();
   void HandleButtons();
   void PrintSelected();

   ClassDef(MyMainFrame, 0)
};

void MyMainFrame::DoSelect()
{
   Printf("Slot DoSelect()");
}

void MyMainFrame::DoExit()
{
   Printf("Slot DoExit()");
   gApplication->Terminate(0);
}

MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h, TTree* tree) :
   TGMainFrame(p, w, h)
{

  // analyze tree
  fTree=tree;
  //tree->Print();
  fitResult* res=0;//new fitResult();
  string branchname("fitResult_v2");
  if(tree->FindBranch(branchname.c_str())==NULL){
    cerr << "Invalid branch ."<<branchname<<endl;
    throw;
   }
  
  tree->SetBranchAddress(branchname.c_str(),&res);
  cout << "Entries: " <<tree->GetEntries() << endl;
  tree->GetEntry(1);
  // get list of waves
 
  

   fwavenames=res->waveNames();
  


   // Create main frame
   
   fListBox = new TGListBox(this, 89);
   fListBox2 = new TGListBox(this, 88);
   fSelected = new TList;
   
   for (int i = 0; i < fwavenames.size(); ++i) {
     cout << fwavenames[i] << endl;
     fListBox->AddEntry(fwavenames[i].c_str(), i);
     fListBox2->AddEntry(fwavenames[i].c_str(), i);
   }
   fListBox->Resize(400,250);
   fListBox2->Resize(400,250);
   
   AddFrame(fListBox, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 5, 5, 5, 5));
   AddFrame(fListBox2, new TGLayoutHints(kLHintsTop | kLHintsRight| kLHintsExpandX, 5, 5, 5, 5));        

   
               
   fCheckMulti = new TGCheckButton(this, "&Mutliple selection", 10);
   AddFrame(fCheckMulti, new TGLayoutHints(kLHintsTop | kLHintsLeft,
                                           5, 5, 5, 5));
   fCheckMulti->Connect("Clicked()", "MyMainFrame", this, "HandleButtons()"); 
   // Create a horizontal frame containing button(s)
   TGHorizontalFrame *hframe = new TGHorizontalFrame(this, 400, 20, kFixedWidth);
   TGTextButton *show = new TGTextButton(hframe, "&Show");
   show->SetToolTipText("Click here to print the selection you made");
   show->Connect("Pressed()", "MyMainFrame", this, "PrintSelected()");
   hframe->AddFrame(show, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));
   TGTextButton *exit = new TGTextButton(hframe, "&Exit ");
   exit->Connect("Pressed()", "MyMainFrame", this, "DoExit()");
   hframe->AddFrame(exit, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));
   AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 2, 2, 5, 1));

   // Set a name to the main frame   
   SetWindowName("List Box");
   MapSubwindows();

   // Initialize the layout algorithm via Resize()
   Resize(GetDefaultSize());

   // Map main frame
   MapWindow();
   fListBox->Select(0);
   fListBox2->Select(0);
}

MyMainFrame::~MyMainFrame()
{
   // Clean up main frame...
   Cleanup();
   if (fSelected) {
      fSelected->Delete();
      delete fSelected;
   }
}

void MyMainFrame::HandleButtons()
{
   // Handle check button.
   Int_t id;
   TGButton *btn = (TGButton *) gTQSender;
   id = btn->WidgetId();

   printf("HandleButton: id = %d\n", id);

   //if (id == 10)  
   //   fListBox->SetMultipleSelections(fCheckMulti->GetState());
}


void MyMainFrame::PrintSelected()
{
   // Writes selected entries in TList if multiselection.
   string w1=fwavenames[fListBox->GetSelected()];
   string w2=fwavenames[fListBox2->GetSelected()];
   cout << w1 << endl;
   cout << w2 << endl;
   // Produce plots
   unsigned int n=fTree->GetEntries();
   
   Int_t colour=1;
   double binwidth=0.060;

   TGraphErrors* gph = new TGraphErrors(n);
   stringstream graphName;
   graphName << "PHI"<<w1<<"---"<<"PHI"<<w2;
   cout << "creating graph   " << graphName.str() << endl;
   gph->SetName (graphName.str().c_str());
   gph->SetTitle(graphName.str().c_str());
   gph->SetMarkerStyle(21);
   gph->SetMarkerSize(0.5);
   gph->SetMarkerColor(colour);
   gph->SetLineColor(colour);

   TGraphErrors* gphP1 = (TGraphErrors*)gph->Clone("gph+1");
   TGraphErrors* gphM1 = (TGraphErrors*)gph->Clone("gph-1");
   gphP1->SetMarkerColor(2);
   gphP1->SetLineColor(2);
   gphM1->SetMarkerColor(3);
   gphM1->SetLineColor(3);

   TGraphErrors* gRe = new TGraphErrors(n);
   graphName.str("");
   graphName.clear();
   graphName << "RE_"<<w1<<"---"<<""<<w2;
   cout << "creating graph   " << graphName.str() << endl;
   gRe->SetName (graphName.str().c_str());
   gRe->SetTitle(graphName.str().c_str());
   gRe->SetMarkerStyle(21);
   gRe->SetMarkerSize(0.5);
   gRe->SetMarkerColor(colour);
   gRe->SetLineColor(colour);

   TGraphErrors* gIm = new TGraphErrors(n);
   graphName.str("");
   graphName.clear();
   graphName << "IM_"<<w1<<"---"<<""<<w2;
   cout << "creating graph   " << graphName.str() << endl;
   gIm->SetName (graphName.str().c_str());
   gIm->SetTitle(graphName.str().c_str());
   gIm->SetMarkerStyle(21);
   gIm->SetMarkerSize(0.5);
   gIm->SetMarkerColor(colour);
   gIm->SetLineColor(colour);
   
   TGraphErrors* g1 = new TGraphErrors(n);
   graphName.str("");
   graphName.clear();
   graphName<<"g"<<w1;
   g1->SetName (graphName.str().c_str());
   g1->SetTitle(graphName.str().c_str());
   g1->SetMarkerStyle(21);
   g1->SetMarkerSize(0.5);
   TGraphErrors* g2 = new TGraphErrors(n);
   graphName.str("");
   graphName.clear();
   graphName<<"g"<<w2;
   g2->SetName (graphName.str().c_str());
   g2->SetTitle(graphName.str().c_str());
   g2->SetMarkerStyle(21);
   g2->SetMarkerSize(0.5);




   fitResult* result=0;
   fTree->SetBranchAddress("fitResult_v2",&result);

   

   for(unsigned int i=0;i<n;++i){
     fTree->GetEntry(i);
     
     if(!result->converged())continue;
     if(!result->hasHessian())continue;
 
     double intensity1=result->intensity(w1.c_str());
     if((numeric_limits<double>::has_infinity && intensity1 == numeric_limits<double>::infinity()) || intensity1!=intensity1)continue;

     g1->SetPoint(i,
		 result->massBinCenter()*0.001,
		 intensity1);
     g1->SetPointError(i,
     		      binwidth*0.5,
     		      result->intensityErr(w1.c_str()));
     
     double intensity2=result->intensity(w2.c_str());
     if((numeric_limits<double>::has_infinity && intensity2 == numeric_limits<double>::infinity()) || intensity2!=intensity2)continue;

     g2->SetPoint(i,
		 result->massBinCenter()*0.001,
		 intensity2);
    
     g2->SetPointError(i,
     		      binwidth*0.5,
     		      result->intensityErr(w2.c_str()));
      

     double ph=result->phase(w1.c_str(),w2.c_str());
     double pherr=result->phaseErr(w1.c_str(),w2.c_str());
     // check if we should make a transformation by 2pi
     // this is needed because of cyclical variable phi
     if(i>11){
       double mpre;
       double phpre;
       gph->GetPoint(i-1,mpre,phpre);
       double diff1=fabs(ph-phpre);
       double diff2=fabs(ph+360-phpre);
       double diff3=fabs(ph-360-phpre);
       if(diff2<diff1 && diff2<diff3)ph+=360;
       else if(diff3<diff1 && diff3<diff2)ph-=360;
     }
     
     
     
     gph->SetPoint(i,
		 result->massBinCenter()*0.001,
		 ph);
     gph->SetPointError(i,
		      binwidth*0.5,
		      pherr);
     
     // add point +- 360 degree
     gphP1->SetPoint(i,
		 result->massBinCenter()*0.001,
		 ph+360);
     gphP1->SetPointError(i,
		      binwidth*0.5,
		      pherr);
     
     gphM1->SetPoint(i,
		 result->massBinCenter()*0.001,
		 ph-360);
     gphM1->SetPointError(i,
			binwidth*0.5,
			pherr);
     
     unsigned int wi1=result->waveIndex(w1);
     unsigned int wi2=result->waveIndex(w2);
     complex<double> rho=result->spinDensityMatrixElem(wi1,wi2);
     TMatrixT<double> rhoCov=result->spinDensityMatrixElemCov(wi1,wi2);
     gRe->SetPoint(i,
		   result->massBinCenter()*0.001,
		   rho.real());
     gRe->SetPointError(i,
			binwidth*0.5,
			sqrt(rhoCov[0][0]));


     gIm->SetPoint(i,
		   result->massBinCenter()*0.001,
		   rho.imag());
     
      gIm->SetPointError(i,
			 binwidth*0.5,
			 sqrt(rhoCov[1][1]));
   }// end loop over bins
   
   
   // plot graphs
   
   TCanvas*c=new TCanvas("c","c",10,10,1200,800);
   c->Divide(2,3);
   c->cd(1);
   gph->Draw("AP");
   gph->GetXaxis()->SetTitle("5#pi mass (GeV/c^2)");
   gph->GetYaxis()->SetTitle("Phase difference");
   gph->GetYaxis()->SetRangeUser(-270,270);

   gphP1->Draw("PSAME");
   gphM1->Draw("PSAME");

   c->cd(3);
   g1->Draw("AP");
   g1->GetXaxis()->SetTitle("5#pi mass (GeV/c^2)");
   g1->GetYaxis()->SetTitle("Intensity");
   c->cd(5);
   g2->Draw("AP");
   g2->GetXaxis()->SetTitle("5#pi mass (GeV/c^2)");
   g2->GetYaxis()->SetTitle("Intensity");
   
   c->cd(2);
   gRe->Draw("AP");
   gRe->GetXaxis()->SetTitle("5#pi mass (GeV/c^2)");
   gRe->GetYaxis()->SetTitle("Re(#rho_{ij})");
   c->cd(4);
   gIm->Draw("AP");
   gIm->GetXaxis()->SetTitle("5#pi mass (GeV/c^2)");
   gIm->GetYaxis()->SetTitle("Im(#rho_{ij})");

}

void plotGui(TString infilename)
{

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
  
   // Popup the GUI...
  new MyMainFrame(gClient->GetRoot(), 20, 20, pwa);
}
