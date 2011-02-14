#include "TFile.h"
#include "TTree.h"
#include "TPDGEntry.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TLegend.h"
#include "TH1D.h"
#include <vector>
using namespace std;


bool cut(TPDGEntry* entry){
  TString excludes[10]={"rho(1570)","f(1)(1510)","eta(1475)","f(2)(1430)","f(1)(1420)","eta(1405)","phi(1020)","f(0)(1500)","f(2)'(1525)"};

  bool excluded=false;
  for(unsigned int i=0;i<9;++i){
    if(entry->name().Contains(excludes[i])){
      excluded=true;
      cout << "Excluding " << entry->name() << endl;
      break;
    }
  }

  return !entry->isLightMeson() || entry->J()<0 || entry->J()==0.5 || entry->mass()>2500 || entry->status()<0 || entry->isExotic() || entry->name().Contains("X") || entry->q()!=0;
}


int 
plotPdg(TString file){
  TFile* infile=TFile::Open(file,"READ");

  TTree* tree=(TTree*)infile->Get("pdg");

  tree->Print();

  TPDGEntry* entry=new TPDGEntry();
  tree->SetBranchAddress("TPDGEntry",&entry);


  TH1D* hm=new TH1D("hm","mass",100,0,6);
  
  
  unsigned int nI0Pp=0;
  unsigned int nI0Pm=0;
  unsigned int nI1Pp=0;
  unsigned int nI1Pm=0;
  unsigned int n=0;

  unsigned int nentries=tree->GetEntries();
  for(unsigned int i=0;i<nentries;++i){
    tree->GetEntry(i);
    if(cut(entry))continue;
    hm->Fill(entry->mass()/1000);

    if(entry->I()==0){
      if(entry->P()==1)++nI0Pp;
      else if(entry->P()==-1)++nI0Pm;

    }
    else if(entry->I()==1){
    if(entry->P()==1)++nI1Pp;
      else if(entry->P()==-1)++nI1Pm;
    }

    n++;

  }

    TGraphErrors* gRegge=new TGraphErrors(n);
  TGraphErrors* gReggeI0Pp=new TGraphErrors(nI0Pp);
  TGraphErrors* gReggeI0Pm=new TGraphErrors(nI0Pm);
  TGraphErrors* gReggeI1Pp=new TGraphErrors(nI1Pp);
  TGraphErrors* gReggeI1Pm=new TGraphErrors(nI1Pm);
  TMultiGraph* mg=new TMultiGraph();
  mg->Add(gReggeI0Pp);
  mg->Add(gReggeI0Pm);
  mg->Add(gReggeI1Pp);
  mg->Add(gReggeI1Pm);
  TLegend* leg=new TLegend(0.58,0.15,0.88,0.4);
  leg->SetFillColor(0);
  leg->AddEntry(gReggeI0Pp,"I=0, P=+","P");
  leg->AddEntry(gReggeI0Pm,"I=0, P=-","P");
  leg->AddEntry(gReggeI1Pp,"I=1, P=+","P");
  leg->AddEntry(gReggeI1Pm,"I=1, P=-","P");
  


  unsigned int counter=0;
  unsigned int cI0Pp=0;
  unsigned int cI0Pm=0;
  unsigned int cI1Pp=0;
  unsigned int cI1Pm=0;
  

  for(unsigned int i=0;i<nentries;++i){
    tree->GetEntry(i);
    if(cut(entry))continue;

    gRegge->SetPoint(counter,entry->J(),entry->mass()*entry->mass()*1E-6);
    gRegge->SetPointError(counter,0,2.*entry->mass()*1E-6*entry->mass_er());
    ++counter;
    entry->Print();

   if(entry->I()==0){
      if(entry->P()==1){
	gReggeI0Pp->SetPoint(cI0Pp,entry->J()-0.15,entry->mass()*entry->mass()*1E-6);
	gReggeI0Pp->SetPointError(cI0Pp,0,2.*entry->mass()*1E-6*entry->mass_er());
	++cI0Pp;
      }
      else if(entry->P()==-1){
	 gReggeI0Pm->SetPoint(cI0Pm,entry->J()-0.05,entry->mass()*entry->mass()*1E-6);
	 gReggeI0Pm->SetPointError(cI0Pm,0,2.*entry->mass()*1E-6*entry->mass_er());
	++cI0Pm;
      }

    }
    else if(entry->I()==1){
      if(entry->P()==1){
	gReggeI1Pp->SetPoint(cI1Pp,entry->J()+0.05,entry->mass()*entry->mass()*1E-6);
	gReggeI1Pp->SetPointError(cI1Pp,0,2.*entry->mass()*1E-6*entry->mass_er());
	++cI1Pp;
      }
      else if(entry->P()==-1){
	gReggeI1Pm->SetPoint(cI1Pm,entry->J()+0.15,entry->mass()*entry->mass()*1E-6);
	gReggeI1Pm->SetPointError(cI1Pm,0,2.*entry->mass()*1E-6*entry->mass_er());
	++cI1Pm;
      }
    }
   
   

  }

  gRegge->SetMarkerStyle(20);
  //gRegge->Draw("AP");

  gReggeI0Pp->SetMarkerStyle(24);
  gReggeI0Pm->SetMarkerStyle(20);
  gReggeI1Pp->SetMarkerStyle(30);
  gReggeI1Pm->SetMarkerStyle(29);gReggeI1Pm->SetMarkerSize(1.5);

  mg->Draw("AP");
  mg->GetXaxis()->SetTitle("J");
  mg->GetYaxis()->SetTitle("m^{2} (GeV^{2}/c^{4})");
  leg->Draw();
  //hm->DrawClone();

  return 0;





}
