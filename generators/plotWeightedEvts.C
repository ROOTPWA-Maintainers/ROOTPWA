#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TLorentzRotation.h"
#include "TH1D.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TMath.h"
#include "TPad.h"
#include "TCanvas.h"
#include <iostream>
#include <vector>
#include "NParticleEvent.h"


using namespace std; 


void plotWeightedEvents(TTree* mctr, TTree* datatr){

gROOT->SetStyle("Plain");

 vector<TH1D*> hM;
 TH1D* hMMC=new TH1D("hMMC","Mass (MC)",100,0.5,2.5);
 hM.push_back(hMMC);
 TH1D* hMData=new TH1D("hMData","Mass (DATA)",100,0.5,2.5);
 hM.push_back(hMData);

 vector<TH1D*> hMIsobar;
 TH1D* hMIsobarMC=new TH1D("hMIsobarMC","Isobar Mass (MC)",100,0.2,2.0);
 hMIsobar.push_back(hMIsobarMC);
 TH1D* hMIsobarData=new TH1D("hMIsobarData","Isobar mass (DATA)",100,0.2,2.0);
 hMIsobar.push_back(hMIsobarData);

 vector<TH1D*> hMIsobar2;
 TH1D* hMIsobar2MC=new TH1D("hMIsobar2MC","Isobar Mass (MC)",100,0.4,2.0);
 hMIsobar2.push_back(hMIsobar2MC);
 TH1D* hMIsobar2Data=new TH1D("hMIsobar2Data","Isobar mass (DATA)",100,0.4,2.0);
 hMIsobar2.push_back(hMIsobar2Data);


vector<TH1D*> hMIsobar3;
 TH1D* hMIsobar3MC=new TH1D("hMIsobar3MC","Isobar Mass (MC)",100,0.2,1.8);
 hMIsobar3.push_back(hMIsobar3MC);
 TH1D* hMIsobar3Data=new TH1D("hMIsobar3Data","Isobar mass (DATA)",100,0.2,1.8);
 hMIsobar3.push_back(hMIsobar3Data);


  vector<TH1D*> hGJ;
 TH1D* hGJMC=new TH1D("hGJMC","Cos Gottfried-Jackson Theta (MC)",20,-1,1);
 hGJ.push_back(hGJMC);
 TH1D* hGJData=new TH1D("hGJData","Cos Gottfried-Jackson Theta (DATA)",20,-1,1); hGJ.push_back(hGJData);

  vector<TH1D*> hGJ2;
 TH1D* hGJ2MC=new TH1D("hGJ2MC","Cos Gottfried-Jackson Theta (MC)",20,-1,1);
 hGJ2.push_back(hGJ2MC);
 TH1D* hGJ2Data=new TH1D("hGJ2Data","Cos Gottfried-Jackson Theta (DATA)",20,-1,1); hGJ2.push_back(hGJ2Data);

  vector<TH1D*> hGJ3;
 TH1D* hGJ3MC=new TH1D("hGJ3MC","Cos Gottfried-Jackson Theta (MC)",20,-1,1);
 hGJ3.push_back(hGJ3MC);
 TH1D* hGJ3Data=new TH1D("hGJ3Data","Cos Gottfried-Jackson Theta (DATA)",20,-1,1); hGJ3.push_back(hGJ3Data);

 
 vector<TH1D*> hTY;
 TH1D* hTYMC=new TH1D("hTYMC","Treiman-Yang Phi (MC)",20,-TMath::Pi(),TMath::Pi());
 TH1D* hTYData=new TH1D("hTYMC","Treiman-Yang Phi (DATA)",20,-TMath::Pi(),TMath::Pi());
 hTY.push_back(hTYMC);
 hTY.push_back(hTYData);


vector<TH1D*> hTY2;
 TH1D* hTY2MC=new TH1D("hTY2MC","Treiman-Yang Phi (MC)",20,-TMath::Pi(),TMath::Pi());
 TH1D* hTY2Data=new TH1D("hTY2MC","Treiman-Yang Phi (DATA)",20,-TMath::Pi(),TMath::Pi());
 hTY2.push_back(hTY2MC);
 hTY2.push_back(hTY2Data);

vector<TH1D*> hTY3;
 TH1D* hTY3MC=new TH1D("hTY3MC","Treiman-Yang Phi (MC)",20,-TMath::Pi(),TMath::Pi());
 TH1D* hTY3Data=new TH1D("hTY3MC","Treiman-Yang Phi (DATA)",20,-TMath::Pi(),TMath::Pi());
 hTY3.push_back(hTY3MC);
 hTY3.push_back(hTY3Data);


 double avweight=1; 

//Loop both over data and mc tree.
 for(unsigned int itree=0;itree<2;++itree){
   TTree* tr= itree==0 ? mctr : datatr;
   if(tr==NULL)continue;

   double weight=1;
   double maxweight=0; 
   TClonesArray* p=new TClonesArray("TLorentzVector");
   TLorentzVector* beam=NULL;
   int qbeam;
   std::vector<int>* q=NULL; 
   if (itree==0)tr->SetBranchAddress("weight",&weight);
   tr->SetBranchAddress("p",&p);
   tr->SetBranchAddress("beam",&beam);
   tr->SetBranchAddress("qbeam",&qbeam);
   tr->SetBranchAddress("q",&q);

 
   TVector3 vertex;

   NParticleEvent event(p,q,beam,&qbeam,&vertex);


   unsigned int nevt=tr->GetEntries();	
   for(unsigned int i=0;i<nevt;++i){
     tr->GetEntry(i);
     if(itree==1)weight=1;
     event.refresh();
     
     if(weight>maxweight)maxweight=weight;
     if(itree==0)avweight+=weight;
     // transform into GJ 
     event.toGJ();
     
     // loop over all states that contain n-1 final state particles
     // and plot GJ angles
     unsigned int npart=event.nParticles();
     unsigned int nstates=event.nStates();
     for(unsigned int is=0;is<nstates;++is){
       const NParticleState& state=event.getState(is);
       if(state.n()==npart){
	 
	 
	 hM[itree]->Fill(state.p().M(),weight);
	 
       }

       if(state.n()==npart-1 && state.q()==0){
	 
	 hGJ[itree]->Fill(state.p().CosTheta(),weight);
	 hTY[itree]->Fill(state.p().Phi(),weight);
	 hMIsobar[itree]->Fill(state.p().M(),weight);
	 
       }
       else if(state.n()==npart-2 && state.q()==-1){
	 hGJ2[itree]->Fill(state.p().CosTheta(),weight);
	 hTY2[itree]->Fill(state.p().Phi(),weight);
	 hMIsobar2[itree]->Fill(state.p().M(),weight);
       }
       else if(state.n()==npart-3 && state.q()==0){
	 hGJ3[itree]->Fill(state.p().CosTheta(),weight);
	 hTY3[itree]->Fill(state.p().Phi(),weight);
	 hMIsobar3[itree]->Fill(state.p().M(),weight);
       }
     }
     
     
     
   }// end loop over events
   if(itree==0)avweight/=(double)nevt;
   cout << "Maxweight=" << maxweight << endl; 
   cout << "Average weight=" << avweight << endl; 
 }
 TCanvas* cm=new TCanvas("PredictM","Weighted Events",20,20,600,800);
 hM[0]->SetLineColor(kRed);
 hM[0]->Draw();
 hM[1]->Draw("same");
 

 TCanvas* c=new TCanvas("Predict","Weighted Events",10,10,600,800);
 c->Divide(3,3);
 c->cd(1);
 hMIsobar[0]->SetLineColor(kRed);
 hMIsobar[0]->Draw();
 double totMC=hMIsobar[0]->Integral();
 double totDATA=hMIsobar[1]->Integral();
 hMIsobar[1]->Scale(totMC/totDATA);
 hMIsobar[1]->Draw("same");
 hMIsobar[0]->GetYaxis()->SetRangeUser(0,hMIsobar[0]->GetMaximum()*1.1);
 gPad->Update();

 c->cd(2);
 hMIsobar2[0]->SetLineColor(kRed);
 hMIsobar2[0]->Draw();
 totMC=hMIsobar2[0]->Integral();
 totDATA=hMIsobar2[1]->Integral();
 hMIsobar2[1]->Scale(totMC/totDATA);
 hMIsobar2[1]->Draw("same");
 hMIsobar2[0]->GetYaxis()->SetRangeUser(0,hMIsobar2[0]->GetMaximum()*1.1);
 gPad->Update();
 
 c->cd(3);
 hMIsobar3[0]->SetLineColor(kRed);
 hMIsobar3[0]->Draw();
 totMC=hMIsobar3[0]->Integral();
 totDATA=hMIsobar3[1]->Integral();
 hMIsobar3[1]->Scale(totMC/totDATA);
 hMIsobar3[1]->Draw("same");
 hMIsobar3[0]->GetYaxis()->SetRangeUser(0,hMIsobar3[0]->GetMaximum()*1.1);
 gPad->Update();
 


 c->cd(4);
 hGJ[0]->SetLineColor(kRed);
 hGJ[0]->Draw();
 
 totMC=hGJ[0]->Integral();
 totDATA=hGJ[1]->Integral();
 hGJ[1]->Sumw2();
 hGJ[1]->Scale(totMC/totDATA);
 hGJ[1]->Draw("same E");
 hGJ[0]->GetYaxis()->SetRangeUser(0,hGJ[0]->GetMaximum()*1.1);
 gPad->Update();

 c->cd(5);
 hGJ2[0]->SetLineColor(kRed);
 hGJ2[0]->Draw();
 totMC=hGJ2[0]->Integral();
 totDATA=hGJ2[1]->Integral();
 hGJ2[1]->Sumw2();
 hGJ2[1]->Scale(totMC/totDATA);
 hGJ2[1]->Draw("same E");
 hGJ2[0]->GetYaxis()->SetRangeUser(0,hGJ2[0]->GetMaximum()*1.1);
 gPad->Update();

 c->cd(6);
 hGJ3[0]->SetLineColor(kRed);
 hGJ3[0]->Draw();
 totMC=hGJ3[0]->Integral();
 totDATA=hGJ3[1]->Integral();
 hGJ3[1]->Sumw2();
 hGJ3[1]->Scale(totMC/totDATA);
 hGJ3[1]->Draw("same E");
 hGJ3[0]->GetYaxis()->SetRangeUser(0,hGJ3[0]->GetMaximum()*1.1);
 gPad->Update();


 c->cd(7);
 hTY[0]->SetLineColor(kRed);
 hTY[0]->Draw();
 
 totMC=hTY[0]->Integral();
 totDATA=hTY[1]->Integral();
 hTY[1]->Sumw2();
 hTY[1]->Scale(totMC/totDATA);
 hTY[1]->Draw("same E");
 hTY[0]->GetYaxis()->SetRangeUser(0,hTY[0]->GetMaximum()*1.1);
 gPad->Update();

 c->cd(8);
 hTY2[0]->SetLineColor(kRed);
 hTY2[0]->Draw();
 totMC=hTY2[0]->Integral();
 totDATA=hTY2[1]->Integral();
 hTY2[1]->Sumw2();
 hTY2[1]->Scale(totMC/totDATA);
 hTY2[1]->Draw("same E");
 hTY2[0]->GetYaxis()->SetRangeUser(0,hTY2[0]->GetMaximum()*1.1);
 gPad->Update();


c->cd(9);
 hTY3[0]->SetLineColor(kRed);
 hTY3[0]->Draw();
 totMC=hTY3[0]->Integral();
 totDATA=hTY3[1]->Integral();
 hTY3[1]->Sumw2();
 hTY3[1]->Scale(totMC/totDATA);
 hTY3[1]->Draw("same E");
 hTY3[0]->GetYaxis()->SetRangeUser(0,hTY3[0]->GetMaximum()*1.1);
 gPad->Update();

 
 
}
