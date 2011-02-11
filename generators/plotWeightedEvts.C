#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TLatex.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TLorentzRotation.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TMath.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TString.h"
#include <iostream>
#include <vector>
#include "NParticleEvent.h"


using namespace std; 


void plotWeightedEvts(TTree* mctr, TTree* datatr, TString outfilename="kineplots.root", TString mass="000"){

gROOT->SetStyle("Plain");

 TString massbin;
 if(mass!="000"){
   massbin=("_m")+mass;massbin.ReplaceAll(" ","");
 }

 int nbninsm=144;
 int nbinsang=80;

 vector<TH1D*> hM;
 TH1D* hMMC=new TH1D("hMMC"+massbin,"Mass (MC)",60,1.0,3.0);
 hM.push_back(hMMC);
 TH1D* hMData=new TH1D("hMData"+massbin,"Mass (DATA)",60,1.0,3.0);
 hM.push_back(hMData);

 vector<TH1D*> hMIsobar;
 TH1D* hMIsobarMC=new TH1D("hMIsobarMC"+massbin,"Isobar Mass (MC) for m_{5#pi}="+mass,nbninsm,0.2,3.8);
 hMIsobar.push_back(hMIsobarMC);
 TH1D* hMIsobarData=new TH1D("hMIsobarData"+massbin,"Isobar mass (DATA)",nbninsm,0.2,3.8);
 hMIsobar.push_back(hMIsobarData);
 TH1D* hDiffMIsobar; //=new TH1D("hDiffMIsobarData","Isobar mass (DATA) Diff",40,0.2,2.2);
TH1D* hDiffMIsobar2;
TH1D* hDiffMIsobar3;
hMIsobar[0]->Sumw2();

 vector<TH1D*> hMIsobar2;
 TH1D* hMIsobar2MC=new TH1D("hMIsobar2MC"+massbin,"Isobar Mass (MC)",nbninsm,0.4,3.6);
 hMIsobar2.push_back(hMIsobar2MC);
 TH1D* hMIsobar2Data=new TH1D("hMIsobar2Data"+massbin,"Isobar mass (DATA)",nbninsm,0.4,3.6);
 hMIsobar2.push_back(hMIsobar2Data);
hMIsobar2[0]->Sumw2();

vector<TH1D*> hMIsobar3;
 TH1D* hMIsobar3MC=new TH1D("hMIsobar3MC"+massbin,"Isobar Mass (MC)",nbninsm,0.2,3.4);
 hMIsobar3.push_back(hMIsobar3MC);
 TH1D* hMIsobar3Data=new TH1D("hMIsobar3Data"+massbin,"Isobar mass (DATA)",nbninsm,0.2,3.4);
 hMIsobar3.push_back(hMIsobar3Data);
hMIsobar3[0]->Sumw2();

  vector<TH1D*> hGJ;
 TH1D* hGJMC=new TH1D("hGJMC"+massbin,"Cos Gottfried-Jackson Theta (MC)",nbinsang,-1,1);
 hGJ.push_back(hGJMC);
 TH1D* hGJData=new TH1D("hGJData"+massbin,"Cos Gottfried-Jackson Theta (DATA)",nbinsang,-1,1); hGJ.push_back(hGJData);

 TH1D* hGJMCraw=new TH1D("hGJMC"+massbin,"Cos Gottfried-Jackson Theta (MC raw)",nbinsang,-1,1);
 
hGJ[0]->Sumw2();


 vector<TH1D*> hThetaLab;
 TH1D* hThetaLabMC=new TH1D("hThetaLabMC"+massbin,"Cos Theta Lab (MC)",nbninsm,0.997,1);
 hThetaLab.push_back(hThetaLabMC);
 TH1D* hThetaLabData=new TH1D("hThetaLabData"+massbin,"Cos Theta Lab (Data)",nbninsm,0.997,1);
 hThetaLab.push_back(hThetaLabData);
 hThetaLab[0]->Sumw2();

 vector<TH2D*> hGJt;
 TH2D* hGJtMC=new TH2D("hGJtMC"+massbin,"Cos GJ Theta vs t' (MC)",nbinsang,-1,1,20,0,0.01);
 hGJt.push_back(hGJtMC);
 TH2D* hGJtData=new TH2D("hGJtData"+massbin,"Cos GJ Theta vs t' (DATA)",nbinsang,-1,1,20,0,0.01); hGJt.push_back(hGJtData);



  vector<TH1D*> hGJ2;
 TH1D* hGJ2MC=new TH1D("hGJ2MC"+massbin,"Cos Gottfried-Jackson Theta (MC)",nbinsang,-1,1);
 hGJ2.push_back(hGJ2MC);
 TH1D* hGJ2Data=new TH1D("hGJ2Data"+massbin,"Cos Gottfried-Jackson Theta (DATA)",nbinsang,-1,1); hGJ2.push_back(hGJ2Data);
hGJ2[0]->Sumw2();

  vector<TH1D*> hGJ3;
 TH1D* hGJ3MC=new TH1D("hGJ3MC"+massbin,"Cos Gottfried-Jackson Theta (MC)",nbinsang,-1,1);
 hGJ3.push_back(hGJ3MC);
 TH1D* hGJ3Data=new TH1D("hGJ3Data"+massbin,"Cos Gottfried-Jackson Theta (DATA)",nbinsang,-1,1); hGJ3.push_back(hGJ3Data);
hGJ3[0]->Sumw2();
 
 vector<TH1D*> hTY;
 TH1D* hTYMC=new TH1D("hTYMC"+massbin,"Treiman-Yang Phi (MC)",nbinsang,-TMath::Pi(),TMath::Pi());
 TH1D* hTYData=new TH1D("hTYMC"+massbin,"Treiman-Yang Phi (DATA)",nbinsang,-TMath::Pi(),TMath::Pi());
 hTY.push_back(hTYMC);
 hTY.push_back(hTYData);
hTY[0]->Sumw2();

vector<TH1D*> hTY2;
 TH1D* hTY2MC=new TH1D("hTY2MC"+massbin,"Treiman-Yang Phi (MC)",nbinsang,-TMath::Pi(),TMath::Pi());
 TH1D* hTY2Data=new TH1D("hTY2MC"+massbin,"Treiman-Yang Phi (DATA)",nbinsang,-TMath::Pi(),TMath::Pi());
 hTY2.push_back(hTY2MC);
 hTY2.push_back(hTY2Data);
hTY2[0]->Sumw2();


vector<TH1D*> hTY3;
 TH1D* hTY3MC=new TH1D("hTY3MC"+massbin,"Treiman-Yang Phi (MC)",nbinsang,-TMath::Pi(),TMath::Pi());
 TH1D* hTY3Data=new TH1D("hTY3MC"+massbin,"Treiman-Yang Phi (DATA)",nbinsang,-TMath::Pi(),TMath::Pi());
 hTY3.push_back(hTY3MC);
 hTY3.push_back(hTY3Data);
hTY3[0]->Sumw2();


 double avweight=1; 

//Loop both over data and mc tree.
 for(unsigned int itree=0;itree<2;++itree){
   TTree* tr= itree==0 ? mctr : datatr;
   if(tr==NULL)continue;

   double weight=1;double impweight=1;
   double maxweight=0; 
   TClonesArray* p=new TClonesArray("TLorentzVector");
   TLorentzVector* beam=NULL;
   int qbeam;
   std::vector<int>* q=NULL; 
   if (itree==0)tr->SetBranchAddress("weight",&weight);
   if (itree==0)tr->SetBranchAddress("impweight",&impweight);
   tr->SetBranchAddress("p",&p);
   tr->SetBranchAddress("beam",&beam);
   tr->SetBranchAddress("qbeam",&qbeam);
   tr->SetBranchAddress("q",&q);

 
   TVector3 vertex;

   NParticleEvent event(p,q,beam,&qbeam,&vertex);


   unsigned int nevt=tr->GetEntries();	
   for(unsigned int i=0;i<nevt;++i){
     tr->GetEntry(i);
     if(itree==1){weight=1;impweight=1;}
     if(weight==0)continue;

     event.refresh();
     
     if(impweight!=0)weight/=impweight;

     if(weight>maxweight)maxweight=weight;
     if(itree==0)avweight+=weight;

     double tprime=event.tprime() ;
     //cout << tprime << endl;
     unsigned int npart=event.nParticles();
     // loop over pions
     for(unsigned int ipart=0;ipart<npart;++ipart){
       TLorentzVector p=event.getParticle(ipart).p();
       hThetaLab[itree]->Fill(p.CosTheta(),weight);
     }


     // transform into GJ 
     event.toGJ();
     
     // loop over all states that contain n-1 final state particles
     // and plot GJ angles
    
     unsigned int nstates=event.nStates();
     for(unsigned int is=0;is<nstates;++is){
       
       const NParticleState& state=event.getState(is);
       if(state.n()==npart){
	 
	 
	 hM[itree]->Fill(state.p().M());
	 
       }

       if(state.n()==npart-1 && state.q()==0){
	 
	 hGJ[itree]->Fill(state.p().CosTheta(),weight);
	 if(itree==0)hGJMCraw->Fill(state.p().CosTheta());

	 hGJt[itree]->Fill(state.p().CosTheta(),tprime,1);
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




 TCanvas* cm=new TCanvas("PredictM"+massbin,"Weighted Events",20,20,600,800);
 if (!cm) cout << " äh " << endl;
 hM[0]->SetLineColor(kRed);
 hM[0]->Draw();
 hM[1]->Draw("same");
 
 TCanvas* gjt=new TCanvas("GJT"+massbin,"Events",40,40,600,800);
 // hGJt[0]->SetLineColor(kRed);
 hGJt[1]->Draw("COLZ");
 //hM[1]->Draw("same");


 TCanvas* c=new TCanvas("KineValidate"+massbin,"Weighted Events",10,10,600,800);
 c->Divide(3,3);
 c->cd(1);
 
 double totMC=hMIsobar[0]->Integral();
 double totDATA=hMIsobar[1]->Integral();
 
 if(totMC!=0)hMIsobar[0]->Scale(totDATA/totMC);
 hMIsobar[0]->SetLineColor(kRed);
 hMIsobar[0]->SetFillColor(kRed);
 hMIsobar[0]->Draw("E4");
 
 hDiffMIsobar=new TH1D(*hMIsobar[0]);
 hDiffMIsobar->Add(hMIsobar[1],-1.);
 hMIsobar[1]->Draw("same");
 hDiffMIsobar->SetLineColor(kOrange-3);
 hDiffMIsobar->Draw("same");

 hMIsobar[0]->GetYaxis()->SetRangeUser(hDiffMIsobar->GetMinimum()*1.5,hMIsobar[0]->GetMaximum()*1.2);
 gPad->Update();

 c->cd(2);

 totMC=hMIsobar2[0]->Integral();
 totDATA=hMIsobar2[1]->Integral();
 //hMIsobar2[0]->Sumw2();
 if(totMC!=0)hMIsobar2[0]->Scale(totDATA/totMC);
 hMIsobar2[0]->SetLineColor(kRed);
 hMIsobar2[0]->SetFillColor(kRed);
 hMIsobar2[0]->Draw("E4");
 
 hMIsobar2[1]->Draw("same");

 hDiffMIsobar2=new TH1D(*hMIsobar2[0]);
 hDiffMIsobar2->Add(hMIsobar2[1],-1.);
 hDiffMIsobar2->SetLineColor(kOrange-3);
 hDiffMIsobar2->Draw("same");


 hMIsobar2[0]->GetYaxis()->SetRangeUser(hDiffMIsobar2->GetMinimum()*1.5,hMIsobar2[0]->GetMaximum()*1.2);
 gPad->Update();
 
 c->cd(3);
 hMIsobar3[0]->Sumw2();
 totMC=hMIsobar3[0]->Integral();
 totDATA=hMIsobar3[1]->Integral();
 if(totMC!=0)hMIsobar3[0]->Scale(totDATA/totMC);
 hMIsobar3[0]->SetLineColor(kRed);
 hMIsobar3[0]->SetFillColor(kRed);
 hMIsobar3[0]->Draw("E4");
 hMIsobar3[1]->Draw("same");

 hDiffMIsobar3=new TH1D(*hMIsobar3[0]);
 hDiffMIsobar3->Add(hMIsobar3[1],-1.);
 hDiffMIsobar3->SetLineColor(kOrange-3);
 hDiffMIsobar3->Draw("same");


 hMIsobar3[0]->GetYaxis()->SetRangeUser(hDiffMIsobar3->GetMinimum()*1.5,hMIsobar3[0]->GetMaximum()*1.2);
 gPad->Update();
 


 c->cd(4);
 totMC=hGJ[0]->Integral();
 totDATA=hGJ[1]->Integral();
 //hGJ[0]->Sumw2();
 if(totMC!=0)hGJ[0]->Scale(totDATA/totMC);
 hGJ[0]->SetLineColor(kRed);
 hGJ[0]->SetFillColor(kRed);
 hGJ[0]->Draw("E4");
  
 hGJ[1]->Draw("same E");

 totMC=hGJMCraw->Integral();
 hGJMCraw->Scale(totDATA/totMC);
 hGJMCraw->Draw("same");

 hGJ[0]->GetYaxis()->SetRangeUser(0,hGJ[0]->GetMaximum()*1.5);
 gPad->Update();

 c->cd(5);
 totMC=hGJ2[0]->Integral();
 totDATA=hGJ2[1]->Integral();
 //hGJ2[0]->Sumw2();
 if(totMC!=0)hGJ2[0]->Scale(totDATA/totMC);
 hGJ2[0]->SetLineColor(kRed);
  hGJ2[0]->SetFillColor(kRed);
 hGJ2[0]->Draw("E4");
 
 hGJ2[1]->Draw("same E");
 hGJ2[0]->GetYaxis()->SetRangeUser(0,hGJ2[0]->GetMaximum()*1.5);
 gPad->Update();

 c->cd(6);
 // hM[0]->SetLineColor(kRed);
//  hM[0]->Draw("");
//  hM[1]->Draw("same");
 if(totMC!=0)hThetaLab[0]->Scale(totDATA/totMC);
 hThetaLab[0]->SetLineColor(kRed);
 hThetaLab[0]->SetMarkerColor(kRed);
 hThetaLab[0]->Draw();
 //gPad->Update();
//gPad->SetLogy();

 hThetaLab[1]->Draw("same");
 gPad->SetLogy();
 gPad->Update();



//  totMC=hGJ3[0]->Integral();
//  totDATA=hGJ3[1]->Integral();
//  //hGJ3[0]->Sumw2();
//  hGJ3[0]->Scale(totDATA/totMC);
//  hGJ3[0]->SetLineColor(kRed);
//  //hGJ3[0]->SetFillColor(kRed);
//  hGJ3[0]->Draw();
 
//  hGJ3[1]->Draw("same E");
//  hGJ3[0]->GetYaxis()->SetRangeUser(0,hGJ3[0]->GetMaximum()*1.5);
//gPad->Update();


 // c->cd(7);
 //totMC=hTY[0]->Integral();
 //totDATA=hTY[1]->Integral();
 //hTY[0]->Sumw2();
 //if(totMC!=0)hTY[0]->Scale(totDATA/totMC);
 //hTY[0]->SetLineColor(kRed);
 //hTY[0]->SetFillColor(kRed);
 //hTY[0]->Draw("E4");
 

 //hTY[1]->Draw("same E");
 //hTY[0]->GetYaxis()->SetRangeUser(0,hTY[0]->GetMaximum()*1.5);
 //gPad->Update();

 //c->cd(8);
 //totMC=hTY2[0]->Integral();
 // totDATA=hTY2[1]->Integral();
//  //hTY2[0]->Sumw2();
//  if(totMC!=0)hTY2[0]->Scale(totDATA/totMC);
//  hTY2[0]->SetLineColor(kRed);
//  hTY2[0]->SetFillColor(kRed);
//  hTY2[0]->Draw("E4");
 
//  hTY2[1]->Draw("same E");
//  hTY2[0]->GetYaxis()->SetRangeUser(0,hTY2[0]->GetMaximum()*1.5);
//  gPad->Update();


// c->cd(9);
//  totMC=hTY3[0]->Integral();
//  totDATA=hTY3[1]->Integral();
//  //hTY3[0]->Sumw2();
//  if(totMC!=0)hTY3[0]->Scale(totDATA/totMC);
//  hTY3[0]->SetLineColor(kRed);
//  hTY3[0]->SetFillColor(kRed);
//  hTY3[0]->Draw("E4");
 
//  hTY3[1]->Draw("same E");
//  hTY3[0]->GetYaxis()->SetRangeUser(0,hTY3[0]->GetMaximum()*1.5);
//  gPad->Update();

 c->cd(1);
 TLatex* Label=new TLatex();
 Label->PaintLatex(2,hMIsobar[0]->GetMaximum()*0.8,0,0.1,massbin.Data());
 
 c->Update();



 TList* Hlist=gDirectory->GetList();
 Hlist->Remove(mctr);
 Hlist->Remove(datatr);
 //Hlist->Remove("hWeights");

 TFile* outfile=TFile::Open(outfilename,"UPDATE");
 TString psFileName=outfilename;
 psFileName.ReplaceAll(".root",massbin);
 psFileName.Append(".ps");
 

 if(totMC!=0){
 // get mass-integrated plots (or create them if not there)
 TH1D* hMIsobarMCGlob=(TH1D*)outfile->Get("hMIsobarMCGlob");
 if(hMIsobarMCGlob==NULL)hMIsobarMCGlob=new TH1D("hMIsobarMCGlob","4#pi Isobar Mass (MC)",nbninsm,0.2,3.8);
 hMIsobarMCGlob->Add(hMIsobar[0]);
 hMIsobarMCGlob->Write();
 TH1D* hMIsobarDataGlob=(TH1D*)outfile->Get("hMIsobarDataGlob");
 if(hMIsobarDataGlob==NULL)hMIsobarDataGlob=new TH1D("hMIsobarDataGlob","4#pi Isobar mass (DATA)",nbninsm,0.2,3.8);
 hMIsobarDataGlob->Add(hMIsobar[1]);
 hMIsobarDataGlob->Write();

TH1D* hMIsobar2MCGlob=(TH1D*)outfile->Get("hMIsobar2MCGlob");
 if(hMIsobar2MCGlob==NULL)hMIsobar2MCGlob=new TH1D("hMIsobar2MCGlob","3#pi Isobar mass (MC)",nbninsm,0.4,3.6);
 hMIsobar2MCGlob->Add(hMIsobar2[0]);
 hMIsobar2MCGlob->Write();

 TH1D* hMIsobar2DataGlob=(TH1D*)outfile->Get("hMIsobar2DataGlob");
 if(hMIsobar2DataGlob==NULL)hMIsobar2DataGlob=new TH1D("hMIsobar2DataGlob","3#pi Isobar mass (DATA)",nbninsm,0.4,3.6);
 hMIsobar2DataGlob->Add(hMIsobar2[1]);
 hMIsobar2DataGlob->Write();

TH1D* hMIsobar3MCGlob=(TH1D*)outfile->Get("hMIsobar3MCGlob");
 if(hMIsobar3MCGlob==NULL)hMIsobar3MCGlob=new TH1D("hMIsobar3MCGlob","2#pi Isobar mass (MC)",nbninsm,0.2,3.4);
 hMIsobar3MCGlob->Add(hMIsobar3[0]);
 hMIsobar3MCGlob->Write();

 TH1D* hMIsobar3DataGlob=(TH1D*)outfile->Get("hMIsobar3DataGlob");
 if(hMIsobar3DataGlob==NULL)hMIsobar3DataGlob=new TH1D("hMIsobar3DataGlob","2#pi Isobar mass (DATA)",nbninsm,0.2,3.4);
 hMIsobar3DataGlob->Add(hMIsobar3[1]);
 hMIsobar3DataGlob->Write();

 
 TH1D* hGJMCGlob=(TH1D*)outfile->Get("hGJMCGlob");
 if(hGJMCGlob==NULL)hGJMCGlob=new TH1D("hGJMCGlob","Gottfried-Jackson Theta (MC)",nbinsang,-1,1);
 hGJMCGlob->Add(hGJ[0]);
 hGJMCGlob->Write();

TH1D* hGJDataGlob=(TH1D*)outfile->Get("hGJDataGlob");
 if(hGJDataGlob==NULL)hGJDataGlob=new TH1D("hGJDataGlob","Gottfried-Jackson Theta (Data)",nbinsang,-1,1);
 hGJDataGlob->Add(hGJ[1]);
 hGJDataGlob->Write();
 
TH1D* hGJ2MCGlob=(TH1D*)outfile->Get("hGJ2MCGlob");
 if(hGJ2MCGlob==NULL)hGJ2MCGlob=new TH1D("hGJ2MCGlob","Gottfried-Jackson Theta (MC)",nbinsang,-1,1);
 hGJ2MCGlob->Add(hGJ2[0]);
 hGJ2MCGlob->Write();

TH1D* hGJ2DataGlob=(TH1D*)outfile->Get("hGJ2DataGlob");
 if(hGJ2DataGlob==NULL)hGJ2DataGlob=new TH1D("hGJ2DataGlob","Gottfried-Jackson Theta (Data)",nbinsang,-1,1);
 hGJ2DataGlob->Add(hGJ2[1]);
 hGJ2DataGlob->Write();
 

 c->cd(7);
 hMIsobarMCGlob->SetLineColor(kRed);
 hMIsobarMCGlob->SetFillColor(kRed);
 hMIsobarMCGlob->Draw("E4");
 hMIsobarDataGlob->Draw("E SAME");

 c->cd(8);
 hMIsobar2MCGlob->SetLineColor(kRed);
 hMIsobar2MCGlob->SetFillColor(kRed);
 hMIsobar2MCGlob->Draw("E4");
 hMIsobar2DataGlob->Draw("E SAME");


 c->cd(9);
 hGJMCGlob->SetLineColor(kRed);
 hGJMCGlob->SetFillColor(kRed);
 hGJMCGlob->Draw("E4");
 hGJDataGlob->Draw("E SAME");
 }

 c->Write();
 TCanvas      dummyCanv("dummy", "dummy");
 c->Print((psFileName));

 int nobj=Hlist->GetEntries();
 std::cout<<"Found "<<nobj<<" Objects in HList"<<std::endl;
 Hlist->Print();
 Hlist->Write();
 outfile->Close();
 
 gjt->Update();

 gROOT->cd();

 
}
