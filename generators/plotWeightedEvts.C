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
#include "NParticleEvent.h"
// #include "plotWeightedEvts.h" maybe you are missing this? not needed since declaration is not needed, too


using namespace std; 


// subroutine to build a graph with errors out of a list of histograms
// the first histo is taken as the value
// the errors are computed from the spread of the values over all histograms

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
    sigma/=(n-1);
    g->SetPointError(ib,dx,sqrt(sigma));
  }
  
  return g;

}

// this moved to the header, this is a temporary solution
// if you need the default values, try to include the header
// I have no clue if it still compiles with ALIC compiler with root
//void plotWeightedEvts(TTree* mctr, TTree* datatr, TString outfilename="kineplots.root", TString mass="000");

void plotWeightedEvts(TTree* mctr, TTree* datatr, TString outfilename, TString mass, unsigned int nsamples=1){

gROOT->SetStyle("Plain");
 Int_t NCont=100;
 
  const Int_t NRGBs = 5;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
 
  gStyle->SetNumberContours(NCont);

 TString massbin;
 if(mass!="000"){
   massbin=("_m")+mass;massbin.ReplaceAll(" ","");
 }

 int nbninsm=144; // 144
 int nbinsang=40; // 80


 gStyle->SetOptStat(0);
 gStyle->SetOptTitle(0);
 gStyle->SetStripDecimals(1);
 TGaxis::SetMaxDigits(4);
 
 Int_t font=132;
 
 gStyle->SetTextFont(font);
 gStyle->SetLabelFont(font,"xyz");
 gStyle->SetTitleFont(font,"xyz");
 
 vector<TH1D*> hM;
 TH1D* hMMC=new TH1D("hMMC"+massbin,"Mass (MC)",60,1.0,3.0);
 hM.push_back(hMMC);
 TH1D* hMData=new TH1D("hMData"+massbin,"Mass (DATA)",60,1.0,3.0);
 hM.push_back(hMData);


 vector<TH1D*> hMIsobar; // 4pi
 vector<TH1D*> hMIsobar2; // 3pi
 vector<TH1D*> hMIsobar3; // 2pi
 vector<TH1D*> hGJ; // 4pi/1pi
 vector<TH1D*> hTY; // 4pi/1pi
 
 vector<TH1D*> hGJ2; // 3pi/2pi
 vector<TH1D*> hHe22; // 2pi2pi helicity frame
 vector<TH1D*> hHe22Phi;
 vector<TH1D*> hHe31; // 3pi1pi helicity frame
 vector<TH1D*> hHe31Phi;
 vector<TH1D*> hHe21; // 2pi1pi helicity frame
 vector<TH1D*> hHe21Phi;


  for(unsigned int isamples=0;isamples<nsamples;++isamples){
    TString blub;blub+=massbin;blub+="_";blub+=isamples;
    TH1D* hMIsobarMC=new TH1D("hMIsobarMC"+blub,"Isobar Mass (MC) for m_{5#pi}="+mass,nbninsm,0.2,3);
    hMIsobar.push_back(hMIsobarMC);
    hMIsobarMC->Sumw2();
    
    TH1D* hMIsobar2MC=new TH1D("hMIsobar2MC"+blub,"Isobar Mass (MC)",nbninsm,0.4,2.7);
    hMIsobar2.push_back(hMIsobar2MC);
    hMIsobar2MC->Sumw2();
    
    TH1D* hMIsobar3MC=new TH1D("hMIsobar3MC"+blub,"Isobar Mass (MC)",nbninsm,0.2,2.2);
    hMIsobar3.push_back(hMIsobar3MC);
    hMIsobar3MC->Sumw2();
    
    TString name="hGJMC";name+=+massbin;name+="_";name+=isamples;
    TH1D* hGJMC=new TH1D(name,"Cos Gottfried-Jackson Theta (MC)",nbinsang,-1,1);
    hGJ.push_back(hGJMC);
    hGJMC->Sumw2();
 
    TH1D* hTYMC=new TH1D("hTYMC"+blub,"Treiman-Yang Phi (MC)",nbinsang,-TMath::Pi(),TMath::Pi());
    hTY.push_back(hTYMC);
    hTYMC->Sumw2();


    TH1D* hGJ2MC=new TH1D("hGJ2MC"+blub,"Cos Gottfried-Jackson Theta (MC)",nbinsang,-1,1);
    hGJ2.push_back(hGJ2MC);
    hGJ2MC->Sumw2();

    TH1D* hHe22MC=new TH1D("hHe22ThMC"+blub,"Helicity Theta 2pi2pi (MC)",nbinsang,-1,1);
    hHe22.push_back(hHe22MC);
    hHe22MC->Sumw2();

    TH1D* hHe22PhiMC=new TH1D("hHe22PhiMC"+blub,"Helicity Phi 2pi2pi (MC)",nbinsang,-TMath::Pi(),TMath::Pi());
    hHe22Phi.push_back(hHe22PhiMC);
    hHe22PhiMC->Sumw2();

 TH1D* hHe31MC=new TH1D("hHe31ThMC"+blub,"Helicity Theta 3pi1pi (MC)",nbinsang,-1,1);
    hHe31.push_back(hHe31MC);
    hHe31MC->Sumw2();

    TH1D* hHe31PhiMC=new TH1D("hHe31PhiMC"+blub,"Helicity Phi 3pi1pi (MC)",nbinsang,-TMath::Pi(),TMath::Pi());
    hHe31Phi.push_back(hHe31PhiMC);
    hHe31PhiMC->Sumw2();


TH1D* hHe21MC=new TH1D("hHe21ThMC"+blub,"Helicity Theta 2pi1pi (MC)",nbinsang,-1,1);
    hHe21.push_back(hHe21MC);
    hHe21MC->Sumw2();

    TH1D* hHe21PhiMC=new TH1D("hHe21PhiMC"+blub,"Helicity Phi 2pi1pi (MC)",nbinsang,-TMath::Pi(),TMath::Pi());
    hHe21Phi.push_back(hHe21PhiMC);
    hHe21PhiMC->Sumw2();
 } // end loop over model samples


 TH1D* hMIsobarData=new TH1D("hMIsobarData"+massbin,"Isobar mass (DATA)",nbninsm,0.2,3);
 hMIsobar.push_back(hMIsobarData);

 TH1D* hDiffMIsobar; //=new TH1D("hDiffMIsobarData","Isobar mass (DATA) Diff",40,0.2,2.2);
TH1D* hDiffMIsobar2;
TH1D* hDiffMIsobar3;


 TH1D* hMIsobar2Data=new TH1D("hMIsobar2Data"+massbin,"Isobar mass (DATA)",nbninsm,0.4,2.7);
 hMIsobar2.push_back(hMIsobar2Data);


 TH1D* hMIsobar3Data=new TH1D("hMIsobar3Data"+massbin,"Isobar mass (DATA)",nbninsm,0.2,2.3);
 hMIsobar3.push_back(hMIsobar3Data);


 TH1D* hGJData=new TH1D("hGJData"+massbin,"Cos Gottfried-Jackson Theta (DATA)",nbinsang,-1,1); hGJ.push_back(hGJData);

 TH1D* hGJMCraw=new TH1D("hGJMCraw"+massbin,"Cos Gottfried-Jackson Theta (MC raw)",nbinsang,-1,1);
 
 TH1D* hHe22Data=new TH1D("hHe22ThData"+massbin,"Helicity Theta 2pi2pi (Data)",nbinsang,-1,1);
    hHe22.push_back(hHe22Data);


    TH1D* hHe22PhiData=new TH1D("hHe22PhiData"+massbin,"Helicity Phi 2pi2pi (Data)",nbinsang,-TMath::Pi(),TMath::Pi());
    hHe22Phi.push_back(hHe22PhiData);

 TH1D* hHe31Data=new TH1D("hHe31ThData"+massbin,"Helicity Theta 3pi1pi (Data)",nbinsang,-1,1);
    hHe31.push_back(hHe31Data);


    TH1D* hHe31PhiData=new TH1D("hHe31PhiData"+massbin,"Helicity Phi 3pi1pi (Data)",nbinsang,-TMath::Pi(),TMath::Pi());
    hHe31Phi.push_back(hHe31PhiData);

 TH1D* hHe21Data=new TH1D("hHe21ThData"+massbin,"Helicity Theta 2pi1pi (Data)",nbinsang,-1,1);
    hHe21.push_back(hHe21Data);


    TH1D* hHe21PhiData=new TH1D("hHe21PhiData"+massbin,"Helicity Phi 2pi1pi (Data)",nbinsang,-TMath::Pi(),TMath::Pi());
    hHe21Phi.push_back(hHe21PhiData);


 vector<TH2D*> hM1vsGJ;
 TH2D* hM1vsGJMC=new TH2D("hM1vsGJMC"+massbin,"M1 Cos Gottfried-Jackson Theta (MC)",nbinsang/4,-1,1,nbninsm/2,0.2,3.8);
 hM1vsGJ.push_back(hM1vsGJMC);
 TH2D* hM1vsGJData=new TH2D("hM1vsGJData"+massbin,"M1 Cos Gottfried-Jackson Theta (Data)",nbinsang/4,-1,1,nbninsm/2,0.2,3.8);
 hM1vsGJ.push_back(hM1vsGJData);


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





 TH1D* hGJ2Data=new TH1D("hGJ2Data"+massbin,"Cos Gottfried-Jackson Theta (DATA)",nbinsang,-1,1); hGJ2.push_back(hGJ2Data);


 //  vector<TH1D*> hGJ3;
 //TH1D* hGJ3MC=new TH1D("hGJ3MC"+massbin,"Cos Gottfried-Jackson Theta (MC)",nbinsang,-1,1);
 //hGJ3.push_back(hGJ3MC);
 //TH1D* hGJ3Data=new TH1D("hGJ3Data"+massbin,"Cos Gottfried-Jackson Theta (DATA)",nbinsang,-1,1); hGJ3.push_back(hGJ3Data);
 //hGJ3[0]->Sumw2();
 

 TH1D* hTYData=new TH1D("hTYData"+massbin,"Treiman-Yang Phi (DATA)",nbinsang,-TMath::Pi(),TMath::Pi());
  hTY.push_back(hTYData);
 

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
   std::vector<double> weights(nsamples);
   if(itree==0){
     for(unsigned int isamples=0;isamples<nsamples;++isamples){
       TString wname("W");wname+=isamples;
       tr->SetBranchAddress(wname.Data(),&weights[isamples]);
     }
   }

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

     double tprime=event.tprime();

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
	 
	 if(itree==0){ // mc-tree
	   hGJMCraw->Fill(state.p().CosTheta());
	   for(unsigned int isamples=0;isamples<nsamples;++isamples){
	     hGJ[isamples]->Fill(state.p().CosTheta(),weights[isamples]);
	     hTY[isamples]->Fill(state.p().Phi(),weights[isamples]);
	     hMIsobar[isamples]->Fill(state.p().M(),weights[isamples]); 
	   } // end loop over sample models
	 }
	 else { // data tree
	   hGJ[nsamples+itree-1]->Fill(state.p().CosTheta(),1);
	   hMIsobar[nsamples+itree-1]->Fill(state.p().M(),1); 
	   //hGJt[itree]->Fill(state.p().CosTheta(),tprime,1);
	   hTY[nsamples+itree-1]->Fill(state.p().Phi(),1);
	   //hM1vsGJ[itree]->Fill(state.p().CosTheta(),state.p().M(),weight); 
	 }
       }
       else if(state.n()==npart-2 && state.q()==-1){
	 if(itree==0){
	   for(unsigned int isamples=0;isamples<nsamples;++isamples){
	     hGJ2[isamples]->Fill(state.p().CosTheta(),weights[isamples]);
	     hMIsobar2[isamples]->Fill(state.p().M(),weights[isamples]);
	   }
	 }
	 else {
	   hGJ2[itree+nsamples-1]->Fill(state.p().CosTheta(),1);
	   hMIsobar2[itree+nsamples-1]->Fill(state.p().M(),1);
	 }
       }
       else if(state.n()==npart-3 && state.q()==0){
	  if(itree==0){
	    for(unsigned int isamples=0;isamples<nsamples;++isamples){
	      //hGJ3[isamples]->Fill(state.p().CosTheta(),weights[isamples]);
	      hMIsobar3[isamples]->Fill(state.p().M(),weights[isamples]);
	    }
	  }
	  else {
	    //hGJ3[isamples]->Fill(state.p().CosTheta(),1);
	    hMIsobar3[itree+nsamples-1]->Fill(state.p().M(),1);
	  }
       }
     } // end loop over states 1
     
     
     
     // loop again and fill helicity frame stuff
     // for this we get rid of GJ transform
     event.refresh(); // we do not need to build states again!
     TLorentzVector pX=event.p();
     // and loooooop
     for(unsigned int is=0;is<nstates;++is){
       
       const NParticleState& state=event.getState(is);
       if(state.n()==4 && state.q()==0){
	 // select sub-systems
	 // (2pi)(2pi) or (3pi)(pi)
	 for(unsigned int js=0;js<nstates;++js){
	   const NParticleState& substate=event.getState(js);
	   // (2pi)(2pi)
	   if(substate.n()==2 && substate.q()==0){
	     //cerr << "FOUND 3PISTATE" << endl;
	     if(substate.isSubstate(&state)){
	       TLorentzVector phe=event.getHelicityFrame(pX,state.p(),substate.p());
	       if(itree==0){ // mc-tree
		 for(unsigned int isamples=0;isamples<nsamples;++isamples){
		   hHe22[isamples]->Fill(phe.CosTheta(),weights[isamples]);
		   hHe22Phi[isamples]->Fill(phe.Phi(),weights[isamples]);
		 } // end loop over sample models
	       }
	       else { // data tree
		 hHe22[nsamples+itree-1]->Fill(phe.CosTheta(),1);
		 hHe22Phi[nsamples+itree-1]->Fill(phe.Phi(),1);
	       }
	       //break; // make sure we do not double-count!
	     }// end (2pi)(2pi) case
	   }
	   else  if(substate.n()==3 && substate.q()==-1){
	     //cerr << "FOUND 3PISTATE" << endl;
	     if(substate.isSubstate(&state)){
	       TLorentzVector phe=event.getHelicityFrame(pX,state.p(),substate.p());
	       if(itree==0){ // mc-tree
		 for(unsigned int isamples=0;isamples<nsamples;++isamples){
		   hHe31[isamples]->Fill(phe.CosTheta(),weights[isamples]);
		   hHe31Phi[isamples]->Fill(phe.Phi(),weights[isamples]);
		 } // end loop over sample models
	       }
	       else { // data tree
		 hHe31[nsamples+itree-1]->Fill(phe.CosTheta(),1);
		 hHe31Phi[nsamples+itree-1]->Fill(phe.Phi(),1);
	       }
	       //break; // make sure we do not double-count!
	     } // end (3pi)(1pi) case
	   }
	 } // end loop over substates js
       } // end 4pi system
       // 3pi system
       if(state.n()==3 && state.q()==-1){
	 // select sub-systems
	 // (2pi)(1pi)
	 for(unsigned int js=0;js<nstates;++js){
	   const NParticleState& substate=event.getState(js);
	   // (2pi)(1pi)
	   if(substate.n()==2 && substate.q()==0){
	     //cerr << "FOUND 3PISTATE" << endl;
	     if(substate.isSubstate(&state)){
	       TLorentzVector phe=event.getHelicityFrame(pX,state.p(),substate.p());
	       if(itree==0){ // mc-tree
		 for(unsigned int isamples=0;isamples<nsamples;++isamples){
		   hHe21[isamples]->Fill(phe.CosTheta(),weights[isamples]);
		   hHe21Phi[isamples]->Fill(phe.Phi(),weights[isamples]);
		 } // end loop over sample models
	       }
	       else { // data tree
		 hHe21[nsamples+itree-1]->Fill(phe.CosTheta(),1);
		 hHe21Phi[nsamples+itree-1]->Fill(phe.Phi(),1);
	       }
	       //break; // make sure we do not double-count!
	     } // end (2pi)(1pi) case
	   }
	 } // end loop over substates ij
       } // end 3pi case
     } // end loop over states is 2
     
     
   }// end loop over events
   if(itree==0)avweight/=(double)nevt;
   cout << "Maxweight=" << maxweight << endl; 
   cout << "Average weight=" << avweight << endl; 
 } // loop over trees




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
 c->Divide(3,4);

 // plot MC predictions
 double totMC=hMIsobar[0]->Integral();
 double totDATA=hMIsobar[nsamples]->Integral();

 for(unsigned int isamples=0;isamples<nsamples;++isamples){ 
   const char* opt= (isamples==0) ? "" : "same"; 
   c->cd(1);
  
   double totMC=hMIsobar[isamples]->Integral();
   double totDATA=hMIsobar[nsamples]->Integral();
   if(totMC!=0)hMIsobar[isamples]->Scale(totDATA/totMC);
   hMIsobar[isamples]->SetLineColor(kRed);
   ///hMIsobar[isamples]->SetFillColor(kRed);
   hMIsobar[isamples]->Draw(opt);
 
   c->cd(2);

   totMC=hMIsobar2[isamples]->Integral();
   totDATA=hMIsobar2[nsamples]->Integral();
   if(totMC!=0)hMIsobar2[isamples]->Scale(totDATA/totMC);
   hMIsobar2[isamples]->SetLineColor(kRed);
   //hMIsobar2[isamples]->SetFillColor(kRed);
   hMIsobar2[isamples]->Draw(opt);
 

 c->cd(3);
 totMC=hMIsobar3[isamples]->Integral();
 totDATA=hMIsobar3[nsamples]->Integral();
 if(totMC!=0)hMIsobar3[isamples]->Scale(totDATA/totMC);
 hMIsobar3[isamples]->SetLineColor(kRed);
 //hMIsobar3[isamples]->SetFillColor(kRed);
 hMIsobar3[isamples]->Draw(opt);

 c->cd(4);
 totDATA=hGJ[nsamples]->Integral();
 totMC=hGJ[isamples]->Integral();
 if(totMC!=0)hGJ[isamples]->Scale(totDATA/totMC);
 hGJ[isamples]->SetLineColor(kRed);
 //hGJ[isamples]->SetFillColor(kRed);
hGJ[isamples]->Draw(opt);
  

c->cd(5);
 totDATA=hTY[nsamples]->Integral();
 totMC=hTY[isamples]->Integral();
 if(totMC!=0)hTY[isamples]->Scale(totDATA/totMC);
 hTY[isamples]->SetLineColor(kRed);
 //hGJ[isamples]->SetFillColor(kRed);
hTY[isamples]->Draw(opt);

 c->cd(6);
 totMC=hGJ2[isamples]->Integral();
 totDATA=hGJ2[nsamples]->Integral();
 //hGJ2[isamples]->Sumw2();
 if(totMC!=0)hGJ2[isamples]->Scale(totDATA/totMC);
 hGJ2[isamples]->SetLineColor(kRed);
  hGJ2[isamples]->SetFillColor(kRed);
 hGJ2[isamples]->Draw(opt);


 c->cd(7);
  totMC=hHe22[isamples]->Integral();
 totDATA=hHe22[nsamples]->Integral();
 //hGJ2[isamples]->Sumw2();
 if(totMC!=0)hHe22[isamples]->Scale(totDATA/totMC);
 hHe22[isamples]->SetLineColor(kRed);
 // hHe[isamples]->SetFillColor(kRed);
 hHe22[isamples]->Draw(opt);

 c->cd(10);
  totMC=hHe22Phi[isamples]->Integral();
 totDATA=hHe22Phi[nsamples]->Integral();
 //hGJ2[isamples]->Sumw2();
 if(totMC!=0)hHe22Phi[isamples]->Scale(totDATA/totMC);
 hHe22Phi[isamples]->SetLineColor(kRed);
 // hHePhi[isamples]->SetFillColor(kRed);
 hHe22Phi[isamples]->Draw(opt);


 c->cd(8);
  totMC=hHe31[isamples]->Integral();
 totDATA=hHe31[nsamples]->Integral();
 //hGJ2[isamples]->Sumw2();
 if(totMC!=0)hHe31[isamples]->Scale(totDATA/totMC);
 hHe31[isamples]->SetLineColor(kRed);
 // hHe[isamples]->SetFillColor(kRed);
 hHe31[isamples]->Draw(opt);

 c->cd(11);
  totMC=hHe31Phi[isamples]->Integral();
 totDATA=hHe31Phi[nsamples]->Integral();
 //hGJ2[isamples]->Sumw2();
 if(totMC!=0)hHe31Phi[isamples]->Scale(totDATA/totMC);
 hHe31Phi[isamples]->SetLineColor(kRed);
 // hHePhi[isamples]->SetFillColor(kRed);
 hHe31Phi[isamples]->Draw(opt);

 c->cd(9);
  totMC=hHe21[isamples]->Integral();
 totDATA=hHe21[nsamples]->Integral();
 //hGJ2[isamples]->Sumw2();
 if(totMC!=0)hHe21[isamples]->Scale(totDATA/totMC);
 hHe21[isamples]->SetLineColor(kRed);
 // hHe[isamples]->SetFillColor(kRed);
 hHe21[isamples]->Draw(opt);

 c->cd(12);
  totMC=hHe21Phi[isamples]->Integral();
 totDATA=hHe21Phi[nsamples]->Integral();
 //hGJ2[isamples]->Sumw2();
 if(totMC!=0)hHe21Phi[isamples]->Scale(totDATA/totMC);
 hHe21Phi[isamples]->SetLineColor(kRed);
 // hHePhi[isamples]->SetFillColor(kRed);
 hHe21Phi[isamples]->Draw(opt);

 }// end loopover samples

 c->cd(6);
 
 hGJ2[nsamples]->Draw("same E");
 hGJ2[0]->GetYaxis()->SetRangeUser(0,hGJ2[0]->GetMaximum()*1.5);
 gPad->Update();

 c->cd(7);
 hHe22[nsamples]->Draw("same E");
 hHe22[0]->GetYaxis()->SetRangeUser(0,hHe22[0]->GetMaximum()*1.5);
 gPad->Update();

c->cd(10);
 hHe22Phi[nsamples]->Draw("same E");
 hHe22Phi[0]->GetYaxis()->SetRangeUser(0,hHe22Phi[0]->GetMaximum()*1.5);
 gPad->Update();

 c->cd(8);
 hHe31[nsamples]->Draw("same E");
 hHe31[0]->GetYaxis()->SetRangeUser(0,hHe31[0]->GetMaximum()*1.5);
 gPad->Update();

c->cd(11);
 hHe31Phi[nsamples]->Draw("same E");
 hHe31Phi[0]->GetYaxis()->SetRangeUser(0,hHe31Phi[0]->GetMaximum()*1.5);
 gPad->Update();

 c->cd(9);
 hHe21[nsamples]->Draw("same E");
 hHe21[0]->GetYaxis()->SetRangeUser(0,hHe21[0]->GetMaximum()*1.5);
 gPad->Update();

c->cd(12);
 hHe21Phi[nsamples]->Draw("same E");
 hHe21Phi[0]->GetYaxis()->SetRangeUser(0,hHe21Phi[0]->GetMaximum()*1.5);
 gPad->Update();

 // hM[isamples]->SetLineColor(kRed);
//  hM[isamples]->Draw("");
//  hM[1]->Draw("same");
 //if(totMC!=0)hThetaLab[itree+]->Scale(totDATA/totMC);
 //hThetaLab[0]->SetLineColor(kRed);
 //hThetaLab[isamples]->SetMarkerColor(kRed);
 //hThetaLab[isamples]->Draw();
 //gPad->Update();
//gPad->SetLogy();

 //hThetaLab[1]->Draw("same");
 //gPad->SetLogy();
 //gPad->Update();



//  totMC=hGJ3[isamples]->Integral();
//  totDATA=hGJ3[1]->Integral();
//  //hGJ3[isamples]->Sumw2();
//  hGJ3[isamples]->Scale(totDATA/totMC);
//  hGJ3[isamples]->SetLineColor(kRed);
//  //hGJ3[isamples]->SetFillColor(kRed);
//  hGJ3[isamples]->Draw();
 
//  hGJ3[1]->Draw("same E");
//  hGJ3[isamples]->GetYaxis()->SetRangeUser(0,hGJ3[isamples]->GetMaximum()*1.5);
//gPad->Update();


 c->cd(7);
 // totMC=hTY[0]->Integral();
//  totDATA=hTY[1]->Integral();
//  hTY[0]->Sumw2();
//  if(totMC!=0)hTY[0]->Scale(totDATA/totMC);
//  hTY[0]->SetLineColor(kRed);
//  hTY[0]->SetFillColor(kRed);
//  hTY[0]->Draw("E");
 

//  hTY[1]->Draw("same E");
//  hTY[0]->GetYaxis()->SetRangeUser(0,hTY[0]->GetMaximum()*1.5);
//  gPad->Update();

//  c->cd(8);
//  totMC=hTY2[0]->Integral();
//  totDATA=hTY2[1]->Integral();
//  //hTY2[0]->Sumw2();
//  if(totMC!=0)hTY2[0]->Scale(totDATA/totMC);
//  hTY2[0]->SetLineColor(kRed);
//  hTY2[0]->SetFillColor(kRed);
//  hTY2[0]->Draw("E");
 
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
//  hTY3[0]->Draw("E");
 
//  hTY3[1]->Draw("same E");
//  hTY3[0]->GetYaxis()->SetRangeUser(0,hTY3[0]->GetMaximum()*1.5);
//  gPad->Update();

 c->cd(1);
 
 hDiffMIsobar=(TH1D*)hMIsobar[0]->Clone("hMIsobarDiff"+massbin);
 hDiffMIsobar->Add(hMIsobar[nsamples],-1.);
 hMIsobar[nsamples]->Draw("same E");
 hDiffMIsobar->SetLineColor(kOrange-3);
 hDiffMIsobar->Draw("same");

 hMIsobar[0]->GetYaxis()->SetRangeUser(hDiffMIsobar->GetMinimum()*1.5,hMIsobar[0]->GetMaximum()*1.2);

 TString lab=massbin;lab.ReplaceAll("_m","m=");
 TLatex* Label=new TLatex(1.5,hMIsobar[0]->GetMaximum()*0.8,lab.Data());
 //Label->PaintLatex(2,hMIsobar[0]->GetMaximum()*0.8,0,0.1,massbin.Data());
 Label->Draw();
 gPad->Update();


c->cd(2);
 hMIsobar2[nsamples]->Draw("same E");
 hDiffMIsobar2=(TH1D*)hMIsobar2[0]->Clone("hMIsobar2Diff"+massbin);
 hDiffMIsobar2->Add(hMIsobar2[nsamples],-1.);
 hDiffMIsobar2->SetLineColor(kOrange-3);
 hDiffMIsobar2->Draw("same E");


hMIsobar2[0]->GetYaxis()->SetRangeUser(hDiffMIsobar2->GetMinimum()*1.5,hMIsobar2[0]->GetMaximum()*1.2);
 gPad->Update();

c->cd(3);
hMIsobar3[nsamples]->Draw("same E");
 hDiffMIsobar3=(TH1D*)hMIsobar3[0]->Clone("hMIsobar3Diff"+massbin);
 hDiffMIsobar3->Add(hMIsobar3[nsamples],-1.);
 hDiffMIsobar3->SetLineColor(kOrange-3);
 hDiffMIsobar3->Draw("same");

 hMIsobar3[0]->GetYaxis()->SetRangeUser(hDiffMIsobar3->GetMinimum()*1.5,hMIsobar3[0]->GetMaximum()*1.2);
 gPad->Update();

 c->cd(4);
 TGraphErrors* gGJ=buildGraph(hGJ,nsamples);
 gGJ->SetLineColor(kMagenta);
 gGJ->SetFillColor(kMagenta);
  gGJ->Draw("same 2");

 hGJ[nsamples]->Draw("same E");

 hGJ[0]->GetYaxis()->SetRangeUser(0,hGJ[0]->GetMaximum()*1.5);
 gPad->Update();

 c->cd(5);
 hTY[nsamples]->Draw("same E");
 hTY[0]->GetYaxis()->SetRangeUser(0,hTY[0]->GetMaximum()*1.5);
 gPad->Update();

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
 hMIsobarDataGlob->Add(hMIsobar[nsamples]);
 hMIsobarDataGlob->Write();

TH1D* hMIsobar2MCGlob=(TH1D*)outfile->Get("hMIsobar2MCGlob");
 if(hMIsobar2MCGlob==NULL)hMIsobar2MCGlob=new TH1D("hMIsobar2MCGlob","3#pi Isobar mass (MC)",nbninsm,0.4,3.6);
 hMIsobar2MCGlob->Add(hMIsobar2[0]);
 hMIsobar2MCGlob->Write();

 TH1D* hMIsobar2DataGlob=(TH1D*)outfile->Get("hMIsobar2DataGlob");
 if(hMIsobar2DataGlob==NULL)hMIsobar2DataGlob=new TH1D("hMIsobar2DataGlob","3#pi Isobar mass (DATA)",nbninsm,0.4,3.6);
 hMIsobar2DataGlob->Add(hMIsobar2[nsamples]);
 hMIsobar2DataGlob->Write();

TH1D* hMIsobar3MCGlob=(TH1D*)outfile->Get("hMIsobar3MCGlob");
 if(hMIsobar3MCGlob==NULL)hMIsobar3MCGlob=new TH1D("hMIsobar3MCGlob","2#pi Isobar mass (MC)",nbninsm,0.2,3.4);
 hMIsobar3MCGlob->Add(hMIsobar3[0]);
 hMIsobar3MCGlob->Write();

 TH1D* hMIsobar3DataGlob=(TH1D*)outfile->Get("hMIsobar3DataGlob");
 if(hMIsobar3DataGlob==NULL)hMIsobar3DataGlob=new TH1D("hMIsobar3DataGlob","2#pi Isobar mass (DATA)",nbninsm,0.2,3.4);
 hMIsobar3DataGlob->Add(hMIsobar3[nsamples]);
 hMIsobar3DataGlob->Write();

 
 TH1D* hGJMCGlob=(TH1D*)outfile->Get("hGJMCGlob");
 if(hGJMCGlob==NULL)hGJMCGlob=new TH1D("hGJMCGlob","Gottfried-Jackson Theta (MC)",nbinsang,-1,1);
 hGJMCGlob->Add(hGJ[0]);
 hGJMCGlob->Write();

TH1D* hGJDataGlob=(TH1D*)outfile->Get("hGJDataGlob");
 if(hGJDataGlob==NULL)hGJDataGlob=new TH1D("hGJDataGlob","Gottfried-Jackson Theta (Data)",nbinsang,-1,1);
 hGJDataGlob->Add(hGJ[nsamples]);
 hGJDataGlob->Write();
 
TH1D* hGJ2MCGlob=(TH1D*)outfile->Get("hGJ2MCGlob");
 if(hGJ2MCGlob==NULL)hGJ2MCGlob=new TH1D("hGJ2MCGlob","Gottfried-Jackson Theta (MC)",nbinsang,-1,1);
 hGJ2MCGlob->Add(hGJ2[0]);
 hGJ2MCGlob->Write();

TH1D* hGJ2DataGlob=(TH1D*)outfile->Get("hGJ2DataGlob");
 if(hGJ2DataGlob==NULL)hGJ2DataGlob=new TH1D("hGJ2DataGlob","Gottfried-Jackson Theta (Data)",nbinsang,-1,1);
 hGJ2DataGlob->Add(hGJ2[nsamples]);
 hGJ2DataGlob->Write();
 

 // c->cd(13);
 // hMIsobarMCGlob->SetLineColor(kRed);
 // hMIsobarMCGlob->SetFillColor(kRed);
 // hMIsobarMCGlob->Draw("E");
 // hMIsobarDataGlob->Draw("E SAME");

 // c->cd(14);
 // hMIsobar2MCGlob->SetLineColor(kRed);
 // hMIsobar2MCGlob->SetFillColor(kRed);
 // hMIsobar2MCGlob->Draw("E");
 // hMIsobar2DataGlob->Draw("E SAME");


 // c->cd(15);
 // hGJMCGlob->SetLineColor(kRed);
 // hGJMCGlob->SetFillColor(kRed);
 // hGJMCGlob->GetYaxis()->SetRangeUser(0,hGJMCGlob->GetMaximum()*1.5);
 // hGJMCGlob->Draw("E");
 // hGJDataGlob->Draw("E SAME");

 // c->cd(10);
 // hM1vsGJMC->Scale(totDATA/totMC);
 // hM1vsGJMC->Draw("COLZ");

 // c->cd(11);
 // hM1vsGJData->Draw("COLZ");

 // c->cd(12);
 // TH2D* hDiff=new TH2D(*hM1vsGJ[0]);
 // hDiff->Divide(hM1vsGJMC,hM1vsGJData,1,1);
 // hDiff->Draw("COLZ");
 // hDiff->GetZaxis()->SetRangeUser(0.7,1.3);

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

//If I'm not wrong this script should compile for root even with a main in it
// that will be simply ignored
int main(int argc, char** argv){

	if (argc != 5){
		cout << " not enough arguments! " << endl;
		return 1;
	}

	TString evtfile(argv[1]);
	TString mcfile(argv[2]);
	TString outfile(argv[3]);
	TString mass(argv[4]);

	TFile* file1=TFile::Open(evtfile,"READ");
	TFile* file2=TFile::Open(mcfile,"READ");
	TTree* data=(TTree*)file1->Get("events");
	TTree* mc=(TTree*)file2->Get("pwevents");

	plotWeightedEvts(mc,data,outfile,mass);

	file1->Close();
	file2->Close();
}
