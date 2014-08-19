#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TLatex.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TLorentzRotation.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TRandom3.h"
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


void plotEvts(TTree* datatr, TString outfilename="kineplots.root", TString mass="000"){

gROOT->SetStyle("Plain");

 TString massbin;
 if(mass!="000"){
   massbin=("_m")+mass;massbin.ReplaceAll(" ","");
 }



 TTree* tr= datatr;
   double impweight=1;
   TClonesArray* p=new TClonesArray("TLorentzVector");
   TLorentzVector* beam=NULL;
   int qbeam;
   std::vector<int>* q=NULL;
   //if (itree==0)tr->SetBranchAddress("weight",&weight);
   //if (itree==0)tr->SetBranchAddress("impweight",&impweight);
   tr->SetBranchAddress("p",&p);
   tr->SetBranchAddress("beam",&beam);
   tr->SetBranchAddress("qbeam",&qbeam);
   tr->SetBranchAddress("q",&q);

   TVector3 vertex;
   NParticleEvent event(p,q,beam,&qbeam,&vertex);

   TF1* acc=new TF1("acc","1-exp([0]*(1-x))",0,2000); // artificial acceptance function

   TH1D* hMass=new TH1D("hmass","Mass",1000,0,20);

   TH1D* hGJ=new TH1D("hGJ","Cos GJ-Theta",20,-1,1);
   TH1D* hGJacc=new TH1D("hGJacc","acc Cos GJ-Theta",20,-1,1);

   TH1D* hExcl=new TH1D("hExcl","Charged Exclusivity",100,0,200);
   TH1D* hT=new TH1D("hT","Charged t",1000,0,2);

   TH1D* hMom=new TH1D("hMOM","events",50,0,50);
   TH1D* hMombad=new TH1D("hMOMbad","events",50,0,50);

   TH2D* hPbad=new TH2D("hPbad","bad events",40,0,0.05,50,0,160);
   TH2D* hP=new TH2D("hP","events",40,0,0.05,50,0,160);
   TH2D* hPXYbad=new TH2D("hPXYbad","bad events",50,-0.05,0.05,50,-0.05,0.05);
   TH2D* hPXY=new TH2D("hPXY","events",50,-0.05,0.05,50,-0.05,0.05);

   acc->SetParameter(0,0.2);

   unsigned int nevt=tr->GetEntries();
   for(unsigned int i=0;i<nevt;++i){
     tr->GetEntry(i);

     event.refresh();

     hMass->Fill(event.p().M());

     if(1){
       TLorentzVector xcharged;
     // save momenta
     vector<TLorentzVector> momenta;
     bool accept=true;
     for(unsigned int ip=0;ip<event.nParticles();++ip){
       if(event.getParticle(ip).q()==0)continue; // only look at charged particles
       momenta.push_back(event.getParticle(ip).p());
       xcharged += momenta[ip];
       double mag=momenta[ip].Vect().Mag();
       hMom->Fill(mag);
       hP->Fill(momenta[ip].Vect().Theta(),
		momenta[ip].Vect().Mag());
       hPXY->Fill(momenta[ip].Vect().X()/mag,
		  momenta[ip].Vect().Y()/mag);
       // cut away events with small momentum particles:
       //double w=acc->Eval(mag);
       //accept &= gRandom->Uniform()<w;

     } // end loop over particles

     hExcl->Fill(xcharged.E());


     event.toGJ();

     // loop over all states that contain n-1 final state particles
     // and plot GJ angles
     unsigned int npart=event.nParticles();
     unsigned int nstates=event.nStates();
     for(unsigned int is=0;is<nstates;++is){
       const NParticleState& state=event.getState(is);
       if(state.qabs()==5 && state.n()==5){
	 double t=state.t();
	 double m=state.p().M();
	 double p2=state.beam().Vect()*state.beam().Vect();
	 double term=m*m-0.019479835;
	 double tmin=term*term*0.25/p2;
	 //cout << "t="<<t<<"     tmin= " << tmin << endl;
	 double tprime=t-tmin;
	 hT->Fill(tprime);
       }

       if(state.n()==npart){




       }

       if(state.n()==npart-1 && state.q()==0){
	 // cut out a bit of the CosGJ distribution
	 hGJ->Fill(state.p().CosTheta());
	 //double w=acc->Eval(state.p().CosTheta());
	 //acc &= gRandom->Uniform()<w;


	 if(accept){
	   hGJacc->Fill(state.p().CosTheta());
	 }
	 else {
	   for(unsigned int ip=0;ip<event.nParticles();++ip){
	     double mag=momenta[ip].Vect().Mag();
	     hMombad->Fill(mag);
	     hPbad->Fill(momenta[ip].Vect().Theta(),
			 momenta[ip].Vect().Mag());
	     hPXYbad->Fill(momenta[ip].Vect().X()/mag,
			 momenta[ip].Vect().Y()/mag);
	   } // end loop over particles

	 }

       }
       else if(state.n()==npart-2 && state.q()==-1){

       }
       else if(state.n()==npart-3 && state.q()==0){

       }
     }
     } // end if(0)


   }// end loop over events






   TCanvas* c=new TCanvas("KineValidate"+massbin,"Events",10,10,1000,800);
   c->Divide(3,3);
   c->cd(1);
   hExcl->Draw();
   c->cd(3);
   hT->Draw();

   c->cd(2);
   hP->Draw("colz");



   c->cd(4);
   hMass->Draw();

   c->cd(5);
hGJ->Draw();
   hGJ->GetYaxis()->SetRangeUser(0,hGJ->GetMaximum()*1.5);
   hGJacc->SetLineColor(kBlue);
   hGJacc->SetFillColor(kBlue);

   hGJacc->Draw("same");

   //hPXYbad->Draw("surf2");
c->cd(6);
   hPXY->Draw("surf2");

   c->cd(7);
 hMombad->Draw();
c->cd(8);
   hMom->Draw();

   TCanvas* cm=new TCanvas("cm","Mass",10,10,1000,800);
   hMass->Draw();
   TF1* mparm=new TF1("fm","(x-[0])*exp((x-[0])*([1]+(x-[0])*[2]))",0.9,3.);
   hMass->Fit(mparm,"","",0.9,3);

}
