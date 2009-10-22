#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TLorentzRotation.h"
#include "TH1D.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TMath.h"
#include "TCanvas.h"


/*@brief Transform collection of final state momenta into Gottfried Jackson Frame
 */
void
toGJ(TLorentzVector* beam, TClonesArray* particles, unsigned int n){
  // calculate center of mass vector
  TLorentzVector p(0,0,0,0);
  for(unsigned int i=0;i<n;++i){
    TLorentzVector* pa=(TLorentzVector*)particles->At(i);
    p+=*pa;
  }

  TLorentzVector tempX=p;
  // rotate event into scattering plane
  // get normal vector
  TVector3 y(0,1,0);
  TVector3 N=beam->Vect().Cross(tempX.Vect());
  TVector3 rot=N.Cross(y);
  TRotation t;
  double a=N.Angle(y);
  t.Rotate(a,rot);
  //t.SetXEulerAngles(N.Phi(),N.Theta()-TMath::Pi()*0.5,-TMath::Pi()*0.5);
  TLorentzRotation T(t);
  TLorentzRotation L1(T);
  tempX*=T;
  //tempX.Vect().Print();
  beam->Transform(T);
  //_beamPi.p().Vect().Print();

  // boost to X rest frame
  TVector3 boost=-tempX.BoostVector();
  TLorentzRotation b;
  b.Boost(boost);
  tempX*=b;
  //tempX.Vect().Print();
  beam->Transform(b);

  // put beam along z-axis
  TVector3 beamdir=beam->Vect();
  //std::cout<<"beamDir before rotation:";beamdir.Print();
  a=beamdir.Angle(TVector3(0,0,1));
  //std::cout<<"angle="<<a<<std::endl;
  TRotation t2;
  t2.Rotate(a,TVector3(0,1,0));
  T=TLorentzRotation(t2);
  beam->Transform(T);

 
  for(unsigned int i=0;i<n;++i){
    TLorentzVector* pa=(TLorentzVector*)particles->At(i);
     pa->Transform(L1);
     pa->Transform(b);
     pa->Transform(T);
  }

}


void plotWeightedEvents(TTree* tr){

gROOT->SetStyle("Plain");

TH1D* h2pi=new TH1D("h2pi","2 pion mass",200,0,1.4);


TH1D* hGJ=new TH1D("hGJ","Cos Gottfried-Jackson Theta",80,-1,1);
TH1D* hTY=new TH1D("hTY","Treiman-Yang Phi",80,-TMath::Pi(),TMath::Pi());

double weight;
  TClonesArray* p=new TClonesArray("TLorentzVector");
  TLorentzVector* beam=NULL;
  tr->SetBranchAddress("weight",&weight);
  tr->SetBranchAddress("p",&p);
  tr->SetBranchAddress("beam",&beam);


unsigned int nevt=tr->GetEntries();	


for(unsigned int i=0;i<nevt;++i){
	tr->GetEntry(i);
	
	// transform into GJ F
	toGJ(beam,p,p->GetEntries());
	
	TLorentzVector pi2;
	TLorentzVector* pi=(TLorentzVector*)p->At(0);
	pi2=*pi;
	pi=(TLorentzVector*)p->At(1);
	pi2+=*pi;

	h2pi->Fill(pi2.M(),weight);

	hGJ->Fill(pi2.CosTheta(),weight);
	hTY->Fill(pi2.Phi(),weight);

}// end loop over events


 TCanvas* c=new TCanvas("Predict","Weighted Events",10,10,500,500);
 c->Divide(2,2);
 c->cd(1);
 h2pi->Draw();
 c->cd(3);
 hGJ->Draw();
  c->cd(4);
 hTY->Draw();
}
