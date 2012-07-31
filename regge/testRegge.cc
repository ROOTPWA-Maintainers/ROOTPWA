
///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "NParticleEvent.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include <iostream>
using namespace std;

#include "reggeprop.h"


int main(int argc, char** argv){

  TFile* infile=TFile::Open(argv[1]);
  if(infile==NULL)return 1;

  TTree* tr= (TTree*)infile->Get("pwevents");
  //. double impweight=1;
   TClonesArray* p=new TClonesArray("TLorentzVector");
   TLorentzVector* beam=NULL;
   int qbeam;
   std::vector<int>* q=NULL; 
   tr->SetBranchAddress("p",&p);
   tr->SetBranchAddress("beam",&beam);
   tr->SetBranchAddress("qbeam",&qbeam);
   tr->SetBranchAddress("q",&q);
 
   TVector3 vertex;
   NParticleEvent event(p,q,beam,&qbeam,&vertex);

  unsigned int n=tr->GetEntries();
  TGraph* piontrajectory=new TGraph(n);piontrajectory->SetName("AlphaPi");
  TGraph* pionpropRe=new TGraph(n);pionpropRe->SetName("PiPropagRe");
  TGraph* pionpropIm=new TGraph(n);pionpropIm->SetName("PiPropagIm");
  TGraph* pionprop=new TGraph(n);pionprop->SetName("PiPropag");

  TGraph* kineS=new TGraph(n);kineS->SetName("S");

 

  reggeprop pionProp;
 
 
 
  // loop over phasespace events
  
  for(unsigned int i=0; i<n; ++i){
    tr->GetEntry(i);
   event.refresh();
    // get kinematic variables 
    // CAREFUL with definitions of momentum transfer t, q2 etc
    // in the regge classes momentum transfer from scattering is always negative


    // setup for classic deck effect (only use first three pions for this exercise)

    double tin=0.019479835; // incoming pion (beam)
    TLorentzVector  pq=event.getBeam() - event.p();
    double tout = pq.M2();

    TLorentzVector p12=event.getParticle(0).p()+event.getParticle(1).p();

    TLorentzVector pRegge = event.getBeam() - p12;
    double t=pRegge.M2();

    double s1=p12.M2(); // mass of upper 2pion system

    double s2=tin; // outgoing scattered pion
    TLorentzVector pCM=event.getBeam()+TLorentzVector(0,0,0,0.93827);
    double s=pCM.M2();

    piontrajectory->SetPoint(i,t,pionProp.alphapi(t));
    
    std::complex<double> amp=pionProp.ampBCP(t,s,tin,tout,s1,s2);

    pionpropRe->SetPoint(i,t,amp.real());
    pionpropIm->SetPoint(i,t,amp.imag());
    pionprop->SetPoint(i,t,norm(amp));
    //kineS->SetPoint(i,t,pionProp.S(tin,t,tout,s1,s,s2));
  }

  TFile* outfile=TFile::Open("reggetest.root","RECREATE");

  piontrajectory->Write();

  pionpropRe->Write();
  pionpropIm->Write();
  pionprop->Write();
  kineS->Write();

  outfile->Close();

  return 0;

}
