
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


#include <iostream>
using namespace std;

#include "reggeprop.h"


int main(int argc, char** argv){




  reggeprop pionProp;


  unsigned int n=2000;
  TGraph* piontrajectory=new TGraph(n);piontrajectory->SetName("AlphaPi");
  TGraph* pionpropRe=new TGraph(n);pionpropRe->SetName("PiPropagRe");
  TGraph* pionpropIm=new TGraph(n);pionpropIm->SetName("PiPropagIm");
  TGraph* pionprop=new TGraph(n);pionprop->SetName("PiPropag");

  TGraph* kineS=new TGraph(n);kineS->SetName("S");

  
  // setup for classic deck effect
  double tin=0.019479835; // incoming pion
  double tout=-0.01; // pomeron
  double s1=0.5; // roughly rho
  double s2=tin; // outgoing scattered pion
  double s=1.5;  // mass squared of 3-pion system
 
  double t=-0.; double tstep=0.0001;
  for(unsigned int i=0; i<n; ++i){
    t-=tstep;
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
