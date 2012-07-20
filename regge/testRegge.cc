
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

  double tin=0.019479835;
  double tout=0.01;
  double s1=0.5; // roughly rhos
  double s2=0.5;
  double s=2.5;
 
  double t=-0; double tstep=0.0001;
  for(unsigned int i=0; i<n; ++i){
    t+=tstep;
    piontrajectory->SetPoint(i,t,pionProp.alphapi(t));
    pionpropRe->SetPoint(i,t,pionProp.amp(t,s,tin,tout,s1,s2).real());
    pionpropIm->SetPoint(i,t,pionProp.amp(t,s,tin,tout,s1,s2).imag());
    pionprop->SetPoint(i,t,norm(pionProp.amp(t,s,tin,tout,s1,s2)));
    kineS->SetPoint(i,t,pionProp.S(tin,t,tout,s1,s,s2));
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
