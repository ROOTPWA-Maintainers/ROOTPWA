
#include <iostream>
#include <string>
#include <complex>
#include <particle.h>

#include <Vec.h>
#include <massDep.h>
#include <particleData.h>

#include "TGraph.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TROOT.h"

using namespace std;

extern particleDataTable PDGtable;

int main(int argc, char** argv){

  PDGtable.initialize();
  particle myp;
  particle pi1;pi1.setMass(PDGtable.get("pi").Mass());pi1.setCharge(1);
  particle pi2;pi2.setMass(PDGtable.get("pi").Mass());pi2.setCharge(-1);


  decay mydec;
  mydec.setL(0);
  mydec.setS(0);
  mydec.addChild(pi1);
  mydec.addChild(pi2);


  myp.setDecay(mydec);
  myp.setWidth(0.6);
  myp.setMass(0.6);

  massDep* dep;
  string opt;
  if(argc>1)opt=(argv[1]);
  if(opt=="AMP")dep=new AMP_M();
  else if(opt=="VES")dep=new AMP_ves();
  else dep=new breitWigner();
  dep->print();

  TApplication app("", 0, 0);
  gROOT->SetStyle("Plain");

  unsigned int n=1000;
  TGraph* intens=new TGraph(n);
  TGraph* argand=new TGraph(n);
  
  double mstart=0.1;
  double mend=2.5;
  double mstep=(mend-mstart)/(double)n;
  
  
  for(unsigned int i=0; i<n; ++i){
    
    double mass=mstart+i*mstep;
    myp.set4P(fourVec(mass,threeVec(0,0,0)));
    complex<double> amp=dep->val(myp);
    //cout << amp << endl;
    intens->SetPoint(i,mass,norm(amp));
    double rho=2*myp.q()/mass;
    
    argand->SetPoint(i,rho*amp.real(),rho*amp.imag());
    
  }
  
  TCanvas* c=new TCanvas("c","Mass Dep",10,10,1200,600);
  c->Divide(2,1);
  c->cd(1);
  intens->SetTitle("Intensity");
  intens->Draw("APC");
  c->cd(2);
  argand->SetTitle("Argand plot");
  argand->Draw("APC");
  c->ForceUpdate();
  c->Flush();
  
  gApplication->SetReturnFromRun(kFALSE);
  gSystem->Run();



  return 0;
}
