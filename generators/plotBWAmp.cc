
#include <iostream>
#include <string>
#include <complex>


#include "TBWProductionAmp.h"

#include "TGraph.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TROOT.h"

using namespace std;


int main(int argc, char** argv){

  TBWProductionAmp bw(atof(argv[1]),atof(argv[2])); 
  TBWProductionAmp bw2(atof(argv[3]),atof(argv[4])); 

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
   
    complex<double> amp=bw.amp(mass)+std::complex<double>(-0.85,-0.1)*bw2.amp(mass);
    //cout << amp << endl;
    intens->SetPoint(i,mass,norm(amp));
    double rho=1;
    
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
