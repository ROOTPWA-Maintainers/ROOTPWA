
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
	string opt = "BW";
	string part1 = "pi";
	string part2 = "pi";
	double mass(0.6);
	double width(0.6);

	if (argc < 2){
		cout << " usage: plotGampMassDep AMP|VES|KACH|BW [inv_mass(def 600MeV)] [width(def 600MeV)] [pi|K (def pi)] [pi|K (def pi)] " << endl;
		//return 0;
	} else {
		//cout << argc << endl;
		if (argc > 1){
			opt=(argv[1]);
		}
		if (argc > 2){
			mass=atoi(argv[2])/1000.;
		}
		if (argc > 3){
			width=atoi(argv[3])/1000.;
		}
		if (argc > 4){
			part1=(argv[4]);
		}
		if (argc > 5){
			part2=(argv[5]);
		}
	}
  cout << " plotting " << opt << " with a mass of " << mass << " and width of " << width << " decaying into " << part1 << " " << part2 << endl;

  PDGtable.initialize();
  particle myp;
  particle pi1;pi1.setMass(PDGtable.get(part1).Mass());pi1.setCharge(1);
  particle pi2;pi2.setMass(PDGtable.get(part2).Mass());pi2.setCharge(-1);


  decay mydec;
  mydec.setL(0);
  mydec.setS(0);
  mydec.addChild(pi1);
  mydec.addChild(pi2);


  myp.setDecay(mydec);
  myp.setWidth(width);
  myp.setMass(mass);

  massDep* dep;
  if(opt=="AMP")dep=new AMP_M();
  else if(opt=="VES")dep=new AMP_ves();
  else if(opt=="KACH")dep=new AMP_kach();
  else if(opt=="BW")dep=new breitWigner();
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
