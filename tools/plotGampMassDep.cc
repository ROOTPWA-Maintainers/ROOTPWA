
#include <iostream>
#include <string>
#include <complex>
#include <particle.h>

#include <Vec.h>
#include <massDep.h>
#include <particleData.h>

#include "TGraph.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TROOT.h"

using namespace std;

extern particleDataTable PDGtable;

void dressGraph(TGraph *g, int font=132){
  

  g->GetXaxis()->SetLabelFont(font);
  g->GetXaxis()->SetTitleFont(font);
  g->GetYaxis()->SetLabelFont(font);
  g->GetYaxis()->SetTitleFont(font);
 
}



int main(int argc, char** argv){
	string opt = "BW";
	string part1 = "pi";
	string part2 = "pi";
	string part  = "";
	double mass(0.6);
	double width(0.6);

	if (argc < 2){
		cout << " usage: plotGampMassDep AMP|VES|KACH|LASS|RPRIME|BW [inv_mass(def 600MeV)] [width(def 600MeV)] [pi|K (def pi)] [pi|K (def pi)] [PDG name (mass and width will be ignored)]" << endl;
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
		if (argc > 6){
		  part = (argv[6]);
		  cout << " using PDG values for " << part << endl; 
		}
	}

  PDGtable.initialize();
  particle myp;
  particle pi1;pi1.setMass(PDGtable.get(part1).Mass());pi1.setCharge(1);
  particle pi2;pi2.setMass(PDGtable.get(part2).Mass());pi2.setCharge(-1);
  if (part != ""){
    mass = PDGtable.get(part).Mass();
    width= PDGtable.get(part).Width();
  }

 cout << " plotting " << opt << " with a mass of " << mass << " and width of " << width << " decaying into " << part1 << " " << part2 << endl;

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
  else if(opt=="LASS")dep=new AMP_LASS();
  else if(opt=="RPRIME")dep=new rhoPrime();
  else if(opt=="BW")dep=new breitWigner();
  else dep=new breitWigner();
  dep->print();

  TApplication app("", 0, 0);
  gROOT->SetStyle("Plain");

  unsigned int n=1000;
  TGraph* intens=new TGraph(n);
  TGraph* argand=new TGraph(n);
  TGraph* phase=new TGraph(n);
  
  double mstart=0.1;
  double mend=2.5;
  double mstep=(mend-mstart)/(double)n;
  
  
  for(unsigned int i=0; i<n; ++i){
    
    double mass=mstart+i*mstep;
    myp.set4P(fourVec(mass,threeVec(0,0,0)));
    complex<double> amp=dep->val(myp);
    //cout << amp << endl;
    intens->SetPoint(i,mass,norm(amp));
    phase->SetPoint(i,mass,arg(amp));
    double rho=2*myp.q()/mass;
    
    argand->SetPoint(i,rho*amp.real(),rho*amp.imag());
    
  }



gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetStripDecimals(1);
  TGaxis::SetMaxDigits(3);

  Int_t font=132;

  gStyle->SetTextFont(font);
  gStyle->SetLabelFont(font);

  

  
  TCanvas* c=new TCanvas("c","Mass Dep",10,10,1200,600);
  c->Divide(2,1);
  c->cd(1);
  intens->SetTitle("Intensity");
  intens->Draw("APC");
  dressGraph(intens);
  intens->GetXaxis()->SetTitle("Mass of (#pi#pi) system (GeV/c^{2})");
  intens->GetYaxis()->SetTitle("Intensity");
  c->cd(2);
  phase->SetTitle("Phase");
  phase->Draw("APC");
  dressGraph(phase);
  phase->GetXaxis()->SetTitle("Mass of (#pi#pi) system (GeV/c^{2})");
    phase->GetYaxis()->SetTitle("Phase");
  /*c->cd(3);
  argand->SetTitle("Argand plot");
  argand->Draw("APC");
  argand->SetMarkerStyle(2);
  */
  c->ForceUpdate();
  c->Flush();
  
  gApplication->SetReturnFromRun(kFALSE);
  gSystem->Run();



  return 0;
}
