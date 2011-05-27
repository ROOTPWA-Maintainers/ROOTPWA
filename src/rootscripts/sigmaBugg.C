// David Buggs sigma parameterization
// From: J.Phys.G34 (2007) 151

#include <complex>
#include <iostream>
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
//#include "SpecFuncMathMore.h"

//using namespace ROOT::Math;
using namespace std;

double const mpi=0.13957018;
double const mpi2=mpi*mpi;
double const mK=0.493677;
double const mK2=mK*mK;
double const mEta=0.547853;
double const mEta2=mEta*mEta;

double const pi=TMath::Pi();

// sigma mass
double const M=0.953;
double const M2=M*M;

double const b1=1.302;
double const b2=0.340;
double const A=2.426;
double const sa=0.5*mpi2;
double const alpha=1.3; 
double const gKK=0.6;
double const gEtaEta=0.2;
double const g4=0.011;

typedef complex<double> comp;

comp const ic(0,1);

// phase space factor analytically continued below threshold
comp rho(double s, double m2){
  if(s>4*m2)return sqrt(1-4*m2/s);
  else return ic*sqrt(4*m2/s-1);
}

comp dampedRho(double s, double m2){
  return exp(-alpha*fabs(s-4*m2))*rho(s,m2);
}

// 4pi rho parametr
comp rho4pi(double s){
  return 1.0/(1 + exp(7.082 -2.845*s));
}


// empirical factor from BES analysis
comp g1(double s){
  return M*(b1+b2*s)*exp(-(s-M2)/A);
}

comp adler(double s){
  return g1(s)*(s-sa)/(M2-s);
}

// pipi width
comp gamma1(double s){
  return adler(s)*rho(s,mpi2);
}



// KK
comp gamma2(double s){
  return gKK*g1(s)*s/M2*dampedRho(s,mK2);
}

// etaeta
comp gamma3(double s){
  return gEtaEta*g1(s)*s/M2*dampedRho(s,mEta2);
}

// 4pi
comp gamma4(double s){
  if(s<16*mpi2)return 0;
  else return M*g4*rho4pi(s)/rho4pi(M2);
}

comp gammaTot(double s){
  return gamma1(s)+gamma2(s)+gamma3(s)+gamma4(s);
}


// dispersion term (only needed for pipi -> r=real
comp myj1(double s){
  double r=rho(s,mpi2).real();
  return 1./pi*(2 + r*TMath::Log((1-r)/(1+r)));
}

comp z1(double s){
  return myj1(s)-myj1(M2);
}




comp D(double s){
  
  return M2 - s - adler(s)*z1(s) - ic*gammaTot(s); 


}


comp tscat(double m){
  double s=m*m;
  return M*gamma1(s)/D(s); 

}

comp tprod(double m){
  double s=m*m;
  return 1./D(s); 

}


void sigmaBugg(){


 unsigned int npoints=400;
  TGraph* gRe=new TGraph(npoints);
  TGraph* gIm=new TGraph(npoints);
  TGraph* gIntense=new TGraph(npoints);
  TGraph* gPhase=new TGraph(npoints);
  gRe->SetTitle("Scattering");
  gIm->SetTitle("Scattering Imag");
  gIntense->SetTitle("Scattering Intensity");
  gPhase->SetTitle("Scattering Phase");

  TGraph* gReProd=new TGraph(npoints);
  TGraph* gImProd=new TGraph(npoints);
  TGraph* gIntenseProd=new TGraph(npoints);
  TGraph* gPhaseProd=new TGraph(npoints);
  gReProd->SetTitle("Production");
  gImProd->SetTitle("Production Imag");
  gIntenseProd->SetTitle("Production Intensity");
  gPhaseProd->SetTitle("Production Phase");

  double step=0.002;
  double M=2*mpi+step;
  for(unsigned int i=0;i<npoints; ++i){
    complex<double> amp=tscat(M); 
    complex<double> pamp=tprod(M);
    double p=sqrt(M*M/4.-mpi2);
   
    gRe->SetPoint(i,M,amp.real()*M/(2*p));
    gIm->SetPoint(i,M,amp.imag()*M/(2*p));
    gIntense->SetPoint(i,M,norm(amp));
    gPhase->SetPoint(i,M,arg(amp));

    gReProd->SetPoint(i,M,pamp.real());
    gImProd->SetPoint(i,M,pamp.imag());
    gIntenseProd->SetPoint(i,M,norm(pamp)*2*p);
    gPhaseProd->SetPoint(i,M,arg(pamp));

    M+=step;
  }
  cout << "finished" << endl;

  TCanvas* c=new TCanvas("c","c",10,10,800,1000);
  c->Divide(3,2);
  c->cd(1);
  gRe->Draw("APL");
  gIm->Draw("PL same");
  c->cd(2);
  gIntense->Draw("APL");
  c->cd(3);
  gPhase->Draw("APL");
  
  c->cd(4);
  gReProd->Draw("APL");
  gImProd->Draw("PL same");
  c->cd(5);
  gIntenseProd->Draw("APL");
   c->cd(6);
  gPhaseProd->Draw("APL");
  

}









