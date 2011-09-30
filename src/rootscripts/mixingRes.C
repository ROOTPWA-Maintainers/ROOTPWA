#include <complex>
#include "TCanvas.h"
#include "TF1.h"

using namespace std;


complex<double> singleMixedBW(double* x, double* par){
  double s=x[0]*x[0];
  double m1=par[0]; // mass of resonance
  double gamma1=par[1]; // width of resonance 
  double m2=par[2];
  double gamma2=par[3];
  complex<double> B(par[4],par[5]); // mixing strength
  
  complex<double> bw1(m1*m1-s,-gamma1*m1);
  complex<double> bw2(m2*m2-s,-gamma2*m2);

  complex<double> D=bw1 - B/bw2;
  return 1./D;
}


double singleMixedBWamp(double* x, double* par){
  return norm(singleMixedBW(x,par));
}

double singleMixedBWphase(double* x, double* par){
  return arg(singleMixedBW(x,par));
}

complex<double> doubleMixedBW(double* x, double* par){
  double* par1=par;
  double* par2=&par[6];
  return singleMixedBW(x,par1)+ 0.5*singleMixedBW(x,par2);
}

double doubleMixedBWamp(double* x, double* par){
  return norm(doubleMixedBW(x,par));
}

double doubleMixedBWphase(double* x, double* par){
  return arg(doubleMixedBW(x,par));
}

void mixingRes(double m1, double m2,double bre, double bim){
  double gamma1=0.3;
  double gamma2=0.3;
  double Bre=bre; double Bim=bim;

  TF1* fAmp1=new TF1("fAmp1",&singleMixedBWamp,1.0,2.5,6);
  TF1* fPhase1=new TF1("fPhase1",&singleMixedBWphase,1.0,2.5,6);
  TF1* fAmp2=new TF1("fAmp2",&singleMixedBWamp,1.0,2.5,6);
  TF1* fPhase2=new TF1("fPhase2",&singleMixedBWphase,1.0,2.5,6);
  TF1* fAmp=new TF1("fAmp",&doubleMixedBWamp,1.0,2.5,12);
  TF1* fPhase=new TF1("fPhase",&doubleMixedBWphase,1.0,2.5,12);

  double par[12]={m1,gamma1,m2,gamma2,Bre,Bim,m2,gamma2,m1,gamma1,Bre,Bim};

  double par1[6]={m1,gamma1,m2,gamma2,Bre,Bim};
  double par2[6]={m2,gamma2,m1,gamma1,Bre,Bim};

  fAmp1->SetParameters(par1);
  fPhase1->SetParameters(par1);
  fAmp2->SetParameters(par2);
  fPhase2->SetParameters(par2);
   fAmp->SetParameters(par);
  fPhase->SetParameters(par);


  TCanvas* c=new TCanvas("Mixing","Mixing",10,10,800,800);
  c->Divide(2,3);
  c->cd(1);fAmp1->Draw();
  c->cd(2);fPhase1->Draw();
  c->cd(3);fAmp2->Draw();
  c->cd(4);fPhase2->Draw();
   c->cd(5);fAmp->Draw();
  c->cd(6);fPhase->Draw();
}
