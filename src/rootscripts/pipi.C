// Following Eef van Beveren & George Rupp
// arXiv:hep-ph/073286v1

#include <complex>
#include <iostream>
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "SpecFuncMathMore.h"

using namespace ROOT::Math;
using namespace std;

double const mpi=0.13957018;
double const mpi2=mpi*mpi;

// confinement spectrum (Equation 1):
double confLevel(unsigned int n){
  double omega=0.19; // strength of confinement potential
  double E01 = 1.30; // zero point energy
  return E01 + 2*n*omega;
}


double ladderSum(double p, unsigned int n=5){
  double E=2*sqrt(p*p+mpi2);
  double sum=0;
  for(unsigned int i=0;i<n;++i){
    double g0l=(i+1)*TMath::Power(4.0,-(double)i);
    sum+= g0l/(confLevel(i)-E);
  }
  return sum;
}

std::complex<double> sph_hankel1( unsigned int l, double x){
  return std::complex<double>(sph_bessel(l,x),sph_neumann(l,x));
}

std::complex<double> tscat(double E){
  double p=sqrt(E*E/4.-mpi2);

  double mu=E/4.; //?????
  double a=2.90; // 1/GeV
  double lambda2=1.29*1.29;
  double sum=ladderSum(p,100);
   
  std::complex<double> Denom(1,0);
  std::complex<double> i(0,1);
  std::complex<double> term=-2.*i*lambda2*mu*p*a*sph_bessel(0,p*a)*sph_hankel1(0,p*a)*sum;

  Denom+=term;

  double Nom=2.*lambda2*mu*a*p*sum*sph_bessel(0,p*a)*sph_bessel(0,p*a);

  return Nom/Denom;
}

std::complex<double> tprod(double E){
  
 double p=sqrt(E*E/4.-mpi2);

  double mu=E/4.; //?????
  double a=2.90; // 1/GeV
  double lambda=1.29;
  double lambda2=lambda*lambda;
  double sum=ladderSum(p,100);

   std::complex<double> Denom(1,0);
  std::complex<double> i(0,1);
  std::complex<double> term=-2.*i*lambda2*mu*p*a*sph_bessel(0,p*a)*sph_hankel1(0,p*a)*sum;

  Denom+=term;

  double Nom=lambda*sph_bessel(0,p*a);

  return Nom/Denom;
}



void pipi(){

  unsigned int npoints=400;
  TGraph* gRe=new TGraph(npoints);
  TGraph* gIm=new TGraph(npoints);
  TGraph* gIntense=new TGraph(npoints);
  gRe->SetTitle("Scattering");
  gIm->SetTitle("Scattering Imag");
  gIntense->SetTitle("Scattering Intensity");

  TGraph* gReProd=new TGraph(npoints);
  TGraph* gImProd=new TGraph(npoints);
  TGraph* gIntenseProd=new TGraph(npoints);
  gReProd->SetTitle("Production");
  gImProd->SetTitle("Production Imag");
  gIntenseProd->SetTitle("Prodcution Intensity");


  double step=0.002;
  double M=2*mpi+step;
  for(unsigned int i=0;i<npoints; ++i){
    complex<double> amp=tscat(M); 
    complex<double> pamp=tprod(M);
    double p=sqrt(M*M/4.-mpi2);
    cout <<  "t("<<M<<")= "<< tscat(M)  << endl;
    gRe->SetPoint(i,M,amp.real()*M/(2*p));
    gIm->SetPoint(i,M,amp.imag()*M/(2*p));
    gIntense->SetPoint(i,M,norm(amp));

    gReProd->SetPoint(i,M,pamp.real()*sqrt(p));
    gImProd->SetPoint(i,M,pamp.imag()*sqrt(p));
    gIntenseProd->SetPoint(i,M,norm(pamp)*2*p);


    M+=step;
  }
  cout << "finished" << endl;

  TCanvas* c=new TCanvas("c","c",10,10,800,800);
  c->Divide(2,2);
  c->cd(1);
  gRe->Draw("APL");
  gIm->Draw("PL same");
  c->cd(3);
  gIntense->Draw("APL");

  c->cd(2);
  gReProd->Draw("APL");
  gImProd->Draw("PL same");
  c->cd(4);
  gIntenseProd->Draw("APL");

}
