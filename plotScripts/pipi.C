// Following Eef van Beveren & George Rupp
// arXiv:hep-ph/073286v1

#include <complex>
#include <iostream>
#include <vector>
#include <map>
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "SpecFuncMathMore.h"

using namespace ROOT::Math;
using namespace std;

double const mpi=0.13957018;
double const mpi2=mpi*mpi;
double const mK=0.493677;
double const meta=0.547853;
double const mK2=mK*mK;

typedef pair<double,double> tchannel;

// confinement spectrum (Equation 1):
double confLevel(unsigned int n){
  double omega=0.19; // strength of confinement potential
  double E01 = 1.30; // zero point energy
  return E01 + 2*n*omega;
}


double ladderSum(double E,  unsigned int n=5){
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

std::complex<double> tscat(double E, const vector<tchannel>& channels,
			   unsigned int i){
  const double a=2.90; // 1/GeV
  const double lambda2=1.29*1.29;
  const complex<double> ic(0,1);

   // loop over channels
  complex<double> sum=0;
  for(unsigned int iCh=0;iCh<channels.size();iCh++){
    double m1=channels[iCh].first;
    double m2=channels[iCh].second;
    //if((m1+m2)>E)continue; // threshold!
    double m22=m1+m2; m22=m22*m22;
    double m12=m1*m2;
    double p=sqrt(m12/m22*fabs((E*E-m22)));
    //if((m1+m2)<E)p=sqrt(m12/m22*(E*E-m22));

    double m212=(m1*m1-m2*m2);
    double mu=(E-m212*m212/(E*E*E))/4;
    double ladder=ladderSum(E,15);
    // cout << iCh << "   "
    // 	 << m1 << "   "
    // 	 << m2 << "   "
    // 	 << E << "   "
    // 	 << m22 << "   "
    // 	 << m12 << "   "
    // 	 << p << "   "
    // 	 << m212 << "   "
    // 	 << mu << "   "
    // 	 << ladder << "   " << endl;

    if(m1+m2<E){
      sum+= - 2.*ic*lambda2*mu*p*a*sph_bessel(0,p*a)*sph_hankel1(0,p*a)*ladder;
    }
    else { // modify to account for imaginary breakup mumentum.
      complex<double> k(0,p);
      double sphBessC=TMath::SinH(p*a)/(p*a); // bessel function analytically cont
      double sphHankelC=sphBessC - TMath::CosH(p*a)/(p*a);
      sum+=  2.*lambda2*mu*a*p*sphBessC*sphHankelC*ladder;
    }
  }

  complex<double> Denom(1,0);
  Denom+=sum;

  // choose input and ouput channel to be the same
  unsigned int j=i;
  double mi1=channels[i].first;
  double mi2=channels[i].second;
  double mi22=mi1+mi2; mi22=mi22*mi22;
  double mi12=mi1*mi2;
  double pi=sqrt(mi12/mi22*(E*E-mi22));
  double mi212=(mi1*mi1-mi2*mi2);
  double mui=(E-mi212*mi212/(E*E*E))/4;
  double mj1=channels[j].first;
  double mj2=channels[j].second;
  double mj22=mj1+mj2; mj22=mj22*mj22;
  double mj12=mj1*mj2;
  double pj=sqrt(mj12/mj22*(E*E-mj22));
  double mj212=(mj1*mj1-mj2*mj2);
  double muj=(E-mj212*mj212/(E*E*E))/4;

  double Nom=2.*lambda2*a*sqrt(mui*muj*pi*pj)*ladderSum(E,15)*sph_bessel(0,pi*a)*sph_bessel(0,pj*a);
  //double Nom=1;

  return Nom/Denom;
}

std::complex<double> tprod(double E,double m1,double m2){
   double m22=m1+m2; m22=m22*m22;
  double m12=m1*m2;
  double p=sqrt(m12/m22*(E*E-m22));


  double mu=E/4.; //?????
  double a=2.90; // 1/GeV
  double lambda=1.29;
  double lambda2=lambda*lambda;
  double sum=ladderSum(E,100);

   std::complex<double> Denom(1,0);
  std::complex<double> i(0,1);
  std::complex<double> term=-2.*i*lambda2*mu*p*a*sph_bessel(0,p*a)*sph_hankel1(0,p*a)*sum;

  Denom+=term;

  double Nom=lambda*sph_bessel(0,p*a);

  return Nom/Denom;
}



void pipi(){

  vector< tchannel > channels;
  channels.push_back(tchannel(mpi,mpi));
  // channels.push_back(tchannel(mpi,mK));
  // channels.push_back(tchannel(mK,mK));
  // channels.push_back(tchannel(mpi,meta));
  // channels.push_back(tchannel(mK,meta));

  unsigned int npoints=400;
  TGraph* gRe=new TGraph(npoints);
  TGraph* gIm=new TGraph(npoints);
  TGraph* gIntense=new TGraph(npoints);
  TGraph* gPhase=new TGraph(npoints);
  gRe->SetTitle("Scattering");
  gIm->SetTitle("Scattering Imag");
  gIntense->SetTitle("Scattering Intensity");

  TGraph* gReProd=new TGraph(npoints);
  TGraph* gImProd=new TGraph(npoints);
  TGraph* gIntenseProd=new TGraph(npoints);
  TGraph* gPhaseProd=new TGraph(npoints);
  TGraph* gArgandProd=new TGraph(npoints);

  gReProd->SetTitle("Production");
  gImProd->SetTitle("Production Imag");
  gIntenseProd->SetTitle("Prodcution Intensity");

  unsigned int channel=0;
  double step=0.002;
  double M=channels[channel].first+channels[channel].second+step;
  for(unsigned int i=0;i<npoints; ++i){
    complex<double> amp=tscat(M,channels,0);
    complex<double> pamp=tprod(M,mpi,mpi);
    double p=sqrt(M*M/4.-mpi2);
    cout <<  "tscat("<<M<<")= "<< tscat(M,channels,0)  << endl;
    gRe->SetPoint(i,M,amp.real()*M/(2*p));
    gIm->SetPoint(i,M,amp.imag()*M/(2*p));
    gIntense->SetPoint(i,M,norm(amp));
    gPhase->SetPoint(i,M,arg(amp));

    gReProd->SetPoint(i,M,pamp.real()*sqrt(p));
    gImProd->SetPoint(i,M,pamp.imag()*sqrt(p));
    gIntenseProd->SetPoint(i,M,norm(pamp)*2*p);
    gPhaseProd->SetPoint(i,M,arg(pamp));
    gArgandProd->SetPoint(i,pamp.real()*sqrt(p),pamp.imag()*sqrt(p));
    M+=step;
  }
  cout << "finished" << endl;

  TCanvas* c=new TCanvas("cnew","c",10,10,800,800);
  c->Divide(3,2);
  c->cd(4);
  gIntense->Draw("APL");
  c->cd(5);
  gPhase->Draw("APL");
  c->cd(6);
 gRe->Draw("APL");
  gIm->Draw("PL same");

  c->cd(1);

  gIntenseProd->Draw("APL");
  c->cd(2);

  gPhaseProd->Draw("APL");
  c->cd(3);
 gReProd->Draw("APL");
 gImProd->Draw("PL same");
 //gArgandProd->Draw("APL");

}
