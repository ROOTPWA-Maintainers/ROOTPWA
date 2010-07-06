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

/** @brief test application for dynamic width resonance parameterizations
 */




#include <complex>
#include <vector>
#include <iostream>
#include <iomanip>
#include "TGraph.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TMath.h"
#include "TApplication.h"
#include "mcPhaseSpace.h"
#include "TROOT.h"
#include "cauchyIntegral.h"

using namespace std;
using namespace rpwa;

typedef complex<double> cd;

const double gChargedPionMass = 0.13957018;  // charged pion rest mass [GeV/c^2] [PDG08]


class massDep 
{
public:
  massDep(double M, double width) 
    : mS(M*M), mM(M), mWidth(width){
  double masses[4] = {gChargedPionMass,gChargedPionMass,
		      gChargedPionMass,gChargedPionMass};
  ps=new mcPhaseSpace(4,masses,0,40,400,500000);
  rho0=ps->eval(M);
}
  complex<double> val(double m);

  mcPhaseSpace* phasespace() const {return ps;}
  double get_rho0()const {return rho0;}
  double get_rho(double m)const {return ps->eval(m);}
  double get_ms(double s) const ;
  double disperse(double* x, double* par) const;

private:
  double mS;
  double mM;
  double mWidth;
  mcPhaseSpace* ps;
  double rho0;

  


};


cd
massDep::val(double m){
  double s=m*m;
  double ms=get_ms(s);
  cd N(mWidth*mM,0);
  cd D(mS-s-ms,-mM*mWidth*ps->eval(m)/rho0);
  return N/D;
}


double
massDep::disperse(double* x, double* par) const{
  // x[0] is s' (which will be integrated over in the cauchyIntegral)
  // par[0] is s (which the dispersion integral depends on)
  
  double nom=mM*mWidth*get_rho(sqrt(x[0]))/rho0;
  double denom=(x[0]-par[0])*(x[0]-mS);
  return nom/denom;

}

double
massDep::get_ms(double s) const {
  if(s==mS)return 0;
  TF1* f=new TF1("ms",this,&massDep::disperse,ps->thres()*ps->thres(),1600,1,"massDep","disperse");
  f->SetParameter(0,s);
  vector<realPole> p;
  double residual_s=mM*mWidth*get_rho(sqrt(s))/rho0/(s-mS);
  double residual_mS=mM*mWidth*get_rho(sqrt(mS))/rho0/(mS-s);
  p.push_back(realPole(s,residual_s));
  p.push_back(realPole(mS,residual_mS));
  double low=ps->thres();low*=low;
  cauchyIntegral cauchy(f,p,low,1600);
  double I=cauchy.eval_Hunter(4);
  //cout << "I= " << setprecision(12);
  //cout << I << endl; 
  return (s-mS)/TMath::Pi()*I;
}


int
main(int argc, char** argv)
{

  TApplication app("", 0, 0);
  gROOT->SetStyle("Plain");

  // double mass=1.0;
  //double q=2.0
    
  
  //  double masses[4] = {gChargedPionMass,gChargedPionMass,
  //	      gChargedPionMass,gChargedPionMass};
  //mcPhaseSpace* ps=new mcPhaseSpace(4,masses,0,4,100,100000);
  //ps->rho(1);

TF1* f=new TF1("f","x*x/(x-10)/(x-6)",-2,20);
//f->SetParameter(0,1);
  vector<realPole> p;
  p.push_back(realPole(10,25));
  p.push_back(realPole(6,-9));

  //p.push_back(realPole(6,1));

 //  for(double k=0;k<=1;k+=1){
//     double res=pow(-1.,k)*TMath::Exp(k-0.5)/TMath::Pi();
//     p.push_back(realPole(k-0.5,res));
    
//   }

  cout << "ROOT-Integral = " << f->Integral(-2,20) << endl;
  cauchyIntegral cauchy(f,p,-2,20);
  cout << "Calculating hunter: " << setprecision(12);
  cout << cauchy.eval_Hunter(4) << endl; 

  

  double mr=1.275;
  //double mr=1.7;
  //double wr=0.185;
  double wr=0.30;

  massDep* myBW=new massDep(mr,wr);

  


  unsigned int nsteps=500;
  double step=0.005;
  double m0=1.00;

  TGraph* intens=new TGraph(nsteps); // intensity
  TGraph* rho=new TGraph(nsteps);    // phase space
  TGraph* argand=new TGraph(nsteps); // argand plot
  TGraph* ms=new TGraph(nsteps); // dispersion

  double rho0=myBW->get_rho0();

  //TF1* f2=new TF1("ms2",myBW,&massDep::disperse,myBW->phasespace()->thres()*myBW->phasespace()->thres(),400,1,"massDep","disperse");

  for(unsigned int i=0; i<nsteps; ++i){
    double m=m0+i*step;
    std::cout << m << ".." << flush;
    complex<double> amp=myBW->val(m);
    double r=myBW->get_rho(m)/rho0;
    intens->SetPoint(i,m,norm(amp)*r);
    
    rho->SetPoint(i,m,r);
    argand->SetPoint(i,amp.real(),amp.imag());
    ms->SetPoint(i,m,myBW->get_ms(m*m));

    // calculate dispersion integral m(s)
    // dispersion part:
   //  TF1* f=new TF1("f",ps,&mcPhaseSpace::Evaluate,0.6,15,0,"mcPhaseSpace","Evaluate");
//     vector<realPole> p;
//     p.push_back(realPole(0,1));
    
//     cauchyIntegral cauchy(f,p,0.6,15);
//     cout << "Calculating hunter: " << setprecision(12);
//     cout << cauchy.eval_Hunter(4) << endl; 
  }

  std::cout << std::endl;

  
  



  TCanvas* c=new TCanvas("c","Mass Dep",10,10,600,600);
  c->Divide(2,2);
  c->cd(1);
  intens->SetTitle("Intensity");
  intens->Draw("APC");
  c->cd(2);
 
  argand->SetTitle("Argand plot");
  argand->Draw("AP");


  c->cd(3);
  //f->Draw("AP");

  c->cd(4);
  ms->Draw("APC");

  c->cd(3);
  myBW->phasespace()->getGraph()->Draw("APC");
  //f2->SetParameter(0,mr*mr);
  //f2->Draw("same");  

  c->ForceUpdate();
  c->Flush();



  gApplication->SetReturnFromRun(kFALSE);
  gSystem->Run();



  return 0;

}
