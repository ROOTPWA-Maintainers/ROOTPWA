#include <complex>
#include <vector>
#include <fstream>
#include <string>


#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

#include "Math/GSLIntegrator.h"
#include "Math/IFunction.h"
#include <complex>
#include "TCanvas.h"
#include "TF1.h"

using namespace std;
 #include <boost/numeric/ublas/vector.hpp>
 #include <boost/numeric/ublas/vector_proxy.hpp>
 #include <boost/numeric/ublas/matrix.hpp>
 #include <boost/numeric/ublas/triangular.hpp>
 #include <boost/numeric/ublas/lu.hpp>
 #include <boost/numeric/ublas/io.hpp>

 namespace ublas = boost::numeric::ublas;
 typedef complex<double> cnum;
 typedef ublas::matrix<cnum> cmatrix;

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


cmatrix D(double* x, double* par){
  double s=x[0]*x[0];
  double m1=par[0]; // mass of resonance
  double gamma1=par[1]; // width of resonance 
  double m2=par[2];
  double gamma2=par[3];
  complex<double> B(par[4],par[5]); // mixing strength
  
  complex<double> bw1(m1*m1-s,-gamma1*m1);
  complex<double> bw2(m2*m2-s,-gamma2*m2);
  
  cmatrix result(2,2);
  
  cnum denom=bw1*bw2-B*B;
  cnum c=1./denom;

  result(0,0)=c*bw2;
  result(0,1)=c*B;
  result(1,0)=c*B;
  result(1,1)=c*bw1;

  return result;
}

cnum amp(double* x, double* par){
  // x[0]=m;
  // par=m1,gamma1,m2,gamma2,B_re/im,
  // production couplings (2xcomplex), decay couplings (2xreal)

  cmatrix gprod(2,1);
  gprod(0,0)=cnum(par[6],par[7]);
  gprod(1,0)=cnum(par[8],par[9]);
  cmatrix gdecay(1,2);
  gdecay(0,0)=cnum(par[10],0);
  gdecay(0,1)=cnum(par[11],0);
  cmatrix gp=prod(D(x,par),gprod);
  cmatrix gpd=prod(gdecay,gp);
  return gpd(0,0);
  
}

double Intens(double* x, double* par){
  return norm(amp(x,par));
}

double Phase(double* x, double* par){
  return arg(amp(x,par));
}


// standard sum of breitwigeners
cnum ampStan(double* x, double* par){
  double s=x[0]*x[0];
  double m1=par[0]; // mass of resonance
  double gamma1=par[1]; // width of resonance 
  double m2=par[2];
  double gamma2=par[3];
  complex<double> B(par[4],par[5]); // mixing strength
  
  complex<double> bw1(m1*m1-s,-gamma1*m1);
  complex<double> bw2(m2*m2-s,-gamma2*m2);

  complex<double> c1(par[6],par[7]);
  complex<double> c2(par[8],par[9]);
    

  return  1./bw1*c1*par[10]+1./bw2*c2*par[11];

}



double IntensStandard(double* x, double* par){
 return norm(ampStan(x,par));
}

double PhaseStandard(double* x, double* par){
  return arg(ampStan(x,par));
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

int 
main(int argc, char** argv){
  
  
  TApplication app("", 0, 0);
  gROOT->SetStyle("Plain");
  gStyle->SetGridStyle(1);

  double m1= 1.8;
  double m2= 1.8;
  double bre=0.4;
  double bim=0.4;
  double gamma1=0.3;
  double gamma2=0.3;
  double Bre=bre; double Bim=bim;

  TF1* fAmp1=new TF1("fAmp1",&singleMixedBWamp,1.0,2.5,6);
  TF1* fPhase1=new TF1("fPhase1",&singleMixedBWphase,1.0,2.5,6);
  TF1* fAmp2=new TF1("fAmp2",&singleMixedBWamp,1.0,2.5,6);
  TF1* fPhase2=new TF1("fPhase2",&singleMixedBWphase,1.0,2.5,6);


  TF1* fAmpS=new TF1("fAmpS",&IntensStandard,1.0,2.5,12);
  TF1* fPhaseS=new TF1("fPhaseS",&PhaseStandard,1.0,2.5,12);

  TF1* fAmp=new TF1("fAmp",&Intens,1.0,2.5,12);
  TF1* fPhase=new TF1("fPhase",&Phase,1.0,2.5,12);

  double par[12]={m1,gamma1,m2,gamma2,Bre,Bim,1,0.,1.0,0,1,1};

  double par1[6]={m1,gamma1,m2,gamma2,Bre,Bim};
  double par2[6]={m2,gamma2,m1,gamma1,Bre,Bim};

  fAmp1->SetParameters(par1);
  fPhase1->SetParameters(par1);
  fAmp2->SetParameters(par2);
  fPhase2->SetParameters(par2);
   fAmp->SetParameters(par);
  fPhase->SetParameters(par);
  fAmpS->SetParameters(par);
  fPhaseS->SetParameters(par);

  TCanvas* c=new TCanvas("Mixing","Mixing",10,10,800,800);
  c->Divide(2,4);
  c->cd(1);fAmp1->Draw();
  c->cd(2);fPhase1->Draw();
  c->cd(3);fAmp2->Draw();
  c->cd(4);fPhase2->Draw();
   c->cd(5);fAmpS->Draw();
   c->cd(6);fPhaseS->Draw();
  c->cd(7);fAmp->Draw();
  c->cd(8);fPhase->Draw();

  gApplication->SetReturnFromRun(kFALSE);
  gSystem->Run();

}
