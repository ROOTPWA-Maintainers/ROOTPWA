#include <complex>
#include <vector>

#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TSystem.h"


#include "TCMatrix.h"

#include "matrix.hpp"
#include "io.hpp"

//typedef TMatrixT<std::complex<double> > TMatrixC;

using namespace std;



 // REMEMBER to update "lu.hpp" header includes from boost-CVS
 #include <vector.hpp>
 #include <vector_proxy.hpp>
 #include <matrix.hpp>
 #include <triangular.hpp>
 #include <lu.hpp>
 #include <io.hpp>

 namespace ublas = boost::numeric::ublas;
 typedef complex<double> cnum;
 typedef ublas::matrix<cnum> cmatrix;

 /* Matrix inversion routine.
    Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
 template<class T>
 bool InvertMatrix (const ublas::matrix<T>& input, ublas::matrix<T>& inverse) {
 	using namespace boost::numeric::ublas;
 	typedef permutation_matrix<std::size_t> pmatrix;
 	// create a working copy of the input
 	matrix<T> A(input);
 	// create a permutation matrix for the LU-factorization
 	pmatrix pm(A.size1());

 	// perform LU-factorization
 	int res = lu_factorize(A,pm);
        if( res != 0 ) return false;

 	// create identity matrix of "inverse"
 	inverse.assign(ublas::identity_matrix<T>(A.size1()));

 	// backsubstitute to get the inverse
 	lu_substitute(A, pm, inverse);

 	return true;
 }




 // simple single channel K function with several poles:
complex<double> K(double s, double f, double s0,
		  const vector<double>& mu2, const vector<double>& g){
  
  unsigned int npoles=mu2.size();
  double k=0;
 // pipi phasespace:
  double rho=sqrt((s-0.0784)/s);
  for(unsigned int i=0;i<npoles;++i){
    double rho0=sqrt((mu2[i]-0.0784)/s);
    k+=g[i]/(mu2[i]-s)*rho/rho0;
  }
  k+=f*(1.0+s0)/(s+s0);
  k*=(s-0.02)/(s+0.2);

  // K-Matrix formula:
  complex<double> denom(1.,-k);
  complex<double> nom(k,0);
  return nom/denom;
}
 // simple breit Wigner sum amplitude:
complex<double> BW(double s, double f, double s0,
		  const vector<double>& mu2, const vector<double>& g){
  
   // pipi phasespace:
  double rho=sqrt((s-0.0784)/s);

  unsigned int npoles=mu2.size();
  complex<double> BW;
  for(unsigned int i=0;i<npoles;++i){
    double gamma=2*g[i]/sqrt(mu2[i]); //gamma=g/2m
    complex<double> nom(gamma,0);
    double rho0=sqrt((mu2[i]-0.0784)/s);
    complex<double> denom(mu2[i]-s,-fabs(gamma)*rho/rho0);
    BW+=nom/denom;
  }
 
  return BW;
}

///////////////////////// T-Matrix a'la Novoseller ///////////////

cmatrix SMatrix(double s){


  cmatrix S(2,2);
  cmatrix SB(2,2); // Background
  cmatrix SR(2,2); // Resonance

  // parameterize Background:
  // SB=(1+iKb)(1-iKb)^-1

  cmatrix Kb(2,2);
  // Kb = \alpha rho O rho
  cmatrix rho(2,2);
  double rho0=sqrt((s-0.0784)/s);
  double rho1=sqrt((s-0.1)/s);
  rho(0,0)=cnum(rho0,0);
  rho(1,1)=cnum(rho1,0);
  cnum alpha(1,0); // could use a polynomial here
  double theta=0.2; // mixing angle between backgrounds
  cnum c1(cos(theta),0);
  cnum s1(sin(theta),0);
  cmatrix O(2,2);
  O(0,0)=c1;O(0,1)=s1;O(1,0)=s1;O(1,1)=-c1;
  cout << O << endl;

  Kb=alpha* prod(cmatrix(prod(rho,O)),rho);

  cmatrix uni(2,2);
  uni(0,0)=cnum(1,0);
  uni(1,1)=cnum(1,0);

  cnum i(0,-1);

  cmatrix A(2,2);
  A=uni-i*Kb;
  cmatrix Ai(2,2);
  InvertMatrix(A,Ai);
  cmatrix B(2,2);
  B=uni+i*Kb;
  
  SB=prod(B,A);
  std::cout<< SB << std::endl;

  SR(0,0)=cnum(1,0);
  SR(0,1)=cnum(2,1);
  SR(1,0)=cnum(3,2);
  SR(1,1)=cnum(-1,0);

 cmatrix SRT=herm(SR);
  
  std::cout << SR << std::endl;
  std::cout << SRT << std::endl;

  cmatrix SRTSB=prod(SRT,SB);
  return prod(SRTSB,SR);

}





///////////////////////// Main ///////////////////////////////////
int
main(int argc, char** argv)
{

  TApplication app("", 0, 0);
  gROOT->SetStyle("Plain");


  std::cout << SMatrix(1.0) << std::endl;

  double mstart=0.9;
  double mstep=0.02;
  unsigned int nsteps=100;
  
  vector<double> mu2; 
  mu2.push_back(2.43*2.43);
  mu2.push_back(1.77*1.77);
  // mu2.push_back(1.97*1.97);
  vector<double> g;
  g.push_back(0.2);
  g.push_back(-0.1);
  //g.push_back(1);

  double s0=4;
  double f=0.2;
  
  TGraph* gI=new TGraph(nsteps);
  TGraph* gPhase=new TGraph(nsteps);
  TGraph* gArgand=new TGraph(nsteps);

  TGraph* gIBW=new TGraph(nsteps);
  TGraph* gPhaseBW=new TGraph(nsteps);
  TGraph* gArgandBW=new TGraph(nsteps);

  for(unsigned int i=0;i<nsteps;++i){
    double m=mstart+(double)i*mstep;
    double s=m*m;
    complex<double> amp=K(s,f,s0,mu2,g);
    gI->SetPoint(i,m,norm(amp));
    double phi=arg(amp);
    if(i>0){
      double lastphi, lastm;gPhase->GetPoint(i-1,lastm,lastphi);
      if(fabs(phi-lastphi)>=fabs(phi+TMath::Pi()-lastphi))phi+=TMath::Pi();
      else if(fabs(phi-lastphi)>=fabs(phi-TMath::Pi()-lastphi))phi-=TMath::Pi();
    }
    gPhase->SetPoint(i,m,phi);
    gArgand->SetPoint(i,amp.real(),amp.imag());
    complex<double> ampBW=BW(s,f,s0,mu2,g);
    gIBW->SetPoint(i,m,norm(ampBW));
    gPhaseBW->SetPoint(i,m,arg(ampBW));
    gArgandBW->SetPoint(i,ampBW.real(),ampBW.imag());
  }
	
  TCanvas* c=new TCanvas("c","c",10,10,1000,600);
  c->Divide(3,2);
  c->cd(1);
  gI->Draw("APC");
  c->cd(2);
  gPhase->Draw("AP");
  c->cd(3);
  gArgand->Draw("AP");
 c->cd(4);
  gIBW->Draw("APC");
  c->cd(5);
  gPhaseBW->Draw("AP");
  c->cd(6);
  gArgandBW->Draw("AP");


  gApplication->SetReturnFromRun(kFALSE);
  gSystem->Run();

  return 0;
}





