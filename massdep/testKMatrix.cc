#include <complex>
#include <vector>

#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TSystem.h"


#include "TCMatrix.h"

//#include "matrix.hpp"
//#include "io.hpp"

//typedef TMatrixT<std::complex<double> > TMatrixC;

using namespace std;



 // REMEMBER to update "lu.hpp" header includes from boost-CVS
 #include <boost/numeric/ublas/vector.hpp>
 #include <boost/numeric/ublas/vector_proxy.hpp>
 #include <boost/numeric/ublas/matrix.hpp>
 #include <boost/numeric/ublas/triangular.hpp>
 #include <boost/numeric/ublas/lu.hpp>
 #include <boost/numeric/ublas/io.hpp>

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

double rhoPiPi(double s){
  return sqrt((s-0.0784)/s);
}
double rho2(double s){
  return sqrt((s-0.2)/s);
}


cmatrix TMatrix(double s){


  cmatrix S(2,2);
  cmatrix SB(2,2); // Background
  cmatrix SR(2,2); // Resonance

  // parameterize Background:
  // SB=(1+iKb)(1-iKb)^-1

  cmatrix Kb(2,2);
  // Kb = \alpha rho O rho
  cmatrix rho(2,2);
  double rho0=rhoPiPi(s);
  double rho1=rho2(s);
  rho(0,0)=cnum(rho0,0);
  rho(1,1)=cnum(rho1,0);
  cnum alpha(0.1,0); // could use a polynomial here
  double theta=3.; // mixing angle between backgrounds
  cnum c1(cos(theta),0);
  cnum s1(sin(theta),0);
  cmatrix O(2,2);
  O(0,0)=c1;O(0,1)=s1;O(1,0)=s1;O(1,1)=-c1;
  

  Kb=alpha * prod(cmatrix(prod(rho,O)),rho);

  cmatrix uni(2,2);
  uni(0,0)=cnum(1,0);
  uni(1,1)=cnum(1,0);
  uni(0,1)=cnum(0,0);
  uni(1,0)=cnum(0,0);


  cnum i(0,-1);

  
  //Kb=uni;
  //Kb(0,0)=cnum(0,0);
  //Kb(1,1)=cnum(0,0);

  cmatrix A(2,2);
  A=uni-i*Kb;
  cmatrix Ai(2,2);
  InvertMatrix(A,Ai);
  cmatrix B(2,2);
  B=uni+i*Kb;
  
  SB=prod(B,Ai);

  std::cout<< "SB=" << SB << std::endl;
  std::cout<< "SB*SBT=" << prod(SB,herm(SB)) << std::endl;

  // two resonances:
  // SR=S1*S2
  // m2>m1
  // Sk=1+[i-xk+sqrt(1+xk^2)]Tk
  // xk=(mk-sqrt(s))/(0.5\Sum_r\gamma_rk)
  // T_k(ij)=(0.5\epsik\epsjk\sqrt(\gamma_ik\gamma_jk))/(mk-sqrt(s)-0.5iSum_r \gamma_rk
  // \gamma_ik=\Gamma_ik\frac{\rho(m)}{\rho(mk)\}

  //channel 1
  double m1=1.8;
  double m2=2.4;
  double G1=0.2;double G11=0.95*G1; double G12=(1.-0.95)*G1;
  double G2=0.2;double G21=0.05*G2; double G22=(1.-0.05)*G2;
  double eps11=1;double eps12=1;double eps21=1;double eps22=1;
  double g11=G11*rhoPiPi(s)/rhoPiPi(m1*m1);
  double g12=G12*rho2(s)/rho2(m1*m1);
  double g21=G21*rhoPiPi(s)/rhoPiPi(m2*m2);
  double g22=G22*rho2(s)/rho2(m2*m2);

  // resonance 1:
  double m=sqrt(s);
  cnum denom1(m1-m,-0.5*(g11+g12));
  cmatrix T1(2,2);
  T1(0,0)=cnum(0.5*eps11*eps11*sqrt(g11*g11),0)/denom1;
  T1(0,1)=cnum(0.5*eps11*eps12*sqrt(g11*g12),0)/denom1;
  T1(1,0)=cnum(0.5*eps12*eps11*sqrt(g12*g11),0)/denom1;
  T1(1,1)=cnum(0.5*eps12*eps12*sqrt(g12*g12),0)/denom1;
  


  double x1=2.*(m1-m)/(g11+g12);
  cmatrix S1=uni+cnum(sqrt(1+x1*x1)-x1,1)*T1;

  std::cout << "S1*S1T=" << prod(S1,herm(S1)) << std::endl;

  // resonance 2:
  cnum denom2(m2-m,-0.5*(g21+g22));
  cmatrix T2(2,2);
  T2(0,0)=cnum(0.5*eps21*eps21*sqrt(g21*g21),0)/denom2;
  T2(0,1)=cnum(0.5*eps21*eps22*sqrt(g21*g22),0)/denom2;
  T2(1,0)=cnum(0.5*eps22*eps21*sqrt(g22*g21),0)/denom2;
  T2(1,1)=cnum(0.5*eps22*eps22*sqrt(g22*g22),0)/denom2;
  
  double x2=2.*(m2-m)/(g21+g22);
  cmatrix S2=uni+cnum(sqrt(1+x2*x2)-x2,1)*T2;

  SR=prod(S1,S2);
  //SR=uni;
  cmatrix SRT=herm(SR);
  std::cout << "SR="<< SR << std::endl;
  std::cout << "SRT="<<SRT << std::endl;
  
  // check unitarity
  std::cout << "SR*SRT=" << prod(SR,SRT) << std::endl;

  //std::cout << SR << std::endl;
  //std::cout << SRT << std::endl;

  cmatrix SRTSB=prod(SRT,SB);
  S=prod(SRTSB,SR);
  cmatrix Sminus1=S-uni;

  // check unitarity
  std::cout << "S="<< S << std::endl;
  std::cout << "S-1="<< Sminus1 << std::endl;

  std::cout << "S*ST=" << prod(S,herm(S)) << std::endl;
 
  cmatrix Tp=(-0.5*i)*Sminus1;
  cmatrix T=Tp;
  std::cout << "T="<< T << std::endl;
  
  cmatrix Tt=herm(T);

  cmatrix Tdiff=T - Tt;
  cmatrix Tprod=(2.*i)*prod(Tt,T);

  
  std::cout << "Tdiff="<< Tdiff << std::endl;

  std::cout << "Tprod="<< Tprod << std::endl;

  cmatrix TUnitarity=Tdiff - Tprod;

  std::cout << "(T-T^t) - 2iT^tT = "<< TUnitarity << std::endl;

  return T;
  
}





///////////////////////// Main ///////////////////////////////////
int
main(int argc, char** argv)
{

  TApplication app("", 0, 0);
  gROOT->SetStyle("Plain");


  double mstart=1.1;
  double mstep=0.02;
  unsigned int nsteps=200;
  
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

  TGraph* gIS=new TGraph(nsteps);
  TGraph* gPhaseS=new TGraph(nsteps);
  TGraph* gArgandS=new TGraph(nsteps);

  for(unsigned int i=0;i<nsteps;++i){
    double m=mstart+(double)i*mstep;
    double s=m*m;
    
    std::cout << "m=" << m << "   ---------------------" <<  std::endl;


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

    cnum ampT=(TMatrix(s))(0,0);
    
    gIS->SetPoint(i,m,norm(ampT));
    gPhaseS->SetPoint(i,m,arg(ampT));
    gArgandS->SetPoint(i,ampT.real(),ampT.imag());
  }
	
  TCanvas* c=new TCanvas("c","c",10,10,1000,1000);
  c->Divide(3,3);
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
 c->cd(7);
  gIS->Draw("APC");
  c->cd(8);
  gPhaseS->Draw("AP");
  c->cd(9);
  gArgandS->Draw("AP");



  gApplication->SetReturnFromRun(kFALSE);
  gSystem->Run();

  return 0;
}





