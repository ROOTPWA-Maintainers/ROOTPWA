//////////////////////////////////////////////////////////////////////
//
//  3pi Rescattering Correction a'la Aitchison & Brehm
//  Phys. Rev. D23 (1981) 1194
//  See that paper for the formulas!
//
/////////////////////////////////////////////////////////////////////


#include <complex>
#include <vector>

#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TSystem.h"
#include "Math/GSLIntegrator.h"
#include "Math/IFunction.h"

#include "matrix.hpp"
#include "io.hpp"

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



//////////////// CONSTANTS ///////////////////////

double const mpi=0.13957018;
double const mpi2=mpi*mpi;
double const mK=0.493677;
double const mK2=mK*mK;
double const mEta=0.547853;
double const mEta2=mEta*mEta;

double const pi=TMath::Pi();





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


//////////////////// Kinematic Functions Eqn 38:  //////////////////////
// s = W^2 = 3body mass (squared)
// t =       2body subenergy as appearing in amplitude
// z =       remaining 2body subenergy integrated over in dispersion relation


cnum R(cnum t, cnum z, cnum s){
  return (t-mpi2)*(z-mpi2)+s*(t+z)-s*s;
}

cnum R23(cnum t, cnum z, cnum s){
  return z*(s-mpi2-t)-2*mpi2*(s-mpi2);
}

cnum R32(cnum t, cnum z, cnum s){
  return t*(s-mpi2-z)-2*mpi2*(s-mpi2);
}

cnum Q(cnum t, cnum z, cnum s){
  return t*z-2*mpi2*(s-mpi2);
}

cnum r23(cnum t, cnum z, cnum s){
  return z+2.*t-s-3*mpi2;
}

cnum r32(cnum t, cnum z, cnum s){
  return t+2.*z-s-3*mpi2;
}

cnum K(cnum z, cnum s){
  cnum W=sqrt(s);
  return sqrt( ((W+mpi)*(W+mpi)-z) * ((W-mpi)*(W-mpi)-z) );
}

cnum phi(cnum t, cnum z, cnum s){
  return -t*z*(t+z-s-3*mpi2)-mpi2*(s-mpi2)*(s-mpi2);
}

cnum Delta(cnum t, cnum z, cnum s){
  cnum K2=K(z,s); cnum K3=K(t,s);
  cnum myR=R(t,z,s);
  cnum rk=(myR-K2*K3)/(myR+K2*K3);
  cnum thetaPhi= phi(t,z,s).real()>0 ? 1 : 0;
  return 1./K3*(log(fabs(rk.real())) + cnum(0,pi)*thetaPhi); 
}


cnum F(cnum t, cnum z, cnum s){
  return (s+mpi2-t)/s;
}

cnum G(cnum t, cnum z, cnum s){
  return (s-mpi2+t)/s;
}

cnum H(cnum t, cnum z, cnum s){
  return ( (z+s-mpi2)*(t+s-mpi2)+4.*s*(t-s+mpi2) )/s;
}





///////////////////////// Isobars ///////////////////////////////////

// analytic function J eqn 39
cnum J(cnum z){
  if(z.real()<4*mpi2)return 0;
  cnum zm=sqrt(z-4*mpi2);
  cnum term1=0.5*zm/sqrt(z);
  cnum num=sqrt(z)+zm;
  cnum denom=sqrt(z)-zm;
  cnum term2=log(num/denom) - cnum(0,pi);
  return term1*term2;
}


// analytic amplitude
cnum ampIso(cnum z, double a, double b, double c, double d, int l){
  cnum A=a+b*z+c/(z-d); // phase shift parameterization eqn 41/42
  cnum bar(1.,0);
  if(l==1)bar=z*(z-4*mpi2); // eqn 40
  cnum denom=A+bar*J(z);
  return 1./denom;
}

// analytic amplitude for epsilon
cnum ampEps(cnum z){
  return ampIso(z,3.12,-4.48,0.32,0.98,0);
}
// analytic amplitude for rho
cnum ampRho(cnum z){
  return ampIso(z,40.5,16.1,116.0,2.89,1);
}


///////// Isobar Matrix //////
// 2x2 matrix containing isobar amplitudes on diagonal /////

cmatrix Iso(cnum z){
  cmatrix result(2,2);
  result(0,0)=ampEps(z);
  result(0,1)=0;result(1,0)=0;
  result(1,1)=ampRho(z);
  return result;
}


////////////// Kernel matrices //////////////////////////////////////

// 1++ Kernel:  rho-epsilon rescattering
// The kernel is a 2x2 matrix 
// representing the rho-rho rho-eps eps-rho and eps-eps rescattering


// See Eqn (44) and (37)
cnum PhiElem(cnum t, cnum z, cnum s, 
	       unsigned int i, unsigned int j){
   cnum K2=K(z,s); cnum K3=K(t,s);
   if(i==0 && j==0) return (R(t,z,s)*Delta(t,z,s)+2.*K2)/(3.*K3*K3);
   else if(i==0 && j==1) return (R23(t,z,s)*Delta(t,z,s)+K2*F(t,z,s))/(K3*K3);
   else if(i==1 && j==0) return R32(t,z,s)*Delta(t,z,s)/3.;
   else if(i==1 && j==1) return 0.5*Q(t,z,s)*Delta(t,z,s);
   else return 0;
}


cmatrix Phi(cnum t, cnum z, cnum s){
  cmatrix result(2,2);
  result(0,0)=PhiElem(t,z,s,0,0);//(R(t,z,s)*Delta(t,z,s)+2.*K2)/(3.*K3*K3);
  result(0,1)=PhiElem(t,z,s,0,1);//(R23(t,z,s)*Delta(t,z,s)+K2*F(t,z,s))/(K3*K3);
  result(1,0)=PhiElem(t,z,s,1,0);//R23(t,z,s)*Delta(t,z,s)/3.;
  result(1,1)=PhiElem(t,z,s,1,1);//0.5*Q(t,z,s)*Delta(t,z,s);
  return result; 
}




cmatrix Kern1(cnum t, cnum z, cnum s){
  cmatrix dPhi=Phi(t,z,s)-Phi(4.*mpi2,z,s);
  return prod(dPhi,Iso(z));
}

cnum Kern1Elem(cnum t, cnum z, cnum s,
		  unsigned int i, unsigned int j){
  cnum dPhi=PhiElem(t,z,s,i,j)-PhiElem(4.*mpi2,z,s,i,j);
  cnum I=(Iso(z))(j,j);
  return dPhi*I;
}




///////////////// Kernel Integrals ///////////////////////////////

/// we need to integrate the 2x2 kernel matrix element wise from 0->tmax
/// tmax is upper boundary of dalitz plot tmax=(W-mpi)^2
/// We also need the integrals convoluted with the basis funtions....

/// using ROOT GSL Integrator to do the integrals

/// Integrand Functor
class kernIntegrand1: public ROOT::Math::IBaseFunctionOneDim {
public:
  kernIntegrand1():_s(0,0),_t(0,0),_i(0),_j(0),_ReIm(0){};
  double DoEval(double x) const;
  ROOT::Math::IBaseFunctionOneDim* Clone() const {return new kernIntegrand1(*this);}

  void setST(cnum s, cnum t){_s=s;_t=t;}
  void setElem(unsigned int i, unsigned int j){_i=i;_j=j;}
  void setReIm(bool flag){_ReIm=flag;} // 0=re; 1=im
  
private:
  cnum _s;
  cnum _t;
  unsigned int _i;
  unsigned int _j;
  bool _ReIm;

};

double
kernIntegrand1::DoEval(double x) const {
  // integration for real z only...
  cnum z(x,0);
  cnum result=Kern1Elem(_t,z,_s,_i,_j);
  if(_ReIm==0)return result.real();
  else return result.imag();
}


/// Integral Function 
cnum kernInt1Elem(cnum t, cnum s, unsigned int i, unsigned int j){
  // Set up Integrator
   ROOT::Math::GSLIntegrator Integrator;
   // set up function
   kernIntegrand1 f;
   f.setST(s,t);
   f.setElem(i,j);
   // do real part first:
   f.setReIm(0);
   Integrator.SetFunction(f);
   double upper=(sqrt(s.real())-mpi);upper*=upper;
   double Re=Integrator.Integral(f,0,upper);
   f.setReIm(1);
   double Im=Integrator.Integral(f,0,upper);
   return cnum(Re,Im);
}



///////////////////////// Main ///////////////////////////////////
int
main(int argc, char** argv)
{

  TApplication app("", 0, 0);
  gROOT->SetStyle("Plain");


  double mstart=0.28;
  double mstep=0.005;
  unsigned int nsteps=200;
  

  TGraph* gIEps=new TGraph(nsteps);
  TGraph* gPhaseEps=new TGraph(nsteps);
  TGraph* gIRho=new TGraph(nsteps);
  TGraph* gPhaseRho=new TGraph(nsteps);

   //double s0=4;
  //double f=0.2;
  
 //  TGraph* gI=new TGraph(nsteps);
//   TGraph* gPhase=new TGraph(nsteps);

//   TGraph* gArgand=new TGraph(nsteps);

//   TGraph* gIBW=new TGraph(nsteps);
//   TGraph* gPhaseBW=new TGraph(nsteps);
//   TGraph* gArgandBW=new TGraph(nsteps);

//   TGraph* gIS=new TGraph(nsteps);
//   TGraph* gPhaseS=new TGraph(nsteps);
//   TGraph* gArgandS=new TGraph(nsteps);

  for(unsigned int i=0;i<nsteps;++i){
    double m=mstart+(double)i*mstep;
    double s=m*m;
    cnum z(s,0);
    std::cout << "m=" << m << "   --------------------- s="<<s <<  std::endl;

    cnum aEps=ampEps(z);
    cnum aRho=ampRho(z);
    gIEps->SetPoint(i,m,norm(aEps));
    gIRho->SetPoint(i,m,norm(aRho));
    gPhaseEps->SetPoint(i,m,arg(aEps));
    gPhaseRho->SetPoint(i,m,arg(aRho));


  }
	
 
  vector<double> masses(5);
  masses[0]=0.8;
  masses[1]=1.0;
  masses[2]=1.2;
  masses[3]=1.4;
  masses[4]=1.6;
  double tstart=0;
  double dt=0.0025;

  vector<TGraph*> gKRe(5);
  vector<TGraph*> gKIm(5);
  TMultiGraph* mgKRe=new TMultiGraph();
  TMultiGraph* mgKIm=new TMultiGraph();


  unsigned int nt=450;
  for(unsigned is=2;is<3;++is){
    gKRe[is]=new TGraph(nt);mgKRe->Add(gKRe[is],"P");
    gKIm[is]=new TGraph(nt);mgKIm->Add(gKIm[is],"P");
    
    cnum s(masses[is]*masses[is],0);
    for(unsigned it=0;it<nt;++it){
      cnum t(it*dt+tstart,0);
      cnum k(0,0);
      if(t.real()<(masses[is]-mpi)*(masses[is]-mpi))k=kernInt1Elem(t,s,0,1);
      gKRe[is]->SetPoint(it,t.real(),k.real());
      gKIm[is]->SetPoint(it,t.real(),k.imag());
    } // end loop over t
  } // end loop over s



 



  TCanvas* c=new TCanvas("c","c",10,10,1000,1000);
  c->Divide(4,3);
  c->cd(1);
  gIEps->Draw("APC");  
  c->cd(2);
  gPhaseEps->Draw("APC");  
  c->cd(3);
  gIRho->Draw("APC");  
  c->cd(4);
  gPhaseRho->Draw("APC");  
  c->cd(5);
  mgKRe->Draw("AP");
  c->cd(6);
  mgKIm->Draw("AP");
 


  gApplication->SetReturnFromRun(kFALSE);
  gSystem->Run();

  return 0;
}





