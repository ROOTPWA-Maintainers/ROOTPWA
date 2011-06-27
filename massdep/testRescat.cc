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
  if(z.real()==0)return cnum(1,0);
  cnum zm=sqrt(z-4*mpi2);
  cnum term1=0.5*zm/sqrt(z);
  cnum num=sqrt(z)+zm;
  cnum denom=sqrt(z)-zm;
  cnum term2=log(num/denom) - cnum(0,pi);
  return term1*term2;
}


// analytic amplitude
cnum ampIso(cnum z, double a, double b, double c, double d, int l){
  //if(z.real()<4*mpi2)return 0;
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

////////// Basis Function Set /////////////////////////////////////

cnum u(unsigned int a, cnum t){
  if(a==0)return J(t);
  cnum tt=t-4*mpi2;
  if(a==1) return sqrt(tt);
  else if(a==2) return tt;
  else if(a==3) return tt*tt;
  return 0;
}


// expansion in u:
cnum fu(cnum t, cmatrix lambda){
  cnum result(0,0);
  for(unsigned int i=0;i<4;++i){
    result += lambda(i,0)*u(i,t);
  }
  return result;
} 


// matrix expansion in u see eqn (49)
cmatrix fuM(cnum t, cmatrix Lambda){
  cmatrix result(2,2);result(0,0)=0;result(0,1)=0;result(1,0)=0;result(1,1)=0;
  for(unsigned int a=0;a<4;++a){
    cmatrix I(2,2);
    cnum ua=u(a,t);
    for(unsigned int i=0;i<2;++i){
      for(unsigned int j=0;j<2;++j){
	I(i,j)=Lambda(a*2+i,j)*ua;
      }
    }
    result+=I; 
  }
  return result;
}


///////////////// Kernel Integrals ///////////////////////////////

/// we need to integrate the 2x2 kernel matrix element wise from 0->tmax
/// tmax is upper boundary of dalitz plot tmax=(W-mpi)^2
/// We also need the integrals convoluted with the basis funtions....

/// using ROOT GSL Integrator to do the integrals

/// Integrand Functor
class kernIntegrand1: public ROOT::Math::IBaseFunctionOneDim {
public:
  kernIntegrand1():_s(0,0),_t(0,0),_i(0),_j(0),_a(0),_ReIm(0){};
  double DoEval(double x) const;
  ROOT::Math::IBaseFunctionOneDim* Clone() const {return new kernIntegrand1(*this);}

  void setST(cnum s, cnum t){_s=s;_t=t;}
  void setElem(unsigned int i, unsigned int j){_i=i;_j=j;}
  void setReIm(bool flag){_ReIm=flag;} // 0=re; 1=im
  void setU(unsigned int a){_a=a;}
  
private:
  cnum _s;
  cnum _t;
  unsigned int _i;
  unsigned int _j;
  unsigned int _a; // basis function a=0 means no basis function
  bool _ReIm;

};

double
kernIntegrand1::DoEval(double x) const {
  // integration for real z only...
  cnum z(x,0);
  cnum result=Kern1Elem(_t,z,_s,_i,_j);
  if(_a!=0)result*=u(_a-1,z);
  if(_ReIm==0)return result.real();
  else return result.imag();
}


/// Integral Function 
cnum kernInt1Elem(cnum t, cnum s, unsigned int i, unsigned int j, unsigned int b=0){
  // Set up Integrator
   ROOT::Math::GSLIntegrator Integrator;
   // set up function
   kernIntegrand1 f;
   f.setST(s,t);
   f.setElem(i,j);
   f.setU(b);
   // do real part first:
   f.setReIm(0);
   Integrator.SetFunction(f);
   double upper=(sqrt(s.real())-mpi);upper*=upper;
   double Re=Integrator.Integral(f,0.,upper);
   f.setReIm(1);
   double Im=Integrator.Integral(f,0.,upper);
   return cnum(Re,Im);
}

cmatrix kernInt1(cnum t, cnum s,unsigned int b=0){
  cmatrix result(2,2);
  for(unsigned int i=0;i<2;++i)
    for(unsigned int j=0;j<2;++j)result(i,j)=kernInt1Elem(t,s,i,j,b);
  return result;
}

/////////// 4-Point Matching of basis functions ////////////////
// returns column vector containing expansion coefficients
cmatrix match(unsigned int i, unsigned int j, // element of kernel
	      unsigned int b, // basis function
	   cnum s, 
	   const vector<double>& tm)  // matching points
{
  // calulate result vector at points;
  cmatrix res(4,1);
  cmatrix um(4,4);
  for(unsigned int k=0;k<4;++k){
    cnum t(tm[k],0);
    cout << "t=" << tm[k] << endl;
     // calulate result vector at points;
    res(k,0)=kernInt1Elem(t,s,i,j,b);
     // calculate basis matrix at points;
     for(unsigned int h=0;h<4;++h){
       um(k,h)=u(h,t);
     }
  }

  cout << um << endl;
  // invert basis matrix to obtain mathing coefficients
  cmatrix umInv(4,4);
  InvertMatrix(um,umInv);
  
  cmatrix result=prod(umInv,res);
  return result;
}


///////////////////////// Main ///////////////////////////////////
int
main(int argc, char** argv)
{

  TApplication app("", 0, 0);
  gROOT->SetStyle("Plain");


  double mstart=0.005;
  double mstep=0.005;
  unsigned int nsteps=200;
  

  TGraph* gIEps=new TGraph(nsteps);
  TGraph* gPhaseEps=new TGraph(nsteps);
  TGraph* gIRho=new TGraph(nsteps);
  TGraph* gPhaseRho=new TGraph(nsteps);

  TGraph* gJRe=new TGraph(nsteps);
  TGraph* gJIm=new TGraph(nsteps);
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
    gJRe->SetPoint(i,m,J(z).real());
    gJIm->SetPoint(i,m,J(z).imag());
 
  }
	
 
  vector<double> masses(5);
  masses[0]=0.7;
  masses[1]=0.9;
  masses[2]=1.2;
  masses[3]=1.4;
  masses[4]=1.6;
  double tstart=0.001;
  double dt=0.01;


  TMultiGraph* mgKeeRe=new TMultiGraph();
  TMultiGraph* mgKerRe=new TMultiGraph();
  TMultiGraph* mgKreRe=new TMultiGraph();
  TMultiGraph* mgKrrRe=new TMultiGraph();
  TMultiGraph* mgKeeIm=new TMultiGraph();
  TMultiGraph* mgKerIm=new TMultiGraph();
  TMultiGraph* mgKreIm=new TMultiGraph();
  TMultiGraph* mgKrrIm=new TMultiGraph();
  
  
  TMultiGraph* mgIepsepsRe=new TMultiGraph();
  TMultiGraph* mgIepsrhoRe=new TMultiGraph();
  TMultiGraph* mgIrhoepsRe=new TMultiGraph();
  TMultiGraph* mgIrhorhoRe=new TMultiGraph();
  TMultiGraph* mgIepsepsIm=new TMultiGraph();
  TMultiGraph* mgIepsrhoIm=new TMultiGraph();
  TMultiGraph* mgIrhoepsIm=new TMultiGraph();
  TMultiGraph* mgIrhorhoIm=new TMultiGraph();


  vector<double> tmatch(4);
  tmatch[0]=0;

 
  for(unsigned is=0;is<5;++is){ // loop over three pion masses s
   
    cnum s(masses[is]*masses[is],0);
    double tHat=(masses[is]-mpi)*(masses[is]-mpi);
    unsigned int nt=(unsigned int)floor(tHat/(double)dt);

    TGraph* gKeeRe=new TGraph(nt);mgKeeRe->Add(gKeeRe,"P");
    TGraph* gKeeIm=new TGraph(nt);mgKeeIm->Add(gKeeIm,"P");
    TGraph* gKeeReFit=new TGraph(nt);mgKeeRe->Add(gKeeReFit,"C");
    gKeeReFit->SetLineColor(kRed);
    TGraph* gKeeImFit=new TGraph(nt);mgKeeIm->Add(gKeeImFit,"C");
    gKeeImFit->SetLineColor(kRed);

    TGraph* gKerRe=new TGraph(nt);mgKerRe->Add(gKerRe,"P");
    TGraph* gKerIm=new TGraph(nt);mgKerIm->Add(gKerIm,"P");
    TGraph* gKerReFit=new TGraph(nt);mgKerRe->Add(gKerReFit,"C");
    gKerReFit->SetLineColor(kRed);
    TGraph* gKerImFit=new TGraph(nt);mgKerIm->Add(gKerImFit,"C");
    gKerImFit->SetLineColor(kRed);

    TGraph* gKreRe=new TGraph(nt);mgKreRe->Add(gKreRe,"P");
    TGraph* gKreIm=new TGraph(nt);mgKreIm->Add(gKreIm,"P");
    TGraph* gKreReFit=new TGraph(nt);mgKreRe->Add(gKreReFit,"C");
    gKreReFit->SetLineColor(kRed);
    TGraph* gKreImFit=new TGraph(nt);mgKreIm->Add(gKreImFit,"C");
    gKreImFit->SetLineColor(kRed);

    TGraph* gKrrRe=new TGraph(nt);mgKrrRe->Add(gKrrRe,"P");
    TGraph* gKrrIm=new TGraph(nt);mgKrrIm->Add(gKrrIm,"P");
    TGraph* gKrrReFit=new TGraph(nt);mgKrrRe->Add(gKrrReFit,"C");
    gKrrReFit->SetLineColor(kRed);
    TGraph* gKrrImFit=new TGraph(nt);mgKrrIm->Add(gKrrImFit,"C");
    gKrrImFit->SetLineColor(kRed);


    TGraph* gIeeRe=new TGraph(nt);mgIepsepsRe->Add(gIeeRe,"C");
    TGraph* gIeeIm=new TGraph(nt);mgIepsepsIm->Add(gIeeIm,"C");
    TGraph* gIerRe=new TGraph(nt);mgIepsrhoRe->Add(gIerRe,"C");
    TGraph* gIerIm=new TGraph(nt);mgIepsrhoIm->Add(gIerIm,"C");
    TGraph* gIreRe=new TGraph(nt);mgIrhoepsRe->Add(gIreRe,"C");
    TGraph* gIreIm=new TGraph(nt);mgIrhoepsIm->Add(gIreIm,"C");
    TGraph* gIrrRe=new TGraph(nt);mgIrhorhoRe->Add(gIrrRe,"C");
    TGraph* gIrrIm=new TGraph(nt);mgIrhorhoIm->Add(gIrrIm,"C");

    // matching points:
    tmatch[1]=(4*mpi2*0.8+tHat*0.2);
    tmatch[2]=0.5*(4*mpi2+tHat);
    tmatch[3]=0.9*tHat;

    // do matching:
  
    cmatrix lambda(8,2);
    cmatrix R(8,8);
    for(unsigned int i=0;i<2;++i){
      for(unsigned int j=0;j<2;++j){
	// I0 -> lambda
	cmatrix lambdaij=match(i,j,0,s,tmatch);
	for(unsigned int a=0;a<4;++a){
	  lambda(i+a*2,j)=lambdaij(a,0);
	}
	//Expand Kernel in Basis functions
        cmatrix LambdaAB(4,1);
	for(unsigned int b=0;b<4;++b){
	  LambdaAB=match(i,j,b+1,s,tmatch);
	  for(unsigned int a=0;a<4;++a){
	    R(i+a*2,j+b*2)=LambdaAB(a,0);
	  }
	}
	
	
      }
    }      
  
    ublas::identity_matrix<cnum> uni(8,8);
    R=uni-R;
    cout << R << endl;

    ///// At this point we have set up everything and are ready to solve ////
    ///// The linear system and thereby the system of integral equations ////
    ////////// Invert R to solve  ////
    cmatrix T(8,8);
    
    InvertMatrix(R,T);

    cout << T << endl;
    
    /// construct expansion coefficients matrix
    cmatrix IA=prod(T,lambda);

    cmatrix lambda00=match(0,0,0,s,tmatch);
    cmatrix lambda01=match(0,1,0,s,tmatch);
    cmatrix lambda10=match(1,0,0,s,tmatch);
    cmatrix lambda11=match(1,1,0,s,tmatch);

    cout << lambda << endl;

    cmatrix la(4,1);
    for(unsigned int a=0;a<4;++a){
      la(a,0)=lambda(0+a*2,1);
    }
    cout << "##### LA: " << endl;
    cout << la << endl;
    cout << lambda00 << endl;
    cout << lambda01 << endl;
    cout << lambda10 << endl;
    cout << lambda11 << endl;

   
    cerr << "Start loop over subenergy" << endl;
    for(unsigned it=0;it<nt;++it){ // loop over subenergy
      cnum t(it*dt+tstart,0);
      if(t.real()<tHat){
	cmatrix k=kernInt1(t,s);
	

	//cnum kfit=fu(t,lambda00);
	cnum kfit=fuM(t,lambda)(0,0);
	
	/// Kernel Integrals (inhomogenous term I0) 
     	gKeeRe->SetPoint(it,t.real(),k(0,0).real());
	gKeeIm->SetPoint(it,t.real(),k(0,0).imag());
	gKeeReFit->SetPoint(it,t.real(),kfit.real());
	gKeeImFit->SetPoint(it,t.real(),kfit.imag());


	//kfit=fu(t,lambda01);
	kfit=fuM(t,lambda)(0,1);

	gKerRe->SetPoint(it,t.real(),k(0,1).real());
	gKerIm->SetPoint(it,t.real(),k(0,1).imag());
	gKerReFit->SetPoint(it,t.real(),kfit.real());
	gKerImFit->SetPoint(it,t.real(),kfit.imag());


	//kfit=fu(t,lambda10);
	kfit=fuM(t,lambda)(1,0);

	gKreRe->SetPoint(it,t.real(),k(1,0).real());
	gKreIm->SetPoint(it,t.real(),k(1,0).imag());
	gKreReFit->SetPoint(it,t.real(),kfit.real());
	gKreImFit->SetPoint(it,t.real(),kfit.imag());


	//kfit=fu(t,lambda11);
	kfit=fuM(t,lambda)(1,1);
	gKrrRe->SetPoint(it,t.real(),k(1,1).real());
	gKrrIm->SetPoint(it,t.real(),k(1,1).imag());
	gKrrReFit->SetPoint(it,t.real(),kfit.real());
	gKrrImFit->SetPoint(it,t.real(),kfit.imag());
	
	/// Full rescattering
	cmatrix I=fuM(t,IA);
	gIeeRe->SetPoint(it,t.real(),I(0,0).real());
	gIeeIm->SetPoint(it,t.real(),I(0,0).imag());
	gIerRe->SetPoint(it,t.real(),I(0,1).real());
	gIerIm->SetPoint(it,t.real(),I(0,1).imag());
	gIreRe->SetPoint(it,t.real(),I(1,0).real());
	gIreIm->SetPoint(it,t.real(),I(1,0).imag());
	gIrrRe->SetPoint(it,t.real(),I(1,1).real());
	gIrrIm->SetPoint(it,t.real(),I(1,1).imag());
      }
    } // end loop over t

  

  } // end loop over s



 



  TCanvas* c=new TCanvas("c","c",10,10,1000,1000);
  c->Divide(4,2);
  c->cd(1);
  gIEps->Draw("APC");  
  c->cd(2);
  gPhaseEps->Draw("APC");  
  c->cd(3);
  gIRho->Draw("APC");  
  c->cd(4);
  gPhaseRho->Draw("APC");  
  c->cd(5);
  gJRe->Draw("APC");
  c->cd(6);
  gJIm->Draw("APC");
 
 TCanvas* cK=new TCanvas("cI0","cI0",10,10,1000,1000);
 cK->Divide(4,2);
 cK->cd(1);
 mgKeeRe->Draw("AC");
 cK->cd(2);
 mgKeeIm->Draw("AC");
 cK->cd(3);
 mgKerRe->Draw("AC");
 cK->cd(4);
 mgKerIm->Draw("AC");
 cK->cd(5);
mgKreRe->Draw("AC");
 cK->cd(6);
 mgKreIm->Draw("AC");
 cK->cd(7);
mgKrrRe->Draw("AC");
 cK->cd(8);
 mgKrrIm->Draw("AC");

 TCanvas* cI=new TCanvas("cI","cI",20,20,1000,1000);
 cI->Divide(4,2);
 cI->cd(1);
 mgIepsepsRe->Draw("AC");
 cI->cd(2);
 mgIepsepsIm->Draw("AC");
 cI->cd(3);
 mgIepsrhoRe->Draw("AC");
 cI->cd(4);
 mgIepsrhoIm->Draw("AC");
 cI->cd(5);
mgIrhoepsRe->Draw("AC");
 cI->cd(6);
 mgIrhoepsIm->Draw("AC");
 cI->cd(7);
mgIrhorhoRe->Draw("AC");
 cI->cd(8);
 mgIrhorhoIm->Draw("AC");

  gApplication->SetReturnFromRun(kFALSE);
  gSystem->Run();

  return 0;
}





