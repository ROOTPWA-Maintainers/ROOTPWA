//-----------------------------------------------------------
//
// Description:
//      Implementation of class reggeprop
//      see reggeprop.hh for details
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


// This Class' Header ------------------
#include "reggeprop.h"


// C/C++ Headers ----------------------
#include "TMath.h"
#include <iostream>

using namespace std;



// Collaborating Class Headers --------


// Class Member definitions -----------


std::complex<double>
reggeprop::ampHarris(double t, double s, double tin, double tout, double s1, double s2 ){


  double a=alphapi(t);
  complex<double> term1(TMath::Pi()*alphapiprime()/TMath::Gamma(1+a),0);
  double sinepia=TMath::Sin(-a*TMath::Pi());
  complex<double> term2(1.+TMath::Cos(-a*TMath::Pi()),sinepia);
  term2/=2.*sinepia;
  complex<double> term3(TMath::Power(S(tin,t,tout,s1,s,s2),a),0);


  cout << t << endl;
  cout << a << endl;
  cout << term1 << endl;
  cout << term2 << endl;
  cout << term3 << endl;
  cout << term1*term2*term3 << endl;
  cout << "----------------------------"<<endl;

  return term1*term3*term2;
}

std::complex<double>
reggeprop::ampPDK(double t, double s, double tin, double tout, double s1, double s2 ){
  return TMath::Exp(_a*t)/(t-spi);
}

std::complex<double>
reggeprop::ampFULL(double t, double s, double tin, double tout, double s1, double s2 ){
  double alpha=alphapi(t);
  std::complex<double> i(0,1);
  std::complex<double> num=TMath::Power(s,alpha)*TMath::Exp(_a*t)*std::exp(-i*TMath::PiOver2()*alpha);
  double denom=1./(t-spi);
  return num*denom;
}

std::complex<double>
reggeprop::ampSMU(double t, double s, double tin, double tout, double s1, double s2 ){

  double smu=s+0.5*(tin-tout-s1)-spi;
 double alpha=alphapi(t);
  std::complex<double> i(0,1);
  std::complex<double> num=TMath::Power(smu,alpha)*TMath::Exp(_a*t)*std::exp(-i*TMath::PiOver2()*alpha);
  double denom=1./(t-spi);
  return num*denom;

}



std::complex<double>
reggeprop::ampBCP(double t, double s, double tin, double tout, double s1, double s2 ){

 double alpha=alphapi(t);
 double ss=S(tin,t,tout,s1,s,s2);
 if(ss<0){
   std::cerr << "Problem in S calculation\n"
	     << " S="<<ss << std::endl
	     << " tin="<<tin << std::endl
	     << " t="<<t << std::endl
	     << " tout="<<tout << std::endl
	     << " s1="<<s1 << std::endl
	     << " s="<<s << std::endl
	     << " s2="<<s2 << std::endl;
   return 0;
 }
  std::complex<double> i(0,1);
  std::complex<double> num=TMath::Power(ss,alpha)*TMath::Exp(_a*t)*std::exp(-i*TMath::PiOver2()*alpha);
  double denom=1./(t-spi);
  return num*denom;

}

double
reggeprop::S(double tin,double t, double tout, double s1, double s, double s2){ // multiparticle kinematic function
  return s-tin-tout+(s1-tin-t)*(s2-tout-t)/(2.*t);
}






// pion regge trajectory - for a reggeon t has to be negative t<0
// here: curved regge trajectory ala Pignotti
double
reggeprop::alphapi(double t) const {
  // spi=m_pi^2
  return (t-spi)/(1.-0.5*(t-spi));
  // factor 0.5 appears in Harris paper
  // but not in Berger paper
}

double
reggeprop::alphapiprime() const{
  return _alphapiprime;
}
