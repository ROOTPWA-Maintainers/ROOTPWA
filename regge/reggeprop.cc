//-----------------------------------------------------------
// File and Version Information:
// $Id$
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
reggeprop::amp(double t, double s, double tin, double tout, double s1, double s2 ){
  

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



double
reggeprop::S(double tin,double t, double tout, double s1, double s, double s2){ // multiparticle kinematic function
  return s-tin-tout+(s1-tin-t)*(s2-tout-t)/(2.*t);
}






double 
reggeprop::alphapi(double t){
  return (t-spi)/(1.-0.5*(t-spi));
}

double 
reggeprop::alphapiprime(){
  return 1.+1./(1.-spi);
}
