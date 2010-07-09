// This Class' Header ------------------
#include "cauchyIntegral.h"

// C/C++ Headers ----------------------
#include "TMath.h"
#include <iostream>

// Collaborating Class Headers --------
#include "gsl/gsl_sf_legendre.h"

// Class Member definitions -----------

using namespace std;


cauchyIntegral::cauchyIntegral(TF1* func, const std::vector<realPole>& poles,
		 double rangelow, double rangeup)
  : _poles(poles),_rlow(rangelow),_rup(rangeup),_range(rangelow+rangeup),_diff(rangeup-rangelow),_func(func)
{
  /// calculated explicitely roots of P_l:
  std::vector<double> r1;
  r1.push_back(0);
  _roots.push_back(r1);
  r1.clear();
  r1.push_back(-1./TMath::Sqrt(3));
  r1.push_back(1./TMath::Sqrt(3));
  _roots.push_back(r1);
  r1.clear();
  r1.push_back(-TMath::Sqrt(3./5.));
  r1.push_back(0);
  r1.push_back(TMath::Sqrt(3./5.));
  _roots.push_back(r1);
  r1.clear();
  r1.push_back(-TMath::Sqrt((3.+2*TMath::Sqrt(6./5.))/7.));
  r1.push_back(-TMath::Sqrt((3.-2*TMath::Sqrt(6./5.))/7.));
  r1.push_back(TMath::Sqrt((3.-2*TMath::Sqrt(6./5.))/7.));
  r1.push_back(TMath::Sqrt((3.+2*TMath::Sqrt(6./5.))/7.));
  _roots.push_back(r1);
  r1.clear();
  r1.push_back(-1./3.*TMath::Sqrt(5.+2*TMath::Sqrt(10./7.)));
  r1.push_back(-1./3.*TMath::Sqrt(5.-2*TMath::Sqrt(10./7.)));
  r1.push_back(0);
  r1.push_back(1./3.*TMath::Sqrt(5.-2*TMath::Sqrt(10./7.)));
  r1.push_back(1./3.*TMath::Sqrt(5.+2*TMath::Sqrt(10./7.)));
  _roots.push_back(r1);
}


double
cauchyIntegral::trafo(double t){
  return 0.5*(_range+_diff*t);
}

double
cauchyIntegral::ofart(double x){
  if(x>_rup || x<_rlow){
    cout << "cauchyIntergral:: x out of range: " << x 
	 << " ("<< _rlow << "," << _rup << ")" << endl;
    throw;
  }
  return (2.* x-_range)/_diff;
}


double 
cauchyIntegral::eval_Hunter(unsigned int degree){
  if(degree>4 || degree < 1) throw;
  // do coordinate transformation
  // x->t x=0.5((a+b)-(b-a)t] dx=0.5(b-a)dt
  // loop over zeros of Pl to calculate G
  double G=0;
  for(unsigned int i=0; i<degree; ++i){
    G+=H(i,degree)*_func->Eval(trafo(_roots[degree-1].at(i)));
  }
  
  //cout << "G = "<< G << endl;

  // loop over signularities
  double R=0;
  for(unsigned int i=0;i<_poles.size();++i){
    double x=ofart(_poles[i].x);
    double res=2./_diff*_poles[i].residual;
    //cout << "Pole" << i << ": " << x << "("<<_poles[i].x<<")" 
    //	 << "  Residual="<<  res << "("<<_poles[i].residual<<")"  <<endl;
    R+=res*gsl_sf_legendre_Ql(degree,x)/gsl_sf_legendre_Pl(degree,x);
  }


  //cout << "R = "<< R << endl;

  // rescale coordinate trafo [-1,1]->[a,b]
  return 0.5*_diff*(G-2.*R);
}



double 
cauchyIntegral::H(unsigned int r,     /// index of zero of legendre polynomial
		  unsigned int degree){
  if(degree>4 || degree < 1) throw;
  if(r>degree) throw;
  
  double pl[degree+2];
  double plderiv[degree+2];
  
  gsl_sf_legendre_Pl_deriv_array(degree+1,_roots[degree-1].at(r),pl,plderiv);
  
  
  double denom=(degree+1)*pl[degree+1]*plderiv[degree];
  

  return -2./denom;
}
