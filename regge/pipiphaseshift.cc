//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of pipiphaseshift classes
//      see pipiphaseshift.h for details
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


// This Class' Header ------------------
#include "pipiphaseshift.h"


// C/C++ Headers ----------------------
#include "TMath.h"
#include <iostream>

using namespace std;



// Collaborating Class Headers --------


// Class Member definitions -----------

std::complex<double>
absPhaseshift::w(std::complex<double> s) const
{
  complex<double> s1=sqrt(_s0-s);
  complex<double> sqs=sqrt(s);
  return (sqs-s1)/(sqs+s1);

}


/////////////////////////////////////////////////////////
s0wave::s0wave(double z02,
	       double B0, double B1, double B2) 
  : _z02(z02),_B0(B0),_B1(B1),_B2(B2),
    _alpha({0.843,0.2}),_beta({1.02,1.33}),
    _mu(1.0), _invmu(1./_mu), _M12(0.888*0.888), _M22(1.327*1.327) {
  _s0=4.0*0.493677*0.493677; /// 4 m_K^2
  _gamma.resize(2,2);
  _gamma(0,0)=3.10;_gamma(0,1)=1.82;
  _gamma(1,0)=1.82;_gamma(1,1)=-7.0;
 
 
}

double
s0wave::cotDelta(std::complex<double> s) const
{
  complex<double> k1=kThr(s,mpic2);
  complex<double> sqs=sqrt(s);
  complex<double> result(0,0);
  // way below the KKbar threshold sqrt(s)<932MeV:
  if(s.real()<0.868624){
    complex<double> myw=w(s);
    complex<double> adler=0.5*sqs/k1 * mpic2/(s-0.5*_z02);
    result = _z02/(mpic*sqs) + _B0 + _B1*myw + _B2*myw*myw;
      result *=adler;
  } // end way below KKbar threshold
  // K-Matrix parameterization
   else{

    //Build K-Matrix
    rmatrix _K(2,2);

    for(unsigned int i=0;i<2;++i)
      for(unsigned int j=0;j<2;++j){
	_K(i,j)=_mu*(_alpha[i]*_alpha[j]/(_M12-s.real()) + _beta[i]*_beta[j]/(_M22-s.real())) +_invmu*_gamma(i,j);
      }
    
    double _detK=_K(0,0)*_K(1,1)-_K(1,0)*_K(0,1);
    
    complex<double> k2=kThr(s,mKc2);
    // intermediate region but still below KKbar threshold = elastic`
    if(s.real()<mKc2){
      result=(1.0+k2*_K(1,1))/(k1*fabs(k2)*_detK + k1*_K(0,0));
    } // end if in intermediate region close but below to KKbar threshold
    // above KKbar 
    else {
      complex<double> k12=k1*k1;
      complex<double> k22=k2*k2;
      complex<double> term1=k12*_K(0,0)*_K(0,0)+ k12*k22*_detK*_detK;
      complex<double> term2=k22*_K(1,1)*_K(1,1);
      complex<double> sum=term1 + term2;
      complex<double> K124=_K(0,1)*_K(0,1);K124*=K124;
      complex<double> tanD=term1 - term2 - 1.0 + sqrt(sum*sum - 4.0*k12*k22*K124);
      result=2.*k1*(_K(0,0)+k22*_K(1,1)*_detK)/tanD;
    }// end if above KKbar threshold
  }// end if sqrt(s)< 932MeV

  return result.real();
}


double
s0wave::eta(std::complex<double> s) const 
{
  //Build K-Matrix
    rmatrix _K(2,2);
    for(unsigned int i=0;i<2;++i)
      for(unsigned int j=0;j<2;++j){
	_K(i,j)=_mu*(_alpha[i]*_alpha[j]/(_M12-s.real()) + _beta[i]*_beta[j]/(_M22-s.real())) +_invmu*_gamma(i,j);
      }

  complex<double> k1=kThr(s,mpic2);
  complex<double> k2=kThr(s,mKc2);
  double detK=_K(0,0)*_K(1,1)-_K(1,0)*_K(0,1);
  complex<double> result=1;
  if(s.real()>0.987354){
    complex<double> term1=1.+k1*k2*detK;term1*=term1;
    complex<double> term2=1.-k1*k2*detK;term2*=term2;
    complex<double> term3=k1*_K(0,0)-k2*_K(1,1);term3*=term3;
    complex<double> term4=k1*_K(0,0)+k2*_K(1,1);term4*=term4;
    result=sqrt( (term1+term3)/(term2+term4) );
  }
  return result.real();
}
