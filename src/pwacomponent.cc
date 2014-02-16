//-----------------------------------------------------------
//
// Description:
//      Implementation of class pwacomponent
//      see pwacomponent.h for details
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

// This Class' Header ------------------
#include "pwacomponent.h"

// C/C++ Headers ----------------------
#include <algorithm>
#include <iostream>

// Collaborating Class Headers --------
#include "TF1.h"

#include "reportingUtils.hpp"

// Class Member definitions -----------

using namespace std;



rpwa::pwachannel::pwachannel(const rpwa::pwachannel& ch) /// cp ctor
  : _C(ch._C),_ps(ch._ps){
}




rpwa::pwacomponent::pwacomponent(const string& name,
				 double m0, double gamma,
				 const map<string,pwachannel >& channels)
  : _name(name), _m0(m0), _m02(m0*m0),_m0min(0),_m0max(5000),_gamma(gamma),_gammamin(0),_gammamax(1000),_fixm(false),_fixgamma(false), _constWidth(false),_channels(channels)
{
  // fill vector with channel for random access
  map<string,pwachannel >::iterator it=_channels.begin();
  while(it!=_channels.end()){
    _vchannels.push_back(&(it->second)); // carefull this is dangerous
    // refactor to vector only use as soon as possible!
    _channelname.push_back(it->first);
    ++it;
  }
}




void
rpwa::pwacomponent::setCouplings(const double* par){
  map<string,pwachannel >::iterator it=_channels.begin();
  unsigned int counter=0;
  while(it!=_channels.end()){
    it->second.setCoupling(complex<double>(par[counter],par[counter+1]));
    counter+=2;
    ++it;
  }
}
void
rpwa::pwacomponent::getCouplings(double* par){
  map<string,pwachannel >::iterator it=_channels.begin();
  unsigned int counter=0;
  while(it!=_channels.end()){
    par[counter]=it->second.C().real();
    par[counter+1]=it->second.C().imag();
    counter+=2;
    ++it;
  }
}


complex<double>
rpwa::pwacomponent::val(double m) const {
  // calculate dynamic width:
  // loop over decay channels
  double gamma=0;
  std::map<std::string,pwachannel >::const_iterator it=_channels.begin();
  double n=(double)numChannels();
  //cout << "NChannels : " << n << endl;
  double ps=1;
  if(!_constWidth){
    ps=0; // no not forget to reset phase space here!
    while(it!=_channels.end()){
      double myps=1.;
      if(it->second.ps()!=NULL){
	double ps0=it->second.ps(_m0);
	myps=(it->second.ps(m))/ps0;
      }
      ps+=myps;
      ++it;
    }
    ps/=n;
  }
  gamma=_gamma*ps;

  //std::cerr << m << "   " << gamma/_gamma << std::endl;
  //std::cerr << _name <<"  compval=" <<gamma*_m0/complex<double>(m*m-_m02,gamma*_m0) << std::endl;
  return _gamma*_m0/complex<double>(_m02-m*m,-gamma*_m0);
}


double
lambda(const double a,
       const double b,
       const double c)
{
  return a * a + b * b + c * c - 2.0 * (a * b + b * c + c * a);
}


complex<double>
q(const double M,
  const double m1,
  const double m2)
{
  double lam = lambda(M * M, m1 * m1, m2 * m2);
  complex<double> ret;
  if (lam < 0)
    return complex<double>(0.0, 0.0);
  return complex<double>(sqrt(lam / (4 * M * M)), 0.0 );
}

complex<double>
rpwa::pwabkg::val(double m) const {
  m-=_m0; // shift baseline mass
  // calculate breakup momentum
  complex<double> p;
  if(m<_m1+_m2)return complex<double>(1,0);
  p=q(m,_m1,_m2);
  //std::cerr << _name <<"  val=" << exp(-_gamma*p) << std::endl;
  return exp(-_gamma*p.real()*p.real());
}


std::ostream& rpwa::operator<< (std::ostream& o,const rpwa::pwacomponent& c){
  o << c.name() << endl
    << "Mass= " << c.m0() << "   Width="<< c.gamma() << "     ConstWidth="<< c. constWidth() << endl
    << "Decay modes: " << endl;
  std::map<std::string,pwachannel >::const_iterator it=c.channels().begin();
  while(it!=c.channels().end()){
    o << it->first << "    C=" << it->second.C() << endl;
      ++it;
  }

  return o;
}
