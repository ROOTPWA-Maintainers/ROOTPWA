//-----------------------------------------------------------
// File and Version Information:
// $Id$
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


// Collaborating Class Headers --------
#include <algorithm>
#include <iostream>

// Class Member definitions -----------

using namespace std;

rpwa::pwacomponent::pwacomponent(const string& name,
				 double m0, double gamma,
				 const map<string,pwachannel >& channels)
  : _name(name), _m0(m0), _m02(m0*m0),_m0min(0),_m0max(5000),_gamma(gamma),_gammamin(0),_gammamax(1000),_fixm(false),_fixgamma(false), _channels(channels)
{}



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
  while(it!=_channels.end()){
    double ps=1;
    if(it->second.ps()!=NULL){
      double ps0=it->second.ps(_m0);
      ps=it->second.ps(m)/ps0;
    }
    gamma+=_gamma*ps/n;
    ++it;
  }

  //std::cerr << _name <<"  compval=" <<gamma*_m0/complex<double>(m*m-_m02,gamma*_m0) << std::endl;
  return gamma*_m0/complex<double>(m*m-_m02,gamma*_m0);
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
    return complex<double>(0.0, sqrt(fabs(lam / (4 * M * M))));
  return complex<double>(sqrt(lam / (4 * M * M)), 0.0 );
}

complex<double>
rpwa::pwabkg::val(double m) const {
  m-=_m0; // shift baseline mass
  // calculate breakup momentum
  complex<double> p=q(m,_m1,_m2);
  //std::cerr << _name <<"  val=" << exp(-_gamma*p) << std::endl;
  return exp(-_gamma*p);
}


vector<string> 
rpwa::pwacompset::wavelist() const {
  vector<string> wl;
  for(unsigned int i=0;i<n();++i){
    const map<string,pwachannel >& channellist=_comp[i]->channels();
    map<string,pwachannel >::const_iterator it=channellist.begin();
    while(it!=channellist.end()){
      if(find(wl.begin(),wl.end(),it->first)==wl.end())
	 wl.push_back(it->first);
      ++it;
    }
  }
  return wl;
}


void
rpwa::pwacompset::setPar(const double* par){ // set parameters
  unsigned int parcount=0;
  for(unsigned int i=0;i<n();++i){
    _comp[i]->setPar(par[parcount],par[parcount+1]);
    parcount+=2;
    _comp[i]->setCouplings(&par[parcount]);
    parcount+=_comp[i]->numChannels()*2; // RE and Im for each channel
  }
}


void 
rpwa::pwacompset::getPar(double* par){       // return parameters 
  unsigned int parcount=0;
  for(unsigned int i=0;i<n();++i){
    par[parcount]=_comp[i]->m0();
    par[parcount+1]=_comp[i]->gamma();
    parcount+=2;
    _comp[i]->getCouplings(&par[parcount]);
    parcount+=_comp[i]->numChannels()*2; // RE and Im for each channel
  }
}


double 
rpwa::pwacompset::intensity(const std::string& wave, double m){
  // loop over all components and pick up those that contribute to this channel
  complex<double> rho(0,0);
  for(unsigned int ic=0;ic<n();++ic){
    if(_comp[ic]->channels().count(wave)==0)continue;
    else {
      rho+=_comp[ic]->val(m)*_comp[ic]->channels().find(wave)->second.C();
    }

  }
  return norm(rho);
}

double 
rpwa::pwacompset::phase(const std::string& wave1,
			double ps1,
			const std::string& wave2,
			double ps2,
			double m){
  // loop over all components and pick up those that contribute to this channel
  complex<double> rho1(0,0);
  complex<double> rho2(0,0);

  for(unsigned int ic=0;ic<n();++ic){
    if(_comp[ic]->channels().count(wave1)!=0){
      rho1+=_comp[ic]->val(m)*_comp[ic]->channels().find(wave1)->second.C();
    }
    if(_comp[ic]->channels().count(wave2)!=0){
      rho2=_comp[ic]->val(m)*_comp[ic]->channels().find(wave2)->second.C();
    }
  }
  rho1*=ps1;
  rho2*=ps2;
  return arg(rho1*conj(rho2));
}


std::ostream& rpwa::operator<< (std::ostream& o,const rpwa::pwacomponent& c){
  o << c.name() << endl
    << "Mass= " << c.m0() << "   Width="<< c.gamma() << endl
    << "Decay modes: " << endl;
  std::map<std::string,pwachannel >::const_iterator it=c.channels().begin();
  while(it!=c.channels().end()){
    o << it->first << "    C=" << it->second.C() << endl;
      ++it;
  }

  return o;
}

std::ostream& rpwa::operator<< (std::ostream& out,const rpwa::pwacompset& cs){
  for(unsigned int i=0;i<cs.n();++i){
    const rpwa::pwacomponent& c =*cs[i];
    out << c << endl;
  }
  return out;
}
