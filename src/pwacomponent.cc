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
#include <algorithm>
#include <iostream>

// Collaborating Class Headers --------
#include "TF1.h"

// Class Member definitions -----------

using namespace std;

rpwa::pwacomponent::pwacomponent(const string& name,
				 double m0, double gamma,
				 const map<string,pwachannel >& channels)
  : _name(name), _m0(m0), _m02(m0*m0),_m0min(0),_m0max(5000),_gamma(gamma),_gammamin(0),_gammamax(1000),_fixm(false),_fixgamma(false), _constWidth(false),_channels(channels)
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
rpwa::pwacompset::setPS(TF1* fPS){
  _phasespace=fPS;
  // check if there are free parameters in the phase space that should be fitted
  unsigned int nparPS=_phasespace->GetNpar();
  // loop over parameters and check limits
  // remember which parameters to let float
  _freePSpar.clear();
  for(unsigned int i=0;i<nparPS;++i){
    double min,max;
    _phasespace->GetParLimits(i,min,max);
    if(min!=max){
      _freePSpar.push_back(i);
      cout << "PS parameter "<< i << " floating in ["
	   << min  << "," << max << "]" << endl;
    }
  }// end loop over parameters
  _numpar+=_freePSpar.size();
}

double 
rpwa::pwacompset::getFreePSPar(unsigned int i){
  if(i<_freePSpar.size())
    return _phasespace->GetParameter(_freePSpar[i]);
  else return 0;
}


void 
rpwa::pwacompset::getFreePSLimits(unsigned int i, double& lower, double& upper){
  if(i<_freePSpar.size()){
    _phasespace->GetParLimits(_freePSpar[i],lower,upper);
  }
}

void
rpwa::pwacompset::setPar(const double* par){ // set parameters
  unsigned int parcount=0;
  // components
  for(unsigned int i=0;i<n();++i){
    _comp[i]->setPar(par[parcount],par[parcount+1]);
    parcount+=2;
    _comp[i]->setCouplings(&par[parcount]);
    parcount+=_comp[i]->numChannels()*2; // RE and Im for each channel
  } // end loop over components
  // phase space
  unsigned int nfreepar=_freePSpar.size();
  for(unsigned int ipar=0;ipar<nfreepar;++ipar){
    _phasespace->SetParameter(_freePSpar[ipar],par[parcount]);
    ++parcount;
  }
}


void 
rpwa::pwacompset::getPar(double* par){       // return parameters 
  unsigned int parcount=0;
  // components
  for(unsigned int i=0;i<n();++i){
    par[parcount]=_comp[i]->m0();
    par[parcount+1]=_comp[i]->gamma();
    parcount+=2;
    _comp[i]->getCouplings(&par[parcount]);
    parcount+=_comp[i]->numChannels()*2; // RE and Im for each channel
  }
 // phase space
  unsigned int nfreepar=_freePSpar.size();
  for(unsigned int ipar=0;ipar<nfreepar;++ipar){
    par[parcount]=_phasespace->GetParameter(_freePSpar[ipar]);
    ++parcount;
  }
}

double 
rpwa::pwacompset::ps(double m){return _phasespace->Eval(m);}

double 
rpwa::pwacompset::intensity(const std::string& wave, double m){
  // loop over all components and pick up those that contribute to this channel
  complex<double> rho(0,0);
  for(unsigned int ic=0;ic<n();++ic){
    if(_comp[ic]->channels().count(wave)==0)continue;
    else {
      rho+=_comp[ic]->val(m)*_comp[ic]->channels().find(wave)->second.C()*sqrt(_comp[ic]->channels().find(wave)->second.ps(m));
    }

  }
  return norm(rho)*_phasespace->Eval(m);
}

double 
rpwa::pwacompset::phase(const std::string& wave, double m){
  // loop over all components and pick up those that contribute to this channel
  complex<double> rho(0,0);
  for(unsigned int ic=0;ic<n();++ic){
    if(_comp[ic]->channels().count(wave)==0)continue;
    else {
      rho+=_comp[ic]->val(m)*_comp[ic]->channels().find(wave)->second.C();
    }

  }
  return arg(rho);
}


double 
rpwa::pwacompset::phase(const std::string& wave1,
			const std::string& wave2,
			double m){
  return arg(overlap(wave1,wave2,m));
}

std::complex<double>
rpwa::pwacompset::overlap(const std::string& wave1,
			  const std::string& wave2,
			  double m){
    // loop over all components and pick up those that contribute to this channel
  complex<double> rho1(0,0);
  complex<double> rho2(0,0);

  for(unsigned int ic=0;ic<n();++ic){
    if(_comp[ic]->channels().count(wave1)!=0){
      rho1+=_comp[ic]->val(m)*_comp[ic]->channels().find(wave1)->second.C()*sqrt(_comp[ic]->channels().find(wave1)->second.ps(m));
    }
    if(_comp[ic]->channels().count(wave2)!=0){
      rho2+=_comp[ic]->val(m)*_comp[ic]->channels().find(wave2)->second.C()*sqrt(_comp[ic]->channels().find(wave2)->second.ps(m));
    }
  }
  return rho1*conj(rho2)*_phasespace->Eval(m);
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
