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


// Class Member definitions -----------

using namespace std;

rpwa::pwacomponent::pwacomponent(const string& name,
				 double m0, double gamma,
				 const map<string,complex<double> >& channels)
  : _name(name), _m0(m0), _m02(m0*m0),_m0min(0),_m0max(5000),_gamma(gamma),_gammamin(0),_gammamax(1000),_fixm(false),_fixgamma(false), _channels(channels)
{}



void
rpwa::pwacomponent::setCouplings(const double* par){
  map<string,complex<double> >::iterator it=_channels.begin();
  unsigned int counter=0;
  while(it!=_channels.end()){
    it->second=complex<double>(par[counter],par[counter+1]);
    counter+=2;
    ++it;
  }
}
void
rpwa::pwacomponent::getCouplings(double* par){
  map<string,complex<double> >::iterator it=_channels.begin();
  unsigned int counter=0;
  while(it!=_channels.end()){
    par[counter]=it->second.real();
    par[counter+1]=it->second.imag();
    counter+=2;
    ++it;
  }
}

    
complex<double>
rpwa::pwacomponent::val(double m) const {
  return _gamma*_m0/complex<double>(m*m-_m02,_gamma*_m0);
}





vector<string> 
rpwa::pwacompset::wavelist() const {
  vector<string> wl;
  for(unsigned int i=0;i<n();++i){
    const map<string,complex<double> >& channellist=_comp[i].channels();
    map<string,complex<double> >::const_iterator it=channellist.begin();
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
    _comp[i].setPar(par[parcount],par[parcount+1]);
    parcount+=2;
    _comp[i].setCouplings(&par[parcount]);
    parcount+=_comp[i].numChannels()*2; // RE and Im for each channel
  }
}


void 
rpwa::pwacompset::getPar(double* par){       // return parameters 
  unsigned int parcount=0;
  for(unsigned int i=0;i<n();++i){
    par[parcount]=_comp[i].m0();
    par[parcount+1]=_comp[i].gamma();
    parcount+=2;
    _comp[i].getCouplings(&par[parcount]);
    parcount+=_comp[i].numChannels()*2; // RE and Im for each channel
  }
}


double 
rpwa::pwacompset::intensity(std::string wave, double m){
  // loop over all components and pick up those that contribute to this channel
  complex<double> rho(0,0);
  for(unsigned int ic=0;ic<n();++ic){
    if(_comp[ic].channels().count(wave)==0)continue;
    else {
      rho+=_comp[ic].val(m)*_comp[ic].channels().find(wave)->second;
    }

  }
  return norm(rho);
}



std::ostream& rpwa::operator<< (std::ostream& o,const rpwa::pwacomponent& c){
  o << c.name() << endl
    << "Mass= " << c.m0() << "   Width="<< c.gamma() << endl
    << "Decay modes: " << endl;
  std::map<std::string,std::complex<double> >::const_iterator it=c.channels().begin();
  while(it!=c.channels().end()){
    o << it->first << "    C=" << it->second << endl;
      ++it;
  }

  return o;
}

std::ostream& rpwa::operator<< (std::ostream& out,const rpwa::pwacompset& cs){
  for(unsigned int i=0;i<cs.n();++i){
    const rpwa::pwacomponent& c =cs[i];
    out << c << endl;
  }
  return out;
}
