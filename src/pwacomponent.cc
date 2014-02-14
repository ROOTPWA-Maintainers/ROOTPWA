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


void
rpwa::pwacompset::setFuncFsmd(TF1* funcFsmd){
  _funcFsmd=funcFsmd;
  // clear list of free parameters
  _freeFsmdPar.clear();

  if (_funcFsmd == NULL)
    return;

  // check if there are free parameters in the phase space that should be fitted
  unsigned int nparFsmd=_funcFsmd->GetNpar();
  // loop over parameters and check limits
  // remember which parameters to let float
  for(unsigned int i=0;i<nparFsmd;++i){
    double min,max;
    _funcFsmd->GetParLimits(i,min,max);
    if(min!=max){
      _freeFsmdPar.push_back(i);
      cout << "final-state mass dependence parameter "<< i << " floating in ["
	   << min  << "," << max << "]" << endl;
    }
  }// end loop over parameters
  _numpar+=_freeFsmdPar.size();
}


// performs mapping from the index of a wave in wavelist() to the components and channels that couple to this wave
bool
rpwa::pwacompset::doMapping() {
	// check that all waves used in a decay channel have been defined
	for(unsigned int i=0; i<n(); ++i) {
		for(map<string, pwachannel>::const_iterator itChan=_comp[i]->channels().begin(); itChan !=_comp[i]->channels().end(); ++itChan) {
			if(find(_waveList.begin(), _waveList.end(), itChan->first) == _waveList.end()) {
				printErr << "wave '" << itChan->first << "' not known in decay of '" << _comp[i]->name() << "'." << endl;
				return false;
			}
		}
	}

	// check that all defined waves are also used
	for(vector<string>::const_iterator itWave=_waveList.begin(); itWave!=_waveList.end(); ++itWave) {
		bool found(false);
		for(unsigned int i=0; i<n(); ++i) {
			for(map<string, pwachannel>::const_iterator itChan=_comp[i]->channels().begin(); itChan !=_comp[i]->channels().end(); ++itChan) {
				if(itChan->first == *itWave) {
					found = true;
					break;
				}
			}
		}

		if(not found) {
			printErr << "wave '" << *itWave << "' defined but not used in any decay." << endl;
			return false;
		}
	}

	for(vector<string>::const_iterator itWave=_waveList.begin(); itWave!=_waveList.end(); ++itWave) {
		// check which components this waves belongs to and which channel they are
		_compChannel.push_back(getCompChannel(*itWave));
	}

	return true;
}


double 
rpwa::pwacompset::getFreeFsmdPar(unsigned int i) const {
  if(i<_freeFsmdPar.size())
    return _funcFsmd->GetParameter(_freeFsmdPar[i]);
  else return 0;
}


void 
rpwa::pwacompset::getFreeFsmdLimits(unsigned int i, double& lower, double& upper) const {
  if(i<_freeFsmdPar.size()){
    _funcFsmd->GetParLimits(_freeFsmdPar[i],lower,upper);
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
  unsigned int nfreepar=_freeFsmdPar.size();
  for(unsigned int ipar=0;ipar<nfreepar;++ipar){
    _funcFsmd->SetParameter(_freeFsmdPar[ipar],par[parcount]);
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
  unsigned int nfreepar=_freeFsmdPar.size();
  for(unsigned int ipar=0;ipar<nfreepar;++ipar){
    par[parcount]=_funcFsmd->GetParameter(_freeFsmdPar[ipar]);
    ++parcount;
  }
}

double 
rpwa::pwacompset::calcFsmd(double m){
  if (_funcFsmd == NULL) {
    return 1.;
  }
  
  return _funcFsmd->Eval(m);
}

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
  return norm(rho*calcFsmd(m));
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
  return (rho1*calcFsmd(m))*conj(rho2*calcFsmd(m));
}

std::complex<double>
rpwa::pwacompset::overlap(unsigned int wave1,
			  unsigned int wave2,
			  double m,
                          const size_t idxMass){
    // loop over all components and pick up those that contribute to this channels
  complex<double> rho1(0,0);
  complex<double> rho2(0,0);

  // get entry from mapping
  vector<pair<unsigned int, unsigned int> > cc1=_compChannel[wave1];
  vector<pair<unsigned int, unsigned int> > cc2=_compChannel[wave2];
  unsigned int nc1=cc1.size();
  unsigned int nc2=cc2.size();

  for(unsigned int ic1=0;ic1<nc1;++ic1){
    unsigned int icomp=cc1[ic1].first;
    unsigned int ichan=cc1[ic1].second;
    rho1+=_comp[icomp]->val(m)*_comp[icomp]->getChannel(ichan).CsqrtPS(m);
  }
 for(unsigned int ic2=0;ic2<nc2;++ic2){
    unsigned int icomp=cc2[ic2].first;
    unsigned int ichan=cc2[ic2].second;
    rho2+=_comp[icomp]->val(m)*_comp[icomp]->getChannel(ichan).CsqrtPS(m);
  }
  return (rho1*calcFsmd(m))*conj(rho2*calcFsmd(m));
}


std::vector<std::pair<unsigned int,unsigned int> >
rpwa::pwacompset::getCompChannel(const std::string& wave) const {
  cerr << "Channel-mapping for wave " << wave << endl;
  std::vector<std::pair<unsigned int,unsigned int> > result;
  for(unsigned int ic=0;ic<n();++ic){
    // loop over channels of component and see if wave is there
    unsigned int nch=_comp[ic]->numChannels();
    for(unsigned int ich=0;ich<nch;++ich){
      if(_comp[ic]->getChannelName(ich)==wave){
	result.push_back(pair<unsigned int,unsigned int>(ic,ich));
	cerr << "     comp("<<ic<<"): " << _comp[ic]->name() << "  ch:"<<ich<<endl;
      }
    } // end loop over channels
  } // end loop over components
  return result;
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






std::ostream& rpwa::operator<< (std::ostream& out,const rpwa::pwacompset& cs){
  for(unsigned int i=0;i<cs.n();++i){
    const rpwa::pwacomponent& c =*cs[i];
    out << c << endl;
  }
  return out;
}
