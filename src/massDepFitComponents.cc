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
#include "massDepFitComponents.h"

// C/C++ Headers ----------------------
#include <algorithm>
#include <iostream>

// Collaborating Class Headers --------
#include "TF1.h"

#include "reportingUtils.hpp"

// Class Member definitions -----------

using namespace std;
using namespace rpwa;


/////////////////////////////
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

rpwa::pwachannel::pwachannel(const rpwa::pwachannel& ch) /// cp ctor
  : _waveName(ch._waveName), _C(ch._C), _phaseSpace(ch._phaseSpace), _ps(ch._ps) {
}
/////////////////////////////


massDepFitComponent::massDepFitComponent(const string& name,
                                         const size_t nrParameters,
                                         const vector<pwachannel>& channels)
	: _name(name),
	  _channels(channels),
	  _nrParameters(nrParameters)
{
}


void
massDepFitComponent::getCouplings(double* par) const
{
	size_t counter=0;
	for(vector<pwachannel>::const_iterator it=_channels.begin(); it!=_channels.end(); ++it, counter+=2) {
		par[counter] = it->C().real();
		par[counter+1] = it->C().imag();
	}
}


void
massDepFitComponent::setCouplings(const double* par)
{
	size_t counter=0;
	for(vector<pwachannel>::iterator it=_channels.begin(); it!=_channels.end(); ++it, counter+=2) {
		it->setCoupling(complex<double>(par[counter], par[counter+1]));
	}
}


ostream&
massDepFitComponent::print(std::ostream& out) const
{
	out << "Decay modes:" << endl;
	for(vector<pwachannel>::const_iterator it=_channels.begin(); it!=_channels.end(); ++it) {
		out << it->getWaveName() << "    C=" << it->C() << endl;
	}
	return out;
}


rpwa::pwacomponent::pwacomponent(const string& name,
                                 const vector<pwachannel >& channels,
                                 const double m0,
                                 const double gamma,
                                 const bool constWidth)
	: massDepFitComponent(name, 2, channels),
	  _m0(m0),
	  _m02(m0*m0),
	  _m0min(0),
	  _m0max(5000),
	  _fixm(false),
	  _gamma(gamma),
	  _gammamin(0),
	  _gammamax(1000),
	  _fixgamma(false), 
	  _constWidth(constWidth)
{
}

    
void
pwacomponent::getParameters(double* par) const
{
	par[0]=_m0;
	par[1]=_gamma;
}


void
pwacomponent::setParameters(const double* par)
{
	_m0=par[0];
	_m02=par[0] * par[0];
	_gamma=par[1];
}


double
pwacomponent::getParameter(const size_t idx) const
{
	if(idx == 0) return _m0;
	if(idx == 1) return _gamma;
	throw;
}


void
pwacomponent::setParameter(const size_t idx, const double value)
{
	if(idx == 0) { _m0 = value; return; }
	if(idx == 1) { _gamma = value; return; }
	throw;
}


bool
pwacomponent::getParameterFixed(const size_t idx) const
{
	if(idx == 0) return _fixm;
	if(idx == 1) return _fixgamma;
	throw;
}


void
pwacomponent::setParameterFixed(const size_t idx, const bool fixed)
{
	if(idx == 0) { _fixm = fixed; return; }
	if(idx == 1) { _fixgamma = fixed; return; }
	throw;
}


pair<double, double>
pwacomponent::getParameterLimits(const size_t idx) const
{
	if(idx == 0) return make_pair(_m0min, _m0max);
	if(idx == 1) return make_pair(_gammamin, _gammamax);
	throw;
}


void
pwacomponent::setParameterLimits(const size_t idx, const pair<double, double>& limits)
{
	if(idx == 0) { _m0min = limits.first; _m0max = limits.second; return; }
	if(idx == 1) { _gammamin = limits.first; _gammamax = limits.second; return; }
	throw;
}


complex<double>
rpwa::pwacomponent::val(const double m) const {
	// calculate dynamic width:
	// loop over decay channels
	double ps = 1.;
	if(!_constWidth){
		ps=0; // no not forget to reset phase space here!
		for(std::vector<pwachannel>::const_iterator itChan=getChannels().begin(); itChan!=getChannels().end(); ++itChan) {
			double myps=1.;
			if(itChan->ps() != NULL){
				double ps0=itChan->ps(_m0);
				myps=(itChan->ps(m))/ps0;
			}
			ps+=myps;
		}
		ps /= getNrChannels();
	}
	const double gamma = _gamma*ps;
  
	return _gamma*_m0 / complex<double>(_m02-m*m, -gamma*_m0);
}


std::ostream&
pwacomponent::print(std::ostream& out) const
{
	out << getName() << endl
	    << "Mass=" << _m0 << "   Width=" << _gamma << "    ConstWidth=" << _constWidth << endl;
	return massDepFitComponent::print(out);
}


pwabkg::pwabkg(const string& name,
               const vector<pwachannel >& channels,
               const double m0,
               const double gamma,
               const double m1,
               const double m2)
	: massDepFitComponent(name, 2, channels),
	  _m0(m0),
	  _m0min(0),
	  _m0max(5000),
	  _fixm(false),
	  _gamma(gamma),
	  _gammamin(0),
	  _gammamax(1000),
	  _fixgamma(false),
	  _m1(m1),
	  _m2(m2)
{
}


void
pwabkg::getParameters(double* par) const
{
	par[0]=_m0;
	par[1]=_gamma;
}


void
pwabkg::setParameters(const double* par)
{
	_m0=par[0];
	_gamma=par[1];
}


double
pwabkg::getParameter(const size_t idx) const
{
	if(idx == 0) return _m0;
	if(idx == 1) return _gamma;
	throw;
}


void
pwabkg::setParameter(const size_t idx, const double value)
{
	if(idx == 0) { _m0 = value; return; }
	if(idx == 1) { _gamma = value; return; }
	throw;
}


bool
pwabkg::getParameterFixed(const size_t idx) const
{
	if(idx == 0) return _fixm;
	if(idx == 1) return _fixgamma;
	throw;
}


void
pwabkg::setParameterFixed(const size_t idx, const bool fixed)
{
	if(idx == 0) { _fixm = fixed; return; }
	if(idx == 1) { _fixgamma = fixed; return; }
	throw;
}


pair<double, double>
pwabkg::getParameterLimits(const size_t idx) const
{
	if(idx == 0) return make_pair(_m0min, _m0max);
	if(idx == 1) return make_pair(_gammamin, _gammamax);
	throw;
}


void
pwabkg::setParameterLimits(const size_t idx, const pair<double, double>& limits)
{
	if(idx == 0) { _m0min = limits.first; _m0max = limits.second; return; }
	if(idx == 1) { _gammamin = limits.first; _gammamax = limits.second; return; }
	throw;
}


complex<double>
pwabkg::val(const double m) const {
	// shift baseline mass
	const double mass = m-_m0;

	// calculate breakup momentum
	if(mass < _m1+_m2) {
		return complex<double>(1,0);
	}
	const complex<double> p = q(mass,_m1,_m2);

	return exp(-_gamma*p.real()*p.real());
}
