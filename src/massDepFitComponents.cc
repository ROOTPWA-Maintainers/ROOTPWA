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

#include <boost/assign/std/vector.hpp>

// Collaborating Class Headers --------
#include "TF1.h"

#include "libConfigUtils.hpp"
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

pwachannel::pwachannel(const pwachannel& ch) /// cp ctor
  : _waveName(ch._waveName), _C(ch._C), _phaseSpace(ch._phaseSpace), _ps(ch._ps) {
}
/////////////////////////////


massDepFitComponent::massDepFitComponent(const string& name)
	: _name(name)
{
}


bool
massDepFitComponent::init(const libconfig::Setting* configComponent,
                          const vector<double>& massBinCenters,
                          const map<string, size_t>& waveIndices,
                          const boost::multi_array<double, 2>& phaseSpaceIntegrals,
                          const bool debug) {
	if(debug) {
		printDebug << "starting initialization of 'massDepFitComponent' for component '" << getName() << "'." << endl;
	}

	const libconfig::Setting* decayChannels = findLibConfigList(*configComponent, "decaychannels");
	if(not decayChannels) {
		printErr << "component '" << getName() << "' has no decay channels." << endl;
		return false;
	}

	const int nrDecayChannels = decayChannels->getLength();
	for(int idxDecayChannel=0; idxDecayChannel<nrDecayChannels; ++idxDecayChannel) {
		const libconfig::Setting* decayChannel = &((*decayChannels)[idxDecayChannel]);

		map<string, libconfig::Setting::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("amp", libconfig::Setting::TypeString)
		                     ("coupling_Re", libconfig::Setting::TypeFloat)
		                     ("coupling_Im", libconfig::Setting::TypeFloat);
		if(not checkIfAllVariablesAreThere(decayChannel, mandatoryArguments)) {
			printErr << "one of the decay channels of the component '" << getName() << "' does not contain all required fields." << endl;
			return false;
		}

		string waveName;
		decayChannel->lookupValue("amp", waveName);

		double couplingReal;
		decayChannel->lookupValue("coupling_Re", couplingReal);

		double couplingImag;
		decayChannel->lookupValue("coupling_Im", couplingImag);

		const std::complex<double> coupling(couplingReal, couplingImag);

		const map<string, size_t>::const_iterator it = waveIndices.find(waveName);
		if(it == waveIndices.end()) {
			printErr << "wave '" << waveName << "' not in fit, but used as decay channel." << endl;
			return false;
		}

		boost::multi_array<double, 2>::const_array_view<1>::type view = phaseSpaceIntegrals[boost::indices[boost::multi_array<double, 2>::index_range()][it->second]];
		_channels.push_back(pwachannel(waveName, coupling, massBinCenters, std::vector<double>(view.begin(), view.end())));

		if(debug) {
			printDebug << "coupling to wave '" << waveName << "' (index " << it->second << ") with coupling " << coupling << "." << endl;
		}
	}

	if(debug) {
		printDebug << "finished initialization of 'massDepFitComponent'." << endl;
	}

	return true;
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
massDepFitComponent::print(ostream& out) const
{
	out << "Decay modes:" << endl;
	for(vector<pwachannel>::const_iterator it=_channels.begin(); it!=_channels.end(); ++it) {
		out << it->getWaveName() << "    C=" << it->C() << endl;
	}
	return out;
}


pwacomponent::pwacomponent(const string& name)
	: massDepFitComponent(name)
{
}


bool
pwacomponent::init(const libconfig::Setting* configComponent,
                   const vector<double>& massBinCenters,
                   const map<string, size_t>& waveIndices,
                   const boost::multi_array<double, 2>& phaseSpaceIntegrals,
                   const bool debug) {
	if(debug) {
		printDebug << "starting initialization of 'pwacomponent' for component '" << getName() << "'." << endl;
	}

	setNrParameters(2);
	if(not massDepFitComponent::init(configComponent, massBinCenters, waveIndices, phaseSpaceIntegrals, debug)) {
		printErr << "error while reading configuration of 'massDepFitComponent' class." << endl;
		return false;
	}

	map<string, libconfig::Setting::Type> mandatoryArguments;
	boost::assign::insert(mandatoryArguments)
	                     ("val", libconfig::Setting::TypeFloat)
	                     ("lower", libconfig::Setting::TypeFloat)
	                     ("upper", libconfig::Setting::TypeFloat)
	                     ("fix", libconfig::Setting::TypeBoolean);

	const libconfig::Setting* configMass = findLibConfigGroup(*configComponent, "mass");
	if(not configMass) {
		printErr << "component '" << getName() << "' has no section 'mass'." << endl;
		return false;
	}

	if(not checkIfAllVariablesAreThere(configMass, mandatoryArguments)) {
		printErr << "the 'mass' section of the component '" << getName() << "' does not contain all required fields." << endl;
		return false;
	}

	configMass->lookupValue("val", _m0);
	configMass->lookupValue("lower", _m0min);
	configMass->lookupValue("upper", _m0max);
	configMass->lookupValue("fix", _fixm);
	_m02 = _m0*_m0;

	const libconfig::Setting* configWidth = findLibConfigGroup(*configComponent, "width");
	if(not configWidth) {
		printErr << "component '" << getName() << "' has no section 'width'." << endl;
		return false;
	}

	if(not checkIfAllVariablesAreThere(configWidth, mandatoryArguments)) {
		printErr << "the 'width' section of the component '" << getName() << "' does not contain all required fields." << endl;
		return false;
	}

	configWidth->lookupValue("val", _gamma);
	configWidth->lookupValue("lower", _gammamin);
	configWidth->lookupValue("upper", _gammamax);
	configWidth->lookupValue("fix", _fixgamma);

	if(not configWidth->lookupValue("dyn", _constWidth)) {
		_constWidth = false;
	}
	_constWidth = !_constWidth;

	printInfo << "component '" << getName() << "':" << endl
	          << "mass: " << _m0 << " MeV/c^2, limits: " << _m0min << "-" << _m0max << " MeV/c^2 " << (_fixm ? "(FIXED)" : "") << endl
	          << "width: " << _gamma << " MeV/c^2, limits: " << _gammamin << "-" << _gammamax << " MeV/c^2 " << (_fixgamma ? "(FIXED) " : "") << (_constWidth ? "-- constant width" : "-- dynamic width") << endl;

	if(debug) {
		printDebug << "finished initialization of 'pwacomponent'." << endl;
	}
	return true;
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
pwacomponent::val(const double m) const {
	// calculate dynamic width:
	// loop over decay channels
	double ps = 1.;
	if(!_constWidth){
		ps=0; // no not forget to reset phase space here!
		for(vector<pwachannel>::const_iterator itChan=getChannels().begin(); itChan!=getChannels().end(); ++itChan) {
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


ostream&
pwacomponent::print(ostream& out) const
{
	out << getName() << endl
	    << "Mass=" << _m0 << "   Width=" << _gamma << "    ConstWidth=" << _constWidth << endl;
	return massDepFitComponent::print(out);
}


pwabkg::pwabkg(const string& name)
	: massDepFitComponent(name)
{
}


bool
pwabkg::init(const libconfig::Setting* configComponent,
             const vector<double>& massBinCenters,
             const map<string, size_t>& waveIndices,
             const boost::multi_array<double, 2>& phaseSpaceIntegrals,
             const bool debug) {
	if(debug) {
		printDebug << "starting initialization of 'pwabkg' for component '" << getName() << "'." << endl;
	}

	setNrParameters(2);
	if(not massDepFitComponent::init(configComponent, massBinCenters, waveIndices, phaseSpaceIntegrals, debug)) {
		printErr << "error while reading configuration of 'massDepFitComponent' class." << endl;
		return false;
	}

	map<string, libconfig::Setting::Type> mandatoryArguments;
	boost::assign::insert(mandatoryArguments)
	                     ("val", libconfig::Setting::TypeFloat)
	                     ("lower", libconfig::Setting::TypeFloat)
	                     ("upper", libconfig::Setting::TypeFloat)
	                     ("fix", libconfig::Setting::TypeBoolean);

	const libconfig::Setting* configMass = findLibConfigGroup(*configComponent, "m0");
	if(not configMass) {
		printErr << "component '" << getName() << "' has no section 'mass'." << endl;
		return false;
	}

	if(not checkIfAllVariablesAreThere(configMass, mandatoryArguments)) {
		printErr << "the 'mass' section of the component '" << getName() << "' does not contain all required fields." << endl;
		return false;
	}

	configMass->lookupValue("val", _m0);
	configMass->lookupValue("lower", _m0min);
	configMass->lookupValue("upper", _m0max);
	configMass->lookupValue("fix", _fixm);

	const libconfig::Setting* configWidth = findLibConfigGroup(*configComponent, "g");
	if(not configWidth) {
		printErr << "component '" << getName() << "' has no section 'width'." << endl;
		return false;
	}

	if(not checkIfAllVariablesAreThere(configWidth, mandatoryArguments)) {
		printErr << "the 'width' section of the component '" << getName() << "' does not contain all required fields." << endl;
		return false;
	}

	configWidth->lookupValue("val", _gamma);
	configWidth->lookupValue("lower", _gammamin);
	configWidth->lookupValue("upper", _gammamax);
	configWidth->lookupValue("fix", _fixgamma);

	map<string, libconfig::Setting::Type> mandatoryArgumentsIsobarMasses;
	boost::assign::insert(mandatoryArgumentsIsobarMasses)
	                     ("mIsobar1", libconfig::Setting::TypeFloat)
	                     ("mIsobar2", libconfig::Setting::TypeFloat);
	if(not checkIfAllVariablesAreThere(configComponent, mandatoryArgumentsIsobarMasses)) {
		printErr << "not all required isobar masses are defined for the component '" << getName() << "'." << endl;
		return false;
	}

	configComponent->lookupValue("mIsobar1", _m1);
	configComponent->lookupValue("mIsobar2", _m2);

	printInfo << "component '" << getName() << "':" << endl
	          << "mass: " << _m0 << " MeV/c^2, limits: " << _m0min << "-" << _m0max << " MeV/c^2 " << (_fixm ? "(FIXED)" : "") << endl
	          << "width: " << _gamma << " MeV/c^2, limits: " << _gammamin << "-" << _gammamax << " MeV/c^2 " << (_fixgamma ? "(FIXED) " : "") << endl
	          << "mass of isobar 1: " << _m1 << " MeV/c^2, mass of isobar 2: " << _m2 << " MeV/c^2" << endl;

	if(debug) {
		printDebug << "finished initialization of 'pwabkg'." << endl;
	}

	return true;
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
