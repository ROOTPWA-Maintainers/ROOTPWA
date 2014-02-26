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
#include <Math/Minimizer.h>
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


std::complex<double>
pwachannel::CsqrtPS(const double mass,
                    const size_t idxMass) const
{
	if(idxMass != numeric_limits<size_t>::max()) {
		return _C * sqrt(_phaseSpace[idxMass]);
	}
  
	return _C*sqrt(_ps->Eval(mass));
}


double
pwachannel::ps(const double mass,
               const size_t idxMass) const
{
	if(idxMass != numeric_limits<size_t>::max()) {
		return _phaseSpace[idxMass];
	}
  
	return _ps->Eval(mass);
}


massDepFitComponent::massDepFitComponent(const string& name,
                                         const size_t nrParameters)
	: _name(name),
	  _nrParameters(nrParameters),
	  _parameters(nrParameters),
	  _parametersFixed(nrParameters),
	  _parametersLimitLower(nrParameters),
	  _parametersLimitedLower(nrParameters),
	  _parametersLimitUpper(nrParameters),
	  _parametersLimitedUpper(nrParameters),
	  _parametersName(nrParameters),
	  _parametersStep(nrParameters)
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

	for(size_t idxParameter=0; idxParameter<_nrParameters; ++idxParameter) {
		if(debug) {
			printDebug << "reading parameter '" << _parametersName[idxParameter] << "'." << endl;
		}

		const libconfig::Setting* configParameter = findLibConfigGroup(*configComponent, _parametersName[idxParameter]);
		if(not configParameter) {
			printErr << "component '" << getName() << "' has no section '" << _parametersName[idxParameter] << "'." << endl;
			return false;
		}

		map<string, libconfig::Setting::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("val", libconfig::Setting::TypeFloat)
		                     ("fix", libconfig::Setting::TypeBoolean);

		if(not checkIfAllVariablesAreThere(configParameter, mandatoryArguments)) {
			printErr << "the '" << _parametersName[idxParameter] << "' section of the component '" << getName() << "' does not contain all required fields." << endl;
			return false;
		}

		configParameter->lookupValue("val", _parameters[idxParameter]);
		bool fixed;
		configParameter->lookupValue("fix", fixed);
		_parametersFixed[idxParameter] = fixed;

		_parametersLimitedLower[idxParameter] = configParameter->lookupValue("lower", _parametersLimitLower[idxParameter]);
		_parametersLimitedUpper[idxParameter] = configParameter->lookupValue("upper", _parametersLimitUpper[idxParameter]);

		double step;
		if(configParameter->lookupValue("step", step)) {
			_parametersStep[idxParameter] = step;
		}
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

		// check that a wave with this wave is not yet in the decay channels
		for(size_t idx=0; idx<_channels.size(); ++idx) {
			if(_channels[idx].getWaveName() == waveName) {
				printErr << "wave '" << waveName << "' defined twice in the decay channels of '" << getName() << "'." << endl;
				return false;
			}
		}

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


bool
massDepFitComponent::update(const libconfig::Setting* configComponent,
                            const ROOT::Math::Minimizer* minimizer,
                            const bool debug) const
{
	if(debug) {
		printDebug << "starting updating of 'massDepFitComponent' for component '" << getName() << "'." << endl;
	}

	for(size_t idxParameter=0; idxParameter<_nrParameters; ++idxParameter) {
		if(debug) {
			printDebug << "updating parameter '" << _parametersName[idxParameter] << "'." << endl;
		}

		libconfig::Setting* configParameter = &((*configComponent)[_parametersName[idxParameter]]);
		if(not configParameter) {
			printErr << "component '" << getName() << "' has no section '" << _parametersName[idxParameter] << "'." << endl;
			return false;
		}

		(*configParameter)["val"] = _parameters[idxParameter];

		if(not configParameter->exists("error")) {
			configParameter->add("error", libconfig::Setting::TypeFloat);
		}
		(*configParameter)["error"] = minimizer->Errors()[minimizer->VariableIndex((getName() + _parametersName[idxParameter]).c_str())];
	}

	const libconfig::Setting* decayChannels = findLibConfigList(*configComponent, "decaychannels");
	if(not decayChannels) {
		printErr << "component '" << getName() << "' has no decay channels." << endl;
		return false;
	}

	const int nrDecayChannels = decayChannels->getLength();
	if(nrDecayChannels != getNrChannels()) {
		printErr << "number of decay channels in configuration file and fit model does not match." << endl;
		return false;
	}

	for(int idxDecayChannel=0; idxDecayChannel<nrDecayChannels; ++idxDecayChannel) {
		const libconfig::Setting* decayChannel = &((*decayChannels)[idxDecayChannel]);

		map<string, libconfig::Setting::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("amp", libconfig::Setting::TypeString);
		if(not checkIfAllVariablesAreThere(decayChannel, mandatoryArguments)) {
			printErr << "one of the decay channels of the component '" << getName() << "' does not contain all required fields." << endl;
			return false;
		}

		string waveName;
		decayChannel->lookupValue("amp", waveName);

		const pwachannel* channel = NULL;
		for(size_t idx=0; idx<nrDecayChannels; ++idx) {
			if(getChannelWaveName(idx) == waveName) {
				channel = &getChannel(idx);
				break;
			}
		}
		if(not channel) {
			printErr << "could not find channel '" << waveName << "' for component '" << getName() << "'." << endl;
			return false;
		}

		(*decayChannel)["coupling_Re"] = channel->C().real();
		(*decayChannel)["coupling_Im"] = channel->C().imag();
	}

	if(debug) {
		printDebug << "finished updating of 'massDepFitComponent'." << endl;
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


void
massDepFitComponent::getParameters(double* par) const
{
	for(size_t idx=0; idx<_nrParameters; ++idx) {
		par[idx] = _parameters[idx];
	}
}


void
massDepFitComponent::setParameters(const double* par)
{
	for(size_t idx=0; idx<_nrParameters; ++idx) {
		_parameters[idx] = par[idx];
	}
}


double
massDepFitComponent::getParameter(const size_t idx) const
{
	return _parameters[idx];
}


bool
massDepFitComponent::getParameterFixed(const size_t idx) const
{
	return _parametersFixed[idx];
}


double
massDepFitComponent::getParameterLimitLower(const size_t idx) const
{
	return _parametersLimitLower[idx];
}


bool
massDepFitComponent::getParameterLimitedLower(const size_t idx) const
{
	return _parametersLimitedLower[idx];
}


double
massDepFitComponent::getParameterLimitUpper(const size_t idx) const
{
	return _parametersLimitUpper[idx];
}


bool
massDepFitComponent::getParameterLimitedUpper(const size_t idx) const
{
	return _parametersLimitedUpper[idx];
}


const std::string&
massDepFitComponent::getParameterName(const size_t idx) const
{
	return _parametersName[idx];
}


double
massDepFitComponent::getParameterStep(const size_t idx) const
{
	return _parametersStep[idx];
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
	: massDepFitComponent(name, 2)
{
	_parametersName[0] = "mass";
	_parametersName[1] = "width";

	_parametersStep[0] = 1.0;
	_parametersStep[1] = 1.0;
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

	if(not massDepFitComponent::init(configComponent, massBinCenters, waveIndices, phaseSpaceIntegrals, debug)) {
		printErr << "error while reading configuration of 'massDepFitComponent' class." << endl;
		return false;
	}

	const libconfig::Setting* configWidth = findLibConfigGroup(*configComponent, "width");
	if(not configWidth) {
		printErr << "component '" << getName() << "' has no section 'width'." << endl;
		return false;
	}

	if(not configWidth->lookupValue("dyn", _constWidth)) {
		_constWidth = false;
	}
	_constWidth = !_constWidth;

	if(debug) {
		ostringstream output;
		output << "    mass: " << _parameters[0] << " MeV/c^2, ";
		if(_parametersLimitedLower[0] && _parametersLimitedUpper[0]) {
			output << "    limits: " << _parametersLimitLower[0] << "-" << _parametersLimitUpper[0] << " MeV/c^2 ";
		} else if(_parametersLimitedLower[0]) {
			output << "    lower limit: " << _parametersLimitLower[0] << " MeV/c^2 ";
		} else if(_parametersLimitedUpper[0]) {
			output << "    upper limit: " << _parametersLimitUpper[0] << " MeV/c^2 ";
		} else {
			output << "    unlimited ";
		}
		output << (_parametersFixed[0] ? "(FIXED)" : "") << endl;

		output << "    width: " << _parameters[1] << " MeV/c^2, ";
		if(_parametersLimitedLower[1] && _parametersLimitedUpper[1]) {
			output << "    limits: " << _parametersLimitLower[1] << "-" << _parametersLimitUpper[1] << " MeV/c^2 ";
		} else if(_parametersLimitedLower[1]) {
			output << "    lower limit: " << _parametersLimitLower[1] << " MeV/c^2 ";
		} else if(_parametersLimitedUpper[1]) {
			output << "    upper limit: " << _parametersLimitUpper[1] << " MeV/c^2 ";
		} else {
			output << "    unlimited ";
		}
		output << (_parametersFixed[1] ? "(FIXED) " : "");
		output << (_constWidth ? "-- constant width" : "-- dynamic width") << endl;

		printDebug << "component '" << getName() << "':" << endl << output.str();

		printDebug << "finished initialization of 'pwacomponent'." << endl;
	}
	return true;
}


complex<double>
pwacomponent::val(const double m) const {
	const double& m0 = _parameters[0];
	const double& gamma0 = _parameters[1];
	// calculate dynamic width:
	// loop over decay channels
	double ps = 1.;
	if(!_constWidth){
		ps=0; // no not forget to reset phase space here!
		for(vector<pwachannel>::const_iterator itChan=getChannels().begin(); itChan!=getChannels().end(); ++itChan) {
			double myps=1.;
			if(itChan->ps() != NULL){
				double ps0=itChan->ps(m0, std::numeric_limits<size_t>::max());
				myps=(itChan->ps(m, std::numeric_limits<size_t>::max()))/ps0;
			}
			ps+=myps;
		}
		ps /= getNrChannels();
	}
	const double gamma = gamma0*ps;
  
	return gamma0*m0 / complex<double>(m0*m0-m*m, -gamma*m0);
}


ostream&
pwacomponent::print(ostream& out) const
{
	out << getName() << endl
	    << "Mass=" << _parameters[0] << "   Width=" << _parameters[1] << "    ConstWidth=" << _constWidth << endl;
	return massDepFitComponent::print(out);
}


pwabkg::pwabkg(const string& name)
	: massDepFitComponent(name, 2)
{
	_parametersName[0] = "m0";
	_parametersName[1] = "g";

	_parametersStep[0] = 1.0;
	_parametersStep[1] = 1.0e-6;
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

	if(not massDepFitComponent::init(configComponent, massBinCenters, waveIndices, phaseSpaceIntegrals, debug)) {
		printErr << "error while reading configuration of 'massDepFitComponent' class." << endl;
		return false;
	}

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

	if(debug) {
		ostringstream output;
		output << "    mass: " << _parameters[0] << " MeV/c^2, ";
		if(_parametersLimitedLower[0] && _parametersLimitedUpper[0]) {
			output << "    limits: " << _parametersLimitLower[0] << "-" << _parametersLimitUpper[0] << " MeV/c^2 ";
		} else if(_parametersLimitedLower[0]) {
			output << "    lower limit: " << _parametersLimitLower[0] << " MeV/c^2 ";
		} else if(_parametersLimitedUpper[0]) {
			output << "    upper limit: " << _parametersLimitUpper[0] << " MeV/c^2 ";
		} else {
			output << "    unlimited ";
		}
		output << (_parametersFixed[0] ? "(FIXED)" : "") << endl;

		output << "    width: " << _parameters[1] << " MeV/c^2, ";
		if(_parametersLimitedLower[1] && _parametersLimitedUpper[1]) {
			output << "    limits: " << _parametersLimitLower[1] << "-" << _parametersLimitUpper[1] << " MeV/c^2 ";
		} else if(_parametersLimitedLower[1]) {
			output << "    lower limit: " << _parametersLimitLower[1] << " MeV/c^2 ";
		} else if(_parametersLimitedUpper[1]) {
			output << "    upper limit: " << _parametersLimitUpper[1] << " MeV/c^2 ";
		} else {
			output << "    unlimited ";
		}
		output << (_parametersFixed[1] ? "(FIXED) " : "") << endl;

		output << "    mass of isobar 1: " << _m1 << " MeV/c^2, mass of isobar 2: " << _m2 << " MeV/c^2" << endl;

		printDebug << "component '" << getName() << "':" << endl << output.str();

		printDebug << "finished initialization of 'pwabkg'." << endl;
	}

	return true;
}


complex<double>
pwabkg::val(const double m) const {
	// shift baseline mass
	const double mass = m - _parameters[0];

	// calculate breakup momentum
	if(mass < _m1+_m2) {
		return complex<double>(1,0);
	}
	const complex<double> p = q(mass,_m1,_m2);

	return exp(-_parameters[1]*p.real()*p.real());
}
