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

#include "massDepFitComponents.h"

#include <boost/assign/std/vector.hpp>

#include <Math/Minimizer.h>

#include "libConfigUtils.hpp"
#include "physUtils.hpp"
#include "reportingUtils.hpp"

using namespace std;


rpwa::massDepFit::channel::channel(const std::string& waveName,
                                   std::complex<double> coupling,
                                   const std::vector<double>& massBinCenters,
                                   const std::vector<double>& phaseSpace)
	: _waveName(waveName),
	  _coupling(coupling),
	  _massBinCenters(massBinCenters),
	  _phaseSpace(phaseSpace)
{
	_interpolator = new ROOT::Math::Interpolator(_massBinCenters, _phaseSpace, ROOT::Math::Interpolation::kLINEAR);
}


rpwa::massDepFit::channel::channel(const rpwa::massDepFit::channel& ch)
	: _waveName(ch._waveName),
	  _coupling(ch._coupling),
	  _massBinCenters(ch._massBinCenters),
	  _phaseSpace(ch._phaseSpace)
{
	_interpolator = new ROOT::Math::Interpolator(_massBinCenters, _phaseSpace, ROOT::Math::Interpolation::kLINEAR);
}


rpwa::massDepFit::channel::~channel() {
	delete _interpolator;
}


rpwa::massDepFit::component::component(const string& name,
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
rpwa::massDepFit::component::init(const libconfig::Setting* configComponent,
                                  const vector<double>& massBinCenters,
                                  const map<string, size_t>& waveIndices,
                                  const boost::multi_array<double, 2>& phaseSpaceIntegrals,
                                  const bool debug) {
	if(debug) {
		printDebug << "starting initialization of 'component' for component '" << getName() << "'." << endl;
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
		_channels.push_back(rpwa::massDepFit::channel(waveName, coupling, massBinCenters, std::vector<double>(view.begin(), view.end())));

		if(debug) {
			printDebug << "coupling to wave '" << waveName << "' (index " << it->second << ") with coupling " << coupling << "." << endl;
		}
	}

	if(debug) {
		printDebug << "finished initialization of 'component'." << endl;
	}

	return true;
}


bool
rpwa::massDepFit::component::update(const libconfig::Setting* configComponent,
                                    const ROOT::Math::Minimizer* minimizer,
                                    const bool debug) const
{
	if(debug) {
		printDebug << "starting updating of 'component' for component '" << getName() << "'." << endl;
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

		const rpwa::massDepFit::channel* channel = NULL;
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

		(*decayChannel)["coupling_Re"] = channel->getCoupling().real();
		(*decayChannel)["coupling_Im"] = channel->getCoupling().imag();
	}

	if(debug) {
		printDebug << "finished updating of 'component'." << endl;
	}

	return true;
}


void
rpwa::massDepFit::component::getCouplings(double* par) const
{
	size_t counter=0;
	for(vector<rpwa::massDepFit::channel>::const_iterator it=_channels.begin(); it!=_channels.end(); ++it, counter+=2) {
		par[counter] = it->getCoupling().real();
		par[counter+1] = it->getCoupling().imag();
	}
}


void
rpwa::massDepFit::component::setCouplings(const double* par)
{
	size_t counter=0;
	for(vector<rpwa::massDepFit::channel>::iterator it=_channels.begin(); it!=_channels.end(); ++it, counter+=2) {
		it->setCoupling(complex<double>(par[counter], par[counter+1]));
	}
}


void
rpwa::massDepFit::component::getParameters(double* par) const
{
	for(size_t idx=0; idx<_nrParameters; ++idx) {
		par[idx] = _parameters[idx];
	}
}


void
rpwa::massDepFit::component::setParameters(const double* par)
{
	for(size_t idx=0; idx<_nrParameters; ++idx) {
		_parameters[idx] = par[idx];
	}
}


ostream&
rpwa::massDepFit::component::print(ostream& out) const
{
	out << "Decay modes:" << endl;
	for(vector<rpwa::massDepFit::channel>::const_iterator it=_channels.begin(); it!=_channels.end(); ++it) {
		out << it->getWaveName() << "    C=" << it->getCoupling() << endl;
	}
	return out;
}


rpwa::massDepFit::pwacomponent::pwacomponent(const string& name)
	: component(name, 2)
{
	_parametersName[0] = "mass";
	_parametersName[1] = "width";

	_parametersStep[0] = 1.0;
	_parametersStep[1] = 1.0;
}


bool
rpwa::massDepFit::pwacomponent::init(const libconfig::Setting* configComponent,
                                     const vector<double>& massBinCenters,
                                     const map<string, size_t>& waveIndices,
                                     const boost::multi_array<double, 2>& phaseSpaceIntegrals,
                                     const bool debug) {
	if(debug) {
		printDebug << "starting initialization of 'pwacomponent' for component '" << getName() << "'." << endl;
	}

	if(not component::init(configComponent, massBinCenters, waveIndices, phaseSpaceIntegrals, debug)) {
		printErr << "error while reading configuration of 'component' class." << endl;
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
rpwa::massDepFit::pwacomponent::val(const double m) const {
	const double& m0 = _parameters[0];
	const double& gamma0 = _parameters[1];
	// calculate dynamic width:
	// loop over decay channels
	double ps = 1.;
	if(!_constWidth){
		ps=0; // no not forget to reset phase space here!
		for(vector<rpwa::massDepFit::channel>::const_iterator itChan=getChannels().begin(); itChan!=getChannels().end(); ++itChan) {
			const double ps0 = itChan->getPhaseSpace(m0, std::numeric_limits<size_t>::max());
			const double ratio = itChan->getPhaseSpace(m, std::numeric_limits<size_t>::max()) / ps0;
			ps += ratio*ratio;
		}
		ps /= getNrChannels();
	}
	const double gamma = gamma0*ps;
  
	return gamma0*m0 / complex<double>(m0*m0-m*m, -gamma*m0);
}


ostream&
rpwa::massDepFit::pwacomponent::print(ostream& out) const
{
	out << getName() << endl
	    << "Mass=" << _parameters[0] << "   Width=" << _parameters[1] << "    ConstWidth=" << _constWidth << endl;
	return component::print(out);
}


rpwa::massDepFit::pwabkg::pwabkg(const string& name)
	: component(name, 2)
{
	_parametersName[0] = "m0";
	_parametersName[1] = "g";

	_parametersStep[0] = 1.0;
	_parametersStep[1] = 1.0e-6;
}


bool
rpwa::massDepFit::pwabkg::init(const libconfig::Setting* configComponent,
                               const vector<double>& massBinCenters,
                               const map<string, size_t>& waveIndices,
                               const boost::multi_array<double, 2>& phaseSpaceIntegrals,
                               const bool debug) {
	if(debug) {
		printDebug << "starting initialization of 'pwabkg' for component '" << getName() << "'." << endl;
	}

	if(not component::init(configComponent, massBinCenters, waveIndices, phaseSpaceIntegrals, debug)) {
		printErr << "error while reading configuration of 'component' class." << endl;
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
rpwa::massDepFit::pwabkg::val(const double m) const {
	// shift baseline mass
	const double mass = m - _parameters[0];

	// calculate breakup momentum
	if(mass < _m1+_m2) {
		return complex<double>(1,0);
	}
	const double q2 = rpwa::breakupMomentumSquared(mass, _m1, _m2);

	return exp(-_parameters[1]*q2);
}
