//-----------------------------------------------------------
//
// Description:
//      Implementation of components and channels
//      see massDepFitComponents.h for details
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


rpwa::massDepFit::channel::channel(const std::string& waveName,
                                   const size_t nrBins,
                                   const std::vector<std::complex<double> >& couplings,
                                   const std::vector<double>& massBinCenters,
                                   const boost::multi_array<double, 2>& phaseSpace)
	: _waveName(waveName),
	  _anchor(false),
	  _nrBins(nrBins),
	  _couplings(couplings),
	  _massBinCenters(massBinCenters),
	  _phaseSpace(phaseSpace)
{
	_interpolator.resize(_nrBins);
	for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
		boost::multi_array<double, 2>::const_array_view<1>::type view = _phaseSpace[boost::indices[idxBin][boost::multi_array<double, 2>::index_range()]];
		_interpolator[idxBin] = new ROOT::Math::Interpolator(_massBinCenters, std::vector<double>(view.begin(), view.end()), ROOT::Math::Interpolation::kLINEAR);
	}
}


rpwa::massDepFit::channel::channel(const rpwa::massDepFit::channel& ch)
	: _waveName(ch._waveName),
	  _anchor(ch._anchor),
	  _nrBins(ch._nrBins),
	  _couplings(ch._couplings),
	  _massBinCenters(ch._massBinCenters),
	  _phaseSpace(ch._phaseSpace)
{
	_interpolator.resize(_nrBins);
	for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
		boost::multi_array<double, 2>::const_array_view<1>::type view = _phaseSpace[boost::indices[idxBin][boost::multi_array<double, 2>::index_range()]];
		_interpolator[idxBin] = new ROOT::Math::Interpolator(_massBinCenters, std::vector<double>(view.begin(), view.end()), ROOT::Math::Interpolation::kLINEAR);
	}
}


rpwa::massDepFit::channel::~channel() {
	for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
		delete _interpolator[idxBin];
	}
}


rpwa::massDepFit::component::component(const std::string& name,
                                       const size_t nrParameters)
	: _name(name),
	  _nrParameters(nrParameters),
	  _nrCouplings(0),
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
                                  const size_t nrBins,
                                  const std::vector<double>& massBinCenters,
                                  const std::map<std::string, size_t>& waveIndices,
                                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                  const bool debug) {
	if(debug) {
		printDebug << "starting initialization of 'component' for component '" << getName() << "'." << std::endl;
	}

	for(size_t idxParameter=0; idxParameter<_nrParameters; ++idxParameter) {
		if(debug) {
			printDebug << "reading parameter '" << _parametersName[idxParameter] << "'." << std::endl;
		}

		const libconfig::Setting* configParameter = findLibConfigGroup(*configComponent, _parametersName[idxParameter]);
		if(not configParameter) {
			printErr << "component '" << getName() << "' has no section '" << _parametersName[idxParameter] << "'." << std::endl;
			return false;
		}

		std::map<std::string, libconfig::Setting::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("val", libconfig::Setting::TypeFloat)
		                     ("fix", libconfig::Setting::TypeBoolean);

		if(not checkIfAllVariablesAreThere(configParameter, mandatoryArguments)) {
			printErr << "the '" << _parametersName[idxParameter] << "' section of the component '" << getName() << "' does not contain all required fields." << std::endl;
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
		printErr << "component '" << getName() << "' has no decay channels." << std::endl;
		return false;
	}

	const int nrDecayChannels = decayChannels->getLength();
	for(int idxDecayChannel=0; idxDecayChannel<nrDecayChannels; ++idxDecayChannel) {
		const libconfig::Setting* decayChannel = &((*decayChannels)[idxDecayChannel]);

		std::map<std::string, libconfig::Setting::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("amp", libconfig::Setting::TypeString);
		if(not checkIfAllVariablesAreThere(decayChannel, mandatoryArguments)) {
			printErr << "one of the decay channels of the component '" << getName() << "' does not contain all required fields." << std::endl;
			return false;
		}

		std::string waveName;
		decayChannel->lookupValue("amp", waveName);

		// check that a wave with this wave is not yet in the decay channels
		for(size_t idx=0; idx<_channels.size(); ++idx) {
			if(_channels[idx].getWaveName() == waveName) {
				printErr << "wave '" << waveName << "' defined twice in the decay channels of '" << getName() << "'." << std::endl;
				return false;
			}
		}

		const std::map<std::string, size_t>::const_iterator it = waveIndices.find(waveName);
		if(it == waveIndices.end()) {
			printErr << "wave '" << waveName << "' not in fit, but used as decay channel." << std::endl;
			return false;
		}

		bool readCouplings = true;

		std::vector<std::complex<double> > couplings;
		if(readCouplings) {
			++_nrCouplings;

			boost::assign::insert(mandatoryArguments)
			                     ("couplings", libconfig::Setting::TypeList);
			if(not checkIfAllVariablesAreThere(decayChannel, mandatoryArguments)) {
				printErr << "one of the decay channels of the component '" << getName() << "' does not contain all required fields." << std::endl;
				return false;
			}

			const libconfig::Setting* configCouplings = findLibConfigList(*decayChannel, "couplings");
			if(not configCouplings) {
				printErr << "decay channel '" << waveName << "' of component '" << getName() << "' has no couplings." << std::endl;
				return false;
			}

			const int nrCouplings = configCouplings->getLength();
			if(nrCouplings < 0 || static_cast<size_t>(nrCouplings) != nrBins) {
				printErr << "decay channel '" << waveName << "' of component '" << getName() << "' has " << nrCouplings << " couplings, not " << nrBins << "." << std::endl;
				return false;
			}

			for(int idxCoupling=0; idxCoupling<nrCouplings; ++idxCoupling) {
				const libconfig::Setting* configCoupling = &((*configCouplings)[idxCoupling]);

				std::map<std::string, libconfig::Setting::Type> mandatoryArguments;
				boost::assign::insert(mandatoryArguments)
				                     ("coupling_Re", libconfig::Setting::TypeFloat)
				                     ("coupling_Im", libconfig::Setting::TypeFloat);
				if(not checkIfAllVariablesAreThere(configCoupling, mandatoryArguments)) {
					printErr << "one of the couplings of the decay channel '" << waveName << "' of the component '" << getName() << "' does not contain all required fields." << std::endl;
					return false;
				}

				double couplingReal;
				configCoupling->lookupValue("coupling_Re", couplingReal);

				double couplingImag;
				configCoupling->lookupValue("coupling_Im", couplingImag);

				const std::complex<double> coupling(couplingReal, couplingImag);
				couplings.push_back(coupling);
			}
		} else {
			for(size_t i=0; i<nrBins; ++i) {
				couplings.push_back(std::complex<double>(1.0, 0.0));
			}
		}

		boost::multi_array<double, 3>::const_array_view<2>::type view = phaseSpaceIntegrals[boost::indices[boost::multi_array<double, 3>::index_range()][boost::multi_array<double, 3>::index_range()][it->second]];
		_channels.push_back(rpwa::massDepFit::channel(waveName, nrBins, couplings, massBinCenters, view));

		if(debug) {
			std::ostringstream output;
			output << "( ";
			for(size_t idxBin=0; idxBin<nrBins; ++idxBin) {
				output << couplings[idxBin] << " ";
			}
			output << ")";
			printDebug << "coupling to wave '" << waveName << "' (index " << it->second << ") with coupling " << output.str() << "." << std::endl;
		}
	}

	if(debug) {
		printDebug << "finished initialization of 'component'." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::component::update(const libconfig::Setting* configComponent,
                                    const ROOT::Math::Minimizer* minimizer,
                                    const bool debug) const
{
	if(debug) {
		printDebug << "starting updating of 'component' for component '" << getName() << "'." << std::endl;
	}

	for(size_t idxParameter=0; idxParameter<_nrParameters; ++idxParameter) {
		if(debug) {
			printDebug << "updating parameter '" << _parametersName[idxParameter] << "'." << std::endl;
		}

		libconfig::Setting* configParameter = &((*configComponent)[_parametersName[idxParameter]]);
		if(not configParameter) {
			printErr << "component '" << getName() << "' has no section '" << _parametersName[idxParameter] << "'." << std::endl;
			return false;
		}

		(*configParameter)["val"] = _parameters[idxParameter];

		if(not configParameter->exists("error")) {
			configParameter->add("error", libconfig::Setting::TypeFloat);
		}

		const std::string varName = getName() + "__" + _parametersName[idxParameter];
		const int varIndex = minimizer->VariableIndex(varName);
		if(varIndex == -1) {
			printErr << "variable '" << varName << "' used to extract the error for the parameter '"
			         << _parametersName[idxParameter] << "' of '" << getName() << "' not known to the minimizer." << std::endl;
			return false;
		}
		(*configParameter)["error"] = minimizer->Errors()[varIndex];
	}

	const libconfig::Setting* decayChannels = findLibConfigList(*configComponent, "decaychannels");
	if(not decayChannels) {
		printErr << "component '" << getName() << "' has no decay channels." << std::endl;
		return false;
	}

	const int nrDecayChannels = decayChannels->getLength();
	if(nrDecayChannels < 0 || static_cast<size_t>(nrDecayChannels) != getNrChannels()) {
		printErr << "number of decay channels in configuration file and fit model does not match." << std::endl;
		return false;
	}

	for(int idxDecayChannel=0; idxDecayChannel<nrDecayChannels; ++idxDecayChannel) {
		const libconfig::Setting* decayChannel = &((*decayChannels)[idxDecayChannel]);

		std::map<std::string, libconfig::Setting::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("amp", libconfig::Setting::TypeString)
		                     ("couplings", libconfig::Setting::TypeList);
		if(not checkIfAllVariablesAreThere(decayChannel, mandatoryArguments)) {
			printErr << "one of the decay channels of the component '" << getName() << "' does not contain all required fields." << std::endl;
			return false;
		}

		std::string waveName;
		decayChannel->lookupValue("amp", waveName);

		const rpwa::massDepFit::channel* channel = NULL;
		for(size_t idx=0; idx<getNrChannels(); ++idx) {
			if(getChannelWaveName(idx) == waveName) {
				channel = &getChannel(idx);
				break;
			}
		}
		if(not channel) {
			printErr << "could not find channel '" << waveName << "' for component '" << getName() << "'." << std::endl;
			return false;
		}

		const libconfig::Setting* configCouplings = findLibConfigList(*decayChannel, "couplings");
		if(not configCouplings) {
			printErr << "decay channel '" << waveName << "' of component '" << getName() << "' has no couplings." << std::endl;
			return false;
		}

		const int nrCouplings = configCouplings->getLength();
		if(nrCouplings < 0 || static_cast<size_t>(nrCouplings) != channel->getNrBins()) {
			printErr << "decay channel '" << waveName << "' of component '" << getName() << "' has only " << nrCouplings << " couplings, not " << channel->getNrBins() << "." << std::endl;
			return false;
		}

		for(int idxCoupling=0; idxCoupling<nrCouplings; ++idxCoupling) {
			const libconfig::Setting* configCoupling = &((*configCouplings)[idxCoupling]);

			(*configCoupling)["coupling_Re"] = channel->getCoupling(idxCoupling).real();
			(*configCoupling)["coupling_Im"] = channel->getCoupling(idxCoupling).imag();
		}
	}

	if(debug) {
		printDebug << "finished updating of 'component'." << std::endl;
	}

	return true;
}


size_t
rpwa::massDepFit::component::getCouplings(double* par) const
{
	size_t counter=0;
	for(std::vector<rpwa::massDepFit::channel>::const_iterator it=_channels.begin(); it!=_channels.end(); ++it) {
		const size_t nrBins = it->getNrBins();
		const std::vector<std::complex<double> >& couplings = it->getCouplings();
		for(size_t idxBin=0; idxBin<nrBins; ++idxBin) {
			par[counter] = couplings[idxBin].real();
			counter += 1;

			if(not it->isAnchor()) {
				par[counter] = couplings[idxBin].imag();
				counter += 1;
			}
		}
	}

	return counter;
}


size_t
rpwa::massDepFit::component::setCouplings(const double* par)
{
	size_t counter=0;
	for(std::vector<rpwa::massDepFit::channel>::iterator it=_channels.begin(); it!=_channels.end(); ++it) {
		std::vector<std::complex<double> > couplings;
		if(it->isAnchor()) {
			const size_t nrBins = it->getNrBins();
			for(size_t idxBin=0; idxBin<nrBins; ++idxBin) {
				couplings.push_back(std::complex<double>(par[counter], 0.));
				counter += 1;
			}
		} else {
			const size_t nrBins = it->getNrBins();
			for(size_t idxBin=0; idxBin<nrBins; ++idxBin) {
				couplings.push_back(std::complex<double>(par[counter], par[counter+1]));
				counter += 2;
			}
		}
		it->setCouplings(couplings);
	}

	return counter;
}


size_t
rpwa::massDepFit::component::getParameters(double* par) const
{
	for(size_t idx=0; idx<_nrParameters; ++idx) {
		par[idx] = _parameters[idx];
	}

	return _nrParameters;
}


size_t
rpwa::massDepFit::component::setParameters(const double* par)
{
	for(size_t idx=0; idx<_nrParameters; ++idx) {
		_parameters[idx] = par[idx];
	}

	return _nrParameters;
}


std::ostream&
rpwa::massDepFit::component::print(std::ostream& out) const
{
	out << "Decay modes:" << std::endl;
	for(std::vector<rpwa::massDepFit::channel>::const_iterator it=_channels.begin(); it!=_channels.end(); ++it) {
		std::ostringstream output;
		output << "( ";
		for(size_t idxBin=0; idxBin<it->getNrBins(); ++idxBin) {
			output << it->getCoupling(idxBin) << " ";
		}
		output << ")";
		out << it->getWaveName() << "    C=" << output.str() << std::endl;
	}
	return out;
}


rpwa::massDepFit::fixedWidthBreitWigner::fixedWidthBreitWigner(const std::string& name)
	: component(name, 2)
{
	_parametersName[0] = "mass";
	_parametersName[1] = "width";

	_parametersStep[0] = 0.001;
	_parametersStep[1] = 0.001;
}


bool
rpwa::massDepFit::fixedWidthBreitWigner::init(const libconfig::Setting* configComponent,
                                              const size_t nrBins,
                                              const std::vector<double>& massBinCenters,
                                              const std::map<std::string, size_t>& waveIndices,
                                              const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                              const bool debug)
{
	if(debug) {
		printDebug << "starting initialization of 'fixedWidthBreitWigner' for component '" << getName() << "'." << std::endl;
	}

	if(not component::init(configComponent, nrBins, massBinCenters, waveIndices, phaseSpaceIntegrals, debug)) {
		printErr << "error while reading configuration of 'component' class." << std::endl;
		return false;
	}

	if(debug) {
		std::ostringstream output;
		output << "    mass: " << _parameters[0] << " GeV/c^2, ";
		if(_parametersLimitedLower[0] && _parametersLimitedUpper[0]) {
			output << "    limits: " << _parametersLimitLower[0] << "-" << _parametersLimitUpper[0] << " GeV/c^2 ";
		} else if(_parametersLimitedLower[0]) {
			output << "    lower limit: " << _parametersLimitLower[0] << " GeV/c^2 ";
		} else if(_parametersLimitedUpper[0]) {
			output << "    upper limit: " << _parametersLimitUpper[0] << " GeV/c^2 ";
		} else {
			output << "    unlimited ";
		}
		output << (_parametersFixed[0] ? "(FIXED)" : "") << std::endl;

		output << "    width: " << _parameters[1] << " GeV/c^2, ";
		if(_parametersLimitedLower[1] && _parametersLimitedUpper[1]) {
			output << "    limits: " << _parametersLimitLower[1] << "-" << _parametersLimitUpper[1] << " GeV/c^2 ";
		} else if(_parametersLimitedLower[1]) {
			output << "    lower limit: " << _parametersLimitLower[1] << " GeV/c^2 ";
		} else if(_parametersLimitedUpper[1]) {
			output << "    upper limit: " << _parametersLimitUpper[1] << " GeV/c^2 ";
		} else {
			output << "    unlimited ";
		}
		output << (_parametersFixed[1] ? "(FIXED)" : "") << std::endl;

		printDebug << "component '" << getName() << "':" << std::endl << output.str();

		printDebug << "finished initialization of 'fixedWidthBreitWigner'." << std::endl;
	}
	return true;
}


std::complex<double>
rpwa::massDepFit::fixedWidthBreitWigner::val(const size_t idxBin,
                                             const double m) const
{
	const double& m0 = _parameters[0];
	const double& gamma0 = _parameters[1];

	return gamma0*m0 / std::complex<double>(m0*m0-m*m, -gamma0*m0);
}


std::ostream&
rpwa::massDepFit::fixedWidthBreitWigner::print(std::ostream& out) const
{
	out << getName() << std::endl
	    << "Mass=" << _parameters[0] << "   Width=" << _parameters[1] << std::endl;
	return component::print(out);
}


rpwa::massDepFit::relativisticBreitWigner::relativisticBreitWigner(const std::string& name)
	: component(name, 2)
{
	_parametersName[0] = "mass";
	_parametersName[1] = "width";

	_parametersStep[0] = 0.001;
	_parametersStep[1] = 0.001;
}


bool
rpwa::massDepFit::relativisticBreitWigner::init(const libconfig::Setting* configComponent,
                                                const size_t nrBins,
                                                const std::vector<double>& massBinCenters,
                                                const std::map<std::string, size_t>& waveIndices,
                                                const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                                const bool debug)
{
	if(debug) {
		printDebug << "starting initialization of 'relativisticBreitWigner' for component '" << getName() << "'." << std::endl;
	}

	if(not component::init(configComponent, nrBins, massBinCenters, waveIndices, phaseSpaceIntegrals, debug)) {
		printErr << "error while reading configuration of 'component' class." << std::endl;
		return false;
	}

	std::map<std::string, libconfig::Setting::Type> mandatoryArguments;
	boost::assign::insert(mandatoryArguments)
	                     ("mIsobar1", libconfig::Setting::TypeArray)
	                     ("mIsobar2", libconfig::Setting::TypeArray)
	                     ("relAngularMom", libconfig::Setting::TypeArray);
	if(not checkIfAllVariablesAreThere(configComponent, mandatoryArguments)) {
		printErr << "not all required arrays are defined for the component '" << getName() << "'." << std::endl;
		return false;
	}

	const libconfig::Setting& configMassIsobar1 = (*configComponent)["mIsobar1"];
	const libconfig::Setting& configMassIsobar2 = (*configComponent)["mIsobar2"];
	const libconfig::Setting& configRelAngularMom = (*configComponent)["relAngularMom"];

	const int nrDecayChannels = configMassIsobar1.getLength();
	if(configMassIsobar2.getLength() != nrDecayChannels || configRelAngularMom.getLength() != nrDecayChannels) {
		printErr << "not all arrays have the same length for component '" << getName() << "'." << std::endl;
		return false;
	}

	if(configMassIsobar1[0].getType() != libconfig::Setting::TypeFloat
	   || configMassIsobar2[0].getType() != libconfig::Setting::TypeFloat
	   || configRelAngularMom[0].getType() != libconfig::Setting::TypeInt) {
		printErr << "not all arrays have the correct type for component '" << getName() << "'." << std::endl;
		return false;
	}

	_ratio.resize(nrDecayChannels);
	_l.resize(nrDecayChannels);
	_m1.resize(nrDecayChannels);
	_m2.resize(nrDecayChannels);

	for(int i=0; i<nrDecayChannels; ++i) {
		_l[i] = configRelAngularMom[i];
		_m1[i] = configMassIsobar1[i];
		_m2[i] = configMassIsobar2[i];
	}

	if(nrDecayChannels > 1) {
		boost::assign::insert(mandatoryArguments)
		                     ("branchingRatio", libconfig::Setting::TypeArray);
		if(not checkIfAllVariablesAreThere(configComponent, mandatoryArguments)) {
			printErr << "not all required arrays are defined for the component '" << getName() << "'." << std::endl;
			return false;
		}

		const libconfig::Setting& configBranchingRatio = (*configComponent)["branchingRatio"];
		if(configBranchingRatio.getLength() != nrDecayChannels) {
			printErr << "not all arrays have the same length for component '" << getName() << "'." << std::endl;
			return false;
		}
		if(configBranchingRatio[0].getType() != libconfig::Setting::TypeFloat) {
			printErr << "not all arrays have the correct type for component '" << getName() << "'." << std::endl;
			return false;
		}

		double sum = 0.;
		for(int i=0; i<nrDecayChannels; ++i) {
			_ratio[i] = configBranchingRatio[i];
			sum += _ratio[i];
		}
		for(size_t i=0; i<_ratio.size(); ++i) {
			_ratio[i] *= 1. / sum;
		}
	} else {
		_ratio[0] = 1.;
	}

	if(debug) {
		std::ostringstream output;
		output << "    mass: " << _parameters[0] << " GeV/c^2, ";
		if(_parametersLimitedLower[0] && _parametersLimitedUpper[0]) {
			output << "    limits: " << _parametersLimitLower[0] << "-" << _parametersLimitUpper[0] << " GeV/c^2 ";
		} else if(_parametersLimitedLower[0]) {
			output << "    lower limit: " << _parametersLimitLower[0] << " GeV/c^2 ";
		} else if(_parametersLimitedUpper[0]) {
			output << "    upper limit: " << _parametersLimitUpper[0] << " GeV/c^2 ";
		} else {
			output << "    unlimited ";
		}
		output << (_parametersFixed[0] ? "(FIXED)" : "") << std::endl;

		output << "    width: " << _parameters[1] << " GeV/c^2, ";
		if(_parametersLimitedLower[1] && _parametersLimitedUpper[1]) {
			output << "    limits: " << _parametersLimitLower[1] << "-" << _parametersLimitUpper[1] << " GeV/c^2 ";
		} else if(_parametersLimitedLower[1]) {
			output << "    lower limit: " << _parametersLimitLower[1] << " GeV/c^2 ";
		} else if(_parametersLimitedUpper[1]) {
			output << "    upper limit: " << _parametersLimitUpper[1] << " GeV/c^2 ";
		} else {
			output << "    unlimited ";
		}
		output << (_parametersFixed[1] ? "(FIXED)" : "") << std::endl;

		output << "    " << _ratio.size() << " decay channels:" << std::endl;
		for(size_t i=0; i<_ratio.size(); ++i) {
			output << "      * decay channel " << i << ", branching ratio: " << _ratio[i] << std::endl;
			output << "        mass of isobar 1: " << _m1[i] << " GeV/c^2, mass of isobar 2: " << _m2[i] << " GeV/c^2" << std::endl;
			output << "        relative orbital angular momentum between isobars: " << _l[i] << " (in units of hbar)" << std::endl;
		}

		printDebug << "component '" << getName() << "':" << std::endl << output.str();

		printDebug << "finished initialization of 'relativisticBreitWigner'." << std::endl;
	}
	return true;
}


std::complex<double>
rpwa::massDepFit::relativisticBreitWigner::val(const size_t idxBin,
                                               const double m) const
{
	const double& m0 = _parameters[0];
	const double& gamma0 = _parameters[1];

	double gamma = 0.;
	for(size_t i=0; i<_ratio.size(); ++i) {
		if(m >= _m1[i] + _m2[i]) {
			// calculate breakup momenta
			const double q = rpwa::breakupMomentum(m, _m1[i], _m2[i]);
			const double q0 = rpwa::breakupMomentum(m0, _m1[i], _m2[i]);

			// calculate barrier factors
			const double f2 = rpwa::barrierFactorSquared(2*_l[i], q);
			const double f20 = rpwa::barrierFactorSquared(2*_l[i], q0);

			gamma += _ratio[i] * q/q0 * f2/f20;
		}
	}
	gamma *= gamma0 * m0/m;

	return gamma0*m0 / std::complex<double>(m0*m0-m*m, -gamma*m0);
}


std::ostream&
rpwa::massDepFit::relativisticBreitWigner::print(std::ostream& out) const
{
	out << getName() << std::endl
	    << "Mass=" << _parameters[0] << "   Width=" << _parameters[1] << std::endl;
	return component::print(out);
}


rpwa::massDepFit::parameterizationA1Bowler::parameterizationA1Bowler(const std::string& name)
	: component(name, 2)
{
	_parametersName[0] = "mass";
	_parametersName[1] = "width";

	_parametersStep[0] = 0.001;
	_parametersStep[1] = 0.001;
}


bool
rpwa::massDepFit::parameterizationA1Bowler::init(const libconfig::Setting* configComponent,
                                                 const size_t nrBins,
                                                 const std::vector<double>& massBinCenters,
                                                 const std::map<std::string, size_t>& waveIndices,
                                                 const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                                 const bool debug)
{
	if(debug) {
		printDebug << "starting initialization of 'parameterizationA1Bowler' for component '" << getName() << "'." << std::endl;
	}

	if(not component::init(configComponent, nrBins, massBinCenters, waveIndices, phaseSpaceIntegrals, debug)) {
		printErr << "error while reading configuration of 'component' class." << std::endl;
		return false;
	}

	// just to be on the safe side, in principle this should be possible,
	// but probably one should then also give branching ratios somewhere.
	if(getNrChannels() != 1) {
		printErr << "component of type 'parameterizationA1Bowler' must have exactly one channel." << std::endl;
		return false;
	}

	if(debug) {
		std::ostringstream output;
		output << "    mass: " << _parameters[0] << " GeV/c^2, ";
		if(_parametersLimitedLower[0] && _parametersLimitedUpper[0]) {
			output << "    limits: " << _parametersLimitLower[0] << "-" << _parametersLimitUpper[0] << " GeV/c^2 ";
		} else if(_parametersLimitedLower[0]) {
			output << "    lower limit: " << _parametersLimitLower[0] << " GeV/c^2 ";
		} else if(_parametersLimitedUpper[0]) {
			output << "    upper limit: " << _parametersLimitUpper[0] << " GeV/c^2 ";
		} else {
			output << "    unlimited ";
		}
		output << (_parametersFixed[0] ? "(FIXED)" : "") << std::endl;

		output << "    width: " << _parameters[1] << " GeV/c^2, ";
		if(_parametersLimitedLower[1] && _parametersLimitedUpper[1]) {
			output << "    limits: " << _parametersLimitLower[1] << "-" << _parametersLimitUpper[1] << " GeV/c^2 ";
		} else if(_parametersLimitedLower[1]) {
			output << "    lower limit: " << _parametersLimitLower[1] << " GeV/c^2 ";
		} else if(_parametersLimitedUpper[1]) {
			output << "    upper limit: " << _parametersLimitUpper[1] << " GeV/c^2 ";
		} else {
			output << "    unlimited ";
		}
		output << (_parametersFixed[1] ? "(FIXED)" : "") << std::endl;

		printDebug << "component '" << getName() << "':" << std::endl << output.str();

		printDebug << "finished initialization of 'parameterizationA1Bowler'." << std::endl;
	}
	return true;
}


std::complex<double>
rpwa::massDepFit::parameterizationA1Bowler::val(const size_t idxBin,
                                                const double m) const
{
	const double& m0 = _parameters[0];
	const double& gamma0 = _parameters[1];

	double gamma = 0.;
	for(size_t i=0; i<getNrChannels(); ++i) {
		const rpwa::massDepFit::channel& channel = getChannel(i);
		const double ps = channel.getPhaseSpace(idxBin, m, std::numeric_limits<size_t>::max());
		const double ps0 = channel.getPhaseSpace(idxBin, m0, std::numeric_limits<size_t>::max());

		gamma += (ps*ps) / (ps0*ps0);
	}
	gamma *= gamma0 * m0/m;

	return gamma0*m0 / std::complex<double>(m0*m0-m*m, -gamma*m0);
}


std::ostream&
rpwa::massDepFit::parameterizationA1Bowler::print(std::ostream& out) const
{
	out << getName() << std::endl
	    << "Mass=" << _parameters[0] << "   Width=" << _parameters[1] << std::endl;
	return component::print(out);
}


rpwa::massDepFit::exponentialBackground::exponentialBackground(const std::string& name)
	: component(name, 2)
{
	_parametersName[0] = "m0";
	_parametersName[1] = "g";

	_parametersStep[0] = 0.001;
	_parametersStep[1] = 1.0;
}


bool
rpwa::massDepFit::exponentialBackground::init(const libconfig::Setting* configComponent,
                                              const size_t nrBins,
                                              const std::vector<double>& massBinCenters,
                                              const std::map<std::string, size_t>& waveIndices,
                                              const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                              const bool debug)
{
	if(debug) {
		printDebug << "starting initialization of 'exponentialBackground' for component '" << getName() << "'." << std::endl;
	}

	if(not component::init(configComponent, nrBins, massBinCenters, waveIndices, phaseSpaceIntegrals, debug)) {
		printErr << "error while reading configuration of 'component' class." << std::endl;
		return false;
	}

	std::map<std::string, libconfig::Setting::Type> mandatoryArgumentsIsobarMasses;
	boost::assign::insert(mandatoryArgumentsIsobarMasses)
	                     ("mIsobar1", libconfig::Setting::TypeFloat)
	                     ("mIsobar2", libconfig::Setting::TypeFloat);
	if(not checkIfAllVariablesAreThere(configComponent, mandatoryArgumentsIsobarMasses)) {
		printErr << "not all required isobar masses are defined for the component '" << getName() << "'." << std::endl;
		return false;
	}

	configComponent->lookupValue("mIsobar1", _m1);
	configComponent->lookupValue("mIsobar2", _m2);

	if(debug) {
		std::ostringstream output;
		output << "    mass: " << _parameters[0] << " GeV/c^2, ";
		if(_parametersLimitedLower[0] && _parametersLimitedUpper[0]) {
			output << "    limits: " << _parametersLimitLower[0] << "-" << _parametersLimitUpper[0] << " GeV/c^2 ";
		} else if(_parametersLimitedLower[0]) {
			output << "    lower limit: " << _parametersLimitLower[0] << " GeV/c^2 ";
		} else if(_parametersLimitedUpper[0]) {
			output << "    upper limit: " << _parametersLimitUpper[0] << " GeV/c^2 ";
		} else {
			output << "    unlimited ";
		}
		output << (_parametersFixed[0] ? "(FIXED)" : "") << std::endl;

		output << "    width: " << _parameters[1] << " GeV/c^2, ";
		if(_parametersLimitedLower[1] && _parametersLimitedUpper[1]) {
			output << "    limits: " << _parametersLimitLower[1] << "-" << _parametersLimitUpper[1] << " GeV/c^2 ";
		} else if(_parametersLimitedLower[1]) {
			output << "    lower limit: " << _parametersLimitLower[1] << " GeV/c^2 ";
		} else if(_parametersLimitedUpper[1]) {
			output << "    upper limit: " << _parametersLimitUpper[1] << " GeV/c^2 ";
		} else {
			output << "    unlimited ";
		}
		output << (_parametersFixed[1] ? "(FIXED) " : "") << std::endl;

		output << "    mass of isobar 1: " << _m1 << " GeV/c^2, mass of isobar 2: " << _m2 << " GeV/c^2" << std::endl;

		printDebug << "component '" << getName() << "':" << std::endl << output.str();

		printDebug << "finished initialization of 'exponentialBackground'." << std::endl;
	}

	return true;
}


std::complex<double>
rpwa::massDepFit::exponentialBackground::val(const size_t idxBin,
                                             const double m) const
{
	// shift baseline mass
	const double mass = m - _parameters[0];

	// calculate breakup momentum
	if(mass < _m1+_m2) {
		return std::complex<double>(1,0);
	}
	const double q2 = rpwa::breakupMomentumSquared(mass, _m1, _m2);

	return exp(-_parameters[1]*q2);
}


std::ostream&
rpwa::massDepFit::exponentialBackground::print(std::ostream& out) const
{
	out << getName() << std::endl
	    << "Mass=" << _parameters[0] << "   Width=" << _parameters[1] << "    Mass isobar 1=" << _m1 << "    Mass isobar 2=" << _m2 << std::endl;
	return component::print(out);
}


rpwa::massDepFit::tPrimeDependentBackground::tPrimeDependentBackground(const std::string& name)
	: component(name, 5)
{
	_parametersName[0] = "m0";
	_parametersName[1] = "c0";
	_parametersName[2] = "c1";
	_parametersName[3] = "c2";
	_parametersName[4] = "c3";

	_parametersStep[0] = 0.001;
	_parametersStep[1] = 1.0e-3;
	_parametersStep[2] = 1.0;
	_parametersStep[3] = 1.0;
	_parametersStep[4] = 1.0;
}


bool
rpwa::massDepFit::tPrimeDependentBackground::setTPrimeMeans(const std::vector<double> tPrimeMeans)
{
	_tPrimeMeans = tPrimeMeans;

	return true;
}


bool
rpwa::massDepFit::tPrimeDependentBackground::init(const libconfig::Setting* configComponent,
                                                  const size_t nrBins,
                                                  const std::vector<double>& massBinCenters,
                                                  const std::map<std::string, size_t>& waveIndices,
                                                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                                  const bool debug)
{
	if(debug) {
		printDebug << "starting initialization of 'tPrimeDependentBackground' for component '" << getName() << "'." << std::endl;
	}

	if(not component::init(configComponent, nrBins, massBinCenters, waveIndices, phaseSpaceIntegrals, debug)) {
		printErr << "error while reading configuration of 'component' class." << std::endl;
		return false;
	}

	std::map<std::string, libconfig::Setting::Type> mandatoryArgumentsIsobarMasses;
	boost::assign::insert(mandatoryArgumentsIsobarMasses)
	                     ("mIsobar1", libconfig::Setting::TypeFloat)
	                     ("mIsobar2", libconfig::Setting::TypeFloat);
	if(not checkIfAllVariablesAreThere(configComponent, mandatoryArgumentsIsobarMasses)) {
		printErr << "not all required isobar masses are defined for the component '" << getName() << "'." << std::endl;
		return false;
	}

	configComponent->lookupValue("mIsobar1", _m1);
	configComponent->lookupValue("mIsobar2", _m2);

	if(_tPrimeMeans.size() != nrBins) {
		printErr << "array of mean t' value in each bin does not contain the correct number of entries (is: " << _tPrimeMeans.size() << ", expected: " << nrBins << ")." << std::endl;
		return false;
	}

	if(debug) {
		std::ostringstream output;
		output << "    mass threshold: " << _parameters[0] << " GeV/c^2, ";
		if(_parametersLimitedLower[0] && _parametersLimitedUpper[0]) {
			output << "    limits: " << _parametersLimitLower[0] << "-" << _parametersLimitUpper[0] << " GeV/c^2 ";
		} else if(_parametersLimitedLower[0]) {
			output << "    lower limit: " << _parametersLimitLower[0] << " GeV/c^2 ";
		} else if(_parametersLimitedUpper[0]) {
			output << "    upper limit: " << _parametersLimitUpper[0] << " GeV/c^2 ";
		} else {
			output << "    unlimited ";
		}
		output << (_parametersFixed[0] ? "(FIXED)" : "") << std::endl;

		for(size_t i=1; i<_parameters.size(); ++i) {
			output << "    c" << i-1 << ": " << _parameters[i] << " , ";
			if(_parametersLimitedLower[i] && _parametersLimitedUpper[i]) {
				output << "    limits: " << _parametersLimitLower[i] << "-" << _parametersLimitUpper[i] << " GeV/c^2 ";
			} else if(_parametersLimitedLower[i]) {
				output << "    lower limit: " << _parametersLimitLower[i] << " GeV/c^2 ";
			} else if(_parametersLimitedUpper[i]) {
				output << "    upper limit: " << _parametersLimitUpper[i] << " GeV/c^2 ";
			} else {
				output << "    unlimited ";
			}
			output << (_parametersFixed[i] ? "(FIXED) " : "") << std::endl;
		}

		output << "    mass of isobar 1: " << _m1 << " GeV/c^2, mass of isobar 2: " << _m2 << " GeV/c^2" << std::endl;

		output << "    for " << nrBins << " bins with mean t' values: " << _tPrimeMeans[0];
		for(size_t i=1; i<_tPrimeMeans.size(); ++i) {
			output << ", " << _tPrimeMeans[i];
		}
		output << std::endl;

		printDebug << "component '" << getName() << "':" << std::endl << output.str();

		printDebug << "finished initialization of 'tPrimeDependentBackground'." << std::endl;
	}

	return true;
}


std::complex<double>
rpwa::massDepFit::tPrimeDependentBackground::val(const size_t idxBin,
                                                 const double m) const
{
	// calculate breakup momentum
	if(m < _m1+_m2) {
		return std::pow(m - _parameters[0], _parameters[1]);
	}
	const double q2 = rpwa::breakupMomentumSquared(m, _m1, _m2);

	// get mean t' value for current bin
	const double tPrime = _tPrimeMeans[idxBin];

	return std::pow(m - _parameters[0], _parameters[1]) * exp(-(_parameters[2] + _parameters[3]*tPrime + _parameters[4]*tPrime*tPrime)*q2);
}


std::ostream&
rpwa::massDepFit::tPrimeDependentBackground::print(std::ostream& out) const
{
	out << getName() << std::endl
	    << "Mass=" << _parameters[0] << "   Width=" << _parameters[1] << "    Mass isobar 1=" << _m1 << "    Mass isobar 2=" << _m2 << std::endl;
	return component::print(out);
}
