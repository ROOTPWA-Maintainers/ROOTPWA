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
                                   const std::vector<double>& massBinCenters,
                                   const boost::multi_array<double, 2>& phaseSpace)
	: _waveName(waveName),
	  _anchor(false),
	  _nrBins(nrBins),
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
	  _massBinCenters(ch._massBinCenters),
	  _phaseSpace(ch._phaseSpace)
{
	_interpolator.resize(_nrBins);
	for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
		boost::multi_array<double, 2>::const_array_view<1>::type view = _phaseSpace[boost::indices[idxBin][boost::multi_array<double, 2>::index_range()]];
		_interpolator[idxBin] = new ROOT::Math::Interpolator(_massBinCenters, std::vector<double>(view.begin(), view.end()), ROOT::Math::Interpolation::kLINEAR);
	}
}


rpwa::massDepFit::channel::~channel()
{
	for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
		delete _interpolator[idxBin];
	}
}


rpwa::massDepFit::component::component(const size_t id,
                                       const std::string& name,
                                       const size_t nrParameters)
	: _id(id),
	  _name(name),
	  _nrParameters(nrParameters),
	  _nrCouplings(0),
	  _nrBranchings(0),
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
                                  rpwa::massDepFit::parameters& fitParameters,
                                  const size_t nrBins,
                                  const std::vector<double>& massBinCenters,
                                  const std::map<std::string, size_t>& waveIndices,
                                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                  const bool useBranchings,
                                  const bool debug)
{
	if(debug) {
		printDebug << "starting initialization of 'component' for component '" << getName() << "'." << std::endl;
	}

	// resize fitParamters without affecting the number of channels first
	// as this is not yet known
	fitParameters.resize(_id+1, 0, _nrParameters, nrBins);

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

		double parameter;
		configParameter->lookupValue("val", parameter);
		fitParameters.setParameter(getId(), idxParameter, parameter);

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

	// resize fitParamters to the correct number of channels
	fitParameters.resize(_id+1, nrDecayChannels, _nrParameters, nrBins);

	std::map<std::string, size_t> couplingsQN;
	std::map<std::string, size_t> branchingDecay;
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
		for(size_t idxChannel=0; idxChannel<_channels.size(); ++idxChannel) {
			if(_channels[idxChannel].getWaveName() == waveName) {
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
		bool readBranching = false;
		size_t couplingIndex = _nrCouplings;
		size_t branchingIndex = _nrBranchings;
		if(useBranchings && nrDecayChannels > 1) {
			const std::string waveQN = waveName.substr(0, 7);
			const std::string waveDecay = waveName.substr(7);
			if(debug) {
				printDebug << "extracted quantum numbers '" << waveQN << "' and decay chain '" << waveDecay << "' from wave name '" << waveName << "'." << std::endl;
			}

			std::map<std::string, size_t>::const_iterator mappingQN = couplingsQN.find(waveQN);
			std::map<std::string, size_t>::const_iterator mappingDecay = branchingDecay.find(waveDecay);
			if(mappingQN == couplingsQN.end() && mappingDecay == branchingDecay.end()) {
				readCouplings = true;
				couplingsQN[waveQN] = couplingIndex;
				readBranching = true;
				branchingDecay[waveDecay] = branchingIndex;
			} else if(mappingQN == couplingsQN.end()) {
				readCouplings = true;
				couplingsQN[waveQN] = couplingIndex;
				readBranching = false;
				branchingIndex = mappingDecay->second;
			} else if(mappingDecay == branchingDecay.end()) {
				readCouplings = false;
				couplingIndex = mappingQN->second;
				readBranching = true;
				branchingDecay[waveDecay] = branchingIndex;
			} else {
				readCouplings = false;
				couplingIndex = mappingQN->second;
				readBranching = false;
				branchingIndex = mappingDecay->second;
			}
		}

		if(readCouplings) {
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
				fitParameters.setCoupling(getId(), _nrCouplings, idxCoupling, coupling);
			}

			_channelsFromCoupling.push_back(_channels.size());
			++_nrCouplings;
		}

		if(readBranching) {
			boost::assign::insert(mandatoryArguments)
			                     ("branching", libconfig::Setting::TypeGroup);
			if(not checkIfAllVariablesAreThere(decayChannel, mandatoryArguments)) {
				printErr << "one of the decay channels of the component '" << getName() << "' does not contain all required fields." << std::endl;
				return false;
			}

			const libconfig::Setting* configBranching = findLibConfigGroup(*decayChannel, "branching");
			if(not configBranching) {
				printErr << "decay channel '" << waveName << "' of component '" << getName() << "' has no branching." << std::endl;
				return false;
			}

			std::map<std::string, libconfig::Setting::Type> mandatoryArguments;
			boost::assign::insert(mandatoryArguments)
			                     ("branching_Re", libconfig::Setting::TypeFloat)
			                     ("branching_Im", libconfig::Setting::TypeFloat);
			if(not checkIfAllVariablesAreThere(configBranching, mandatoryArguments)) {
				printErr << "branching of the decay channel '" << waveName << "' of the component '" << getName() << "' does not contain all required fields." << std::endl;
				return false;
			}

			double branchingReal;
			configBranching->lookupValue("branching_Re", branchingReal);

			double branchingImag;
			configBranching->lookupValue("branching_Im", branchingImag);

			if(branchingIndex == 0) {
				// the first branching should always be 1.
				if(branchingReal != 1.0 || branchingImag != 0.0) {
					printWarn << "branching of the decay channel '" << waveName << "' of the component '" << getName() << "' forced to 1." << std::endl;
					branchingReal = 1.0;
					branchingImag = 0.0;
				}
			}

			const std::complex<double> branching(branchingReal, branchingImag);
			fitParameters.setBranching(getId(), _nrBranchings, branching);
			++_nrBranchings;
		} else {
			assert(branchingIndex == 0);
			assert(_nrBranchings == 0);

			fitParameters.setBranching(getId(), _nrBranchings, std::complex<double>(1., 0.));
		}

		boost::multi_array<double, 3>::const_array_view<2>::type view = phaseSpaceIntegrals[boost::indices[boost::multi_array<double, 3>::index_range()][boost::multi_array<double, 3>::index_range()][it->second]];
		_channels.push_back(rpwa::massDepFit::channel(waveName, nrBins, massBinCenters, view));
		_channelsCoupling.push_back(couplingIndex);
		_channelsBranching.push_back(branchingIndex);
	}

	if(debug) {
		printDebug << "finished initialization of 'component'." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::component::update(const libconfig::Setting* configComponent,
                                    const rpwa::massDepFit::parameters& fitParameters,
                                    const ROOT::Math::Minimizer* minimizer,
                                    const bool useBranchings,
                                    const bool debug) const
{
	if(debug) {
		printDebug << "starting updating of 'component' for component '" << getName() << "'." << std::endl;
		print(printDebug);
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

		(*configParameter)["val"] = fitParameters.getParameter(getId(), idxParameter);

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
		                     ("amp", libconfig::Setting::TypeString);
		if(not checkIfAllVariablesAreThere(decayChannel, mandatoryArguments)) {
			printErr << "one of the decay channels of the component '" << getName() << "' does not contain all required fields." << std::endl;
			return false;
		}

		std::string waveName;
		decayChannel->lookupValue("amp", waveName);

		const rpwa::massDepFit::channel* channel = NULL;
		size_t channelIdx = 0;
		for(size_t idxChannel=0; idxChannel<getNrChannels(); ++idxChannel) {
			if(getChannelWaveName(idxChannel) == waveName) {
				channel = &getChannel(idxChannel);
				channelIdx = idxChannel;
				break;
			}
		}
		if(not channel) {
			printErr << "could not find channel '" << waveName << "' for component '" << getName() << "'." << std::endl;
			return false;
		}

		if(channelIdx == _channelsCoupling[channelIdx]) {
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
			if(nrCouplings < 0 || static_cast<size_t>(nrCouplings) != channel->getNrBins()) {
				printErr << "decay channel '" << waveName << "' of component '" << getName() << "' has only " << nrCouplings << " couplings, not " << channel->getNrBins() << "." << std::endl;
				return false;
			}

			for(int idxCoupling=0; idxCoupling<nrCouplings; ++idxCoupling) {
				const libconfig::Setting* configCoupling = &((*configCouplings)[idxCoupling]);

				(*configCoupling)["coupling_Re"] = fitParameters.getCoupling(getId(), channelIdx, idxCoupling).real();
				(*configCoupling)["coupling_Im"] = fitParameters.getCoupling(getId(), channelIdx, idxCoupling).imag();
			}
		}

		if(useBranchings && nrDecayChannels > 1 && channelIdx == _channelsBranching[channelIdx]) {
			boost::assign::insert(mandatoryArguments)
			                     ("branching", libconfig::Setting::TypeGroup);
			if(not checkIfAllVariablesAreThere(decayChannel, mandatoryArguments)) {
				printErr << "one of the decay channels of the component '" << getName() << "' does not contain all required fields." << std::endl;
				return false;
			}

			const libconfig::Setting* configBranching = findLibConfigGroup(*decayChannel, "branching");
			if(not configBranching) {
				printErr << "decay channel '" << waveName << "' of component '" << getName() << "' has no branching." << std::endl;
				return false;
			}

			(*configBranching)["branching_Re"] = fitParameters.getBranching(getId(), _channelsBranching[channelIdx]).real();
			(*configBranching)["branching_Im"] = fitParameters.getBranching(getId(), _channelsBranching[channelIdx]).imag();
		}
	}

	if(debug) {
		printDebug << "finished updating of 'component'." << std::endl;
	}

	return true;
}


size_t
rpwa::massDepFit::component::importCouplings(const double* par,
                                             rpwa::massDepFit::parameters& fitParameters)
{
	size_t counter=0;
	for(size_t idxCoupling=0; idxCoupling<_nrCouplings; ++idxCoupling) {
		rpwa::massDepFit::channel& channel = _channels[_channelsFromCoupling[idxCoupling]];
		if(channel.isAnchor()) {
			const size_t nrBins = channel.getNrBins();
			for(size_t idxBin=0; idxBin<nrBins; ++idxBin) {
				fitParameters.setCoupling(getId(), idxCoupling, idxBin, std::complex<double>(par[counter], 0.));
				counter += 1;
			}
		} else {
			const size_t nrBins = channel.getNrBins();
			for(size_t idxBin=0; idxBin<nrBins; ++idxBin) {
				fitParameters.setCoupling(getId(), idxCoupling, idxBin, std::complex<double>(par[counter], par[counter+1]));
				counter += 2;
			}
		}
	}

	return counter;
}


size_t
rpwa::massDepFit::component::importBranchings(const double* par,
                                              rpwa::massDepFit::parameters& fitParameters)
{
	// branching with idx 0 is always real and fixed to 1
	fitParameters.setBranching(getId(), 0, std::complex<double>(1., 0.));

	size_t counter=0;
	for(size_t idxBranching=1; idxBranching<_nrBranchings; ++idxBranching) {
		fitParameters.setBranching(getId(), idxBranching, std::complex<double>(par[counter], par[counter+1]));
		counter += 2;
	}

	return counter;
}


size_t
rpwa::massDepFit::component::importParameters(const double* par,
                                              rpwa::massDepFit::parameters& fitParameters)
{
	for(size_t idxParameter=0; idxParameter<_nrParameters; ++idxParameter) {
		fitParameters.setParameter(getId(), idxParameter, par[idxParameter]);
	}

	return _nrParameters;
}


std::ostream&
rpwa::massDepFit::component::print(std::ostream& out) const
{
	out << "Decay modes:" << std::endl;
	for(size_t idxChannel=0; idxChannel<_channels.size(); ++idxChannel) {
		const rpwa::massDepFit::channel& channel = _channels[idxChannel];
		out << "    " << idxChannel << ", '" << channel.getWaveName() << "'" << std::endl;
	}
	return out;
}


rpwa::massDepFit::fixedWidthBreitWigner::fixedWidthBreitWigner(const size_t id,
                                                               const std::string& name)
	: component(id, name, 2)
{
	_parametersName[0] = "mass";
	_parametersName[1] = "width";

	_parametersStep[0] = 0.001;
	_parametersStep[1] = 0.001;
}


bool
rpwa::massDepFit::fixedWidthBreitWigner::init(const libconfig::Setting* configComponent,
                                              rpwa::massDepFit::parameters& fitParameters,
                                              const size_t nrBins,
                                              const std::vector<double>& massBinCenters,
                                              const std::map<std::string, size_t>& waveIndices,
                                              const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                              const bool useBranchings,
                                              const bool debug)
{
	if(debug) {
		printDebug << "starting initialization of 'fixedWidthBreitWigner' for component '" << getName() << "'." << std::endl;
	}

	if(not component::init(configComponent, fitParameters, nrBins, massBinCenters, waveIndices, phaseSpaceIntegrals, useBranchings, debug)) {
		printErr << "error while reading configuration of 'component' class." << std::endl;
		return false;
	}

	if(debug) {
		print(printDebug);
		printDebug << "finished initialization of 'fixedWidthBreitWigner'." << std::endl;
	}

	return true;
}


std::complex<double>
rpwa::massDepFit::fixedWidthBreitWigner::val(const rpwa::massDepFit::parameters& fitParameters,
                                             const size_t idxBin,
                                             const double m) const
{
	const double& m0 = fitParameters.getParameter(getId(), 0);
	const double& gamma0 = fitParameters.getParameter(getId(), 1);

	return gamma0*m0 / std::complex<double>(m0*m0-m*m, -gamma0*m0);
}


std::ostream&
rpwa::massDepFit::fixedWidthBreitWigner::print(std::ostream& out) const
{
	out << "component " << getId() << " '" << getName() << "' (fixedWidthBreitWigner):" << std::endl;

	out << "    mass ";
	if(_parametersLimitedLower[0] && _parametersLimitedUpper[0]) {
		out << "limits: " << _parametersLimitLower[0] << "-" << _parametersLimitUpper[0] << " GeV/c^2";
	} else if(_parametersLimitedLower[0]) {
		out << "lower limit: " << _parametersLimitLower[0] << " GeV/c^2";
	} else if(_parametersLimitedUpper[0]) {
		out << "upper limit: " << _parametersLimitUpper[0] << " GeV/c^2";
	} else {
		out << "unlimited";
	}
	out << (_parametersFixed[0] ? " (FIXED)" : "") << std::endl;

	out << "    width ";
	if(_parametersLimitedLower[1] && _parametersLimitedUpper[1]) {
		out << "limits: " << _parametersLimitLower[1] << "-" << _parametersLimitUpper[1] << " GeV/c^2";
	} else if(_parametersLimitedLower[1]) {
		out << "lower limit: " << _parametersLimitLower[1] << " GeV/c^2";
	} else if(_parametersLimitedUpper[1]) {
		out << "upper limit: " << _parametersLimitUpper[1] << " GeV/c^2";
	} else {
		out << "unlimited";
	}
	out << (_parametersFixed[1] ? " (FIXED)" : "") << std::endl;

	return component::print(out);
}


rpwa::massDepFit::dynamicWidthBreitWigner::dynamicWidthBreitWigner(const size_t id,
                                                                   const std::string& name)
	: component(id, name, 2)
{
	_parametersName[0] = "mass";
	_parametersName[1] = "width";

	_parametersStep[0] = 0.001;
	_parametersStep[1] = 0.001;
}


bool
rpwa::massDepFit::dynamicWidthBreitWigner::init(const libconfig::Setting* configComponent,
                                                rpwa::massDepFit::parameters& fitParameters,
                                                const size_t nrBins,
                                                const std::vector<double>& massBinCenters,
                                                const std::map<std::string, size_t>& waveIndices,
                                                const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                                const bool useBranchings,
                                                const bool debug)
{
	if(debug) {
		printDebug << "starting initialization of 'dynamicWidthBreitWigner' for component '" << getName() << "'." << std::endl;
	}

	if(not component::init(configComponent, fitParameters, nrBins, massBinCenters, waveIndices, phaseSpaceIntegrals, useBranchings, debug)) {
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
		print(printDebug);
		printDebug << "finished initialization of 'dynamicWidthBreitWigner'." << std::endl;
	}
	return true;
}


std::complex<double>
rpwa::massDepFit::dynamicWidthBreitWigner::val(const rpwa::massDepFit::parameters& fitParameters,
                                               const size_t idxBin,
                                               const double m) const
{
	const double& m0 = fitParameters.getParameter(getId(), 0);
	const double& gamma0 = fitParameters.getParameter(getId(), 1);

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
rpwa::massDepFit::dynamicWidthBreitWigner::print(std::ostream& out) const
{
	out << "component " << getId() << " '" << getName() << "' (dynamicWidthBreitWigner):" << std::endl;

	out << "    mass ";
	if(_parametersLimitedLower[0] && _parametersLimitedUpper[0]) {
		out << "limits: " << _parametersLimitLower[0] << "-" << _parametersLimitUpper[0] << " GeV/c^2";
	} else if(_parametersLimitedLower[0]) {
		out << "lower limit: " << _parametersLimitLower[0] << " GeV/c^2";
	} else if(_parametersLimitedUpper[0]) {
		out << "upper limit: " << _parametersLimitUpper[0] << " GeV/c^2";
	} else {
		out << "unlimited";
	}
	out << (_parametersFixed[0] ? " (FIXED)" : "") << std::endl;

	out << "    width ";
	if(_parametersLimitedLower[1] && _parametersLimitedUpper[1]) {
		out << "limits: " << _parametersLimitLower[1] << "-" << _parametersLimitUpper[1] << " GeV/c^2";
	} else if(_parametersLimitedLower[1]) {
		out << "lower limit: " << _parametersLimitLower[1] << " GeV/c^2";
	} else if(_parametersLimitedUpper[1]) {
		out << "upper limit: " << _parametersLimitUpper[1] << " GeV/c^2";
	} else {
		out << "unlimited";
	}
	out << (_parametersFixed[1] ? " (FIXED)" : "") << std::endl;

	out << "    " << _ratio.size() << " decay channels:" << std::endl;
	for(size_t i=0; i<_ratio.size(); ++i) {
		out << "      * decay channel " << i << ", branching ratio: " << _ratio[i] << std::endl;
		out << "        mass of isobar 1: " << _m1[i] << " GeV/c^2, mass of isobar 2: " << _m2[i] << " GeV/c^2" << std::endl;
		out << "        relative orbital angular momentum between isobars: " << _l[i] << " (in units of hbar)" << std::endl;
	}

	return component::print(out);
}


rpwa::massDepFit::parameterizationA1Bowler::parameterizationA1Bowler(const size_t id,
                                                                     const std::string& name)
	: component(id, name, 2),
	  _interpolator(NULL)
{
	_parametersName[0] = "mass";
	_parametersName[1] = "width";

	_parametersStep[0] = 0.001;
	_parametersStep[1] = 0.001;
}


rpwa::massDepFit::parameterizationA1Bowler::~parameterizationA1Bowler()
{
	if (_interpolator != NULL) {
		delete _interpolator;
	}
}


bool
rpwa::massDepFit::parameterizationA1Bowler::init(const libconfig::Setting* configComponent,
                                                 rpwa::massDepFit::parameters& fitParameters,
                                                 const size_t nrBins,
                                                 const std::vector<double>& massBinCenters,
                                                 const std::map<std::string, size_t>& waveIndices,
                                                 const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                                 const bool useBranchings,
                                                 const bool debug)
{
	if(debug) {
		printDebug << "starting initialization of 'parameterizationA1Bowler' for component '" << getName() << "'." << std::endl;
	}

	if(not component::init(configComponent, fitParameters, nrBins, massBinCenters, waveIndices, phaseSpaceIntegrals, useBranchings, debug)) {
		printErr << "error while reading configuration of 'component' class." << std::endl;
		return false;
	}

	// just to be on the safe side, in principle this should be possible,
	// one should then use an array of interpolators (one for each channel),
	// and probably would have to calculate the involved integrals for each
	// channel
	if(getNrChannels() != 1) {
		printErr << "component of type 'parameterizationA1Bowler' must have exactly one channel." << std::endl;
		return false;
	}

	const libconfig::Setting* integrals = findLibConfigList(*configComponent, "integral");
	if(not integrals) {
		printErr << "component '" << getName() << "' has no phase-space integrals required for the Bowler parameterisation." << std::endl;
		return false;
	}

	const int nrValues = integrals->getLength();

	std::vector<double> masses;
	std::vector<double> values;

	for(int idx=0; idx<nrValues; ++idx) {
		const libconfig::Setting* integral = &((*integrals)[idx]);
		if (integral->getType() != libconfig::Setting::TypeArray) {
			printErr << "phase-space integrals of component '" << getName() << "' has to consist of arrays of two floating point numbers." << std::endl;
			return false;
		}
		if (integral->getLength() != 2) {
			printErr << "phase-space integrals of component '" << getName() << "' has to consist of arrays of two floating point numbers." << std::endl;
			return false;
		}
		if ((*integral)[0].getType() != libconfig::Setting::TypeFloat) {
			printErr << "phase-space integrals of component '" << getName() << "' has to consist of arrays of two floating point numbers." << std::endl;
			return false;
		}

		if (masses.size() > 0 && masses.back() > (double)((*integral)[0])) {
			printErr << "masses of phase-space integrals of component '" << getName() << "' have to be strictly ordered." << std::endl;
			return false;
		}

		masses.push_back((*integral)[0]);
		values.push_back((*integral)[1]);
	}

	assert(_interpolator == NULL); // huh, init called twice?
	_interpolator = new ROOT::Math::Interpolator(masses, values, ROOT::Math::Interpolation::kLINEAR);

	if(debug) {
		print(printDebug);
		printDebug << "finished initialization of 'parameterizationA1Bowler'." << std::endl;
	}
	return true;
}


std::complex<double>
rpwa::massDepFit::parameterizationA1Bowler::val(const rpwa::massDepFit::parameters& fitParameters,
                                                const size_t idxBin,
                                                const double m) const
{
	const double& m0 = fitParameters.getParameter(getId(), 0);
	const double& gamma0 = fitParameters.getParameter(getId(), 1);

	double gamma = 0.;
	for(size_t i=0; i<getNrChannels(); ++i) {
		const double ps = _interpolator->Eval(m);
		const double ps0 = _interpolator->Eval(m0);

		gamma += ps / ps0;
	}
	gamma *= gamma0 * m0/m;

	return gamma0*m0 / std::complex<double>(m0*m0-m*m, -gamma*m0);
}


std::ostream&
rpwa::massDepFit::parameterizationA1Bowler::print(std::ostream& out) const
{
	out << "component " << getId() << " '" << getName() << "' (parameterizationA1Bowler):" << std::endl;

	out << "    mass ";
	if(_parametersLimitedLower[0] && _parametersLimitedUpper[0]) {
		out << "limits: " << _parametersLimitLower[0] << "-" << _parametersLimitUpper[0] << " GeV/c^2";
	} else if(_parametersLimitedLower[0]) {
		out << "lower limit: " << _parametersLimitLower[0] << " GeV/c^2";
	} else if(_parametersLimitedUpper[0]) {
		out << "upper limit: " << _parametersLimitUpper[0] << " GeV/c^2";
	} else {
		out << "unlimited";
	}
	out << (_parametersFixed[0] ? " (FIXED)" : "") << std::endl;

	out << "    width ";
	if(_parametersLimitedLower[1] && _parametersLimitedUpper[1]) {
		out << "limits: " << _parametersLimitLower[1] << "-" << _parametersLimitUpper[1] << " GeV/c^2";
	} else if(_parametersLimitedLower[1]) {
		out << "lower limit: " << _parametersLimitLower[1] << " GeV/c^2";
	} else if(_parametersLimitedUpper[1]) {
		out << "upper limit: " << _parametersLimitUpper[1] << " GeV/c^2";
	} else {
		out << "unlimited";
	}
	out << (_parametersFixed[1] ? " (FIXED)" : "") << std::endl;

	return component::print(out);
}


rpwa::massDepFit::exponentialBackground::exponentialBackground(const size_t id,
                                                               const std::string& name)
	: component(id, name, 2)
{
	_parametersName[0] = "m0";
	_parametersName[1] = "g";

	_parametersStep[0] = 0.001;
	_parametersStep[1] = 1.0;
}


bool
rpwa::massDepFit::exponentialBackground::init(const libconfig::Setting* configComponent,
                                              rpwa::massDepFit::parameters& fitParameters,
                                              const size_t nrBins,
                                              const std::vector<double>& massBinCenters,
                                              const std::map<std::string, size_t>& waveIndices,
                                              const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                              const bool useBranchings,
                                              const bool debug)
{
	if(debug) {
		printDebug << "starting initialization of 'exponentialBackground' for component '" << getName() << "'." << std::endl;
	}

	if(not component::init(configComponent, fitParameters, nrBins, massBinCenters, waveIndices, phaseSpaceIntegrals, useBranchings, debug)) {
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
		print(printDebug);
		printDebug << "finished initialization of 'exponentialBackground'." << std::endl;
	}

	return true;
}


std::complex<double>
rpwa::massDepFit::exponentialBackground::val(const rpwa::massDepFit::parameters& fitParameters,
                                             const size_t idxBin,
                                             const double m) const
{
	// shift baseline mass
	const double mass = m - fitParameters.getParameter(getId(), 0);

	// calculate breakup momentum
	if(mass < _m1+_m2) {
		return std::complex<double>(1,0);
	}
	const double q2 = rpwa::breakupMomentumSquared(mass, _m1, _m2);

	return exp(-fitParameters.getParameter(getId(), 1)*q2);
}


std::ostream&
rpwa::massDepFit::exponentialBackground::print(std::ostream& out) const
{
	out << "component " << getId() << " '" << getName() << "' (exponentialBackground):" << std::endl;

	out << "    mass threshold ";
	if(_parametersLimitedLower[0] && _parametersLimitedUpper[0]) {
		out << "limits: " << _parametersLimitLower[0] << "-" << _parametersLimitUpper[0] << " GeV/c^2";
	} else if(_parametersLimitedLower[0]) {
		out << "lower limit: " << _parametersLimitLower[0] << " GeV/c^2";
	} else if(_parametersLimitedUpper[0]) {
		out << "upper limit: " << _parametersLimitUpper[0] << " GeV/c^2";
	} else {
		out << "unlimited";
	}
	out << (_parametersFixed[0] ? " (FIXED)" : "") << std::endl;

	out << "    width ";
	if(_parametersLimitedLower[1] && _parametersLimitedUpper[1]) {
		out << "limits: " << _parametersLimitLower[1] << "-" << _parametersLimitUpper[1] << " GeV/c^2";
	} else if(_parametersLimitedLower[1]) {
		out << "lower limit: " << _parametersLimitLower[1] << " GeV/c^2";
	} else if(_parametersLimitedUpper[1]) {
		out << "upper limit: " << _parametersLimitUpper[1] << " GeV/c^2";
	} else {
		out << "unlimited";
	}
	out << (_parametersFixed[1] ? " (FIXED) " : "") << std::endl;

	out << "    mass of isobar 1: " << _m1 << " GeV/c^2, mass of isobar 2: " << _m2 << " GeV/c^2" << std::endl;

	return component::print(out);
}


rpwa::massDepFit::tPrimeDependentBackground::tPrimeDependentBackground(const size_t id,
                                                                       const std::string& name)
	: component(id, name, 5)
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
                                                  rpwa::massDepFit::parameters& fitParameters,
                                                  const size_t nrBins,
                                                  const std::vector<double>& massBinCenters,
                                                  const std::map<std::string, size_t>& waveIndices,
                                                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                                  const bool useBranchings,
                                                  const bool debug)
{
	if(debug) {
		printDebug << "starting initialization of 'tPrimeDependentBackground' for component '" << getName() << "'." << std::endl;
	}

	if(not component::init(configComponent, fitParameters, nrBins, massBinCenters, waveIndices, phaseSpaceIntegrals, useBranchings, debug)) {
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
		print(printDebug);
		printDebug << "finished initialization of 'tPrimeDependentBackground'." << std::endl;
	}

	return true;
}


std::complex<double>
rpwa::massDepFit::tPrimeDependentBackground::val(const rpwa::massDepFit::parameters& fitParameters,
                                                 const size_t idxBin,
                                                 const double m) const
{
	// calculate breakup momentum
	if(m < _m1+_m2) {
		return std::pow(m - fitParameters.getParameter(getId(), 0), fitParameters.getParameter(getId(), 1));
	}
	const double q2 = rpwa::breakupMomentumSquared(m, _m1, _m2);

	// get mean t' value for current bin
	const double tPrime = _tPrimeMeans[idxBin];

	return std::pow(m - fitParameters.getParameter(getId(), 0), fitParameters.getParameter(getId(), 1)) * exp(-(fitParameters.getParameter(getId(), 2) + fitParameters.getParameter(getId(), 3)*tPrime + fitParameters.getParameter(getId(), 4)*tPrime*tPrime)*q2);
}


std::ostream&
rpwa::massDepFit::tPrimeDependentBackground::print(std::ostream& out) const
{
	out << "component " << getId() << " '" << getName() << "' (tPrimeDependentBackground):" << std::endl;

	out << "    mass threshold ";
	if(_parametersLimitedLower[0] && _parametersLimitedUpper[0]) {
		out << "limits: " << _parametersLimitLower[0] << "-" << _parametersLimitUpper[0] << " GeV/c^2";
	} else if(_parametersLimitedLower[0]) {
		out << "lower limit: " << _parametersLimitLower[0] << " GeV/c^2";
	} else if(_parametersLimitedUpper[0]) {
		out << "upper limit: " << _parametersLimitUpper[0] << " GeV/c^2";
	} else {
		out << "unlimited";
	}
	out << (_parametersFixed[0] ? " (FIXED)" : "") << std::endl;

	for(size_t i=1; i<5; ++i) {
		out << "    c" << i-1 << " ";
		if(_parametersLimitedLower[i] && _parametersLimitedUpper[i]) {
			out << "limits: " << _parametersLimitLower[i] << "-" << _parametersLimitUpper[i] << " GeV/c^2";
		} else if(_parametersLimitedLower[i]) {
			out << "lower limit: " << _parametersLimitLower[i] << " GeV/c^2";
		} else if(_parametersLimitedUpper[i]) {
			out << "upper limit: " << _parametersLimitUpper[i] << " GeV/c^2";
		} else {
			out << "unlimited";
		}
		out << (_parametersFixed[i] ? " (FIXED) " : "") << std::endl;
	}

	out << "    mass of isobar 1: " << _m1 << " GeV/c^2, mass of isobar 2: " << _m2 << " GeV/c^2" << std::endl;

	out << "    for " << _tPrimeMeans.size() << " bins with mean t' values: " << _tPrimeMeans[0];
	for(size_t i=1; i<_tPrimeMeans.size(); ++i) {
		out << ", " << _tPrimeMeans[i];
	}
	out << std::endl;

	return component::print(out);
}
