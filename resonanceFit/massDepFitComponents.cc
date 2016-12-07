///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010-2012 Sebastian Neubert (TUM)
//    Copyright 2014,2015 Sebastian Uhl (TUM)
//
//    This file is part of ROOTPWA
//
//    ROOTPWA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ROOTPWA is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ROOTPWA.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      implementation of components and channels for the resonance fit
//
//-------------------------------------------------------------------------


#include "massDepFitComponents.h"

#include <boost/assign/std/vector.hpp>

#include <yaml-cpp/yaml.h>

#include "physUtils.hpp"
#include "reportingUtils.hpp"
#include "yamlCppUtils.hpp"


rpwa::massDepFit::channel::channel(const size_t waveIdx,
                                   const std::string& waveName,
                                   const size_t nrBins,
                                   const std::vector<double>& massBinCenters,
                                   const boost::multi_array<double, 2>& phaseSpace)
	: _waveIdx(waveIdx),
	  _waveName(waveName),
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
	: _waveIdx(ch._waveIdx),
	  _waveName(ch._waveName),
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
                                       const std::string& type,
                                       const size_t nrParameters)
	: _id(id),
	  _name(name),
	  _type(type),
	  _nrParameters(nrParameters),
	  _nrCouplings(0),
	  _nrBranchings(0),
	  _parametersStart(nrParameters),
	  _parametersError(nrParameters),
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
rpwa::massDepFit::component::init(const YAML::Node& configComponent,
                                  rpwa::massDepFit::parameters& fitParameters,
                                  rpwa::massDepFit::parameters& fitParametersError,
                                  const size_t nrBins,
                                  const std::vector<double>& massBinCenters,
                                  const std::map<std::string, size_t>& waveIndices,
                                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                  const bool useBranchings,
                                  const bool debug)
{
	if(debug) {
		printDebug << "start initializing 'component' for component '" << getName() << "'." << std::endl;
	}

	// resize fitParamters without affecting the number of channels first
	// as this is not yet known
	fitParameters.resize(_id+1, 0, _nrParameters, nrBins);
	fitParametersError.resize(_id+1, 0, _nrParameters, nrBins);

	for(size_t idxParameter=0; idxParameter<_nrParameters; ++idxParameter) {
		if(debug) {
			printDebug << "reading parameter '" << _parametersName[idxParameter] << "'." << std::endl;
		}

		const YAML::Node& configParameter = configComponent[_parametersName[idxParameter]];
		if(not configParameter) {
			printErr << "component '" << getName() << "' does not define parameter '" << _parametersName[idxParameter] << "'." << std::endl;
			return false;
		}

		std::map<std::string, YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("val", YamlCppUtils::TypeFloat)
		                     ("fix", YamlCppUtils::TypeBoolean);
		if(not checkIfAllVariablesAreThere(configParameter, mandatoryArguments)) {
			printErr << "'" << _parametersName[idxParameter] << "' of component '" << getName() << "' does not contain all required variables." << std::endl;
			return false;
		}

		const double parameter = configParameter["val"].as<double>();
		_parametersStart[idxParameter] = parameter;
		fitParameters.setParameter(getId(), idxParameter, parameter);

		if(configParameter["error"]) {
			if(checkVariableType(configParameter["error"], YamlCppUtils::TypeFloat)) {
				const double error = configParameter["error"].as<double>();
				_parametersError[idxParameter] = error;
				fitParametersError.setParameter(getId(), idxParameter, error);
			} else {
				printErr << "variable 'error' for parameter '" << _parametersName[idxParameter] << "' of component '" << getName() << "' defined, but not a floating point number." << std::endl;
				return false;
			}
		}

		const bool fixed = configParameter["fix"].as<bool>();
		_parametersFixed[idxParameter] = fixed;

		if(configParameter["lower"]) {
			if(checkVariableType(configParameter["lower"], YamlCppUtils::TypeFloat)) {
				_parametersLimitedLower[idxParameter] = true;
				_parametersLimitLower[idxParameter] = configParameter["lower"].as<double>();
			} else {
				printErr << "variable 'lower' for parameter '" << _parametersName[idxParameter] << "' of component '" << getName() << "' defined, but not a floating point number." << std::endl;
				return false;
			}
		} else {
			_parametersLimitedLower[idxParameter] = false;
		}
		if(configParameter["upper"]) {
			if(checkVariableType(configParameter["upper"], YamlCppUtils::TypeFloat)) {
				_parametersLimitedUpper[idxParameter] = true;
				_parametersLimitUpper[idxParameter] = configParameter["upper"].as<double>();
			} else {
				printErr << "variable 'upper' for parameter '" << _parametersName[idxParameter] << "' of component '" << getName() << "' defined, but not a floating point number." << std::endl;
				return false;
			}
		} else {
			_parametersLimitedUpper[idxParameter] = false;
		}

		if(configParameter["step"]) {
			if(checkVariableType(configParameter["step"], YamlCppUtils::TypeFloat)) {
				_parametersStep[idxParameter] = configParameter["step"].as<double>();
			} else {
				printErr << "variable 'step' for parameter '" << _parametersName[idxParameter] << "' of component '" << getName() << "' defined, but not a floating point number." << std::endl;
				return false;
			}
		}
	}

	const YAML::Node& decayChannels = configComponent["decaychannels"];
	if(not decayChannels) {
		printErr << "component '" << getName() << "' has no decay channels." << std::endl;
		return false;
	}

	const size_t nrDecayChannels = decayChannels.size();

	// resize fitParamters to the correct number of channels
	fitParameters.resize(_id+1, nrDecayChannels, _nrParameters, nrBins);
	fitParametersError.resize(_id+1, nrDecayChannels, _nrParameters, nrBins);

	std::map<std::string, size_t> couplingsQN;
	std::map<std::string, size_t> branchingDecay;
	for(size_t idxDecayChannel=0; idxDecayChannel<nrDecayChannels; ++idxDecayChannel) {
		const YAML::Node& decayChannel = decayChannels[idxDecayChannel];

		std::map<std::string, YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("amp", YamlCppUtils::TypeString);
		if(not checkIfAllVariablesAreThere(decayChannel, mandatoryArguments)) {
			printErr << "one of the decay channels of the component '" << getName() << "' does not contain all required variables." << std::endl;
			return false;
		}

		const std::string waveName = decayChannel["amp"].as<std::string>();

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
			// new coupling, so this should be the last coupling up to now
			assert(couplingIndex == _nrCouplings);

			boost::assign::insert(mandatoryArguments)
			                     ("couplings", YamlCppUtils::TypeSequence);
			if(not checkIfAllVariablesAreThere(decayChannel, mandatoryArguments)) {
				printErr << "one of the decay channels of the component '" << getName() << "' does not contain all required variables." << std::endl;
				return false;
			}

			const YAML::Node& configCouplings = decayChannel["couplings"];
			if(not configCouplings) {
				printErr << "decay channel '" << waveName << "' of component '" << getName() << "' has no couplings." << std::endl;
				return false;
			}

			const size_t nrCouplings = configCouplings.size();
			if(nrCouplings != nrBins) {
				printErr << "decay channel '" << waveName << "' of component '" << getName() << "' has " << nrCouplings << " couplings, not " << nrBins << "." << std::endl;
				return false;
			}

			for(size_t idxCoupling=0; idxCoupling<nrCouplings; ++idxCoupling) {
				const YAML::Node& configCoupling = configCouplings[idxCoupling];
				if(not checkVariableType(configCoupling, YamlCppUtils::TypeSequence)) {
					printErr << "one of the couplings of the decay channel '" << waveName << "' of the component '" << getName() << "' is not a YAML sequence." << std::endl;
					return false;
				}
				if(configCoupling.size() != 2) {
					printErr << "one of the couplings of the decay channel '" << waveName << "' of the component '" << getName() << "' does not contain exactly two entries." << std::endl;
					return false;
				}

				if(not checkVariableType(configCoupling[0], YamlCppUtils::TypeFloat)) {
					printErr << "real part of one of the couplings of the decay channel '" << waveName << "' of the component '" << getName() << "' is not a floating point number." << std::endl;
					return false;
				}
				if(not checkVariableType(configCoupling[1], YamlCppUtils::TypeFloat)) {
					printErr << "imaginary part of one of the couplings of the decay channel '" << waveName << "' of the component '" << getName() << "' is not a floating point number." << std::endl;
					return false;
				}

				const double couplingReal = configCoupling[0].as<double>();
				const double couplingImag = configCoupling[1].as<double>();

				const std::complex<double> coupling(couplingReal, couplingImag);
				fitParameters.setCoupling(getId(), couplingIndex, idxCoupling, coupling);
			}

			_channelsFromCoupling.push_back(_channels.size());
			++_nrCouplings;
		}

		if(readBranching) {
			// new branching, so this should be the last branching up to now
			assert(branchingIndex == _nrBranchings);

			boost::assign::insert(mandatoryArguments)
			                     ("branching", YamlCppUtils::TypeSequence);
			if(not checkIfAllVariablesAreThere(decayChannel, mandatoryArguments)) {
				printErr << "one of the decay channels of the component '" << getName() << "' does not contain all required variables." << std::endl;
				return false;
			}

			const YAML::Node& configBranching = decayChannel["branching"];
			if(configBranching.size() != 2) {
				printErr << "branching of the decay channel '" << waveName << "' of the component '" << getName() << "' does not contain exactly two entries." << std::endl;
				return false;
			}

				if(not checkVariableType(configBranching[0], YamlCppUtils::TypeFloat)) {
					printErr << "real part of branching of the decay channel '" << waveName << "' of the component '" << getName() << "' is not a floating point number." << std::endl;
					return false;
				}
				if(not checkVariableType(configBranching[1], YamlCppUtils::TypeFloat)) {
					printErr << "imaginary part of branching of the decay channel '" << waveName << "' of the component '" << getName() << "' is not a floating point number." << std::endl;
					return false;
				}

			double branchingReal = configBranching[0].as<double>();
			double branchingImag = configBranching[1].as<double>();

			if(branchingIndex == 0) {
				// the first branching should always be 1.
				if(branchingReal != 1.0 || branchingImag != 0.0) {
					printWarn << "branching of the decay channel '" << waveName << "' of the component '" << getName() << "' forced to 1." << std::endl;
					branchingReal = 1.0;
					branchingImag = 0.0;
				}
			}

			const std::complex<double> branching(branchingReal, branchingImag);
			fitParameters.setBranching(getId(), branchingIndex, branching);

			_channelsFromBranching.push_back(_channels.size());
			++_nrBranchings;
		}

		boost::multi_array<double, 3>::const_array_view<2>::type view = phaseSpaceIntegrals[boost::indices[boost::multi_array<double, 3>::index_range()][boost::multi_array<double, 3>::index_range()][it->second]];
		_channels.push_back(rpwa::massDepFit::channel(it->second, waveName, nrBins, massBinCenters, view));
		_channelsCoupling.push_back(couplingIndex);
		_channelsBranching.push_back(branchingIndex);

		if (not readDecayChannel(decayChannel, idxDecayChannel, debug)) {
			printErr << "error while reading extra information for decay channel at index " << idxDecayChannel << " for component '" << getName() << "'." << std::endl;
			return false;
		}
	}

	// if no branchings are read at all make sure that there is one equal to 1.
	if (_nrBranchings == 0) {
		fitParameters.setBranching(getId(), 0, std::complex<double>(1., 0.));
	}

	if(debug) {
		printDebug << "finished initializing 'component'." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::component::readDecayChannel(const YAML::Node& /*decayChannel*/,
                                              const size_t /*idxDecayChannel*/,
                                              const bool /*debug*/)
{
	return true;
}


bool
rpwa::massDepFit::component::write(YAML::Emitter& yamlOutput,
                                   const rpwa::massDepFit::parameters& fitParameters,
                                   const rpwa::massDepFit::parameters& fitParametersError,
                                   const bool useBranchings,
                                   const bool debug) const
{
	if(debug) {
		printDebug << "start writing 'component' for component '" << getName() << "'." << std::endl;
	}

	yamlOutput << YAML::Key << "name";
	yamlOutput << YAML::Value << getName();

	yamlOutput << YAML::Key << "type";
	yamlOutput << YAML::Value << getType();

	for(size_t idxParameter=0; idxParameter<_nrParameters; ++idxParameter) {
		if(debug) {
			printDebug << "writing parameter '" << _parametersName[idxParameter] << "'." << std::endl;
		}

		yamlOutput << YAML::Key << _parametersName[idxParameter];
		yamlOutput << YAML::Value;

		yamlOutput << YAML::BeginMap;

		yamlOutput << YAML::Key << "val";
		yamlOutput << YAML::Value << fitParameters.getParameter(getId(), idxParameter);

		yamlOutput << YAML::Key << "error";
		yamlOutput << YAML::Value << fitParametersError.getParameter(getId(), idxParameter);

		if(_parametersLimitedLower[idxParameter]) {
			yamlOutput << YAML::Key << "lower";
			yamlOutput << YAML::Value << _parametersLimitLower[idxParameter];
		}

		if(_parametersLimitedUpper[idxParameter]) {
			yamlOutput << YAML::Key << "upper";
			yamlOutput << YAML::Value << _parametersLimitUpper[idxParameter];
		}

		yamlOutput << YAML::Key << "step";
		yamlOutput << YAML::Value << _parametersStep[idxParameter];

		yamlOutput << YAML::Key << "fix";
		yamlOutput << YAML::Value << _parametersFixed[idxParameter];

		yamlOutput << YAML::EndMap;
	}

	yamlOutput << YAML::Key << "decaychannels";
	yamlOutput << YAML::Value;
	yamlOutput << YAML::BeginSeq;

	const size_t nrDecayChannels = getNrChannels();
	for(size_t idxDecayChannel=0; idxDecayChannel<nrDecayChannels; ++idxDecayChannel) {
		yamlOutput << YAML::BeginMap;

		yamlOutput << YAML::Key << "amp";
		yamlOutput << YAML::Value << getChannel(idxDecayChannel).getWaveName();

		if(idxDecayChannel == _channelsFromCoupling[_channelsCoupling[idxDecayChannel]]) {
			yamlOutput << YAML::Key << "couplings";
			yamlOutput << YAML::Value;
			yamlOutput << YAML::BeginSeq;

			const size_t nrBins = getChannel(idxDecayChannel).getNrBins();
			for(size_t idxBin=0; idxBin<nrBins; ++idxBin) {
				yamlOutput << YAML::Flow;
				yamlOutput << YAML::BeginSeq;

				yamlOutput << fitParameters.getCoupling(getId(), _channelsCoupling[idxDecayChannel], idxBin).real();
				yamlOutput << fitParameters.getCoupling(getId(), _channelsCoupling[idxDecayChannel], idxBin).imag();

				yamlOutput << YAML::EndSeq;
				yamlOutput << YAML::Block;
			}

			yamlOutput << YAML::EndSeq;
		}

		if(useBranchings && nrDecayChannels > 1) {
			if(idxDecayChannel == _channelsFromBranching[_channelsBranching[idxDecayChannel]]) {
				yamlOutput << YAML::Key << "branching";
				yamlOutput << YAML::Value;

				yamlOutput << YAML::Flow;
				yamlOutput << YAML::BeginSeq;

				yamlOutput << fitParameters.getBranching(getId(), _channelsBranching[idxDecayChannel]).real();
				yamlOutput << fitParameters.getBranching(getId(), _channelsBranching[idxDecayChannel]).imag();

				yamlOutput << YAML::EndSeq;
				yamlOutput << YAML::Block;

			}
		}

		writeDecayChannel(yamlOutput, idxDecayChannel, debug);

		yamlOutput << YAML::EndMap;
	}

	yamlOutput << YAML::EndSeq;

	if(debug) {
		printDebug << "finished writing 'component'." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::component::writeDecayChannel(YAML::Emitter& /*yamlOutput*/,
                                               const size_t /*idxDecayChannel*/,
                                               const bool /*debug*/) const
{
	return true;
}


size_t
rpwa::massDepFit::component::importCouplings(const double* par,
                                             rpwa::massDepFit::parameters& fitParameters,
                                             rpwa::massDepFit::cache& cache)
{
	size_t counter=0;
	for(size_t idxCoupling=0; idxCoupling<_nrCouplings; ++idxCoupling) {
		rpwa::massDepFit::channel& channel = _channels[_channelsFromCoupling[idxCoupling]];
		const size_t nrBins = channel.getNrBins();
		for(size_t idxBin=0; idxBin<nrBins; ++idxBin) {
			std::complex<double> coupling;
			if(channel.isAnchor()) {
				coupling = std::complex<double>(par[counter], 0.);
				counter += 1;
			} else {
				coupling = std::complex<double>(par[counter], par[counter+1]);
				counter += 2;
			}

			if (fitParameters.getCoupling(getId(), idxCoupling, idxBin) != coupling) {
				fitParameters.setCoupling(getId(), idxCoupling, idxBin, coupling);

				// invalidate the cache
				for(size_t idxChannel=0; idxChannel<_channels.size(); ++idxChannel) {
					if (_channelsCoupling[idxChannel] == idxCoupling) {
						cache.setCoupling(getId(), idxChannel, idxBin, std::numeric_limits<size_t>::max(), 0.);
						cache.setProdAmp(_channels[idxChannel].getWaveIdx(), idxBin, std::numeric_limits<size_t>::max(), 0.);
					}
				}
			}
		}
	}

	return counter;
}


size_t
rpwa::massDepFit::component::importBranchings(const double* par,
                                              rpwa::massDepFit::parameters& fitParameters,
                                              rpwa::massDepFit::cache& cache)
{
	// branching with idx 0 is always real and fixed to 1
	std::complex<double> branching(1., 0.);
	if (fitParameters.getBranching(getId(), 0) != branching) {
		fitParameters.setBranching(getId(), 0, branching);

		// invalidate the cache
		cache.setCoupling(getId(), 0, std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), 0.);
		cache.setProdAmp(_channels[0].getWaveIdx(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), 0.);
	}

	size_t counter=0;
	for(size_t idxBranching=1; idxBranching<_nrBranchings; ++idxBranching) {
		branching = std::complex<double>(par[counter], par[counter+1]);
		counter += 2;

		if (fitParameters.getBranching(getId(), idxBranching) != branching) {
			fitParameters.setBranching(getId(), idxBranching, branching);

			// invalidate the cache
			for(size_t idxChannel=0; idxChannel<_channels.size(); ++idxChannel) {
				if (_channelsBranching[idxChannel] == idxBranching) {
					cache.setCoupling(getId(), idxChannel, std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), 0.);
					cache.setProdAmp(_channels[idxChannel].getWaveIdx(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), 0.);
				}
			}
		}
	}

	return counter;
}


size_t
rpwa::massDepFit::component::importParameters(const double* par,
                                              rpwa::massDepFit::parameters& fitParameters,
                                              rpwa::massDepFit::cache& cache)
{
	bool invalidateCache = false;
	for(size_t idxParameter=0; idxParameter<_nrParameters; ++idxParameter) {
		if (fitParameters.getParameter(getId(), idxParameter) != par[idxParameter]) {
			fitParameters.setParameter(getId(), idxParameter, par[idxParameter]);
			invalidateCache = true;
		}
	}

	if (invalidateCache) {
		cache.setComponent(getId(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), 0.);
		for(size_t idxChannel=0; idxChannel<_channels.size(); ++idxChannel) {
			cache.setProdAmp(_channels[idxChannel].getWaveIdx(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), 0.);
		}
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
	: component(id, name, "fixedWidthBreitWigner", 2)
{
	_parametersName[0] = "mass";
	_parametersName[1] = "width";

	_parametersStep[0] = 0.001;
	_parametersStep[1] = 0.001;
}


bool
rpwa::massDepFit::fixedWidthBreitWigner::init(const YAML::Node& configComponent,
                                              rpwa::massDepFit::parameters& fitParameters,
                                              rpwa::massDepFit::parameters& fitParametersError,
                                              const size_t nrBins,
                                              const std::vector<double>& massBinCenters,
                                              const std::map<std::string, size_t>& waveIndices,
                                              const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                              const bool useBranchings,
                                              const bool debug)
{
	if(debug) {
		printDebug << "start initializing 'fixedWidthBreitWigner' for component '" << getName() << "'." << std::endl;
	}

	if(not component::init(configComponent, fitParameters, fitParametersError, nrBins, massBinCenters, waveIndices, phaseSpaceIntegrals, useBranchings, debug)) {
		printErr << "error while reading configuration of 'component' class." << std::endl;
		return false;
	}

	if(debug) {
		print(printDebug);
		printDebug << "finished initializing 'fixedWidthBreitWigner'." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::fixedWidthBreitWigner::write(YAML::Emitter& yamlOutput,
                                               const rpwa::massDepFit::parameters& fitParameters,
                                               const rpwa::massDepFit::parameters& fitParametersError,
                                               const bool useBranchings,
                                               const bool debug) const
{
	if(debug) {
		printDebug << "start writing 'fixedWidthBreitWigner' for component '" << getName() << "'." << std::endl;
		print(printDebug);
	}

	yamlOutput << YAML::BeginMap;

	if(not component::write(yamlOutput, fitParameters, fitParametersError, useBranchings, debug)) {
		printErr << "error while writing 'component' part of 'fixedWidthBreitWigner'." << std::endl;
		return false;
	}

	yamlOutput << YAML::EndMap;

	if(debug) {
		printDebug << "finished writing 'fixedWidthBreitWigner'." << std::endl;
	}

	return true;
}


std::complex<double>
rpwa::massDepFit::fixedWidthBreitWigner::val(const rpwa::massDepFit::parameters& fitParameters,
                                             rpwa::massDepFit::cache& cache,
                                             const size_t idxBin,
                                             const double m,
                                             const size_t idxMass) const
{
	if (idxMass != std::numeric_limits<size_t>::max()) {
		const std::complex<double> component = cache.getComponent(getId(), idxBin, idxMass);
		if (component != 0.) {
			return component;
		}
	}

	const double& m0 = fitParameters.getParameter(getId(), 0);
	const double& gamma0 = fitParameters.getParameter(getId(), 1);

	const std::complex<double> component = gamma0*m0 / std::complex<double>(m0*m0-m*m, -gamma0*m0);

	if (idxMass != std::numeric_limits<size_t>::max()) {
		cache.setComponent(getId(), std::numeric_limits<size_t>::max(), idxMass, component);
	}

	return component;
}


std::ostream&
rpwa::massDepFit::fixedWidthBreitWigner::print(std::ostream& out) const
{
	out << "component " << getId() << " '" << getName() << "' (fixedWidthBreitWigner):" << std::endl;

	out << "    mass ";
	out << "start value: " << _parametersStart[0] << " +/- " << _parametersError[0] << " ";
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
	out << "start value: " << _parametersStart[1] << " +/- " << _parametersError[1] << " ";
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
	: component(id, name, "dynamicWidthBreitWigner", 2)
{
	_parametersName[0] = "mass";
	_parametersName[1] = "width";

	_parametersStep[0] = 0.001;
	_parametersStep[1] = 0.001;
}


bool
rpwa::massDepFit::dynamicWidthBreitWigner::init(const YAML::Node& configComponent,
                                                rpwa::massDepFit::parameters& fitParameters,
                                                rpwa::massDepFit::parameters& fitParametersError,
                                                const size_t nrBins,
                                                const std::vector<double>& massBinCenters,
                                                const std::map<std::string, size_t>& waveIndices,
                                                const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                                const bool useBranchings,
                                                const bool debug)
{
	if(debug) {
		printDebug << "start initializing 'dynamicWidthBreitWigner' for component '" << getName() << "'." << std::endl;
	}

	if(not component::init(configComponent, fitParameters, fitParametersError, nrBins, massBinCenters, waveIndices, phaseSpaceIntegrals, useBranchings, debug)) {
		printErr << "error while reading configuration of 'component' class." << std::endl;
		return false;
	}

	const size_t nrDecayChannels = getNrChannels();

	const YAML::Node& extraDecayChannels = configComponent["extradecaychannels"];
	const size_t nrExtraDecayChannels = extraDecayChannels ? extraDecayChannels.size() : 0;
	for(size_t idxExtraDecayChannel=0; idxExtraDecayChannel<nrExtraDecayChannels; ++idxExtraDecayChannel) {
		const YAML::Node& decayChannel = extraDecayChannels[idxExtraDecayChannel];
		if (not readDecayChannel(decayChannel, nrDecayChannels + idxExtraDecayChannel, debug)) {
			printErr << "error while reading extra information for extra decay channel at index " << idxExtraDecayChannel << " for component '" << getName() << "'." << std::endl;
			return false;
		}
	}

	assert(_ratio.size() == _m1.size() && _ratio.size() == _m2.size() && _ratio.size() == _l.size());
	if (_ratio.size() != nrDecayChannels + nrExtraDecayChannels) {
		printErr << "inconsistent vector size for decay channel information of component '" << getName() << "'." << std::endl;
		return false;
	}

	if(nrDecayChannels + nrExtraDecayChannels > 1) {
		double sum = 0.;
		for(size_t i=0; i<_ratio.size(); ++i) {
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
		printDebug << "finished initializing 'dynamicWidthBreitWigner'." << std::endl;
	}
	return true;
}


bool
rpwa::massDepFit::dynamicWidthBreitWigner::readDecayChannel(const YAML::Node& decayChannel,
                                                            const size_t idxDecayChannel,
                                                            const bool debug)
{
	if(debug) {
		printDebug << "start initializing decay channel at index " << idxDecayChannel << " in 'dynamicWidthBreitWigner' for component '" << getName() << "'." << std::endl;
	}

	// make sure there are no holes in the vectors
	assert(_ratio.size() == _m1.size() && _ratio.size() == _m2.size() && _ratio.size() == _l.size());
	if (_ratio.size() != idxDecayChannel) {
		printErr << "wrong index for decay channel to read for component '" << getName() << "', got " << idxDecayChannel << ", expected " << _ratio.size() << "." << std::endl;
		return false;
	}

	std::map<std::string, YamlCppUtils::Type> mandatoryArguments;
	boost::assign::insert(mandatoryArguments)
	                     ("mIsobar1", YamlCppUtils::TypeFloat)
	                     ("mIsobar2", YamlCppUtils::TypeFloat)
	                     ("relAngularMom", YamlCppUtils::TypeInt)
	                     ("branchingRatio", YamlCppUtils::TypeFloat);
	if(not checkIfAllVariablesAreThere(decayChannel, mandatoryArguments)) {
		printErr << "decay channel at index " << idxDecayChannel << " of component '" << getName() << "' does not contain all required variables." << std::endl;
		return false;
	}

	_m1.push_back(decayChannel["mIsobar1"].as<double>());
	_m2.push_back(decayChannel["mIsobar2"].as<double>());
	_l.push_back(decayChannel["relAngularMom"].as<int>());
	_ratio.push_back(decayChannel["branchingRatio"].as<double>());

	return true;
}


bool
rpwa::massDepFit::dynamicWidthBreitWigner::write(YAML::Emitter& yamlOutput,
                                                 const rpwa::massDepFit::parameters& fitParameters,
                                                 const rpwa::massDepFit::parameters& fitParametersError,
                                                 const bool useBranchings,
                                                 const bool debug) const
{
	if(debug) {
		printDebug << "start writing 'dynamicWidthBreitWigner' for component '" << getName() << "'." << std::endl;
		print(printDebug);
	}

	yamlOutput << YAML::BeginMap;

	if(not component::write(yamlOutput, fitParameters, fitParametersError, useBranchings, debug)) {
		printErr << "error while writing 'component' part of 'dynamicWidthBreitWigner'." << std::endl;
		return false;
	}

	if (_ratio.size() > getNrChannels()) {
		yamlOutput << YAML::Key << "extradecaychannels";
		yamlOutput << YAML::Value;
		yamlOutput << YAML::BeginSeq;

		const size_t nrDecayChannels = _ratio.size();
		for (size_t idxDecayChannel=getNrChannels(); idxDecayChannel<nrDecayChannels; ++idxDecayChannel) {
			yamlOutput << YAML::BeginMap;
			writeDecayChannel(yamlOutput, idxDecayChannel, debug);
			yamlOutput << YAML::EndMap;
		}

		yamlOutput << YAML::EndSeq;
	}

	yamlOutput << YAML::EndMap;

	if(debug) {
		printDebug << "finished writing 'dynamicWidthBreitWigner'." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::dynamicWidthBreitWigner::writeDecayChannel(YAML::Emitter& yamlOutput,
                                                             const size_t idxDecayChannel,
                                                             const bool debug) const
{
	if(debug) {
		printDebug << "start writing decay channel at index " << idxDecayChannel << " in 'dynamicWidthBreitWigner' for component '" << getName() << "'." << std::endl;
	}

	yamlOutput << YAML::Key << "mIsobar1";
	yamlOutput << YAML::Value << _m1[idxDecayChannel];

	yamlOutput << YAML::Key << "mIsobar2";
	yamlOutput << YAML::Value <<_m2[idxDecayChannel];

	yamlOutput << YAML::Key << "relAngularMom";
	yamlOutput << YAML::Value <<_l[idxDecayChannel];

	yamlOutput << YAML::Key << "branchingRatio";
	yamlOutput << YAML::Value <<_ratio[idxDecayChannel];

	return true;
}


std::complex<double>
rpwa::massDepFit::dynamicWidthBreitWigner::val(const rpwa::massDepFit::parameters& fitParameters,
                                               rpwa::massDepFit::cache& cache,
                                               const size_t idxBin,
                                               const double m,
                                               const size_t idxMass) const
{
	if (idxMass != std::numeric_limits<size_t>::max()) {
		const std::complex<double> component = cache.getComponent(getId(), idxBin, idxMass);
		if (component != 0.) {
			return component;
		}
	}

	const double& m0 = fitParameters.getParameter(getId(), 0);
	const double& gamma0 = fitParameters.getParameter(getId(), 1);

	double gamma = 0.;
	for(size_t i=0; i<_ratio.size(); ++i) {
		// calculation of the break-up momentum will fail if the
		// sum of the two isobar masses is heavier than the mother
		// mass, which is probably a good thing as some more advanced
		// stuff should be used for at- or sub-threshold decays.
		// but if the branching ratio for this channel is 0, then
		// simply ignore this channel
		if (_ratio[i] == 0.) {
			continue;
		}

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

	const std::complex<double> component = gamma0*m0 / std::complex<double>(m0*m0-m*m, -gamma*m0);

	if (idxMass != std::numeric_limits<size_t>::max()) {
		cache.setComponent(getId(), std::numeric_limits<size_t>::max(), idxMass, component);
	}

	return component;
}


std::ostream&
rpwa::massDepFit::dynamicWidthBreitWigner::print(std::ostream& out) const
{
	out << "component " << getId() << " '" << getName() << "' (dynamicWidthBreitWigner):" << std::endl;

	out << "    mass ";
	out << "start value: " << _parametersStart[0] << " +/- " << _parametersError[0] << " ";
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
	out << "start value: " << _parametersStart[1] << " +/- " << _parametersError[1] << " ";
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


rpwa::massDepFit::integralWidthBreitWigner::integralWidthBreitWigner(const size_t id,
                                                                     const std::string& name)
	: component(id, name, "integralWidthBreitWigner", 2)
{
	_parametersName[0] = "mass";
	_parametersName[1] = "width";

	_parametersStep[0] = 0.001;
	_parametersStep[1] = 0.001;
}


rpwa::massDepFit::integralWidthBreitWigner::~integralWidthBreitWigner()
{
	for (size_t i=0; i<_interpolator.size(); ++i) {
		if (_interpolator[i] != NULL) {
			delete _interpolator[i];
		}
	}
	_interpolator.clear();
}


bool
rpwa::massDepFit::integralWidthBreitWigner::init(const YAML::Node& configComponent,
                                                 rpwa::massDepFit::parameters& fitParameters,
                                                 rpwa::massDepFit::parameters& fitParametersError,
                                                 const size_t nrBins,
                                                 const std::vector<double>& massBinCenters,
                                                 const std::map<std::string, size_t>& waveIndices,
                                                 const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                                 const bool useBranchings,
                                                 const bool debug)
{
	if(debug) {
		printDebug << "start initializing 'integralWidthBreitWigner' for component '" << getName() << "'." << std::endl;
	}

	if(not component::init(configComponent, fitParameters, fitParametersError, nrBins, massBinCenters, waveIndices, phaseSpaceIntegrals, useBranchings, debug)) {
		printErr << "error while reading configuration of 'component' class." << std::endl;
		return false;
	}

	const size_t nrDecayChannels = getNrChannels();

	const YAML::Node& extraDecayChannels = configComponent["extradecaychannels"];
	const size_t nrExtraDecayChannels = extraDecayChannels ? extraDecayChannels.size() : 0;
	for(size_t idxExtraDecayChannel=0; idxExtraDecayChannel<nrExtraDecayChannels; ++idxExtraDecayChannel) {
		const YAML::Node& decayChannel = extraDecayChannels[idxExtraDecayChannel];
		if (not readDecayChannel(decayChannel, nrDecayChannels + idxExtraDecayChannel, debug)) {
			printErr << "error while reading extra information for extra decay channel at index " << idxExtraDecayChannel << " for component '" << getName() << "'." << std::endl;
			return false;
		}
	}

	assert(_ratio.size() == _masses.size() && _ratio.size() == _values.size() && _ratio.size() == _interpolator.size());
	if (_ratio.size() != nrDecayChannels + nrExtraDecayChannels) {
		printErr << "inconsistent vector size for decay channel information of component '" << getName() << "'." << std::endl;
		return false;
	}

	if(nrDecayChannels + nrExtraDecayChannels > 1) {
		double sum = 0.;
		for(size_t i=0; i<_ratio.size(); ++i) {
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
		printDebug << "finished initializing 'integralWidthBreitWigner'." << std::endl;
	}
	return true;
}


bool
rpwa::massDepFit::integralWidthBreitWigner::readDecayChannel(const YAML::Node& decayChannel,
                                                             const size_t idxDecayChannel,
                                                             const bool debug)
{
	if(debug) {
		printDebug << "start initializing decay channel at index " << idxDecayChannel << " in 'integralWidthBreitWigner' for component '" << getName() << "'." << std::endl;
	}

	// make sure there are no holes in the vectors
	assert(_ratio.size() == _masses.size() && _ratio.size() == _values.size() && _ratio.size() == _interpolator.size());
	if (_ratio.size() != idxDecayChannel) {
		printErr << "wrong index for decay channel to read for component '" << getName() << "', got " << idxDecayChannel << ", expected " << _ratio.size() << "." << std::endl;
		return false;
	}

	std::map<std::string, YamlCppUtils::Type> mandatoryArguments;
	boost::assign::insert(mandatoryArguments)
	                     ("integral", YamlCppUtils::TypeSequence)
	                     ("branchingRatio", YamlCppUtils::TypeFloat);
	if(not checkIfAllVariablesAreThere(decayChannel, mandatoryArguments)) {
		printErr << "decay channel at index " << idxDecayChannel << " of component '" << getName() << "' does not contain all required variables." << std::endl;
		return false;
	}

	const YAML::Node& integrals = decayChannel["integral"];

	const size_t nrValues = integrals.size();
	if(nrValues < 2) {
		printErr << "phase-space integrals of component '" << getName() << "' has to contain at least two points." << std::endl;
		return false;
	}

	std::vector<double> masses;
	std::vector<double> values;

	for(size_t idx=0; idx<nrValues; ++idx) {
		const YAML::Node& integral = integrals[idx];
		if(not integral.IsSequence()) {
			printErr << "phase-space integrals of component '" << getName() << "' has to consist of arrays of two floating point numbers." << std::endl;
			return false;
		}
		if(integral.size() != 2) {
			printErr << "phase-space integrals of component '" << getName() << "' has to consist of arrays of two floating point numbers." << std::endl;
			return false;
		}
		if(not checkVariableType(integral[0], YamlCppUtils::TypeFloat)) {
			printErr << "phase-space integrals of component '" << getName() << "' has to consist of arrays of two floating point numbers." << std::endl;
			return false;
		}
		if(not checkVariableType(integral[1], YamlCppUtils::TypeFloat)) {
			printErr << "phase-space integrals of component '" << getName() << "' has to consist of arrays of two floating point numbers." << std::endl;
			return false;
		}

		const double mass = integral[0].as<double>();
		const double val = integral[1].as<double>();

		if (masses.size() > 0 && masses.back() > mass) {
			printErr << "masses of phase-space integrals of component '" << getName() << "' have to be strictly ordered." << std::endl;
			return false;
		}

		masses.push_back(mass);
		values.push_back(val);
	}

	_masses.push_back(masses);
	_values.push_back(values);
	_interpolator.push_back(new ROOT::Math::Interpolator(masses, values, ROOT::Math::Interpolation::kLINEAR));
	_ratio.push_back(decayChannel["branchingRatio"].as<double>());

	return true;
}


bool
rpwa::massDepFit::integralWidthBreitWigner::write(YAML::Emitter& yamlOutput,
                                                  const rpwa::massDepFit::parameters& fitParameters,
                                                  const rpwa::massDepFit::parameters& fitParametersError,
                                                  const bool useBranchings,
                                                  const bool debug) const
{
	if(debug) {
		printDebug << "start writing 'integralWidthBreitWigner' for component '" << getName() << "'." << std::endl;
		print(printDebug);
	}

	yamlOutput << YAML::BeginMap;

	if(not component::write(yamlOutput, fitParameters, fitParametersError, useBranchings, debug)) {
		printErr << "error while writing 'component' part of 'integralWidthBreitWigner'." << std::endl;
		return false;
	}

	if (_ratio.size() > getNrChannels()) {
		yamlOutput << YAML::Key << "extradecaychannels";
		yamlOutput << YAML::Value;
		yamlOutput << YAML::BeginSeq;

		const size_t nrDecayChannels = _ratio.size();
		for (size_t idxDecayChannel=getNrChannels(); idxDecayChannel<nrDecayChannels; ++idxDecayChannel) {
			yamlOutput << YAML::BeginMap;
			writeDecayChannel(yamlOutput, idxDecayChannel, debug);
			yamlOutput << YAML::EndMap;
		}

		yamlOutput << YAML::EndSeq;
	}

	yamlOutput << YAML::EndMap;

	if(debug) {
		printDebug << "finished writing 'integralWidthBreitWigner'." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::integralWidthBreitWigner::writeDecayChannel(YAML::Emitter& yamlOutput,
                                                              const size_t idxDecayChannel,
                                                              const bool debug) const
{
	if(debug) {
		printDebug << "start writing decay channel at index " << idxDecayChannel << " in 'integralWidthBreitWigner' for component '" << getName() << "'." << std::endl;
	}

	yamlOutput << YAML::Key << "integral";
	yamlOutput << YAML::Value;

	yamlOutput << YAML::Flow;
	yamlOutput << YAML::BeginSeq;
	for (size_t idx=0; idx<_masses[idxDecayChannel].size(); ++idx) {
		yamlOutput << YAML::BeginSeq;
		yamlOutput << _masses[idxDecayChannel][idx];
		yamlOutput << _values[idxDecayChannel][idx];
		yamlOutput << YAML::EndSeq;
	}
	yamlOutput << YAML::EndSeq;
	yamlOutput << YAML::Block;

	yamlOutput << YAML::Key << "branchingRatio";
	yamlOutput << YAML::Value <<_ratio[idxDecayChannel];

	return true;
}


std::complex<double>
rpwa::massDepFit::integralWidthBreitWigner::val(const rpwa::massDepFit::parameters& fitParameters,
                                                rpwa::massDepFit::cache& cache,
                                                const size_t idxBin,
                                                const double m,
                                                const size_t idxMass) const
{
	if (idxMass != std::numeric_limits<size_t>::max()) {
		const std::complex<double> component = cache.getComponent(getId(), idxBin, idxMass);
		if (component != 0.) {
			return component;
		}
	}

	const double& m0 = fitParameters.getParameter(getId(), 0);
	const double& gamma0 = fitParameters.getParameter(getId(), 1);

	double gamma = 0.;
	for(size_t i=0; i<_ratio.size(); ++i) {
		// save some time not calculating stuff that is ignored
		if (_ratio[i] == 0.) {
			continue;
		}

		const double ps = _interpolator[i]->Eval(m);
		const double ps0 = _interpolator[i]->Eval(m0);

		gamma += _ratio[i] * ps / ps0;
	}
	gamma *= gamma0;

	const std::complex<double> component = gamma0*m0 / std::complex<double>(m0*m0-m*m, -gamma*m0);

	if (idxMass != std::numeric_limits<size_t>::max()) {
		cache.setComponent(getId(), std::numeric_limits<size_t>::max(), idxMass, component);
	}

	return component;
}


std::ostream&
rpwa::massDepFit::integralWidthBreitWigner::print(std::ostream& out) const
{
	out << "component " << getId() << " '" << getName() << "' (integralWidthBreitWigner):" << std::endl;

	out << "    mass ";
	out << "start value: " << _parametersStart[0] << " +/- " << _parametersError[0] << " ";
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
	out << "start value: " << _parametersStart[1] << " +/- " << _parametersError[1] << " ";
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


rpwa::massDepFit::constantBackground::constantBackground(const size_t id,
                                                         const std::string& name)
	: component(id, name, "constantBackground", 0)
{
}


bool
rpwa::massDepFit::constantBackground::init(const YAML::Node& configComponent,
                                           rpwa::massDepFit::parameters& fitParameters,
                                           rpwa::massDepFit::parameters& fitParametersError,
                                           const size_t nrBins,
                                           const std::vector<double>& massBinCenters,
                                           const std::map<std::string, size_t>& waveIndices,
                                           const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                           const bool useBranchings,
                                           const bool debug)
{
	if(debug) {
		printDebug << "start initializing 'constantBackground' for component '" << getName() << "'." << std::endl;
	}

	if(not component::init(configComponent, fitParameters, fitParametersError, nrBins, massBinCenters, waveIndices, phaseSpaceIntegrals, useBranchings, debug)) {
		printErr << "error while reading configuration of 'component' class." << std::endl;
		return false;
	}

	if(getNrChannels() != 1) {
		printErr << "component '" << getName() << "' of type 'constantBackground' needs to have exactly one decay channel." << std::endl;
		return false;
	}

	if(debug) {
		print(printDebug);
		printDebug << "finished initializing 'constantBackground'." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::constantBackground::write(YAML::Emitter& yamlOutput,
                                            const rpwa::massDepFit::parameters& fitParameters,
                                            const rpwa::massDepFit::parameters& fitParametersError,
                                            const bool useBranchings,
                                            const bool debug) const
{
	if(debug) {
		printDebug << "start writing 'constantBackground' for component '" << getName() << "'." << std::endl;
		print(printDebug);
	}

	yamlOutput << YAML::BeginMap;

	if(not component::write(yamlOutput, fitParameters, fitParametersError, useBranchings, debug)) {
		printErr << "error while writing 'component' part of 'constantBackground'." << std::endl;
		return false;
	}

	yamlOutput << YAML::EndMap;

	if(debug) {
		printDebug << "finished writing 'constantBackground'." << std::endl;
	}

	return true;
}


std::complex<double>
rpwa::massDepFit::constantBackground::val(const rpwa::massDepFit::parameters& /*fitParameters*/,
                                          rpwa::massDepFit::cache& /*cache*/,
                                          const size_t /*idxBin*/,
                                          const double /*m*/,
                                          const size_t /*idxMass*/) const
{
	return 1.;
}


std::ostream&
rpwa::massDepFit::constantBackground::print(std::ostream& out) const
{
	out << "component " << getId() << " '" << getName() << "' (constantBackground):" << std::endl;

	return component::print(out);
}


rpwa::massDepFit::exponentialBackground::exponentialBackground(const size_t id,
                                                               const std::string& name)
	: component(id, name, "exponentialBackground", 1)
{
	_parametersName[0] = "g";

	_parametersStep[0] = 1.0;
}


bool
rpwa::massDepFit::exponentialBackground::init(const YAML::Node& configComponent,
                                              rpwa::massDepFit::parameters& fitParameters,
                                              rpwa::massDepFit::parameters& fitParametersError,
                                              const size_t nrBins,
                                              const std::vector<double>& massBinCenters,
                                              const std::map<std::string, size_t>& waveIndices,
                                              const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                              const bool useBranchings,
                                              const bool debug)
{
	if(debug) {
		printDebug << "start initializing 'exponentialBackground' for component '" << getName() << "'." << std::endl;
	}

	if(not component::init(configComponent, fitParameters, fitParametersError, nrBins, massBinCenters, waveIndices, phaseSpaceIntegrals, useBranchings, debug)) {
		printErr << "error while reading configuration of 'component' class." << std::endl;
		return false;
	}

	if(getNrChannels() != 1) {
		printErr << "component '" << getName() << "' of type 'exponentialBackground' needs to have exactly one decay channel." << std::endl;
		return false;
	}

	std::map<std::string, YamlCppUtils::Type> mandatoryArgumentsIsobarMasses;
	boost::assign::insert(mandatoryArgumentsIsobarMasses)
	                     ("mIsobar1", YamlCppUtils::TypeFloat)
	                     ("mIsobar2", YamlCppUtils::TypeFloat);
	if(not checkIfAllVariablesAreThere(configComponent, mandatoryArgumentsIsobarMasses)) {
		printErr << "not all required isobar masses are defined for the component '" << getName() << "'." << std::endl;
		return false;
	}

	_m1 = configComponent["mIsobar1"].as<double>();
	_m2 = configComponent["mIsobar2"].as<double>();

	if (configComponent["relAngularMom"]) {
		if (checkVariableType(configComponent["relAngularMom"], YamlCppUtils::TypeInt)) {
			_l = configComponent["relAngularMom"].as<int>();
		} else {
			printErr << "variable 'relAngularMom' for component '" << getName() << "' defined, but not an integer." << std::endl;
			return false;
		}
	} else {
		printInfo << "variable 'relAngularMom' for component '" << getName() << "' not defined, using default value 0." << std::endl;
		_l = 0;
	}

	if (configComponent["exponent"]) {
		if (checkVariableType(configComponent["exponent"], YamlCppUtils::TypeFloat)) {
			_exponent = configComponent["exponent"].as<double>();
		} else {
			printErr << "variable 'exponent' for component '" << getName() << "' defined, but not a floating point number." << std::endl;
			return false;
		}
	} else {
		printInfo << "variable 'exponent' for component '" << getName() << "' not defined, using default value 2.0." << std::endl;
		_exponent = 2.0;
	}

	const double q = rpwa::breakupMomentum(massBinCenters.back(), _m1, _m2);
	const double f2 = rpwa::barrierFactorSquared(2*_l, q);
	_norm = 1. / (q*f2);

	if(debug) {
		print(printDebug);
		printDebug << "finished initializing 'exponentialBackground'." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::exponentialBackground::write(YAML::Emitter& yamlOutput,
                                               const rpwa::massDepFit::parameters& fitParameters,
                                               const rpwa::massDepFit::parameters& fitParametersError,
                                               const bool useBranchings,
                                               const bool debug) const
{
	if(debug) {
		printDebug << "start writing 'exponentialBackground' for component '" << getName() << "'." << std::endl;
		print(printDebug);
	}

	yamlOutput << YAML::BeginMap;

	if(not component::write(yamlOutput, fitParameters, fitParametersError, useBranchings, debug)) {
		printErr << "error while writing 'component' part of 'exponentialBackground'." << std::endl;
		return false;
	}

	yamlOutput << YAML::Key << "mIsobar1";
	yamlOutput << YAML::Value << _m1;

	yamlOutput << YAML::Key << "mIsobar2";
	yamlOutput << YAML::Value << _m2;

	yamlOutput << YAML::Key << "relAngularMom";
	yamlOutput << YAML::Value <<_l;

	yamlOutput << YAML::Key << "exponent";
	yamlOutput << YAML::Value <<_exponent;

	yamlOutput << YAML::EndMap;

	if(debug) {
		printDebug << "finished writing 'exponentialBackground'." << std::endl;
	}

	return true;
}


std::complex<double>
rpwa::massDepFit::exponentialBackground::val(const rpwa::massDepFit::parameters& fitParameters,
                                             rpwa::massDepFit::cache& cache,
                                             const size_t idxBin,
                                             const double m,
                                             const size_t idxMass) const
{
	if (idxMass != std::numeric_limits<size_t>::max()) {
		const std::complex<double> component = cache.getComponent(getId(), idxBin, idxMass);
		if (component != 0.) {
			return component;
		}
	}

	// calculate breakup momentum
	if(m < _m1+_m2) {
		return std::complex<double>(1,0);
	}
	const double q = rpwa::breakupMomentum(m, _m1, _m2);
	const double f2 = rpwa::barrierFactorSquared(2*_l, q);
	const double c = std::pow(q*f2 * _norm, _exponent);

	const std::complex<double> component = exp(-fitParameters.getParameter(getId(), 0)*c);

	if (idxMass != std::numeric_limits<size_t>::max()) {
		cache.setComponent(getId(), std::numeric_limits<size_t>::max(), idxMass, component);
	}

	return component;
}


std::ostream&
rpwa::massDepFit::exponentialBackground::print(std::ostream& out) const
{
	out << "component " << getId() << " '" << getName() << "' (exponentialBackground):" << std::endl;

	out << "    width ";
	out << "start value: " << _parametersStart[0] << " +/- " << _parametersError[0] << " ";
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

	out << "    mass of isobar 1: " << _m1 << " GeV/c^2, mass of isobar 2: " << _m2 << " GeV/c^2" << std::endl;
	out << "    relative orbital angular momentum between isobars: " << _l << " (in units of hbar)" << std::endl;
	out << "    exponent of break-up momentum times barrier-factor squared: " << _exponent << std::endl;

	return component::print(out);
}


rpwa::massDepFit::tPrimeDependentBackground::tPrimeDependentBackground(const size_t id,
                                                                       const std::string& name)
	: component(id, name, "tPrimeDependentBackground", 5)
{
	_parametersName[0] = "m0";
	_parametersName[1] = "c0";
	_parametersName[2] = "c1";
	_parametersName[3] = "c2";
	_parametersName[4] = "c3";

	_parametersStep[0] = 0.001;
	_parametersStep[1] = 0.001;
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
rpwa::massDepFit::tPrimeDependentBackground::init(const YAML::Node& configComponent,
                                                  rpwa::massDepFit::parameters& fitParameters,
                                                  rpwa::massDepFit::parameters& fitParametersError,
                                                  const size_t nrBins,
                                                  const std::vector<double>& massBinCenters,
                                                  const std::map<std::string, size_t>& waveIndices,
                                                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                                  const bool useBranchings,
                                                  const bool debug)
{
	if(debug) {
		printDebug << "start initializing 'tPrimeDependentBackground' for component '" << getName() << "'." << std::endl;
	}

	if(not component::init(configComponent, fitParameters, fitParametersError, nrBins, massBinCenters, waveIndices, phaseSpaceIntegrals, useBranchings, debug)) {
		printErr << "error while reading configuration of 'component' class." << std::endl;
		return false;
	}

	if(getNrChannels() != 1) {
		printErr << "component '" << getName() << "' of type 'tPrimeDependentBackground' needs to have exactly one decay channel." << std::endl;
		return false;
	}

	std::map<std::string, YamlCppUtils::Type> mandatoryArgumentsIsobarMasses;
	boost::assign::insert(mandatoryArgumentsIsobarMasses)
	                     ("mIsobar1", YamlCppUtils::TypeFloat)
	                     ("mIsobar2", YamlCppUtils::TypeFloat);
	if(not checkIfAllVariablesAreThere(configComponent, mandatoryArgumentsIsobarMasses)) {
		printErr << "not all required isobar masses are defined for the component '" << getName() << "'." << std::endl;
		return false;
	}

	_m1 = configComponent["mIsobar1"].as<double>();
	_m2 = configComponent["mIsobar2"].as<double>();

	if (configComponent["relAngularMom"]) {
		if (checkVariableType(configComponent["relAngularMom"], YamlCppUtils::TypeInt)) {
			_l = configComponent["relAngularMom"].as<int>();
		} else {
			printErr << "variable 'relAngularMom' for component '" << getName() << "' defined, but not an integer." << std::endl;
			return false;
		}
	} else {
		printInfo << "variable 'relAngularMom' for component '" << getName() << "' not defined, using default value 0." << std::endl;
		_l = 0;
	}

	if (configComponent["exponent"]) {
		if (checkVariableType(configComponent["exponent"], YamlCppUtils::TypeFloat)) {
			_exponent = configComponent["exponent"].as<double>();
		} else {
			printErr << "variable 'exponent' for component '" << getName() << "' defined, but not a floating point number." << std::endl;
			return false;
		}
	} else {
		printInfo << "variable 'exponent' for component '" << getName() << "' not defined, using default value 2.0." << std::endl;
		_exponent = 2.0;
	}

	const double q = rpwa::breakupMomentum(massBinCenters.back(), _m1, _m2);
	const double f2 = rpwa::barrierFactorSquared(2*_l, q);
	_norm = 1. / (q*f2);

	if(_tPrimeMeans.size() != nrBins) {
		printErr << "array of mean t' value in each bin does not contain the correct number of entries (is: " << _tPrimeMeans.size() << ", expected: " << nrBins << ")." << std::endl;
		return false;
	}

	if(debug) {
		print(printDebug);
		printDebug << "finished initializing 'tPrimeDependentBackground'." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::tPrimeDependentBackground::write(YAML::Emitter& yamlOutput,
                                                   const rpwa::massDepFit::parameters& fitParameters,
                                                   const rpwa::massDepFit::parameters& fitParametersError,
                                                   const bool useBranchings,
                                                   const bool debug) const
{
	if(debug) {
		printDebug << "start writing 'tPrimeDependentBackground' for component '" << getName() << "'." << std::endl;
		print(printDebug);
	}

	yamlOutput << YAML::BeginMap;

	if(not component::write(yamlOutput, fitParameters, fitParametersError, useBranchings, debug)) {
		printErr << "error while writing 'component' part of 'tPrimeDependentBackground'." << std::endl;
		return false;
	}

	yamlOutput << YAML::Key << "mIsobar1";
	yamlOutput << YAML::Value << _m1;

	yamlOutput << YAML::Key << "mIsobar2";
	yamlOutput << YAML::Value << _m2;

	yamlOutput << YAML::Key << "relAngularMom";
	yamlOutput << YAML::Value <<_l;

	yamlOutput << YAML::Key << "exponent";
	yamlOutput << YAML::Value <<_exponent;

	yamlOutput << YAML::EndMap;

	if(debug) {
		printDebug << "finished writing 'tPrimeDependentBackground'." << std::endl;
	}

	return true;
}


std::complex<double>
rpwa::massDepFit::tPrimeDependentBackground::val(const rpwa::massDepFit::parameters& fitParameters,
                                                 rpwa::massDepFit::cache& cache,
                                                 const size_t idxBin,
                                                 const double m,
                                                 const size_t idxMass) const
{
	if (idxMass != std::numeric_limits<size_t>::max()) {
		const std::complex<double> component = cache.getComponent(getId(), idxBin, idxMass);
		if (component != 0.) {
			return component;
		}
	}

	// calculate breakup momentum
	if(m < _m1+_m2) {
		return std::pow(m - fitParameters.getParameter(getId(), 0), fitParameters.getParameter(getId(), 1));
	}
	const double q = rpwa::breakupMomentum(m, _m1, _m2);
	const double f2 = rpwa::barrierFactorSquared(2*_l, q);
	const double c = std::pow(q*f2 * _norm, _exponent);

	// get mean t' value for current bin
	const double tPrime = _tPrimeMeans[idxBin];
	const double tPrimePol = fitParameters.getParameter(getId(), 2) + fitParameters.getParameter(getId(), 3)*tPrime + fitParameters.getParameter(getId(), 4)*tPrime*tPrime;

	const double mPre = std::pow(m - fitParameters.getParameter(getId(), 0), fitParameters.getParameter(getId(), 1));

	const std::complex<double> component = mPre * exp(-tPrimePol*c);

	if (idxMass != std::numeric_limits<size_t>::max()) {
		cache.setComponent(getId(), idxBin, idxMass, component);
	}

	return component;
}


std::ostream&
rpwa::massDepFit::tPrimeDependentBackground::print(std::ostream& out) const
{
	out << "component " << getId() << " '" << getName() << "' (tPrimeDependentBackground):" << std::endl;

	out << "    mass threshold ";
	out << "start value: " << _parametersStart[0] << " +/- " << _parametersError[0] << " ";
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
		out << "start value: " << _parametersStart[i] << " +/- " << _parametersError[i] << " ";
		if(_parametersLimitedLower[i] && _parametersLimitedUpper[i]) {
			out << "limits: " << _parametersLimitLower[i] << "-" << _parametersLimitUpper[i] << " GeV/c^2";
		} else if(_parametersLimitedLower[i]) {
			out << "lower limit: " << _parametersLimitLower[i] << " GeV/c^2";
		} else if(_parametersLimitedUpper[i]) {
			out << "upper limit: " << _parametersLimitUpper[i] << " GeV/c^2";
		} else {
			out << "unlimited";
		}
		out << (_parametersFixed[i] ? " (FIXED)" : "") << std::endl;
	}

	out << "    mass of isobar 1: " << _m1 << " GeV/c^2, mass of isobar 2: " << _m2 << " GeV/c^2" << std::endl;
	out << "    relative orbital angular momentum between isobars: " << _l << " (in units of hbar)" << std::endl;
	out << "    exponent of break-up momentum times barrier-factor squared: " << _exponent << std::endl;

	out << "    for " << _tPrimeMeans.size() << " bin" << ((_tPrimeMeans.size()>1)?"s":"") << " with mean t' value" << ((_tPrimeMeans.size()>1)?"s":"") << ": " << _tPrimeMeans[0];
	for(size_t i=1; i<_tPrimeMeans.size(); ++i) {
		out << ", " << _tPrimeMeans[i];
	}
	out << std::endl;

	return component::print(out);
}


rpwa::massDepFit::exponentialBackgroundIntegral::exponentialBackgroundIntegral(const size_t id,
                                                                               const std::string& name)
	: component(id, name, "exponentialBackgroundIntegral", 1),
	  _interpolator(NULL)
{
	_parametersName[0] = "g";

	_parametersStep[0] = 1.0;
}


rpwa::massDepFit::exponentialBackgroundIntegral::~exponentialBackgroundIntegral()
{
	if (_interpolator != NULL) {
		delete _interpolator;
	}
}


bool
rpwa::massDepFit::exponentialBackgroundIntegral::init(const YAML::Node& configComponent,
                                                      rpwa::massDepFit::parameters& fitParameters,
                                                      rpwa::massDepFit::parameters& fitParametersError,
                                                      const size_t nrBins,
                                                      const std::vector<double>& massBinCenters,
                                                      const std::map<std::string, size_t>& waveIndices,
                                                      const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                                      const bool useBranchings,
                                                      const bool debug)
{
	if(debug) {
		printDebug << "start initializing 'exponentialBackgroundIntegral' for component '" << getName() << "'." << std::endl;
	}

	if(not component::init(configComponent, fitParameters, fitParametersError, nrBins, massBinCenters, waveIndices, phaseSpaceIntegrals, useBranchings, debug)) {
		printErr << "error while reading configuration of 'component' class." << std::endl;
		return false;
	}

	if(getNrChannels() != 1) {
		printErr << "component '" << getName() << "' of type 'exponentialBackgroundIntegral' needs to have exactly one decay channel." << std::endl;
		return false;
	}

	std::map<std::string, YamlCppUtils::Type> mandatoryArguments;
	boost::assign::insert(mandatoryArguments)
	                     ("integral", YamlCppUtils::TypeSequence)
	                     ("exponent", YamlCppUtils::TypeFloat);
	if(not checkIfAllVariablesAreThere(configComponent, mandatoryArguments)) {
		printErr << "not all required variables are defined for the component '" << getName() << "'." << std::endl;
		return false;
	}

	const YAML::Node& integrals = configComponent["integral"];

	const size_t nrValues = integrals.size();
	if(nrValues < 2) {
		printErr << "phase-space integral of component '" << getName() << "' has to contain at least two points." << std::endl;
		return false;
	}

	_masses.clear();
	_values.clear();

	for(size_t idx=0; idx<nrValues; ++idx) {
		const YAML::Node& integral = integrals[idx];
		if(not integral.IsSequence()) {
			printErr << "phase-space integrals of component '" << getName() << "' has to consist of arrays of two floating point numbers." << std::endl;
			return false;
		}
		if(integral.size() != 2) {
			printErr << "phase-space integrals of component '" << getName() << "' has to consist of arrays of two floating point numbers." << std::endl;
			return false;
		}
		if(not checkVariableType(integral[0], YamlCppUtils::TypeFloat)) {
			printErr << "phase-space integrals of component '" << getName() << "' has to consist of arrays of two floating point numbers." << std::endl;
			return false;
		}
		if(not checkVariableType(integral[1], YamlCppUtils::TypeFloat)) {
			printErr << "phase-space integrals of component '" << getName() << "' has to consist of arrays of two floating point numbers." << std::endl;
			return false;
		}

		const double mass = integral[0].as<double>();
		const double val = integral[1].as<double>();

		if (_masses.size() > 0 && _masses.back() > mass) {
			printErr << "masses of phase-space integrals of component '" << getName() << "' have to be strictly ordered." << std::endl;
			return false;
		}

		_masses.push_back(mass);
		_values.push_back(val);
	}

	assert(_interpolator == NULL);
	_interpolator = new ROOT::Math::Interpolator(_masses, _values, ROOT::Math::Interpolation::kLINEAR);

	_exponent = configComponent["exponent"].as<double>();

	_norm = 1. / (massBinCenters.back() * _interpolator->Eval(massBinCenters.back()));

	if(debug) {
		print(printDebug);
		printDebug << "finished initializing 'exponentialBackgroundIntegral'." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::exponentialBackgroundIntegral::write(YAML::Emitter& yamlOutput,
                                                       const rpwa::massDepFit::parameters& fitParameters,
                                                       const rpwa::massDepFit::parameters& fitParametersError,
                                                       const bool useBranchings,
                                                       const bool debug) const
{
	if(debug) {
		printDebug << "start writing 'exponentialBackgroundIntegral' for component '" << getName() << "'." << std::endl;
		print(printDebug);
	}

	yamlOutput << YAML::BeginMap;

	if(not component::write(yamlOutput, fitParameters, fitParametersError, useBranchings, debug)) {
		printErr << "error while writing 'component' part of 'exponentialBackgroundIntegral'." << std::endl;
		return false;
	}

	yamlOutput << YAML::Key << "integral";
	yamlOutput << YAML::Value;

	yamlOutput << YAML::Flow;
	yamlOutput << YAML::BeginSeq;
	for (size_t idx=0; idx<_masses.size(); ++idx) {
		yamlOutput << YAML::BeginSeq;
		yamlOutput << _masses[idx];
		yamlOutput << _values[idx];
		yamlOutput << YAML::EndSeq;
	}
	yamlOutput << YAML::EndSeq;
	yamlOutput << YAML::Block;

	yamlOutput << YAML::Key << "exponent";
	yamlOutput << YAML::Value <<_exponent;

	yamlOutput << YAML::EndMap;

	if(debug) {
		printDebug << "finished writing 'exponentialBackgroundIntegral'." << std::endl;
	}

	return true;
}


std::complex<double>
rpwa::massDepFit::exponentialBackgroundIntegral::val(const rpwa::massDepFit::parameters& fitParameters,
                                                     rpwa::massDepFit::cache& cache,
                                                     const size_t idxBin,
                                                     const double m,
                                                     const size_t idxMass) const
{
	if (idxMass != std::numeric_limits<size_t>::max()) {
		const std::complex<double> component = cache.getComponent(getId(), idxBin, idxMass);
		if (component != 0.) {
			return component;
		}
	}

	const double ps = _interpolator->Eval(m);
	const double c = std::pow(m * ps * _norm, _exponent);

	const std::complex<double> component = exp(-fitParameters.getParameter(getId(), 0)*c);

	if (idxMass != std::numeric_limits<size_t>::max()) {
		cache.setComponent(getId(), std::numeric_limits<size_t>::max(), idxMass, component);
	}

	return component;
}


std::ostream&
rpwa::massDepFit::exponentialBackgroundIntegral::print(std::ostream& out) const
{
	out << "component " << getId() << " '" << getName() << "' (exponentialBackgroundIntegral):" << std::endl;

	out << "    width ";
	out << "start value: " << _parametersStart[0] << " +/- " << _parametersError[0] << " ";
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

	out << "    exponent of phase-space integral: " << _exponent << std::endl;

	return component::print(out);
}


rpwa::massDepFit::tPrimeDependentBackgroundIntegral::tPrimeDependentBackgroundIntegral(const size_t id,
                                                                                       const std::string& name)
	: component(id, name, "tPrimeDependentBackgroundIntegral", 5),
	  _interpolator(NULL)
{
	_parametersName[0] = "m0";
	_parametersName[1] = "c0";
	_parametersName[2] = "c1";
	_parametersName[3] = "c2";
	_parametersName[4] = "c3";

	_parametersStep[0] = 0.001;
	_parametersStep[1] = 0.001;
	_parametersStep[2] = 1.0;
	_parametersStep[3] = 1.0;
	_parametersStep[4] = 1.0;
}


rpwa::massDepFit::tPrimeDependentBackgroundIntegral::~tPrimeDependentBackgroundIntegral()
{
	if (_interpolator != NULL) {
		delete _interpolator;
	}
}


bool
rpwa::massDepFit::tPrimeDependentBackgroundIntegral::setTPrimeMeans(const std::vector<double> tPrimeMeans)
{
	_tPrimeMeans = tPrimeMeans;

	return true;
}


bool
rpwa::massDepFit::tPrimeDependentBackgroundIntegral::init(const YAML::Node& configComponent,
                                                          rpwa::massDepFit::parameters& fitParameters,
                                                          rpwa::massDepFit::parameters& fitParametersError,
                                                          const size_t nrBins,
                                                          const std::vector<double>& massBinCenters,
                                                          const std::map<std::string, size_t>& waveIndices,
                                                          const boost::multi_array<double, 3>& phaseSpaceIntegrals,
                                                          const bool useBranchings,
                                                          const bool debug)
{
	if(debug) {
		printDebug << "start initializing 'tPrimeDependentBackgroundIntegral' for component '" << getName() << "'." << std::endl;
	}

	if(not component::init(configComponent, fitParameters, fitParametersError, nrBins, massBinCenters, waveIndices, phaseSpaceIntegrals, useBranchings, debug)) {
		printErr << "error while reading configuration of 'component' class." << std::endl;
		return false;
	}

	if(getNrChannels() != 1) {
		printErr << "component '" << getName() << "' of type 'tPrimeDependentBackgroundIntegral' needs to have exactly one decay channel." << std::endl;
		return false;
	}

	std::map<std::string, YamlCppUtils::Type> mandatoryArguments;
	boost::assign::insert(mandatoryArguments)
	                     ("integral", YamlCppUtils::TypeSequence)
	                     ("exponent", YamlCppUtils::TypeFloat);
	if(not checkIfAllVariablesAreThere(configComponent, mandatoryArguments)) {
		printErr << "not all required variables are defined for the component '" << getName() << "'." << std::endl;
		return false;
	}

	const YAML::Node& integrals = configComponent["integral"];

	const size_t nrValues = integrals.size();
	if(nrValues < 2) {
		printErr << "phase-space integral of component '" << getName() << "' has to contain at least two points." << std::endl;
		return false;
	}

	_masses.clear();
	_values.clear();

	for(size_t idx=0; idx<nrValues; ++idx) {
		const YAML::Node& integral = integrals[idx];
		if(not integral.IsSequence()) {
			printErr << "phase-space integrals of component '" << getName() << "' has to consist of arrays of two floating point numbers." << std::endl;
			return false;
		}
		if(integral.size() != 2) {
			printErr << "phase-space integrals of component '" << getName() << "' has to consist of arrays of two floating point numbers." << std::endl;
			return false;
		}
		if(not checkVariableType(integral[0], YamlCppUtils::TypeFloat)) {
			printErr << "phase-space integrals of component '" << getName() << "' has to consist of arrays of two floating point numbers." << std::endl;
			return false;
		}
		if(not checkVariableType(integral[1], YamlCppUtils::TypeFloat)) {
			printErr << "phase-space integrals of component '" << getName() << "' has to consist of arrays of two floating point numbers." << std::endl;
			return false;
		}

		const double mass = integral[0].as<double>();
		const double val = integral[1].as<double>();

		if (_masses.size() > 0 && _masses.back() > mass) {
			printErr << "masses of phase-space integrals of component '" << getName() << "' have to be strictly ordered." << std::endl;
			return false;
		}

		_masses.push_back(mass);
		_values.push_back(val);
	}

	assert(_interpolator == NULL);
	_interpolator = new ROOT::Math::Interpolator(_masses, _values, ROOT::Math::Interpolation::kLINEAR);

	_exponent = configComponent["exponent"].as<double>();

	_norm = 1. / (massBinCenters.back() * _interpolator->Eval(massBinCenters.back()));

	if(_tPrimeMeans.size() != nrBins) {
		printErr << "array of mean t' value in each bin does not contain the correct number of entries (is: " << _tPrimeMeans.size() << ", expected: " << nrBins << ")." << std::endl;
		return false;
	}

	if(debug) {
		print(printDebug);
		printDebug << "finished initializing 'tPrimeDependentBackgroundIntegral'." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::tPrimeDependentBackgroundIntegral::write(YAML::Emitter& yamlOutput,
                                                           const rpwa::massDepFit::parameters& fitParameters,
                                                           const rpwa::massDepFit::parameters& fitParametersError,
                                                           const bool useBranchings,
                                                           const bool debug) const
{
	if(debug) {
		printDebug << "start writing 'tPrimeDependentBackgroundIntegral' for component '" << getName() << "'." << std::endl;
		print(printDebug);
	}

	yamlOutput << YAML::BeginMap;

	if(not component::write(yamlOutput, fitParameters, fitParametersError, useBranchings, debug)) {
		printErr << "error while writing 'component' part of 'tPrimeDependentBackgroundIntegral'." << std::endl;
		return false;
	}

	yamlOutput << YAML::Key << "integral";
	yamlOutput << YAML::Value;

	yamlOutput << YAML::Flow;
	yamlOutput << YAML::BeginSeq;
	for (size_t idx=0; idx<_masses.size(); ++idx) {
		yamlOutput << YAML::BeginSeq;
		yamlOutput << _masses[idx];
		yamlOutput << _values[idx];
		yamlOutput << YAML::EndSeq;
	}
	yamlOutput << YAML::EndSeq;
	yamlOutput << YAML::Block;

	yamlOutput << YAML::Key << "exponent";
	yamlOutput << YAML::Value <<_exponent;

	yamlOutput << YAML::EndMap;

	if(debug) {
		printDebug << "finished writing 'tPrimeDependentBackgroundIntegral'." << std::endl;
	}

	return true;
}


std::complex<double>
rpwa::massDepFit::tPrimeDependentBackgroundIntegral::val(const rpwa::massDepFit::parameters& fitParameters,
                                                         rpwa::massDepFit::cache& cache,
                                                         const size_t idxBin,
                                                         const double m,
                                                         const size_t idxMass) const
{
	if (idxMass != std::numeric_limits<size_t>::max()) {
		const std::complex<double> component = cache.getComponent(getId(), idxBin, idxMass);
		if (component != 0.) {
			return component;
		}
	}

	const double ps = _interpolator->Eval(m);
	const double c = std::pow(m * ps * _norm, _exponent);

	// get mean t' value for current bin
	const double tPrime = _tPrimeMeans[idxBin];
	const double tPrimePol = fitParameters.getParameter(getId(), 2) + fitParameters.getParameter(getId(), 3)*tPrime + fitParameters.getParameter(getId(), 4)*tPrime*tPrime;

	const double mPre = std::pow(m - fitParameters.getParameter(getId(), 0), fitParameters.getParameter(getId(), 1));

	const std::complex<double> component = mPre * exp(-tPrimePol*c);

	if (idxMass != std::numeric_limits<size_t>::max()) {
		cache.setComponent(getId(), idxBin, idxMass, component);
	}

	return component;
}


std::ostream&
rpwa::massDepFit::tPrimeDependentBackgroundIntegral::print(std::ostream& out) const
{
	out << "component " << getId() << " '" << getName() << "' (tPrimeDependentBackgroundIntegral):" << std::endl;

	out << "    mass threshold ";
	out << "start value: " << _parametersStart[0] << " +/- " << _parametersError[0] << " ";
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
		out << "start value: " << _parametersStart[i] << " +/- " << _parametersError[i] << " ";
		if(_parametersLimitedLower[i] && _parametersLimitedUpper[i]) {
			out << "limits: " << _parametersLimitLower[i] << "-" << _parametersLimitUpper[i] << " GeV/c^2";
		} else if(_parametersLimitedLower[i]) {
			out << "lower limit: " << _parametersLimitLower[i] << " GeV/c^2";
		} else if(_parametersLimitedUpper[i]) {
			out << "upper limit: " << _parametersLimitUpper[i] << " GeV/c^2";
		} else {
			out << "unlimited";
		}
		out << (_parametersFixed[i] ? " (FIXED)" : "") << std::endl;
	}

	out << "    exponent of phase-space integral: " << _exponent << std::endl;

	out << "    for " << _tPrimeMeans.size() << " bin" << ((_tPrimeMeans.size()>1)?"s":"") << " with mean t' value" << ((_tPrimeMeans.size()>1)?"s":"") << ": " << _tPrimeMeans[0];
	for(size_t i=1; i<_tPrimeMeans.size(); ++i) {
		out << ", " << _tPrimeMeans[i];
	}
	out << std::endl;

	return component::print(out);
}
