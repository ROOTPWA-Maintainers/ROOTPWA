///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010-2012 Sebastian Neubert (TUM)
//    Copyright 2014-2016 Sebastian Uhl (TUM)
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


#include "components.h"

#include <set>

#include <physUtils.hpp>
#include <reportingUtils.hpp>

#include "resonanceFitHelper.h"


namespace {


	// initialize the bins that can be used for the caching, each value is
	// only cached for the given bin in val()
	std::vector<std::vector<size_t> >
	initBinsEqualValuesForSingleBins(const size_t nrBins)
	{
		std::vector<std::vector<size_t> > binsEqualValues(nrBins, std::vector<size_t>(1, 0));
		for(size_t idxBin = 0; idxBin < nrBins; ++idxBin) {
			binsEqualValues[idxBin][0] = idxBin;
		}

		return binsEqualValues;
	}


	// initialize the bins that can be used for the caching, each value is
	// used in all bins that have the same mass binning
	std::vector<std::vector<size_t> >
	initBinsEqualValuesFromMassBinning(const std::vector<size_t>& nrMassBins,
	                                   const boost::multi_array<double, 2>& massBinCenters)
	{
		const size_t nrBins = nrMassBins.size();

		std::vector<std::vector<size_t> > binsEqualValues(nrBins, std::vector<size_t>());
		for(size_t idxBin = 0; idxBin < nrBins; ++idxBin) {
			if(binsEqualValues[idxBin].size() != 0) {
				continue;
			}

			std::vector<size_t> bins(1, idxBin);
			for(size_t jdxBin = idxBin+1; jdxBin < nrBins; ++jdxBin) {
				if(nrMassBins[idxBin] != nrMassBins[jdxBin]) {
					continue;
				}

				bool sameBinning = true;
				for(size_t idxMass = 0; idxMass < nrMassBins[idxBin]; ++idxMass) {
					if(massBinCenters[idxBin][idxMass] != massBinCenters[jdxBin][idxMass]) {
						sameBinning = false;
						break;
					}
				}
				if(sameBinning) {
					bins.push_back(jdxBin);
				}
			}

			for(size_t jdxBin = 0; jdxBin < bins.size(); ++jdxBin) {
				assert(binsEqualValues[bins[jdxBin]].size() == 0);

				binsEqualValues[bins[jdxBin]] = bins;
			}
		}

		return binsEqualValues;
	}


	template<typename T>
	void
	checkParameters(const T& component)
	{
		const std::vector<rpwa::resonanceFit::parameter> defaultParameters = T::getDefaultParameters();

		// check that all parameters are there
		if(component.getNrParameters() != defaultParameters.size()) {
			printErr << "inconsistent number of parameters for component '" << component.getName() << "', "
			         << "expected " << defaultParameters.size() << ", found " << component.getNrParameters() << ". "
			         << "Aborting..." << std::endl;
			throw;
		}

		// check that the correct parameters are at the correct position
		for(size_t idxParameter = 0; idxParameter < component.getNrParameters(); ++idxParameter) {
			if(component.getParameter(idxParameter).name() != defaultParameters[idxParameter].name()) {
				printErr << "inconsistent naming of parameters for component '" << component.getName() << "', "
				         << "expected '" << defaultParameters[idxParameter].name() << "', "
				         << "found '" << component.getParameter(idxParameter).name() << "' for parameter at index " << idxParameter << ". "
				         << "Aborting..." << std::endl;
				throw;
			}
		}
	}


}


rpwa::resonanceFit::component::channel::channel(const std::string& waveName,
                                                const std::vector<size_t>& waveIndices,
                                                const std::vector<size_t>& nrMassBins,
                                                const boost::multi_array<double, 2>& massBinCenters,
                                                const boost::multi_array<double, 3>& phaseSpaceIntegrals)
	: _waveName(waveName),
	  _waveIndices(waveIndices)
{
	// get dimensions from one array and make sure that all other arrays
	// have the same dimensions
	const size_t nrBins = nrMassBins.size();
	if(nrBins == 0) {
		printErr << "number of bins is zero, cannot perform the fit. Aborting..." << std::endl;
		throw;
	}

	const size_t maxNrMassBins = *std::max_element(nrMassBins.begin(), nrMassBins.end());
	if(maxNrMassBins == 0) {
		printErr << "maximal number of mass bins is zero, cannot perform the fit. Aborting..." << std::endl;
		throw;
	}

	const size_t maxNrWaves = *(phaseSpaceIntegrals.shape()+2);
	if(maxNrWaves == 0) {
		printErr << "maximal number of waves is zero, cannot perform the fit. Aborting..." << std::endl;
		throw;
	}

	checkSize(waveIndices,
	          nrBins, "number of bins is not correct for wave indices.");
	checkSize(nrMassBins,
	          nrBins, "number of bins is not correct for number of mass bins.");
	checkSize(massBinCenters,
	          nrBins, "number of bins is not correct for centers of mass bins.",
	          maxNrMassBins, "maximal number of mass bins is not correct for centers of mass bins.");
	checkSize(phaseSpaceIntegrals,
	          nrBins, "number of bins is not correct for phase-space integrals.",
	          maxNrMassBins, "maximal number of mass bins is not correct for phase-space integrals.",
	          maxNrWaves, "maximal number of waves is not correct for phase-space integrals.");

	_anchors.resize(nrBins, false),

	_interpolators.resize(nrBins);
	_phaseSpaceIntegrals.resize(boost::extents[nrBins][maxNrMassBins]);
	for(size_t idxBin = 0; idxBin < nrBins; ++idxBin) {
		const size_t idxWave = waveIndices[idxBin];

		// skip bins that do not contain the current wave
		if(idxWave == std::numeric_limits<size_t>::max()) {
			continue;
		}

		_bins.push_back(idxBin);

		const boost::multi_array<double, 1>& viewM = massBinCenters[boost::indices[idxBin][boost::multi_array<double, 2>::index_range(0, nrMassBins[idxBin])]];
		const boost::multi_array<double, 1>& viewInt = phaseSpaceIntegrals[boost::indices[idxBin][boost::multi_array<double, 3>::index_range(0, nrMassBins[idxBin])][idxWave]];

		// copy values at mass bin centers
		adjustSizeAndSet(_phaseSpaceIntegrals, idxBin, viewInt);

		// set up interpolator for when the requested mass is not a
		// mass bin center
		_interpolators[idxBin] = std::make_shared<ROOT::Math::Interpolator>(std::vector<double>(viewM.begin(), viewM.end()), std::vector<double>(viewInt.begin(), viewInt.end()), ROOT::Math::Interpolation::kLINEAR);
	}
}


rpwa::resonanceFit::component::component(const size_t id,
                                         const std::string& name,
                                         const std::vector<rpwa::resonanceFit::parameter>& parameters,
                                         const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
                                         const std::vector<size_t>& nrMassBins,
                                         const boost::multi_array<double, 2>& massBinCenters,
                                         const bool useBranchings)
	: _id(id),
	  _name(name),
	  _channels(decayChannels),
	  _nrCouplings(0),
	  _nrBranchings(0),
	  _parameters(parameters)
{
	// get dimensions from one array and make sure that all other arrays
	// have the same dimensions
	const size_t nrBins = nrMassBins.size();
	if(nrBins == 0) {
		printErr << "number of bins is zero, cannot perform the fit. Aborting..." << std::endl;
		throw;
	}

	const size_t maxNrMassBins = *std::max_element(nrMassBins.begin(), nrMassBins.end());
	if(maxNrMassBins == 0) {
		printErr << "maximal number of mass bins is zero, cannot perform the fit. Aborting..." << std::endl;
		throw;
	}

	checkSize(nrMassBins,
	          nrBins, "number of bins is not correct for number of mass bins.");
	checkSize(massBinCenters,
	          nrBins, "number of bins is not correct for centers of mass bins.",
	          maxNrMassBins, "maximal number of mass bins is not correct for centers of mass bins.");

	// initialize the bins that can be used for the caching, the function
	// value might be the same in several bins, but by default each value
	// is only cached for the given bin in val()
	_binsEqualValues = initBinsEqualValuesForSingleBins(nrBins);

	std::vector<std::set<std::string>> waveNames(nrBins, std::set<std::string>());
	std::map<std::pair<std::string, size_t>, size_t> couplingsQN;
	for(size_t idxDecayChannel = 0; idxDecayChannel < _channels.size(); ++idxDecayChannel) {
		const std::string waveName = _channels[idxDecayChannel].getWaveName();

		// get list of bins this wave is defined in
		const std::vector<size_t>& binsForWave = _channels[idxDecayChannel].getBins();

		// check that each wave only appears once in each bin
		for(size_t i = 0; i < binsForWave.size(); ++i) {
			const size_t idxBin = binsForWave[i];

			if(waveNames[idxBin].count(waveName) != 0) {
				printErr << "wave '" << waveName << "' used for more than one decay channel in component '" << name << "' (bin " << idxBin << "). Aborting..." << std::endl;
				throw;
			}

			waveNames[idxBin].insert(waveName);
		}

		bool readCouplings = true;
		bool readBranching = false;
		bool fixedBranching = true;
		size_t couplingIndex = _nrCouplings;
		size_t branchingIndex = _nrBranchings;
		if(useBranchings and _channels.size() > 1) {
			const std::string waveQN = waveName.substr(0, waveName.find("="));
			const std::string waveDecay = waveName.substr(waveName.find("=")+1);

			for(size_t i = 0; i < binsForWave.size(); ++i) {
				const size_t idxBin = binsForWave[i];

				std::map<std::pair<std::string, size_t>, size_t>::const_iterator mappingQN = couplingsQN.find(std::make_pair(waveQN, idxBin));
				if(mappingQN == couplingsQN.end()) {
					couplingsQN[std::make_pair(waveQN, idxBin)] = couplingIndex;

					if(i > 1 and (not readCouplings or not fixedBranching)) {
						printErr << "inconsistent status of couplings and branchings across bins." << std::endl;
						throw;
					}

					readCouplings = true;
					fixedBranching = true;
				} else {
					if(i > 1 and (readCouplings or fixedBranching or couplingIndex != mappingQN->second)) {
						printErr << "inconsistent status of couplings and branchings across bins." << std::endl;
						throw;
					}

					readCouplings = false;
					fixedBranching = false;
					couplingIndex = mappingQN->second;
				}

				if(i > 1 and (not readBranching)) {
					printErr << "inconsistent status of couplings and branchings across bins." << std::endl;
					throw;
				}

				readBranching = true;
			}
		}

		if(readCouplings) {
			// new coupling, so this should be the last coupling up to now
			assert(couplingIndex == _nrCouplings);

			_channelsFromCoupling.push_back(idxDecayChannel);
			++_nrCouplings;
		}

		if(readBranching) {
			// new branching, so this should be the last branching up to now
			assert(branchingIndex == _nrBranchings);

			_channelsFromBranching.push_back(idxDecayChannel);
			_branchingsFixed.push_back(fixedBranching);
			++_nrBranchings;
		}

		_channelsCoupling.push_back(couplingIndex);
		_channelsBranching.push_back(branchingIndex);
	}

	// if no branchings are read at all make sure that there is one equal to 1.
	if(_nrBranchings == 0) {
		_channelsFromBranching.push_back(0);
		_branchingsFixed.push_back(true);
		_nrBranchings++;
	}
}


void
rpwa::resonanceFit::component::setChannelAnchor(const size_t idxBin,
                                                const size_t idxDecayChannel)
{
	_channels[idxDecayChannel].setAnchor(idxBin);
}


size_t
rpwa::resonanceFit::component::importCouplings(const double* par,
                                               rpwa::resonanceFit::parameters& fitParameters,
                                               rpwa::resonanceFit::cache& cache) const
{
	size_t counter = 0;
	for(size_t idxCoupling = 0; idxCoupling < _nrCouplings; ++idxCoupling) {
		const rpwa::resonanceFit::component::channel& channel = _channels[_channelsFromCoupling[idxCoupling]];
		const std::vector<size_t>& bins = channel.getBins();
		for(size_t i = 0; i < bins.size(); ++i) {
			const size_t idxBin = bins[i];
			std::complex<double> coupling;
			if(channel.isAnchor(idxBin)) {
				coupling = std::complex<double>(par[counter], 0.);
				counter += 1;
			} else {
				coupling = std::complex<double>(par[counter], par[counter+1]);
				counter += 2;
			}

			if(fitParameters.getCoupling(getId(), idxCoupling, idxBin) != coupling) {
				fitParameters.setCoupling(getId(), idxCoupling, idxBin, coupling);

				// invalidate the cache
				for(size_t idxChannel = 0; idxChannel < _channels.size(); ++idxChannel) {
					if(_channelsCoupling[idxChannel] == idxCoupling) {
						cache.setCoupling(getId(), idxChannel, idxBin, std::numeric_limits<size_t>::max(), 0.);
						cache.setProdAmp(_channels[idxChannel].getWaveIndices()[idxBin], idxBin, std::numeric_limits<size_t>::max(), 0.);
					}
				}
			}
		}
	}

	return counter;
}


size_t
rpwa::resonanceFit::component::importBranchings(const double* par,
                                                rpwa::resonanceFit::parameters& fitParameters,
                                                rpwa::resonanceFit::cache& cache) const
{
	size_t counter = 0;
	for(size_t idxBranching = 0; idxBranching < _nrBranchings; ++idxBranching) {
		std::complex<double> branching(1., 0.);
		if(not _branchingsFixed[idxBranching]) {
			branching = std::complex<double>(par[counter], par[counter+1]);
			counter += 2;
		}

		if(fitParameters.getBranching(getId(), idxBranching) != branching) {
			fitParameters.setBranching(getId(), idxBranching, branching);

			// invalidate the cache
			for(size_t idxChannel = 0; idxChannel < _channels.size(); ++idxChannel) {
				if(_channelsBranching[idxChannel] == idxBranching) {
					cache.setCoupling(getId(), idxChannel, std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), 0.);
					for(std::vector<size_t>::const_iterator idxBin = _channels[idxChannel].getBins().begin(); idxBin != _channels[idxChannel].getBins().end(); ++idxBin) {
						cache.setProdAmp(_channels[idxChannel].getWaveIndices()[*idxBin], *idxBin, std::numeric_limits<size_t>::max(), 0.);
					}
				}
			}
		}
	}

	return counter;
}


size_t
rpwa::resonanceFit::component::importParameters(const double* par,
                                                rpwa::resonanceFit::parameters& fitParameters,
                                                rpwa::resonanceFit::cache& cache) const
{
	bool invalidateCache = false;
	for(size_t idxParameter = 0; idxParameter < _parameters.size(); ++idxParameter) {
		if(fitParameters.getParameter(getId(), idxParameter) != par[idxParameter]) {
			fitParameters.setParameter(getId(), idxParameter, par[idxParameter]);
			invalidateCache = true;
		}
	}

	if(invalidateCache) {
		cache.setComponent(getId(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), 0.);
		for(size_t idxChannel = 0; idxChannel < _channels.size(); ++idxChannel) {
			for(std::vector<size_t>::const_iterator idxBin = _channels[idxChannel].getBins().begin(); idxBin != _channels[idxChannel].getBins().end(); ++idxBin) {
				cache.setProdAmp(_channels[idxChannel].getWaveIndices()[*idxBin], *idxBin, std::numeric_limits<size_t>::max(), 0.);
			}
		}
	}

	return _parameters.size();
}


std::complex<double>
rpwa::resonanceFit::component::val(const rpwa::resonanceFit::parameters& fitParameters,
                                   rpwa::resonanceFit::cache& cache,
                                   const size_t idxBin,
                                   const double mass,
                                   const size_t idxMass) const
{
	if(idxMass != std::numeric_limits<size_t>::max()) {
		const std::complex<double> component = cache.getComponent(getId(), idxBin, idxMass);
		if(component != 0.) {
			return component;
		}
	}

	const std::complex<double> component = val(fitParameters, idxBin, mass);

	if(idxMass != std::numeric_limits<size_t>::max()) {
		if(_binsEqualValues[idxBin].size() == _binsEqualValues.size()) {
			// the value is the same for all bins
			cache.setComponent(getId(), std::numeric_limits<size_t>::max(), idxMass, component);
		} else {
			for(std::vector<size_t>::const_iterator it = _binsEqualValues[idxBin].begin(); it != _binsEqualValues[idxBin].end(); ++it) {
				cache.setComponent(getId(), *it, idxMass, component);
			}
		}
	}

	return component;
}


std::ostream&
rpwa::resonanceFit::component::print(std::ostream& out, const bool newLine) const
{
	out << "component " << getId() << " '" << getName() << "' (" << getType() << "):" << std::endl;

	for(size_t idxParameter = 0; idxParameter < getNrParameters(); ++idxParameter) {
		const rpwa::resonanceFit::parameter& parameter = getParameter(idxParameter);

		out << "    * '" << parameter.name() << "' ";
		out << "start value: " << parameter.startValue() << " +/- " << parameter.startError() << " ";
		if(parameter.limitedLower() and parameter.limitedUpper()) {
			out << "limits: " << parameter.limitLower() << "-" << parameter.limitUpper() << " GeV/c^2";
		} else if(parameter.limitedLower()) {
			out << "lower limit: " << parameter.limitLower() << " GeV/c^2";
		} else if(parameter.limitedUpper()) {
			out << "upper limit: " << parameter.limitUpper() << " GeV/c^2";
		} else {
			out << "unlimited";
		}
		out << (parameter.fixed() ? " (FIXED)" : "") << std::endl;
	}

	out << "    use equal values for each mass bin in group of bins: ";
	std::set<size_t> printed;
	for(size_t idxBin = 0; idxBin < _binsEqualValues.size(); ++idxBin) {
		if(printed.count(idxBin) > 0) {
			continue;
		}
		out << ((printed.size() > 0) ? ", " : "") << "{";
		for(size_t idx = 0; idx < _binsEqualValues[idxBin].size(); ++idx) {
			out << ((idx > 0) ? ", " : "") << _binsEqualValues[idxBin][idx];
			printed.insert(_binsEqualValues[idxBin][idx]);
		}
		out << "}";
	}
	out << std::endl;

	out << "Decay modes:";
	for(size_t idxChannel = 0; idxChannel < _channels.size(); ++idxChannel) {
		const rpwa::resonanceFit::component::channel& channel = _channels[idxChannel];
		out << std::endl << "    " << idxChannel << ": '" << channel.getWaveName() << "'";
	}

	if(newLine) {
		out << std::endl;
	}

	return out;
}


std::vector<rpwa::resonanceFit::parameter>
rpwa::resonanceFit::fixedWidthBreitWigner::getDefaultParameters()
{
	std::vector<rpwa::resonanceFit::parameter> defaultParameters(2);
	defaultParameters[0].setName("mass");
	defaultParameters[0].setStep(0.001);
	defaultParameters[1].setName("width");
	defaultParameters[1].setStep(0.001);

	return defaultParameters;
}


rpwa::resonanceFit::fixedWidthBreitWigner::fixedWidthBreitWigner(const size_t id,
                                                                 const std::string& name,
                                                                 const std::vector<rpwa::resonanceFit::parameter>& parameters,
                                                                 const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
                                                                 const std::vector<size_t>& nrMassBins,
                                                                 const boost::multi_array<double, 2>& massBinCenters,
                                                                 const bool useBranchings)
	: component(id,
	            name,
	            parameters,
	            decayChannels,
	            nrMassBins,
	            massBinCenters,
	            useBranchings)
{
	// initialize the bins that can be used for the caching, the function
	// value for this component type is the same in bins with an equal
	// mass binning
	_binsEqualValues = initBinsEqualValuesFromMassBinning(nrMassBins, massBinCenters);

	// check presence of all parameters
	checkParameters(*this);
}


std::complex<double>
rpwa::resonanceFit::fixedWidthBreitWigner::val(const rpwa::resonanceFit::parameters& fitParameters,
                                               const size_t /*idxBin*/,
                                               const double mass) const
{
	const double& m0 = fitParameters.getParameter(getId(), 0);
	const double& gamma0 = fitParameters.getParameter(getId(), 1);

	const std::complex<double> component = gamma0*m0 / std::complex<double>(m0*m0-mass*mass, -gamma0*m0);

	return component;
}


std::vector<rpwa::resonanceFit::parameter>
rpwa::resonanceFit::dynamicWidthBreitWigner::getDefaultParameters()
{
	std::vector<rpwa::resonanceFit::parameter> defaultParameters(2);
	defaultParameters[0].setName("mass");
	defaultParameters[0].setStep(0.001);
	defaultParameters[1].setName("width");
	defaultParameters[1].setStep(0.001);

	return defaultParameters;
}


rpwa::resonanceFit::dynamicWidthBreitWigner::dynamicWidthBreitWigner(const size_t id,
                                                                     const std::string& name,
                                                                     const std::vector<rpwa::resonanceFit::parameter>& parameters,
                                                                     const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
                                                                     const std::vector<size_t>& nrMassBins,
                                                                     const boost::multi_array<double, 2>& massBinCenters,
                                                                     const bool useBranchings,
                                                                     const std::vector<double>& branchingRatio,
                                                                     const std::vector<int>& relAngularMom,
                                                                     const std::vector<double>& mIsobar1,
                                                                     const std::vector<double>& mIsobar2)
	: component(id,
	            name,
	            parameters,
	            decayChannels,
	            nrMassBins,
	            massBinCenters,
	            useBranchings),
	  _ratio(branchingRatio),
	  _l(relAngularMom),
	  _m1(mIsobar1),
	  _m2(mIsobar2)
{
	// initialize the bins that can be used for the caching, the function
	// value for this component type is the same in bins with an equal
	// mass binning
	_binsEqualValues = initBinsEqualValuesFromMassBinning(nrMassBins, massBinCenters);

	// check presence of all parameters
	checkParameters(*this);

	const size_t totalDecayChannels = _ratio.size();
	if(totalDecayChannels < getNrChannels()) {
		printErr << "total number of decay channels is smaller than number of proper decay channels." << std::endl;
		throw;
	}
	checkSize(_ratio,
	          totalDecayChannels, "total number of decay channels is not correct for branching ratios.");
	checkSize(_l,
	          totalDecayChannels, "total number of decay channels is not correct for relative orbital angular momenta.");
	checkSize(_m1,
	          totalDecayChannels, "total number of decay channels is not correct for mass of isobar 1.");
	checkSize(_m2,
	          totalDecayChannels, "total number of decay channels is not correct for mass of isobar 2.");

	const double sum = std::accumulate(_ratio.begin(), _ratio.end(), 0.0);
	std::transform(_ratio.begin(), _ratio.end(), _ratio.begin(), [&sum](const double v){ return v/sum; });
}


std::complex<double>
rpwa::resonanceFit::dynamicWidthBreitWigner::val(const rpwa::resonanceFit::parameters& fitParameters,
                                                 const size_t /*idxBin*/,
                                                 const double mass) const
{
	const double& m0 = fitParameters.getParameter(getId(), 0);
	const double& gamma0 = fitParameters.getParameter(getId(), 1);

	double gamma = 0.;
	for(size_t i = 0; i < _ratio.size(); ++i) {
		// calculation of the break-up momentum will fail if the
		// sum of the two isobar masses is heavier than the mother
		// mass, which is probably a good thing as some more advanced
		// stuff should be used for at- or sub-threshold decays.
		// but if the branching ratio for this channel is 0, then
		// simply ignore this channel
		if(_ratio[i] == 0.) {
			continue;
		}

		if(mass >= _m1[i] + _m2[i]) {
			// calculate breakup momenta
			const double q = rpwa::breakupMomentum(mass, _m1[i], _m2[i]);
			const double q0 = rpwa::breakupMomentum(m0, _m1[i], _m2[i]);

			// calculate barrier factors
			const double f2 = rpwa::barrierFactorSquared(2*_l[i], q);
			const double f20 = rpwa::barrierFactorSquared(2*_l[i], q0);

			gamma += _ratio[i] * q/q0 * f2/f20;
		}
	}
	gamma *= gamma0 * m0/mass;

	const std::complex<double> component = gamma0*m0 / std::complex<double>(m0*m0-mass*mass, -gamma*m0);

	return component;
}


std::ostream&
rpwa::resonanceFit::dynamicWidthBreitWigner::print(std::ostream& out, const bool newLine) const
{
	component::print(out, true);

	out << "Decay modes considered for dynamic width (" << _ratio.size() << "):";
	for(size_t i = 0; i < _ratio.size(); ++i) {
		out << std::endl
		    << "    * decay mode " << i << ", branching ratio: " << _ratio[i] << std::endl
		    << "      mass of isobar 1: " << _m1[i] << " GeV/c^2, mass of isobar 2: " << _m2[i] << " GeV/c^2" << std::endl
		    << "      relative orbital angular momentum between isobars: " << _l[i] << " (in units of hbar)";
	}

	if(newLine) {
		out << std::endl;
	}

	return out;
}


std::vector<rpwa::resonanceFit::parameter>
rpwa::resonanceFit::integralWidthBreitWigner::getDefaultParameters()
{
	std::vector<rpwa::resonanceFit::parameter> defaultParameters(2);
	defaultParameters[0].setName("mass");
	defaultParameters[0].setStep(0.001);
	defaultParameters[1].setName("width");
	defaultParameters[1].setStep(0.001);

	return defaultParameters;
}


rpwa::resonanceFit::integralWidthBreitWigner::integralWidthBreitWigner(const size_t id,
                                                                       const std::string& name,
                                                                       const std::vector<rpwa::resonanceFit::parameter>& parameters,
                                                                       const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
                                                                       const std::vector<size_t>& nrMassBins,
                                                                       const boost::multi_array<double, 2>& massBinCenters,
                                                                       const bool useBranchings,
                                                                       const std::vector<double>& branchingRatio,
                                                                       const std::vector<std::vector<double> >& masses,
                                                                       const std::vector<std::vector<double> >& values)
	: component(id,
	            name,
	            parameters,
	            decayChannels,
	            nrMassBins,
	            massBinCenters,
	            useBranchings),
	  _ratio(branchingRatio),
	  _masses(masses),
	  _values(values)
{
	// initialize the bins that can be used for the caching, the function
	// value for this component type is the same in bins with an equal
	// mass binning
	_binsEqualValues = initBinsEqualValuesFromMassBinning(nrMassBins, massBinCenters);

	// check presence of all parameters
	checkParameters(*this);

	const size_t totalDecayChannels = _ratio.size();
	if(totalDecayChannels < getNrChannels()) {
		printErr << "total number of decay channels is smaller than number of proper decay channels." << std::endl;
		throw;
	}
	checkSize(_ratio,
	          totalDecayChannels, "total number of decay channels is not correct for branching ratios.");
	checkSize(_masses,
	          totalDecayChannels, "total number of decay channels is not correct for mass points of phase-space integrals.");
	checkSize(_values,
	          totalDecayChannels, "total number of decay channels is not correct for values of phase-space integrals.");

	for(size_t idxDecayChannel = 0; idxDecayChannel < totalDecayChannels; ++idxDecayChannel) {
		_interpolator.push_back(std::make_shared<ROOT::Math::Interpolator>(_masses[idxDecayChannel], _values[idxDecayChannel], ROOT::Math::Interpolation::kLINEAR));
	}

	const double sum = std::accumulate(_ratio.begin(), _ratio.end(), 0.0);
	std::transform(_ratio.begin(), _ratio.end(), _ratio.begin(), [&sum](const double v){ return v/sum; });
}


std::complex<double>
rpwa::resonanceFit::integralWidthBreitWigner::val(const rpwa::resonanceFit::parameters& fitParameters,
                                                  const size_t /*idxBin*/,
                                                  const double mass) const
{
	const double& m0 = fitParameters.getParameter(getId(), 0);
	const double& gamma0 = fitParameters.getParameter(getId(), 1);

	double gamma = 0.;
	for(size_t i = 0; i < _ratio.size(); ++i) {
		// save some time not calculating stuff that is ignored
		if(_ratio[i] == 0.) {
			continue;
		}

		const double ps = _interpolator[i]->Eval(mass);
		const double ps0 = _interpolator[i]->Eval(m0);

		gamma += _ratio[i] * ps / ps0;
	}
	gamma *= gamma0;

	const std::complex<double> component = gamma0*m0 / std::complex<double>(m0*m0-mass*mass, -gamma*m0);

	return component;
}


std::vector<rpwa::resonanceFit::parameter>
rpwa::resonanceFit::constantBackground::getDefaultParameters()
{
	std::vector<rpwa::resonanceFit::parameter> defaultParameters(0);

	return defaultParameters;
}


rpwa::resonanceFit::constantBackground::constantBackground(const size_t id,
                                                           const std::string& name,
                                                           const std::vector<rpwa::resonanceFit::parameter>& parameters,
                                                           const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
                                                           const std::vector<size_t>& nrMassBins,
                                                           const boost::multi_array<double, 2>& massBinCenters,
                                                           const bool useBranchings)
	: component(id,
	            name,
	            parameters,
	            decayChannels,
	            nrMassBins,
	            massBinCenters,
	            useBranchings)
{
	if(getNrChannels() != 1) {
		printErr << "component '" << getName() << "' of type 'constantBackground' needs to have exactly one decay channel." << std::endl;
		throw;
	}

	// initialize the bins that can be used for the caching, the function
	// value for this component type is the same in bins with an equal
	// mass binning
	_binsEqualValues = initBinsEqualValuesFromMassBinning(nrMassBins, massBinCenters);

	// check presence of all parameters
	checkParameters(*this);
}


std::complex<double>
rpwa::resonanceFit::constantBackground::val(const rpwa::resonanceFit::parameters& /*fitParameters*/,
                                            const size_t /*idxBin*/,
                                            const double /*m*/) const
{
	return 1.;
}


std::vector<rpwa::resonanceFit::parameter>
rpwa::resonanceFit::exponentialBackground::getDefaultParameters()
{
	std::vector<rpwa::resonanceFit::parameter> defaultParameters(1);
	defaultParameters[0].setName("g");
	defaultParameters[0].setStep(1.0);

	return defaultParameters;
}


rpwa::resonanceFit::exponentialBackground::exponentialBackground(const size_t id,
                                                                 const std::string& name,
                                                                 const std::vector<rpwa::resonanceFit::parameter>& parameters,
                                                                 const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
                                                                 const std::vector<size_t>& nrMassBins,
                                                                 const boost::multi_array<double, 2>& massBinCenters,
                                                                 const bool useBranchings,
                                                                 const int relAngularMom,
                                                                 const double mIsobar1,
                                                                 const double mIsobar2,
                                                                 const double exponent)
	: component(id,
	            name,
	            parameters,
	            decayChannels,
	            nrMassBins,
	            massBinCenters,
	            useBranchings),
	  _l(relAngularMom),
	  _m1(mIsobar1),
	  _m2(mIsobar2),
	  _exponent(exponent)
{
	if(getNrChannels() != 1) {
		printErr << "component '" << getName() << "' of type 'exponentialBackground' needs to have exactly one decay channel." << std::endl;
		throw;
	}

	// initialize the bins that can be used for the caching, the function
	// value for this component type is the same in bins with an equal
	// mass binning
	_binsEqualValues = initBinsEqualValuesFromMassBinning(nrMassBins, massBinCenters);

	// check presence of all parameters
	checkParameters(*this);

	// select the maximum of the mass only from the used bins
	const std::vector<size_t>& bins = getChannel(0).getBins();
	std::vector<double> maxMasses(bins.size());
	for(size_t i = 0; i < bins.size(); ++i) {
		const size_t idxBin = bins[i];
		maxMasses[i] = massBinCenters[idxBin][nrMassBins[idxBin] - 1];
	}
	const double maxMass = *std::max_element(maxMasses.begin(), maxMasses.end());
	const double q = rpwa::breakupMomentum(maxMass, _m1, _m2);
	const double f2 = rpwa::barrierFactorSquared(2*_l, q);
	_norm = 1. / (q*f2);
}


std::complex<double>
rpwa::resonanceFit::exponentialBackground::val(const rpwa::resonanceFit::parameters& fitParameters,
                                               const size_t /*idxBin*/,
                                               const double mass) const
{
	// calculate breakup momentum
	if(mass < _m1+_m2) {
		return std::complex<double>(1,0);
	}
	const double q = rpwa::breakupMomentum(mass, _m1, _m2);
	const double f2 = rpwa::barrierFactorSquared(2*_l, q);
	const double c = std::pow(q*f2 * _norm, _exponent);

	const std::complex<double> component = exp(-fitParameters.getParameter(getId(), 0)*c);

	return component;
}


std::ostream&
rpwa::resonanceFit::exponentialBackground::print(std::ostream& out, const bool newLine) const
{
	component::print(out, true);

	out << "    mass of isobar 1: " << _m1 << " GeV/c^2, mass of isobar 2: " << _m2 << " GeV/c^2" << std::endl;
	out << "    relative orbital angular momentum between isobars: " << _l << " (in units of hbar)" << std::endl;
	out << "    exponent of break-up momentum times barrier-factor squared: " << _exponent;

	if(newLine) {
		out << std::endl;
	}

	return out;
}


std::vector<rpwa::resonanceFit::parameter>
rpwa::resonanceFit::tPrimeDependentBackground::getDefaultParameters()
{
	std::vector<rpwa::resonanceFit::parameter> defaultParameters(5);
	defaultParameters[0].setName("m0");
	defaultParameters[0].setStep(0.001);
	defaultParameters[1].setName("c0");
	defaultParameters[1].setStep(0.001);
	defaultParameters[2].setName("c1");
	defaultParameters[2].setStep(1.0);
	defaultParameters[3].setName("c2");
	defaultParameters[3].setStep(1.0);
	defaultParameters[4].setName("c3");
	defaultParameters[4].setStep(1.0);

	return defaultParameters;
}


rpwa::resonanceFit::tPrimeDependentBackground::tPrimeDependentBackground(const size_t id,
                                                                         const std::string& name,
                                                                         const std::vector<rpwa::resonanceFit::parameter>& parameters,
                                                                         const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
                                                                         const std::vector<size_t>& nrMassBins,
                                                                         const boost::multi_array<double, 2>& massBinCenters,
                                                                         const bool useBranchings,
                                                                         const std::vector<double>& tPrimeMeans,
                                                                         const int relAngularMom,
                                                                         const double mIsobar1,
                                                                         const double mIsobar2,
                                                                         const double exponent)
	: component(id,
	            name,
	            parameters,
	            decayChannels,
	            nrMassBins,
	            massBinCenters,
	            useBranchings),
	  _tPrimeMeans(tPrimeMeans),
	  _l(relAngularMom),
	  _m1(mIsobar1),
	  _m2(mIsobar2),
	  _exponent(exponent)
{
	if(getNrChannels() != 1) {
		printErr << "component '" << getName() << "' of type 'tPrimeDependentBackground' needs to have exactly one decay channel." << std::endl;
		throw;
	}

	// check presence of all parameters
	checkParameters(*this);

	const size_t nrBins = nrMassBins.size();
	checkSize(_tPrimeMeans,
	          nrBins, "number of bins is not correct for mean t' value per bin.");

	// select the maximum of the mass only from the used bins
	const std::vector<size_t>& bins = getChannel(0).getBins();
	std::vector<double> maxMasses(bins.size());
	for(size_t i = 0; i < bins.size(); ++i) {
		const size_t idxBin = bins[i];
		maxMasses[i] = massBinCenters[idxBin][nrMassBins[idxBin] - 1];
	}
	const double maxMass = *std::max_element(maxMasses.begin(), maxMasses.end());
	const double q = rpwa::breakupMomentum(maxMass, _m1, _m2);
	const double f2 = rpwa::barrierFactorSquared(2*_l, q);
	_norm = 1. / (q*f2);
}


std::complex<double>
rpwa::resonanceFit::tPrimeDependentBackground::val(const rpwa::resonanceFit::parameters& fitParameters,
                                                   const size_t idxBin,
                                                   const double mass) const
{
	// calculate breakup momentum
	if(mass < _m1+_m2) {
		return std::pow(mass - fitParameters.getParameter(getId(), 0), fitParameters.getParameter(getId(), 1));
	}
	const double q = rpwa::breakupMomentum(mass, _m1, _m2);
	const double f2 = rpwa::barrierFactorSquared(2*_l, q);
	const double c = std::pow(q*f2 * _norm, _exponent);

	// get mean t' value for current bin
	const double tPrime = _tPrimeMeans[idxBin];
	const double tPrimePol = fitParameters.getParameter(getId(), 2) + fitParameters.getParameter(getId(), 3)*tPrime + fitParameters.getParameter(getId(), 4)*tPrime*tPrime;

	const double mPre = std::pow(mass - fitParameters.getParameter(getId(), 0), fitParameters.getParameter(getId(), 1));

	const std::complex<double> component = mPre * exp(-tPrimePol*c);

	return component;
}


std::ostream&
rpwa::resonanceFit::tPrimeDependentBackground::print(std::ostream& out, const bool newLine) const
{
	component::print(out, true);

	out << "    mass of isobar 1: " << _m1 << " GeV/c^2, mass of isobar 2: " << _m2 << " GeV/c^2" << std::endl;
	out << "    relative orbital angular momentum between isobars: " << _l << " (in units of hbar)" << std::endl;
	out << "    exponent of break-up momentum times barrier-factor squared: " << _exponent << std::endl;

	out << "    for " << _tPrimeMeans.size() << " bin" << ((_tPrimeMeans.size()>1)?"s":"") << " with mean t' value" << ((_tPrimeMeans.size()>1)?"s":"") << ": " << _tPrimeMeans[0];
	for(size_t i = 1; i < _tPrimeMeans.size(); ++i) {
		out << ", " << _tPrimeMeans[i];
	}

	if(newLine) {
		out << std::endl;
	}

	return out;
}


std::vector<rpwa::resonanceFit::parameter>
rpwa::resonanceFit::exponentialBackgroundIntegral::getDefaultParameters()
{
	std::vector<rpwa::resonanceFit::parameter> defaultParameters(1);
	defaultParameters[0].setName("g");
	defaultParameters[0].setStep(1.0);

	return defaultParameters;
}


rpwa::resonanceFit::exponentialBackgroundIntegral::exponentialBackgroundIntegral(const size_t id,
                                                                                 const std::string& name,
                                                                                 const std::vector<rpwa::resonanceFit::parameter>& parameters,
                                                                                 const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
                                                                                 const std::vector<size_t>& nrMassBins,
                                                                                 const boost::multi_array<double, 2>& massBinCenters,
                                                                                 const bool useBranchings,
                                                                                 const std::vector<double>& masses,
                                                                                 const std::vector<double>& values,
                                                                                 const double exponent)
	: component(id,
	            name,
	            parameters,
	            decayChannels,
	            nrMassBins,
	            massBinCenters,
	            useBranchings),
	  _masses(masses),
	  _values(values),
	  _exponent(exponent)
{
	if(getNrChannels() != 1) {
		printErr << "component '" << getName() << "' of type 'exponentialBackgroundIntegral' needs to have exactly one decay channel." << std::endl;
		throw;
	}

	// initialize the bins that can be used for the caching, the function
	// value for this component type is the same in bins with an equal
	// mass binning
	_binsEqualValues = initBinsEqualValuesFromMassBinning(nrMassBins, massBinCenters);

	// check presence of all parameters
	checkParameters(*this);

	assert(not _interpolator);
	_interpolator = std::make_shared<ROOT::Math::Interpolator>(_masses, _values, ROOT::Math::Interpolation::kLINEAR);

	// select the maximum of the mass only from the used bins
	const std::vector<size_t>& bins = getChannel(0).getBins();
	std::vector<double> maxMasses(bins.size());
	for(size_t i = 0; i < bins.size(); ++i) {
		const size_t idxBin = bins[i];
		maxMasses[i] = massBinCenters[idxBin][nrMassBins[idxBin] - 1];
	}
	const double maxMass = *std::max_element(maxMasses.begin(), maxMasses.end());
	_norm = 1. / (maxMass * _interpolator->Eval(maxMass));
}


std::complex<double>
rpwa::resonanceFit::exponentialBackgroundIntegral::val(const rpwa::resonanceFit::parameters& fitParameters,
                                                       const size_t /*idxBin*/,
                                                       const double mass) const
{
	const double ps = _interpolator->Eval(mass);
	const double c = std::pow(mass * ps * _norm, _exponent);

	const std::complex<double> component = exp(-fitParameters.getParameter(getId(), 0)*c);

	return component;
}


std::ostream&
rpwa::resonanceFit::exponentialBackgroundIntegral::print(std::ostream& out, const bool newLine) const
{
	component::print(out, true);

	out << "    exponent of phase-space integral: " << _exponent;

	if(newLine) {
		out << std::endl;
	}

	return out;
}


std::vector<rpwa::resonanceFit::parameter>
rpwa::resonanceFit::tPrimeDependentBackgroundIntegral::getDefaultParameters()
{
	std::vector<rpwa::resonanceFit::parameter> defaultParameters(5);
	defaultParameters[0].setName("m0");
	defaultParameters[0].setStep(0.001);
	defaultParameters[1].setName("c0");
	defaultParameters[1].setStep(0.001);
	defaultParameters[2].setName("c1");
	defaultParameters[2].setStep(1.0);
	defaultParameters[3].setName("c2");
	defaultParameters[3].setStep(1.0);
	defaultParameters[4].setName("c3");
	defaultParameters[4].setStep(1.0);

	return defaultParameters;
}


rpwa::resonanceFit::tPrimeDependentBackgroundIntegral::tPrimeDependentBackgroundIntegral(const size_t id,
                                                                                         const std::string& name,
                                                                                         const std::vector<rpwa::resonanceFit::parameter>& parameters,
                                                                                         const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
                                                                                         const std::vector<size_t>& nrMassBins,
                                                                                         const boost::multi_array<double, 2>& massBinCenters,
                                                                                         const bool useBranchings,
                                                                                         const std::vector<double>& tPrimeMeans,
                                                                                         const std::vector<double>& masses,
                                                                                         const std::vector<double>& values,
                                                                                         const double exponent)
	: component(id,
	            name,
	            parameters,
	            decayChannels,
	            nrMassBins,
	            massBinCenters,
	            useBranchings),
	  _tPrimeMeans(tPrimeMeans),
	  _masses(masses),
	  _values(values),
	  _exponent(exponent)
{
	if(getNrChannels() != 1) {
		printErr << "component '" << getName() << "' of type 'tPrimeDependentBackgroundIntegral' needs to have exactly one decay channel." << std::endl;
		throw;
	}

	// check presence of all parameters
	checkParameters(*this);

	const size_t nrBins = nrMassBins.size();
	checkSize(_tPrimeMeans,
	          nrBins, "number of bins is not correct for mean t' value per bin.");

	assert(not _interpolator);
	_interpolator = std::make_shared<ROOT::Math::Interpolator>(_masses, _values, ROOT::Math::Interpolation::kLINEAR);

	// select the maximum of the mass only from the used bins
	const std::vector<size_t>& bins = getChannel(0).getBins();
	std::vector<double> maxMasses(bins.size());
	for(size_t i = 0; i < bins.size(); ++i) {
		const size_t idxBin = bins[i];
		maxMasses[i] = massBinCenters[idxBin][nrMassBins[idxBin] - 1];
	}
	const double maxMass = *std::max_element(maxMasses.begin(), maxMasses.end());
	_norm = 1. / (maxMass * _interpolator->Eval(maxMass));
}


std::complex<double>
rpwa::resonanceFit::tPrimeDependentBackgroundIntegral::val(const rpwa::resonanceFit::parameters& fitParameters,
                                                           const size_t idxBin,
                                                           const double mass) const
{
	const double ps = _interpolator->Eval(mass);
	const double c = std::pow(mass * ps * _norm, _exponent);

	// get mean t' value for current bin
	const double tPrime = _tPrimeMeans[idxBin];
	const double tPrimePol = fitParameters.getParameter(getId(), 2) + fitParameters.getParameter(getId(), 3)*tPrime + fitParameters.getParameter(getId(), 4)*tPrime*tPrime;

	const double mPre = std::pow(mass - fitParameters.getParameter(getId(), 0), fitParameters.getParameter(getId(), 1));

	const std::complex<double> component = mPre * exp(-tPrimePol*c);

	return component;
}


std::ostream&
rpwa::resonanceFit::tPrimeDependentBackgroundIntegral::print(std::ostream& out, const bool newLine) const
{
	component::print(out, true);

	out << "    exponent of phase-space integral: " << _exponent << std::endl;

	out << "    for " << _tPrimeMeans.size() << " bin" << ((_tPrimeMeans.size()>1)?"s":"") << " with mean t' value" << ((_tPrimeMeans.size()>1)?"s":"") << ": " << _tPrimeMeans[0];
	for(size_t i = 1; i < _tPrimeMeans.size(); ++i) {
		out << ", " << _tPrimeMeans[i];
	}

	if(newLine) {
		out << std::endl;
	}

	return out;
}
