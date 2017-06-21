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
//      implementation of fit model for resonance fit
//
//-------------------------------------------------------------------------


#include "model.h"

#include <set>

#include <reportingUtils.hpp>

#include "cache.h"
#include "components.h"
#include "fsmd.h"
#include "input.h"
#include "resonanceFitHelper.h"


rpwa::resonanceFit::model::model(const rpwa::resonanceFit::inputConstPtr& fitInput,
                                 const std::vector<rpwa::resonanceFit::componentPtr>& comp,
                                 const rpwa::resonanceFit::fsmdPtr& fsmd,
                                 const std::vector<std::string>& anchorWaveNames,
                                 const std::vector<std::string>& anchorComponentNames)
	: _mappingEqualInAllBins(false),
	  _nrParameters(0),
	  _maxChannelsInComponent(0),
	  _maxParametersInComponent(0),
	  _anchorWaveNames(anchorWaveNames),
	  _anchorComponentNames(anchorComponentNames),
	  _anchorWaveIndices(anchorWaveNames.size(), std::numeric_limits<size_t>::max()),
	  _anchorComponentIndices(anchorComponentNames.size(), std::numeric_limits<size_t>::max())
{
	// check size of same parameters
	const size_t nrBins = fitInput->nrBins();
	if(nrBins == 0) {
		printErr << "number of bins is zero, cannot perform the fit. Aborting..." << std::endl;
		throw;
	}

	checkSize(_anchorWaveNames,
	          nrBins, "number of bins is not correct for anchor wave names.");
	checkSize(_anchorComponentNames,
	          nrBins, "number of bins is not correct for anchor component names.");

	for(size_t idxComponent = 0; idxComponent < comp.size(); ++idxComponent) {
		_components.push_back(comp[idxComponent]);

		// number of resonance parameters
		_nrParameters += comp[idxComponent]->getNrParameters();
		_maxParametersInComponent = std::max(_maxParametersInComponent, comp[idxComponent]->getNrParameters());

		// number of coupling parameters
		for(size_t idxCoupling = 0; idxCoupling < comp[idxComponent]->getNrCouplings(); ++idxCoupling) {
			const rpwa::resonanceFit::component::channel& channel = comp[idxComponent]->getChannelFromCouplingIdx(idxCoupling);
			_nrParameters += 2 * channel.getBins().size();
		}

		// number of branching parameters (some branchings are always real and fixed to 1)
		for(size_t idxBranching = 0; idxBranching < comp[idxComponent]->getNrBranchings(); ++idxBranching) {
			if(not comp[idxComponent]->isBranchingFixed(idxBranching)) {
				_nrParameters += 2;
			}
		}

		// maximum number of channels in one component
		_maxChannelsInComponent = std::max(_maxChannelsInComponent, comp[idxComponent]->getNrChannels());
	}

	_fsmd = fsmd;
	if(_fsmd) {
		size_t sumNrParameters = 0;
		const size_t maxNrBins = _fsmd->isSameFunctionForAllBins() ? 1 : _fsmd->getNrBins();
		for(size_t idxBin = 0; idxBin < maxNrBins; ++idxBin) {
			sumNrParameters += _fsmd->getNrParameters(idxBin);
		}
		_maxParametersInComponent = std::max(_maxParametersInComponent, sumNrParameters);
		_nrParameters += sumNrParameters;
	}

	if(not initMapping(fitInput)) {
		printErr << "error while mapping the waves to the decay channels and components." << std::endl;
		throw;
	}
}


// performs mapping from the index of a wave in wavelist() to the components and channels that couple to this wave
bool
rpwa::resonanceFit::model::initMapping(const rpwa::resonanceFit::inputConstPtr& fitInput)
{
	// check that all waves used in a decay channel have been defined
	const size_t nrComponents = _components.size();
	std::set<std::string> componentNames;
	for(size_t idxComponent=0; idxComponent<nrComponents; ++idxComponent) {
		const componentConstPtr& component = _components[idxComponent];

		// check that each component is only defined once
		if(componentNames.count(component->getName()) != 0) {
			printErr << "component '" << component->getName() << "' defined more than once." << std::endl;
			return false;
		}
		componentNames.insert(component->getName());

		for(size_t idxChannel = 0; idxChannel < component->getNrChannels(); ++idxChannel) {
			const rpwa::resonanceFit::component::channel& channel = component->getChannel(idxChannel);

			const std::vector<size_t>& bins = channel.getBins();
			for(size_t i = 0; i < bins.size(); ++i) {
				const size_t idxBin = bins[i];
				const rpwa::resonanceFit::input::bin& fitInputBin = fitInput->getBin(idxBin);

				bool found = false;
				for(size_t idxWave = 0; idxWave < fitInputBin.nrWaves(); ++idxWave) {
					const rpwa::resonanceFit::input::bin::wave& wave = fitInputBin.getWave(idxWave);

					if(wave.waveName() == channel.getWaveName()) {
						found = true;
					}
				}

				if(not found) {
					printErr << "wave '" << channel.getWaveName() << "' not known in decay of '" << component->getName() << "' in bin " << idxBin << "." << std::endl;
					return false;
				}
			}
		}
	}

	// check that all defined waves are also used
	for(size_t idxBin = 0; idxBin < fitInput->nrBins(); ++idxBin) {
		const rpwa::resonanceFit::input::bin& fitInputBin = fitInput->getBin(idxBin);

		for(size_t idxWave = 0; idxWave < fitInputBin.nrWaves(); ++idxWave) {
			const rpwa::resonanceFit::input::bin::wave& wave = fitInputBin.getWave(idxWave);

			bool found(false);
			for(size_t idxComponent = 0; idxComponent < nrComponents; ++idxComponent) {
				const componentConstPtr& component = _components[idxComponent];
				for(size_t idxChannel = 0; idxChannel < component->getNrChannels(); ++idxChannel) {
					const rpwa::resonanceFit::component::channel& channel = component->getChannel(idxChannel);

					if(channel.getWaveName() == wave.waveName()) {
						found = true;
						break;
					}
				}
			}

			if(not found) {
				printErr << "wave '" << wave.waveName() << "' in bin " << idxBin << " defined but not used in any decay." << std::endl;
				return false;
			}
		}
	}

	// set up mapping from wave to component/channel and find anchor wave for each bin
	_anchorChannelIndices.resize(fitInput->nrBins(), std::numeric_limits<size_t>::max());
	_waveComponentChannel.resize(boost::extents[fitInput->nrBins()][std::max_element(fitInput->bins().begin(), fitInput->bins().end(), [](const rpwa::resonanceFit::input::bin& binMax, const rpwa::resonanceFit::input::bin& bin){ return binMax.nrWaves() < bin.nrWaves(); })->nrWaves()]);
	for(size_t idxBin = 0; idxBin < fitInput->nrBins(); ++idxBin) {
		const rpwa::resonanceFit::input::bin& fitInputBin = fitInput->getBin(idxBin);

		for(size_t idxWave = 0; idxWave < fitInputBin.nrWaves(); ++idxWave) {
			const rpwa::resonanceFit::input::bin::wave& wave = fitInputBin.getWave(idxWave);

			// check which components this waves belongs to and which channel they are
			for(size_t idxComponent = 0; idxComponent < nrComponents; ++idxComponent) {
				const componentConstPtr& component = _components[idxComponent];
				// loop over channels of component and see if wave is there
				for(size_t idxChannel = 0; idxChannel < component->getNrChannels(); ++idxChannel) {
					const rpwa::resonanceFit::component::channel& channel = component->getChannel(idxChannel);

					if(channel.getWaveName() == wave.waveName()) {
						_waveComponentChannel[idxBin][idxWave].push_back(std::pair<size_t, size_t>(idxComponent,idxChannel));

						if(_anchorWaveNames[idxBin] == wave.waveName() and _anchorComponentNames[idxBin] == component->getName()) {
							if(_anchorWaveIndices[idxBin] != std::numeric_limits<size_t>::max()) {
								printErr << "second anchor wave found." << std::endl;
								return false;
							}
							_anchorWaveIndices[idxBin] = idxWave;

							if(_anchorComponentIndices[idxBin] != std::numeric_limits<size_t>::max()) {
								printErr << "second anchor component found." << std::endl;
								return false;
							}
							_anchorComponentIndices[idxBin] = idxComponent;

							if(_anchorChannelIndices[idxBin] != std::numeric_limits<size_t>::max()) {
								printErr << "second anchor channel found." << std::endl;
								return false;
							}
							_anchorChannelIndices[idxBin] = idxChannel;
						}
					}
				} // end loop over channels
			} // end loop over components
		}

		// test that anchor wave is present
		if(_anchorWaveIndices[idxBin] == std::numeric_limits<size_t>::max() or
		   _anchorComponentIndices[idxBin] == std::numeric_limits<size_t>::max() or
		   _anchorChannelIndices[idxBin] == std::numeric_limits<size_t>::max()) {
			printErr << "anchor wave '" << _anchorWaveNames[idxBin] << "' in component '" << _anchorComponentNames[idxBin] << "' not found." << std::endl;
			return false;
		}

		const componentPtr& anchorComponent = _components[_anchorComponentIndices[idxBin]];

		size_t firstChannel = 0;
		// loop over channels of component and see if wave is there
		for(size_t idxChannel = 0; idxChannel < anchorComponent->getNrChannels(); ++idxChannel) {
			const rpwa::resonanceFit::component::channel& channel = anchorComponent->getChannel(idxChannel);
			const std::vector<size_t>& bins = channel.getBins();
			if(std::find(bins.begin(), bins.end(), idxBin) != bins.end()) {
				firstChannel = idxChannel;
				break;
			}
		}

		if(_anchorChannelIndices[idxBin] != firstChannel) {
			printErr << "anchor wave '" << _anchorWaveNames[idxBin] << "' has to be channel number " << firstChannel << " in anchor component '" << _anchorComponentNames[idxBin] << "' in bin " << idxBin << "." << std::endl;
			return false;
		}

		anchorComponent->setChannelAnchor(idxBin, _anchorChannelIndices[idxBin]);
	}
	_nrParameters -= fitInput->nrBins();

	// can we simplify stuff by assuming all bins have the same mapping?
	_mappingEqualInAllBins = true;
	for(size_t idxBin = 1; idxBin < fitInput->nrBins(); ++idxBin) {
		if(_anchorWaveIndices[idxBin] != _anchorWaveIndices[0]) {
			_mappingEqualInAllBins = false;
			break;
		}
		if(_anchorComponentIndices[idxBin] != _anchorComponentIndices[0]) {
			_mappingEqualInAllBins = false;
			break;
		}
		if(_anchorChannelIndices[idxBin] != _anchorChannelIndices[0]) {
			_mappingEqualInAllBins = false;
			break;
		}
		for(size_t idxWave = 0; idxWave < fitInput->getBin(idxBin).nrWaves(); ++idxWave) {
			if(_waveComponentChannel[idxBin][idxWave] != _waveComponentChannel[0][idxWave]) {
				_mappingEqualInAllBins = false;
				break;
			}
		}
		if(not _mappingEqualInAllBins) {
			break;
		}
	}

	return true;
}


void
rpwa::resonanceFit::model::importParameters(const double* par,
                                            rpwa::resonanceFit::parameters& parameters,
                                            rpwa::resonanceFit::cache& cache) const
{
	size_t parcount=0;

	// couplings
	for(size_t idxComponent=0; idxComponent<_components.size(); ++idxComponent){
		parcount += _components[idxComponent]->importCouplings(&par[parcount], parameters, cache);
	}

	// branchings
	for(size_t idxComponent=0; idxComponent<_components.size(); ++idxComponent){
		parcount += _components[idxComponent]->importBranchings(&par[parcount], parameters, cache);
	}

	// parameters
	for(size_t idxComponent=0; idxComponent<_components.size(); ++idxComponent){
		parcount += _components[idxComponent]->importParameters(&par[parcount], parameters, cache);
	}

	// final-state mass-dependence
	if(_fsmd) {
		parcount += _fsmd->importParameters(&par[parcount], parameters, cache);
	}

	assert(_nrParameters == parcount);
}


std::complex<double>
rpwa::resonanceFit::model::productionAmplitude(const rpwa::resonanceFit::parameters& fitParameters,
                                               rpwa::resonanceFit::cache& cache,
                                               const size_t idxWave,
                                               const size_t idxBin,
                                               const double mass,
                                               const size_t idxMass) const
{
	if (idxMass != std::numeric_limits<size_t>::max()) {
		const std::complex<double> prodAmp = cache.getProdAmp(idxWave, idxBin, idxMass);
		if (prodAmp != 0.) {
			return prodAmp;
		}
	}

	// loop over all components and pick up those that contribute to this channels
	std::complex<double> prodAmp(0., 0.);

	// get entry from mapping
	const std::vector<std::pair<size_t, size_t> >& components = _waveComponentChannel[idxBin][idxWave];
	const size_t nrComponents = components.size();

	for(size_t idxComponents = 0; idxComponents < nrComponents; ++idxComponents) {
		const size_t idxComponent = components[idxComponents].first;
		const size_t idxChannel = components[idxComponents].second;
		prodAmp += _components[idxComponent]->val(fitParameters, cache, idxBin, mass, idxMass) * _components[idxComponent]->getCouplingPhaseSpace(fitParameters, cache, idxChannel, idxBin, mass, idxMass);
	}

	if(_fsmd) {
		prodAmp *= _fsmd->val(fitParameters, cache, idxBin, mass, idxMass);
	}

	if (idxMass != std::numeric_limits<size_t>::max()) {
		cache.setProdAmp(idxWave, idxBin, idxMass, prodAmp);
	}

	return prodAmp;
}


double
rpwa::resonanceFit::model::intensity(const rpwa::resonanceFit::parameters& fitParameters,
                                     rpwa::resonanceFit::cache& cache,
                                     const size_t idxWave,
                                     const size_t idxBin,
                                     const double mass,
                                     const size_t idxMass) const
{
	const std::complex<double> prodAmp = productionAmplitude(fitParameters, cache, idxWave, idxBin, mass, idxMass);

	return norm(prodAmp);
}


double
rpwa::resonanceFit::model::phaseAbsolute(const rpwa::resonanceFit::parameters& fitParameters,
                                         rpwa::resonanceFit::cache& cache,
                                         const size_t idxWave,
                                         const size_t idxBin,
                                         const double mass,
                                         const size_t idxMass) const
{
	const std::complex<double> prodAmp = productionAmplitude(fitParameters, cache, idxWave, idxBin, mass, idxMass);

	return arg(prodAmp);
}


std::complex<double>
rpwa::resonanceFit::model::spinDensityMatrix(const rpwa::resonanceFit::parameters& fitParameters,
                                             rpwa::resonanceFit::cache& cache,
                                             const size_t idxWave,
                                             const size_t jdxWave,
                                             const size_t idxBin,
                                             const double mass,
                                             const size_t idxMass) const
{
	const std::complex<double> prodAmpI = productionAmplitude(fitParameters, cache, idxWave, idxBin, mass, idxMass);
	const std::complex<double> prodAmpJ = productionAmplitude(fitParameters, cache, jdxWave, idxBin, mass, idxMass);

	return prodAmpI * conj(prodAmpJ);
}


double
rpwa::resonanceFit::model::phase(const rpwa::resonanceFit::parameters& fitParameters,
                                 rpwa::resonanceFit::cache& cache,
                                 const size_t idxWave,
                                 const size_t jdxWave,
                                 const size_t idxBin,
                                 const double mass,
                                 const size_t idxMass) const
{
	return arg(spinDensityMatrix(fitParameters, cache, idxWave, jdxWave, idxBin, mass, idxMass));
}


std::ostream&
rpwa::resonanceFit::model::print(std::ostream& out, const bool newLine) const
{
	for(unsigned int i=0;i<_components.size();++i){
		const rpwa::resonanceFit::component& c = *_components[i];
		c.print(out, newLine or (i != (_components.size()-1)));
	}

	if(_fsmd) {
		if(not newLine) {
			out << std::endl;
		}
		_fsmd->print(out, newLine);
	}

	return out;
}
