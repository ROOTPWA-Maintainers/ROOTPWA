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
//      implementation of fit model for resonance fit
//
//-------------------------------------------------------------------------


#include "massDepFitModel.h"

#include "massDepFitCache.h"
#include "massDepFitComponents.h"
#include "massDepFitFsmd.h"
#include "reportingUtils.hpp"


rpwa::massDepFit::model::model(const bool useBranchings)
	: _nrParameters(0),
	  _fsmd(NULL),
	  _maxChannelsInComponent(0),
	  _maxParametersInComponent(0),
	  _useBranchings(useBranchings),
	  _idxAnchorWave(std::numeric_limits<size_t>::max()),
	  _idxAnchorComponent(std::numeric_limits<size_t>::max()),
	  _idxAnchorChannel(std::numeric_limits<size_t>::max())
{
}


rpwa::massDepFit::model::~model()
{
	for(std::vector<rpwa::massDepFit::component*>::iterator it=_components.begin(); it!=_components.end(); ++it) {
		delete *it;
	}
	_components.clear();

	if(_fsmd != NULL) {
		delete _fsmd;
		_fsmd = NULL;
	}
}


bool
rpwa::massDepFit::model::init(const std::vector<std::string>& waveNames,
                              const std::string& anchorWaveName,
                              const std::string& anchorComponentName)
{
	_waveNames = waveNames;

	if(not initMapping(anchorWaveName, anchorComponentName)) {
		printErr << "error while mapping the waves to the decay channels and components." << std::endl;
		return false;
	}

	return true;
}


void
rpwa::massDepFit::model::add(rpwa::massDepFit::component* comp)
{
	_components.push_back(comp);

	// number of resonance parameters
	_nrParameters += comp->getNrParameters();
	_maxParametersInComponent = std::max(_maxParametersInComponent, comp->getNrParameters());

	// number of coupling parameters
	_nrParameters += 2 * comp->getNrCouplings() * comp->getChannel(0).getNrBins();

	// number of branching parameters (first branching is always real and fixed to 1)
	if(_useBranchings && comp->getNrChannels() > 1) {
		_nrParameters += 2 * comp->getNrBranchings() - 2;
	}

	// maximum number of channels in one component
	_maxChannelsInComponent = std::max(_maxChannelsInComponent, comp->getNrChannels());
}


void
rpwa::massDepFit::model::setFsmd(rpwa::massDepFit::fsmd* fsmd)
{
	if(_fsmd != NULL) {
		_nrParameters -= _fsmd->getNrParameters();
		delete _fsmd;
	}

	_fsmd = fsmd;

	if(_fsmd != NULL) {
		_nrParameters += _fsmd->getNrParameters();
		_maxParametersInComponent = std::max(_maxParametersInComponent, _fsmd->getNrParameters());
	}
}


// performs mapping from the index of a wave in wavelist() to the components and channels that couple to this wave
bool
rpwa::massDepFit::model::initMapping(const std::string& anchorWaveName,
                                     const std::string& anchorComponentName)
{
	// check that all waves used in a decay channel have been defined
	const size_t nrComponents = _components.size();
	for(size_t idxComponent=0; idxComponent<nrComponents; ++idxComponent) {
		for(std::vector<rpwa::massDepFit::channel>::const_iterator itChan=_components[idxComponent]->getChannels().begin(); itChan !=_components[idxComponent]->getChannels().end(); ++itChan) {
			if(find(_waveNames.begin(), _waveNames.end(), itChan->getWaveName()) == _waveNames.end()) {
				printErr << "wave '" << itChan->getWaveName() << "' not known in decay of '" << _components[idxComponent]->getName() << "'." << std::endl;
				return false;
			}
		}
	}

	// check that all defined waves are also used
	for(std::vector<std::string>::const_iterator itWave=_waveNames.begin(); itWave!=_waveNames.end(); ++itWave) {
		bool found(false);
		for(size_t idxComponent=0; idxComponent<nrComponents; ++idxComponent) {
			for(std::vector<rpwa::massDepFit::channel>::const_iterator itChan=_components[idxComponent]->getChannels().begin(); itChan !=_components[idxComponent]->getChannels().end(); ++itChan) {
				if(itChan->getWaveName() == *itWave) {
					found = true;
					break;
				}
			}
		}

		if(not found) {
			printErr << "wave '" << *itWave << "' defined but not used in any decay." << std::endl;
			return false;
		}
	}

	_waveComponentChannel.resize(_waveNames.size());
	for(size_t idxWave=0; idxWave<_waveNames.size(); ++idxWave) {
		// check which components this waves belongs to and which channel they are
		for(size_t idxComponent=0; idxComponent<nrComponents; ++idxComponent) {
			// loop over channels of component and see if wave is there
			const size_t nrChannels = _components[idxComponent]->getNrChannels();
			for(size_t idxChannel=0; idxChannel<nrChannels; ++idxChannel) {
				if(_components[idxComponent]->getChannelWaveName(idxChannel) == _waveNames[idxWave]) {
					_waveComponentChannel[idxWave].push_back(std::pair<size_t, size_t>(idxComponent,idxChannel));

					if(anchorWaveName == _waveNames[idxWave] && anchorComponentName == _components[idxComponent]->getName()) {
						_idxAnchorWave = idxWave;
						_idxAnchorComponent = idxComponent;
						_idxAnchorChannel = idxChannel;
					}
				}
			} // end loop over channels
		} // end loop over components
	}

	// test that anchor wave is present
	if(_idxAnchorWave == std::numeric_limits<size_t>::max() || _idxAnchorComponent == std::numeric_limits<size_t>::max() || _idxAnchorChannel == std::numeric_limits<size_t>::max()) {
		printErr << "anchor wave '" << anchorWaveName << "' in component '" << anchorComponentName << "' not found." << std::endl;
		return false;
	}
	if(_idxAnchorChannel != 0) {
		printErr << "anchor wave '" << anchorWaveName << "' has to be first channel in anchor component '" << anchorComponentName << "'." << std::endl;
		return false;
	}
	_components[_idxAnchorComponent]->setChannelAnchor(_idxAnchorChannel, true);
	_nrParameters -= _components[_idxAnchorComponent]->getChannel(_idxAnchorChannel).getNrBins();

	std::ostringstream output;
	for(size_t idxWave=0; idxWave<_waveNames.size(); idxWave++) {
		output << "    wave '" << _waveNames[idxWave] << "' (index " << idxWave << ") used in" << std::endl;
		for(size_t idxComponents=0; idxComponents<_waveComponentChannel[idxWave].size(); idxComponents++) {
			const size_t idxComponent = _waveComponentChannel[idxWave][idxComponents].first;
			output << "        component " << idxComponents << ": " << idxComponent << " '" << _components[idxComponent]->getName() << "'";
			if(_idxAnchorWave == idxWave && _idxAnchorComponent == idxComponent) {
				output << " (anchor)";
			}
			output << std::endl;
		}
	}
	for(size_t idxComponent=0; idxComponent<nrComponents; ++idxComponent) {
		output << "    component '" << _components[idxComponent]->getName() << "' (index " << idxComponent << ") used in" << std::endl;

		const size_t nrChannels = _components[idxComponent]->getNrChannels();
		for(size_t idxChannel=0; idxChannel<nrChannels; ++idxChannel) {
			output << "        channel " << idxChannel << ": " << _components[idxComponent]->getChannelWaveName(idxChannel)
			       << (_components[idxComponent]->getChannel(idxChannel).isAnchor() ? " (anchor)" : "") << std::endl;
		}
	}
	printInfo << _waveNames.size() << " waves and " << _components.size() << " components in fit model:" << std::endl
	          << output.str();

	return true;
}


void
rpwa::massDepFit::model::importParameters(const double* par,
                                          rpwa::massDepFit::parameters& parameters,
                                          rpwa::massDepFit::cache& cache) const
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
	if(_fsmd != NULL) {
		parcount += _fsmd->importParameters(&par[parcount], parameters);
	}
}


std::complex<double>
rpwa::massDepFit::model::productionAmplitude(const rpwa::massDepFit::parameters& fitParameters,
                                             rpwa::massDepFit::cache& cache,
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
	const std::vector<std::pair<size_t, size_t> >& components = _waveComponentChannel[idxWave];
	size_t nrComponents = components.size();

	for(unsigned int idxComponents=0; idxComponents<nrComponents; ++idxComponents) {
		size_t idxComponent = components[idxComponents].first;
		size_t idxChannel = components[idxComponents].second;
		prodAmp += _components[idxComponent]->val(fitParameters, cache, idxBin, mass, idxMass) * _components[idxComponent]->getCouplingPhaseSpace(fitParameters, idxChannel, idxBin, mass, idxMass);
	}

	if(_fsmd != NULL) {
		prodAmp *= _fsmd->val(fitParameters, mass, idxMass);
	}

	if (idxMass != std::numeric_limits<size_t>::max()) {
		cache.setProdAmp(idxWave, idxBin, idxMass, prodAmp);
	}

	return prodAmp;
}


double
rpwa::massDepFit::model::intensity(const rpwa::massDepFit::parameters& fitParameters,
                                   rpwa::massDepFit::cache& cache,
                                   const size_t idxWave,
                                   const size_t idxBin,
                                   const double mass,
                                   const size_t idxMass) const
{
	const std::complex<double> prodAmp = productionAmplitude(fitParameters, cache, idxWave, idxBin, mass, idxMass);

	return norm(prodAmp);
}


double
rpwa::massDepFit::model::phaseAbsolute(const rpwa::massDepFit::parameters& fitParameters,
                                       rpwa::massDepFit::cache& cache,
                                       const size_t idxWave,
                                       const size_t idxBin,
                                       const double mass,
                                       const size_t idxMass) const
{
	const std::complex<double> prodAmp = productionAmplitude(fitParameters, cache, idxWave, idxBin, mass, idxMass);

	return arg(prodAmp);
}


std::complex<double>
rpwa::massDepFit::model::spinDensityMatrix(const rpwa::massDepFit::parameters& fitParameters,
                                           rpwa::massDepFit::cache& cache,
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
rpwa::massDepFit::model::phase(const rpwa::massDepFit::parameters& fitParameters,
                               rpwa::massDepFit::cache& cache,
                               const size_t idxWave,
                               const size_t jdxWave,
                               const size_t idxBin,
                               const double mass,
                               const size_t idxMass) const
{
	return arg(spinDensityMatrix(fitParameters, cache, idxWave, jdxWave, idxBin, mass, idxMass));
}


std::ostream&
rpwa::massDepFit::model::print(std::ostream& out) const
{
	for(unsigned int i=0;i<_components.size();++i){
		const rpwa::massDepFit::component& c = *_components[i];
		c.print(out);
	}
	return out;
}
