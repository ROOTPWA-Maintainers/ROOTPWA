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

#include "massDepFitModel.h"

#include "massDepFitComponents.h"
#include "massDepFitFsmd.h"
#include "reportingUtils.hpp"


rpwa::massDepFit::model::model()
	: _nrParameters(0),
	  _fsmd(NULL),
	  _useBranchings(false),
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

	// number of coupling parameters
	_nrParameters += 2 * comp->getNrCouplings() * comp->getChannel(0).getNrBins();

	// number of branching parameters (first branching is always real and fixed to 1)
	if(_useBranchings && comp->getNrChannels() > 1) {
		_nrParameters += 2 * comp->getNrBranchings() - 2;
	}
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
rpwa::massDepFit::model::getParameters(double* par) const
{
	size_t parcount=0;

	// couplings
	for(size_t idxComponent=0; idxComponent<_components.size(); ++idxComponent) {
		parcount += _components[idxComponent]->getCouplings(&par[parcount]);
	}

	// branchings
	if(_useBranchings) {
		for(size_t idxComponent=0; idxComponent<_components.size(); ++idxComponent){
			parcount += _components[idxComponent]->getBranchings(&par[parcount]);
		}
	}

	// parameters
	for(size_t idxComponent=0; idxComponent<_components.size(); ++idxComponent) {
		parcount += _components[idxComponent]->getParameters(&par[parcount]);
	}

	// final-state mass dependence
	if(_fsmd != NULL) {
		parcount += _fsmd->getParameters(&par[parcount]);
	}
}


void
rpwa::massDepFit::model::setParameters(const double* par)
{
	size_t parcount=0;

	// couplings
	for(size_t idxComponent=0; idxComponent<_components.size(); ++idxComponent){
		parcount += _components[idxComponent]->setCouplings(&par[parcount]);
	}

	// branchings
	if(_useBranchings) {
		for(size_t idxComponent=0; idxComponent<_components.size(); ++idxComponent){
			parcount += _components[idxComponent]->setBranchings(&par[parcount]);
		}
	}

	// parameters
	for(size_t idxComponent=0; idxComponent<_components.size(); ++idxComponent){
		parcount += _components[idxComponent]->setParameters(&par[parcount]);
	}

	// final-state mass-dependence
	if(_fsmd != NULL) {
		parcount += _fsmd->setParameters(&par[parcount]);
	}
}


std::complex<double>
rpwa::massDepFit::model::productionAmplitude(const size_t idxWave,
                                             const size_t idxBin,
                                             const double mass,
                                             const size_t idxMass) const
{
	// loop over all components and pick up those that contribute to this channels
	std::complex<double> prodAmp(0., 0.);

	// get entry from mapping
	const std::vector<std::pair<size_t, size_t> >& components = _waveComponentChannel[idxWave];
	size_t nrComponents = components.size();

	for(unsigned int idxComponents=0; idxComponents<nrComponents; ++idxComponents) {
		size_t idxComponent = components[idxComponents].first;
		size_t idxChannel = components[idxComponents].second;
		prodAmp += _components[idxComponent]->val(idxBin, mass) * _components[idxComponent]->getCouplingPhaseSpace(idxChannel, idxBin, mass, idxMass);
	}

	if(_fsmd != NULL) {
		prodAmp *= _fsmd->val(mass, idxMass);
	}

	return prodAmp;
}


double
rpwa::massDepFit::model::intensity(const size_t idxWave,
                                   const size_t idxBin,
                                   const double mass,
                                   const size_t idxMass) const
{
	const std::complex<double> prodAmp = productionAmplitude(idxWave, idxBin, mass, idxMass);

	return norm(prodAmp);
}


double
rpwa::massDepFit::model::phaseAbsolute(const size_t idxWave,
                                       const size_t idxBin,
                                       const double mass,
                                       const size_t idxMass) const
{
	const std::complex<double> prodAmp = productionAmplitude(idxWave, idxBin, mass, idxMass);

	return arg(prodAmp);
}


std::complex<double>
rpwa::massDepFit::model::spinDensityMatrix(const size_t idxWave,
                                           const size_t jdxWave,
                                           const size_t idxBin,
                                           const double mass,
                                           const size_t idxMass) const
{
	const std::complex<double> prodAmpI = productionAmplitude(idxWave, idxBin, mass, idxMass);
	const std::complex<double> prodAmpJ = productionAmplitude(jdxWave, idxBin, mass, idxMass);

	return prodAmpI * conj(prodAmpJ);
}


double
rpwa::massDepFit::model::phase(const size_t idxWave,
                               const size_t jdxWave,
                               const size_t idxBin,
                               const double mass,
                               const size_t idxMass) const
{
	return arg(spinDensityMatrix(idxWave, jdxWave, idxBin, mass, idxMass));
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
