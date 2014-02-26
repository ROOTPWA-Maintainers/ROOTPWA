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

#include <TF1.h>

#include "massDepFitComponents.h"
#include "reportingUtils.hpp"

using namespace std;
using namespace rpwa;


rpwa::massDepFit::model::model()
	: _nrParameters(0),
	  _idxAnchorWave(numeric_limits<size_t>::max()),
	  _idxAnchorComponent(numeric_limits<size_t>::max()),
	  _idxAnchorChannel(numeric_limits<size_t>::max()),
	  _fsmdFunction(NULL),
	  _fsmdFixed(false)
{
}


bool
rpwa::massDepFit::model::init(const vector<string>& waveNames,
                              const vector<double>& massBinCenters,
                              const std::string& anchorWaveName,
                              const std::string& anchorComponentName)
{
	_waveNames = waveNames;

	if(not initMapping(anchorWaveName, anchorComponentName)) {
		printErr << "error while mapping the waves to the decay channels and components." << endl;
		return false;
	}

	if(not initFsmd(massBinCenters)) {
		printErr << "error precalculating the final-state mass dependence." << endl;
		return false;
	}

	return true;
}


void
rpwa::massDepFit::model::add(rpwa::massDepFit::component* comp)
{
	_components.push_back(comp);
	_nrParameters += comp->getNrParameters() + 2*comp->getNrChannels();
}


void
rpwa::massDepFit::model::setFsmdFunction(TF1* fsmdFunction)
{
	_fsmdFunction = fsmdFunction;

	// clear list of free parameters
	_fsmdFreeParameterIndices.clear();

	if (_fsmdFunction == NULL) {
		_fsmdFixed = true;
		return;
	}

	// check if there are free parameters in the phase space that should be fitted
	unsigned int nparFsmd = _fsmdFunction->GetNpar();
	// loop over parameters and check limits
	// remember which parameters to let float
	for(unsigned int i=0;i<nparFsmd;++i) {
		double min, max;
		_fsmdFunction->GetParLimits(i,min,max);
		if(min!=max) {
			_fsmdFreeParameterIndices.push_back(i);
		}
	}// end loop over parameters

	_nrParameters += _fsmdFreeParameterIndices.size();
	if(_fsmdFreeParameterIndices.size() == 0) {
		_fsmdFixed = true;
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
		for(vector<rpwa::massDepFit::channel>::const_iterator itChan=_components[idxComponent]->getChannels().begin(); itChan !=_components[idxComponent]->getChannels().end(); ++itChan) {
			if(find(_waveNames.begin(), _waveNames.end(), itChan->getWaveName()) == _waveNames.end()) {
				printErr << "wave '" << itChan->getWaveName() << "' not known in decay of '" << _components[idxComponent]->getName() << "'." << endl;
				return false;
			}
		}
	}

	// check that all defined waves are also used
	for(vector<string>::const_iterator itWave=_waveNames.begin(); itWave!=_waveNames.end(); ++itWave) {
		bool found(false);
		for(size_t idxComponent=0; idxComponent<nrComponents; ++idxComponent) {
			for(vector<rpwa::massDepFit::channel>::const_iterator itChan=_components[idxComponent]->getChannels().begin(); itChan !=_components[idxComponent]->getChannels().end(); ++itChan) {
				if(itChan->getWaveName() == *itWave) {
					found = true;
					break;
				}
			}
		}

		if(not found) {
			printErr << "wave '" << *itWave << "' defined but not used in any decay." << endl;
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
					_waveComponentChannel[idxWave].push_back(pair<size_t, size_t>(idxComponent,idxChannel));

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
	if(_idxAnchorWave == numeric_limits<size_t>::max() || _idxAnchorComponent == numeric_limits<size_t>::max() || _idxAnchorChannel == numeric_limits<size_t>::max()) {
		printErr << "anchor wave '" << anchorWaveName << "' in component '" << anchorComponentName << "' not found." << endl;
		return false;
	}

	ostringstream output;
	for(size_t idxWave=0; idxWave<_waveNames.size(); idxWave++) {
		output << "    wave '" << _waveNames[idxWave] << "' (index " << idxWave << ") used in" << endl;
		for(size_t idxComponents=0; idxComponents<_waveComponentChannel[idxWave].size(); idxComponents++) {
			const size_t idxComponent = _waveComponentChannel[idxWave][idxComponents].first;
			output << "        component " << idxComponents << ": " << idxComponent << " '" << _components[idxComponent]->getName() << "'";
			if(_idxAnchorWave == idxWave && _idxAnchorComponent == idxComponent) {
				output << " (anchor)";
			}
			output << endl;
		}
	}
	for(size_t idxComponent=0; idxComponent<nrComponents; ++idxComponent) {
		output << "    component '" << _components[idxComponent]->getName() << "' (index " << idxComponent << ") used in" << endl;

		const size_t nrChannels = _components[idxComponent]->getNrChannels();
		for(size_t idxChannel=0; idxChannel<nrChannels; ++idxChannel) {
			output << "        channel " << idxChannel << ": " << _components[idxComponent]->getChannelWaveName(idxChannel);
			if(_idxAnchorComponent == idxComponent && _idxAnchorChannel == idxChannel) {
				output << " (anchor)";
			}
			output << endl;
		}
	}
	printInfo << _waveNames.size() << " waves and " << _components.size() << " components in fit model:" << endl
	          << output.str();

	return true;
}


bool
rpwa::massDepFit::model::initFsmd(const vector<double>& massBinCenters)
{
	const size_t nrMassBins = massBinCenters.size();

	_fsmdValues.resize(nrMassBins);

	for(size_t idxMassBin=0; idxMassBin<nrMassBins; ++idxMassBin) {
		_fsmdValues[idxMassBin] = calcFsmd(massBinCenters[idxMassBin], numeric_limits<size_t>::max());
	}

	return true;
}


double 
rpwa::massDepFit::model::getFsmdParameter(const size_t idx) const
{
	return _fsmdFunction->GetParameter(_fsmdFreeParameterIndices[idx]);
}


void 
rpwa::massDepFit::model::getFsmdParameterLimits(const size_t idx, double& lower, double& upper) const
{
	return _fsmdFunction->GetParLimits(_fsmdFreeParameterIndices[idx], lower, upper);
}


void 
rpwa::massDepFit::model::getParameters(double* par) const
{
	size_t parcount=0;
	// components
	for(size_t idxComponent=0; idxComponent<_components.size(); ++idxComponent) {
		_components[idxComponent]->getParameters(&par[parcount]);
		parcount += _components[idxComponent]->getNrParameters();

		_components[idxComponent]->getCouplings(&par[parcount]);
		parcount += 2*_components[idxComponent]->getNrChannels();
	}

	// final-state mass dependence
	unsigned int nfreepar=_fsmdFreeParameterIndices.size();
	for(unsigned int ipar=0; ipar<nfreepar; ++ipar) {
		par[parcount] = _fsmdFunction->GetParameter(_fsmdFreeParameterIndices[ipar]);
		++parcount;
	}
}


void
rpwa::massDepFit::model::setParameters(const double* par)
{
	size_t parcount=0;
	// components
	for(size_t idxComponent=0; idxComponent<_components.size(); ++idxComponent){
		_components[idxComponent]->setParameters(&par[parcount]);
		parcount += _components[idxComponent]->getNrParameters();

		_components[idxComponent]->setCouplings(&par[parcount]);
		parcount += 2*_components[idxComponent]->getNrChannels();
	} // end loop over components

	// final-state mass-dependence
	unsigned int nfreepar=_fsmdFreeParameterIndices.size();
	for(unsigned int ipar=0;ipar<nfreepar;++ipar){
		_fsmdFunction->SetParameter(_fsmdFreeParameterIndices[ipar], par[parcount]);
		++parcount;
	}
}


double 
rpwa::massDepFit::model::calcFsmd(const double mass,
                                  const size_t idxMass) const
{
	if(_fsmdFixed && idxMass != numeric_limits<size_t>::max()) {
		return _fsmdValues[idxMass];
	}

	if(not _fsmdFunction) {
		return 1.;
	}

	return _fsmdFunction->Eval(mass);
}


complex<double>
rpwa::massDepFit::model::productionAmplitude(const size_t idxWave,
                                             const double mass,
                                             const size_t idxMass) const
{
	// loop over all components and pick up those that contribute to this channels
	complex<double> prodAmp(0., 0.);

	// get entry from mapping
	const vector<pair<size_t, size_t> >& components = _waveComponentChannel[idxWave];
	size_t nrComponents = components.size();

	for(unsigned int idxComponents=0; idxComponents<nrComponents; ++idxComponents) {
		size_t idxComponent = components[idxComponents].first;
		size_t idxChannel = components[idxComponents].second;
		prodAmp += _components[idxComponent]->val(mass) * _components[idxComponent]->getChannel(idxChannel).getCouplingPhaseSpace(mass, idxMass);
	}

	prodAmp *= calcFsmd(mass, idxMass);

	return prodAmp;
}


double 
rpwa::massDepFit::model::intensity(const size_t idxWave,
                                   const double mass,
                                   const size_t idxMass) const
{
	const complex<double> prodAmp = productionAmplitude(idxWave, mass, idxMass);

	return norm(prodAmp);
}


double 
rpwa::massDepFit::model::phaseAbsolute(const size_t idxWave,
                                       const double mass,
                                       const size_t idxMass) const
{
	const complex<double> prodAmp = productionAmplitude(idxWave, mass, idxMass);

	return arg(prodAmp);
}


complex<double>
rpwa::massDepFit::model::spinDensityMatrix(const size_t idxWave,
                                           const size_t jdxWave,
                                           const double mass,
                                           const size_t idxMass) const
{
	const complex<double> prodAmpI = productionAmplitude(idxWave, mass, idxMass);
	const complex<double> prodAmpJ = productionAmplitude(jdxWave, mass, idxMass);

	return prodAmpI * conj(prodAmpJ);
}


double 
rpwa::massDepFit::model::phase(const size_t idxWave,
                               const size_t jdxWave,
                               const double mass,
                               const size_t idxMass) const
{
	return arg(spinDensityMatrix(idxWave, jdxWave, mass, idxMass));
}


ostream&
rpwa::massDepFit::model::print(ostream& out) const
{
	for(unsigned int i=0;i<_components.size();++i){
		const rpwa::massDepFit::component& c = *_components[i];
		c.print(out);
	}
	return out;
}
