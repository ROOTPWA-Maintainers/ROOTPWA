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


massDepFitModel::massDepFitModel()
	: _numpar(0),
	  _idxAnchorWave(numeric_limits<size_t>::max()),
	  _idxAnchorComponent(numeric_limits<size_t>::max()),
	  _idxAnchorChannel(numeric_limits<size_t>::max()),
	  _fsmdFunction(NULL),
	  _fsmdFixed(false)
{
}


bool
massDepFitModel::init(const vector<string>& waveNames,
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
massDepFitModel::add(massDepFitComponent* comp)
{
    _comp.push_back(comp);
    _numpar+=comp->getNrParameters() + 2*comp->getNrChannels();
}


void
massDepFitModel::setFsmdFunction(TF1* fsmdFunction)
{
  _fsmdFunction=fsmdFunction;
  // clear list of free parameters
  _fsmdFreeParameters.clear();

  if (_fsmdFunction == NULL) {
      _fsmdFixed = true;
    return;
  }

  // check if there are free parameters in the phase space that should be fitted
  unsigned int nparFsmd=_fsmdFunction->GetNpar();
  // loop over parameters and check limits
  // remember which parameters to let float
  for(unsigned int i=0;i<nparFsmd;++i){
    double min,max;
    _fsmdFunction->GetParLimits(i,min,max);
    if(min!=max){
      _fsmdFreeParameters.push_back(i);
      cout << "final-state mass dependence parameter "<< i << " floating in ["
	   << min  << "," << max << "]" << endl;
    }
  }// end loop over parameters
  _numpar+=_fsmdFreeParameters.size();
  if(_fsmdFreeParameters.size() == 0) _fsmdFixed = true;
}


// performs mapping from the index of a wave in wavelist() to the components and channels that couple to this wave
bool
massDepFitModel::initMapping(const std::string& anchorWaveName,
                             const std::string& anchorComponentName)
{
	// check that all waves used in a decay channel have been defined
	for(unsigned int i=0; i<n(); ++i) {
		for(vector<pwachannel>::const_iterator itChan=_comp[i]->getChannels().begin(); itChan !=_comp[i]->getChannels().end(); ++itChan) {
			if(find(_waveNames.begin(), _waveNames.end(), itChan->getWaveName()) == _waveNames.end()) {
				printErr << "wave '" << itChan->getWaveName() << "' not known in decay of '" << _comp[i]->getName() << "'." << endl;
				return false;
			}
		}
	}

	// check that all defined waves are also used
	for(vector<string>::const_iterator itWave=_waveNames.begin(); itWave!=_waveNames.end(); ++itWave) {
		bool found(false);
		for(unsigned int i=0; i<n(); ++i) {
			for(vector<pwachannel>::const_iterator itChan=_comp[i]->getChannels().begin(); itChan !=_comp[i]->getChannels().end(); ++itChan) {
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

	_compChannel.resize(_waveNames.size());
	const size_t nrComponents = n();
	for(size_t idxWave=0; idxWave<_waveNames.size(); ++idxWave) {
		// check which components this waves belongs to and which channel they are
		for(size_t idxComponent=0; idxComponent<nrComponents; ++idxComponent) {
			// loop over channels of component and see if wave is there
			const size_t nrChannels = _comp[idxComponent]->getNrChannels();
			for(size_t idxChannel=0; idxChannel<nrChannels; ++idxChannel) {
				if(_comp[idxComponent]->getChannelWaveName(idxChannel) == _waveNames[idxWave]) {
					_compChannel[idxWave].push_back(pair<size_t, size_t>(idxComponent,idxChannel));

					if(anchorWaveName == _waveNames[idxWave] && anchorComponentName == _comp[idxComponent]->getName()) {
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
		for(size_t idxComponents=0; idxComponents<_compChannel[idxWave].size(); idxComponents++) {
			const size_t idxComponent = _compChannel[idxWave][idxComponents].first;
			output << "        component " << idxComponents << ": " << idxComponent << " '" << _comp[idxComponent]->getName() << "'";
			if(_idxAnchorWave == idxWave && _idxAnchorComponent == idxComponent) {
				output << " (anchor)";
			}
			output << endl;
		}
	}
	for(size_t idxComponent=0; idxComponent<nrComponents; ++idxComponent) {
		output << "    component '" << _comp[idxComponent]->getName() << "' (index " << idxComponent << ") used in" << endl;

		const size_t nrChannels = _comp[idxComponent]->getNrChannels();
		for(size_t idxChannel=0; idxChannel<nrChannels; ++idxChannel) {
			output << "        channel " << idxChannel << ": " << _comp[idxComponent]->getChannelWaveName(idxChannel);
			if(_idxAnchorComponent == idxComponent && _idxAnchorChannel == idxChannel) {
				output << " (anchor)";
			}
			output << endl;
		}
	}
	printInfo << _waveNames.size() << " waves and " << n() << " components in fit model:" << endl
	          << output.str();

	return true;
}


bool
massDepFitModel::initFsmd(const vector<double>& massBinCenters)
{
	const size_t nrMassBins = massBinCenters.size();

	_fsmdValues.resize(nrMassBins);

	for(size_t idxMassBin=0; idxMassBin<nrMassBins; ++idxMassBin) {
		_fsmdValues[idxMassBin] = calcFsmd(massBinCenters[idxMassBin]);
	}

	return true;
}


double 
massDepFitModel::getFreeFsmdPar(unsigned int i) const
{
  if(i<_fsmdFreeParameters.size())
    return _fsmdFunction->GetParameter(_fsmdFreeParameters[i]);
  else return 0;
}


void 
massDepFitModel::getFreeFsmdLimits(unsigned int i, double& lower, double& upper) const
{
  if(i<_fsmdFreeParameters.size()){
    _fsmdFunction->GetParLimits(_fsmdFreeParameters[i],lower,upper);
  }
}


void
massDepFitModel::setPar(const double* par)
{
  size_t parcount=0;
  // components
  for(unsigned int i=0;i<n();++i){
    _comp[i]->setParameters(&par[parcount]);
    parcount+=_comp[i]->getNrParameters();
    _comp[i]->setCouplings(&par[parcount]);
    parcount+=_comp[i]->getNrChannels()*2; // RE and Im for each channel
  } // end loop over components
  // phase space
  unsigned int nfreepar=_fsmdFreeParameters.size();
  for(unsigned int ipar=0;ipar<nfreepar;++ipar){
    _fsmdFunction->SetParameter(_fsmdFreeParameters[ipar],par[parcount]);
    ++parcount;
  }
}


void 
massDepFitModel::getPar(double* par)
{
  size_t parcount=0;
  // components
  for(unsigned int i=0;i<n();++i){
    _comp[i]->getParameters(&par[parcount]);
    parcount+=_comp[i]->getNrParameters();
    _comp[i]->getCouplings(&par[parcount]);
    parcount+=_comp[i]->getNrChannels()*2; // RE and Im for each channel
  }
 // phase space
  unsigned int nfreepar=_fsmdFreeParameters.size();
  for(unsigned int ipar=0;ipar<nfreepar;++ipar){
    par[parcount]=_fsmdFunction->GetParameter(_fsmdFreeParameters[ipar]);
    ++parcount;
  }
}


double 
massDepFitModel::calcFsmd(const double mass,
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
massDepFitModel::productionAmplitude(const size_t idxWave,
                                     const double mass,
                                     const size_t idxMass) const
{
	// loop over all components and pick up those that contribute to this channels
	complex<double> prodAmp(0., 0.);

	// get entry from mapping
	const vector<pair<size_t, size_t> >& components = _compChannel[idxWave];
	size_t nrComponents = components.size();

	for(unsigned int idxComponents=0; idxComponents<nrComponents; ++idxComponents) {
		size_t idxComponent = components[idxComponents].first;
		size_t idxChannel = components[idxComponents].second;
		prodAmp += _comp[idxComponent]->val(mass) * _comp[idxComponent]->getChannel(idxChannel).CsqrtPS(mass);
	}

	prodAmp *= calcFsmd(mass, idxMass);

	return prodAmp;
}


double 
massDepFitModel::intensity(const size_t idxWave,
                           const double mass,
                           const size_t idxMass) const
{
	const complex<double> prodAmp = productionAmplitude(idxWave, mass, idxMass);

	return norm(prodAmp);
}


double 
massDepFitModel::phaseAbsolute(const size_t idxWave,
                               const double mass,
                               const size_t idxMass) const
{
	const complex<double> prodAmp = productionAmplitude(idxWave, mass, idxMass);

	return arg(prodAmp);
}


complex<double>
massDepFitModel::spinDensityMatrix(const size_t idxWave,
                                   const size_t jdxWave,
                                   const double mass,
                                   const size_t idxMass) const
{
	const complex<double> prodAmpI = productionAmplitude(idxWave, mass, idxMass);
	const complex<double> prodAmpJ = productionAmplitude(jdxWave, mass, idxMass);

	return prodAmpI * conj(prodAmpJ);
}


double 
massDepFitModel::phase(const size_t idxWave,
                       const size_t jdxWave,
                       const double mass,
                       const size_t idxMass) const
{
	return arg(spinDensityMatrix(idxWave, jdxWave, mass, idxMass));
}


ostream&
massDepFitModel::print(ostream& out) const
{
	for(unsigned int i=0;i<_comp.size();++i){
		const massDepFitComponent& c = *_comp[i];
		c.print(out);
	}
	return out;
}
