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

#include "pwacomponent.h"
#include "reportingUtils.hpp"


using namespace std;
using namespace rpwa;


bool
massDepFitModel::init(const vector<string>& waveNames,
                      const vector<double>& massBinCenters)
{
	_waveNames = waveNames;

	if(not initMapping()) {
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
massDepFitModel::add(pwacomponent* comp)
{
    _comp.push_back(comp);
    _numpar+=comp->numPar();
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
massDepFitModel::initMapping()
{
	// check that all waves used in a decay channel have been defined
	for(unsigned int i=0; i<n(); ++i) {
		for(map<string, pwachannel>::const_iterator itChan=_comp[i]->channels().begin(); itChan !=_comp[i]->channels().end(); ++itChan) {
			if(find(_waveNames.begin(), _waveNames.end(), itChan->first) == _waveNames.end()) {
				printErr << "wave '" << itChan->first << "' not known in decay of '" << _comp[i]->name() << "'." << endl;
				return false;
			}
		}
	}

	// check that all defined waves are also used
	for(vector<string>::const_iterator itWave=_waveNames.begin(); itWave!=_waveNames.end(); ++itWave) {
		bool found(false);
		for(unsigned int i=0; i<n(); ++i) {
			for(map<string, pwachannel>::const_iterator itChan=_comp[i]->channels().begin(); itChan !=_comp[i]->channels().end(); ++itChan) {
				if(itChan->first == *itWave) {
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

	for(vector<string>::const_iterator itWave=_waveNames.begin(); itWave!=_waveNames.end(); ++itWave) {
		// check which components this waves belongs to and which channel they are
		_compChannel.push_back(getCompChannel(*itWave));
	}

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
  unsigned int parcount=0;
  // components
  for(unsigned int i=0;i<n();++i){
    _comp[i]->setPar(par[parcount],par[parcount+1]);
    parcount+=2;
    _comp[i]->setCouplings(&par[parcount]);
    parcount+=_comp[i]->numChannels()*2; // RE and Im for each channel
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
  unsigned int parcount=0;
  // components
  for(unsigned int i=0;i<n();++i){
    par[parcount]=_comp[i]->m0();
    par[parcount+1]=_comp[i]->gamma();
    parcount+=2;
    _comp[i]->getCouplings(&par[parcount]);
    parcount+=_comp[i]->numChannels()*2; // RE and Im for each channel
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
	const vector<pair<unsigned int, unsigned int> >& components = _compChannel[idxWave];
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


vector<pair<unsigned int,unsigned int> >
massDepFitModel::getCompChannel(const string& wave) const
{
  cerr << "Channel-mapping for wave " << wave << endl;
  vector<pair<unsigned int,unsigned int> > result;
  for(unsigned int ic=0;ic<n();++ic){
    // loop over channels of component and see if wave is there
    unsigned int nch=_comp[ic]->numChannels();
    for(unsigned int ich=0;ich<nch;++ich){
      if(_comp[ic]->getChannelName(ich)==wave){
	result.push_back(pair<unsigned int,unsigned int>(ic,ich));
	cerr << "     comp("<<ic<<"): " << _comp[ic]->name() << "  ch:"<<ich<<endl;
      }
    } // end loop over channels
  } // end loop over components
  return result;
}


ostream&
massDepFitModel::print(ostream& out) const
{
	for(unsigned int i=0;i<_comp.size();++i){
		const pwacomponent& c =*_comp[i];
		out << c << endl;
	}
	return out;
}
