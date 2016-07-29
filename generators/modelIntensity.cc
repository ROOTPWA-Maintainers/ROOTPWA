#include"modelIntensity.h"
#include"ampIntegralMatrix.h"
#include"waveDescription.h"


rpwa::modelIntensity::modelIntensity(fitResultPtr fitResult)
	: _fitResult(fitResult),
	  _decayAmplitudesInitialized(false),
	  _decayAmplitudes(fitResult->nmbWaves()),
	  _refls(fitResult->nmbWaves(), 0),
	  _phaseSpaceIntegralsLoaded(false),
	  _phaseSpaceIntegrals(fitResult->nmbWaves()),
	  _decayAmplitudesFromXDecay(false)
{
	// use the phase-space integrals from fit result if available
	if (fitResult->phaseSpaceIntegralVector().size() == fitResult->nmbWaves()) {
		_phaseSpaceIntegralsLoaded = true;
		_phaseSpaceIntegrals       = fitResult->phaseSpaceIntegralVector();
	}

	_waveIndicesWithoutFlat = fitResult->waveIndicesMatchingPattern("^(?!flat$).*$");
}


bool
rpwa::modelIntensity::addDecayAmplitude(rpwa::isobarAmplitudePtr decayAmplitude)
{
	const std::string waveName  = waveDescription::waveNameFromTopology(*(decayAmplitude->decayTopology()));
	const int         waveIndex = _fitResult->waveIndex(waveName);
	if (waveIndex == -1) {
		printWarn << "cannot find wave '" << waveName << "' in fit result provided." << std::endl;
		return false;
	}
	if (_decayAmplitudes[waveIndex]) {
		printWarn << "decay amplitude for wave '" << waveName << "' is added a second time." << std::endl;
		return false;
	}

	_decayAmplitudesInitialized = false;
	_decayAmplitudes[waveIndex] = decayAmplitude;
	_refls          [waveIndex] = decayAmplitude->decayTopology()->XIsobarDecayVertex()->parent()->reflectivity();
	_allRefls.insert(decayAmplitude->decayTopology()->XIsobarDecayVertex()->parent()->reflectivity());
	return true;
}


bool
rpwa::modelIntensity::loadPhaseSpaceIntegral(const rpwa::ampIntegralMatrix& integralMatrix)
{
	_phaseSpaceIntegralsLoaded = false;
	_phaseSpaceIntegrals.resize(_fitResult->nmbWaves());

	for (size_t wave = 0; wave < _fitResult->nmbWaves(); ++wave) {
		const std::string& waveName = _fitResult->waveName(wave);
		if (waveName == "flat") {
			_phaseSpaceIntegrals[wave] = 1.;
		} else {
			if (not integralMatrix.containsWave(waveName)) {
				printWarn << "wave '" << waveName << "' not in the integral matrix." << std::endl;
				return false;
			}
			const std::complex<double> integral = integralMatrix.element(waveName, waveName);
			_phaseSpaceIntegrals[wave] = std::sqrt(integral.real());
		}
	}

	_phaseSpaceIntegralsLoaded = true;
	return true;
}


bool
rpwa::modelIntensity::initDecayAmplitudes(const std::vector<std::string>& decayKinParticleNames)
{
	std::vector<std::string> prodKinParticleNames(1);
	return initDecayAmplitudes(prodKinParticleNames, decayKinParticleNames, true);
}


bool
rpwa::modelIntensity::initDecayAmplitudes(const std::vector<std::string>& prodKinParticleNames,
                                          const std::vector<std::string>& decayKinParticleNames)
{
	return initDecayAmplitudes(prodKinParticleNames, decayKinParticleNames, false);
}


// if 'initDecayAmplitudes' is called with a 'const std::vector' as an
// argument we have to assume that the 'prodKinParticleNames' is already
// set up, and should not try to set it again. if it is called with a
// non-const vector, which is the case if it is called from the
// 'initDecayAmplitudes' used for decays starting from the X, then the
// name of the production particle should be set.
namespace {

	template<typename T>
	void setProdKinParticleNames(T& prodKinParticleNames, const std::string& name);


	template<>
	void setProdKinParticleNames<std::vector<std::string> >(std::vector<std::string>& prodKinParticleNames, const std::string& name)
	{
		prodKinParticleNames.resize(1);
		prodKinParticleNames[0] = name;
	}


	template<>
	void setProdKinParticleNames<const std::vector<std::string> >(const std::vector<std::string>&, const std::string&)
	{
		// do nothing
	}

}


template<typename T>
bool
rpwa::modelIntensity::initDecayAmplitudes(T&                              prodKinParticleNames,
                                          const std::vector<std::string>& decayKinParticleNames,
                                          const bool                      fromXDecay)
{
	_decayAmplitudesInitialized = false;
	_decayAmplitudesFromXDecay  = fromXDecay;

	//initialize the decay amplitudes
	for (size_t wave = 0; wave < _fitResult->nmbWaves(); ++wave) {
		const std::string& waveName = _fitResult->waveName(wave);
		// no decay amplitude for 'flat' wave
		if (waveName == "flat") {
			if (_decayAmplitudes[wave]) {
				printWarn << "decay amplitude for flat wave was added." << std::endl;
				return false;
			}
			continue;
		}
		if (not _decayAmplitudes[wave]) {
			printWarn << "decay amplitude for wave '" << waveName << "' was not added." << std::endl;
			return false;
		}

		// the name of the production particle is not yet known, it has
		// to be obtained from the X decay vertex
		if (prodKinParticleNames.size() == 0 or prodKinParticleNames[0] == "") {
			const isobarDecayTopologyPtr decay  = _decayAmplitudes[wave]->decayTopology();
			const isobarDecayVertexPtr   vertex = decay->XIsobarDecayVertex();
			setProdKinParticleNames(prodKinParticleNames, vertex->parent()->name());
		}

		// start decay amplitude from X decay
		if (_decayAmplitudesFromXDecay) {
			const isobarDecayTopologyPtr decay  = _decayAmplitudes[wave]->decayTopology();
			const isobarDecayVertexPtr   vertex = decay->XIsobarDecayVertex();
			_decayAmplitudes[wave]->setDecayTopology(rpwa::createIsobarDecayTopology(decay->subDecayConsistent(vertex)));
		}

		_decayAmplitudes[wave]->init();
		if (not _decayAmplitudes[wave]->decayTopology()->initKinematicsData(prodKinParticleNames, decayKinParticleNames)) {
			printWarn << "could not initialize kinematics data for decay amplitude of wave '" << waveName << "'." << std::endl;
			return false;
		}
	}

	_decayAmplitudesInitialized = true;
	return true;
}


double
rpwa::modelIntensity::getIntensity(const std::vector<unsigned int>& waveIndices,
                                   const std::vector<TVector3>&     prodKinMomenta,
                                   const std::vector<TVector3>&     decayKinMomenta) const
{
	const std::vector<std::complex<double> > decayAmplitudes = getDecayAmplitudes(prodKinMomenta, decayKinMomenta);

	double intensity = 0;
	for (std::set<int>::const_iterator it=_allRefls.begin(); it!=_allRefls.end(); ++it) {
		std::complex<double> amp = 0;
		for (size_t i=0; i<waveIndices.size(); ++i) {
			const unsigned int waveIndex = waveIndices[i];
			if (_refls[waveIndex] == *it) {
				amp += _fitResult->prodAmp(waveIndex) * decayAmplitudes[waveIndex];
			}
		}
		intensity += std::norm(amp);
	}

	return intensity;
}


std::vector<std::complex<double> >
rpwa::modelIntensity::getDecayAmplitudes(const std::vector<TVector3>& prodKinMomenta,
                                         const std::vector<TVector3>& decayKinMomenta) const
{
	if (not _decayAmplitudesInitialized) {
		printErr << "decay amplitudes not initialized, cannot evaluate model. Aborting..." << std::endl;
		throw;
	}
	if (not _phaseSpaceIntegralsLoaded) {
		printErr << "integrals not loaded, cannot evaluate model. Aborting..." << std::endl;
		throw;
	}

	std::vector<std::complex<double> > decayAmplitudes(_decayAmplitudes.size());
	for (size_t wave = 0; wave < _fitResult->nmbWaves(); ++wave) {
		// 'flat' wave
		if (not _decayAmplitudes[wave]) {
			decayAmplitudes[wave] = 1. / _phaseSpaceIntegrals[wave];
			continue;
		}

		if (not _decayAmplitudes[wave]->decayTopology()->readKinematicsData(prodKinMomenta, decayKinMomenta)) {
			printErr << "could not read kinematics data for wave '" << _fitResult->waveName(wave) << "'. Aborting..." << std::endl;
			throw;
		}
		decayAmplitudes[wave] = _decayAmplitudes[wave]->amplitude() / _phaseSpaceIntegrals[wave];
	}
	return decayAmplitudes;
}


std::ostream&
rpwa::modelIntensity::print(std::ostream& out) const
{
	out << "status of model intensity object:" << std::endl
	    << "    number of waves ................. " << _fitResult->nmbWaves() << std::endl
	    << "    wave names" << std::endl;
	for (size_t wave = 0; wave < _fitResult->nmbWaves(); ++wave) {
		out << "                                        " << _fitResult->waveName(wave) << std::endl;
	}
	out << "    decay amplitudes initialized .... " << yesNo(_decayAmplitudesInitialized) << std::endl
	    << "    integrals loaded ................ " << yesNo(_phaseSpaceIntegralsLoaded) << std::endl
	    << "    decay amplitudes from X decay ... " << yesNo(_decayAmplitudesFromXDecay) << std::endl;
	return out;
}
