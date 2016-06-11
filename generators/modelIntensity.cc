#include"modelIntensity.h"
#include"ampIntegralMatrix.h"
#include"waveDescription.h"


rpwa::modelIntensity::modelIntensity(fitResultPtr fitResult)
	: _fitResult(fitResult),
	  _amplitudesInitialized(false),
	  _amplitudes(fitResult->nmbWaves()),
	  _refls(fitResult->nmbWaves(), 0),
	  _integralsLoaded(false),
	  _integrals(fitResult->nmbWaves()),
	  _amplitudesFromXDecay(false)
{
	// use the phase-space integrals from fit result if available
	if (fitResult->phaseSpaceIntegralVector().size() == fitResult->nmbWaves()) {
		_integralsLoaded = true;
		_integrals       = fitResult->phaseSpaceIntegralVector();
	}
}


bool
rpwa::modelIntensity::addAmplitude(rpwa::isobarAmplitudePtr amplitude)
{
	const std::string waveName  = waveDescription::waveNameFromTopology(*(amplitude->decayTopology()));
	const int         waveIndex = _fitResult->waveIndex(waveName);
	if (waveIndex == -1) {
		printErr << "wave '" << waveName << "' not used in fit result provided. Aborting..." << std::endl;
		return false;
	}
	if (_amplitudes[waveIndex]) {
		printErr << "amplitude for wave '" << waveName << "' is added a second time. Aborting..." << std::endl;
		return false;
	}

	_amplitudesInitialized = false;
	_amplitudes[waveIndex] = amplitude;
	_refls     [waveIndex] = amplitude->decayTopology()->XIsobarDecayVertex()->parent()->reflectivity();
	_allRefls.insert(amplitude->decayTopology()->XIsobarDecayVertex()->parent()->reflectivity());
	return true;
}


bool
rpwa::modelIntensity::addIntegral(const rpwa::ampIntegralMatrix& integralMatrix)
{
	_integralsLoaded = false;
	_integrals.resize(_fitResult->nmbWaves());

	for (size_t wave = 0; wave < _fitResult->nmbWaves(); ++wave) {
		const std::string& waveName = _fitResult->waveName(wave);
		if (waveName == "flat") {
			_integrals[wave] = 1.;
		} else {
			if (not integralMatrix.containsWave(waveName)) {
				std::cout << "wave '" << waveName << "' not in the integral matrix. Aborting..." << std::endl;
				return false;
			}
			const std::complex<double> integral = integralMatrix.element(waveName, waveName);
			_integrals[wave] = std::sqrt(integral.real());
		}
	}

	_integralsLoaded = true;
	return true;
}


bool
rpwa::modelIntensity::initAmplitudes(const std::vector<std::string>& decayKinParticleNames)
{
	// get X name from first amplitude
	std::vector<std::string> prodKinParticleNames(1);
	for (size_t wave = 0; wave < _fitResult->nmbWaves(); ++wave) {
		const std::string& waveName = _fitResult->waveName(wave);
		// no amplitude for 'flat' wave
		if (waveName == "flat") {
			if (_amplitudes[wave]) {
				printErr << "amplitude for flat wave was added. Aborting..." << std::endl;
				return false;
			}
			continue;
		}
		if (not _amplitudes[wave]) {
			printErr << "amplitude for wave '" << waveName << "' was not added. Aborting..." << std::endl;
			return false;
		}

		const isobarDecayTopologyPtr decay  = _amplitudes[wave]->decayTopology();
		const isobarDecayVertexPtr   vertex = decay->XIsobarDecayVertex();
		prodKinParticleNames[0] = vertex->parent()->name();
		break;
	}

	return initAmplitudes(prodKinParticleNames, decayKinParticleNames, true);
}


bool
rpwa::modelIntensity::initAmplitudes(const std::vector<std::string>& prodKinParticleNames,
                                     const std::vector<std::string>& decayKinParticleNames,
                                     const bool                      fromXDecay)
{
	_amplitudesInitialized = false;
	_amplitudesFromXDecay  = fromXDecay;

	for (size_t wave = 0; wave < _fitResult->nmbWaves(); ++wave) {
		const std::string& waveName = _fitResult->waveName(wave);
		// no amplitude for 'flat' wave
		if (waveName == "flat") {
			if (_amplitudes[wave]) {
				printErr << "amplitude for flat wave was added. Aborting..." << std::endl;
				return false;
			}
			continue;
		}
		if (not _amplitudes[wave]) {
			printErr << "amplitude for wave '" << waveName << "' was not added. Aborting..." << std::endl;
			return false;
		}

		// start amplitude from X decay
		if (_amplitudesFromXDecay) {
			const isobarDecayTopologyPtr decay  = _amplitudes[wave]->decayTopology();
			const isobarDecayVertexPtr   vertex = decay->XIsobarDecayVertex();
			_amplitudes[wave]->setDecayTopology(rpwa::createIsobarDecayTopology(decay->subDecayConsistent(vertex)));
		}
		_amplitudes[wave]->init();
		if (not _amplitudes[wave]->decayTopology()->initKinematicsData(prodKinParticleNames, decayKinParticleNames)) {
			printErr << "could not initialize kinematics data for amplitude of wave '" << waveName << "'. Aborting..." << std::endl;
			return false;
		}
	}

	_amplitudesInitialized = true;
	return true;
}


double
rpwa::modelIntensity::getIntensity(const std::vector<TVector3>& decayKinMomenta) const
{
	if (not _amplitudesFromXDecay) {
		printErr << "amplitudes not starting from X decay, but not production kinematics provided. Aborting..." << std::endl;
		throw;
	}

	return getIntensity(std::vector<TVector3>(1), decayKinMomenta);
}


double
rpwa::modelIntensity::getIntensity(const std::vector<TVector3>& prodKinMomenta,
                                   const std::vector<TVector3>& decayKinMomenta) const
{
	const std::vector<std::complex<double> > amplitudes = getAmplitudes(prodKinMomenta, decayKinMomenta);

	std::vector<unsigned int> waveIndices; for (unsigned int i=0; i<_amplitudes.size()-1; ++i) waveIndices.push_back(i);
	double intensity = 0;
	for (std::set<int>::const_iterator it=_allRefls.begin(); it!=_allRefls.end(); ++it) {
		std::complex<double> amp = 0;
		for (size_t i=0; i<waveIndices.size(); ++i) {
			const unsigned int waveIndex = waveIndices[i];
			if (_refls[waveIndex] == *it) {
				amp += _fitResult->prodAmp(waveIndex) * amplitudes[waveIndex];
			}
		}
		intensity += std::norm(amp);
	}

	return intensity;
}


std::vector<std::complex<double> >
rpwa::modelIntensity::getAmplitudes(const std::vector<TVector3>& prodKinMomenta,
                                    const std::vector<TVector3>& decayKinMomenta) const
{
	if (not _amplitudesInitialized) {
		printErr << "amplitudes not initialized, cannot evaluate model. Aborting..." << std::endl;
		throw;
	}
	if (not _integralsLoaded) {
		printErr << "integrals not loaded, cannot evaluate model. Aborting..." << std::endl;
		throw;
	}

	std::vector<std::complex<double> > amplitudes(_amplitudes.size());
	for (size_t wave = 0; wave < _fitResult->nmbWaves(); ++wave) {
		// 'flat' wave
		if (not _amplitudes[wave]) {
			amplitudes[wave] = 1. / _integrals[wave];
			continue;
		}

		if (not _amplitudes[wave]->decayTopology()->readKinematicsData(prodKinMomenta, decayKinMomenta)) {
			printErr << "could not read kinematics data for wave '" << _fitResult->waveName(wave) << "'. Aborting..." << std::endl;
			throw;
		}
		amplitudes[wave] = _amplitudes[wave]->amplitude() / _integrals[wave];
	}
	return amplitudes;
}


std::ostream&
rpwa::modelIntensity::print(std::ostream& out) const
{
	out << "status of model intensity object:" << std::endl
	    << "    number of waves ........... " << _fitResult->nmbWaves() << std::endl
	    << "    wave names" << std::endl;
	for (size_t wave = 0; wave < _fitResult->nmbWaves(); ++wave) {
		out << "                                " << _fitResult->waveName(wave) << std::endl;
	}
	out << "    amplitudes initialized .... " << yesNo(_amplitudesInitialized) << std::endl
	    << "    integrals loaded .......... " << yesNo(_integralsLoaded) << std::endl
	    << "    amplitudes from X decay ... " << yesNo(_amplitudesFromXDecay) << std::endl;
	return out;
}
