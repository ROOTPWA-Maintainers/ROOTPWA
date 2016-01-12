#include"modelIntensity.h"
#include"ampIntegralMatrix.h"
#include"waveDescription.h"


rpwa::modelIntensity::modelIntensity()
	: _allIntsLoaded(true),
	  _amplitudesInitialized(false),
	  _mBeam(0.),
	  _mTarget(0.),
	  _target("p+"),
	  _waveNames(),
	  _initialStateParticles(),
	  _finalStateParticles(),
	  _finalStateMasses(),
	  _prodKinMomenta(std::vector<TVector3>(1, TVector3())),
	  _reflectivities(),
	  _integrals(),
	  _transitionAmplitudes(),
	  _amplitudes()
{
}


double
rpwa::modelIntensity::getIntensity(const std::vector<TVector3>& prodKinMomenta,
                                   const std::vector<TVector3>& decayKinMomenta) const
{
	std::vector<std::complex<double> > amplitudes = getAmplitudes(prodKinMomenta, decayKinMomenta);
	double retVal = 0.;
	for (size_t incoherentSubAmplitude = 0; incoherentSubAmplitude < amplitudes.size(); ++incoherentSubAmplitude) {
		retVal += std::norm(amplitudes[incoherentSubAmplitude]);
	}
	if (std::isnan(retVal)) {
		printErr << "NaN intensity encountered" << std::endl;
		throw;
	}
	return retVal;
}


std::vector<std::complex<double> >
rpwa::modelIntensity::getAmplitudes(const std::vector<TVector3>& prodKinMomenta,
                                    const std::vector<TVector3>& decayKinMomenta) const
{
	std::vector<std::complex<double> > retVal(2, std::complex<double>(0., 0.)); // at the moment use only positive and negative reflectivity as incoherernt amplitudes (means rank = 1)
	if (!_allIntsLoaded) {
		printErr << "integrals not loaded, can't evaluate model" << std::endl;
		throw;
	}
	if (!_amplitudesInitialized) {
		printErr << "amplitudes not initialized, can't evaluate model" << std::endl;
		throw;
	}
	for (size_t amp = 0; amp < _amplitudes.size(); ++amp) {
		if (!_amplitudes[amp]->decayTopology()->readKinematicsData(prodKinMomenta, decayKinMomenta)) {
			printErr << "could not readKinematicsData() for wave: " << _waveNames[amp] << std::endl;
			throw;
		}
		const std::complex<double> amplitude = (_amplitudes[amp])->amplitude() * _transitionAmplitudes[amp] / _integrals[amp];
		if (_reflectivities[amp] > 0) {
			retVal[0] += amplitude;
		} else {
			retVal[1] += amplitude;
		}
	}
	return retVal;
}


bool
rpwa::modelIntensity::setParticles(const std::vector<std::string>& initialState,
                                   const std::vector<std::string>& finalState)
{
	if (_finalStateMasses.size() != 0 && _finalStateMasses.size() != finalState.size()) {
		printErr << "number of final state particles does not match number of final state masses" << std::endl;
		return false;
	}
	_initialStateParticles = initialState;
	_finalStateParticles   = finalState;
	return true;
}


bool
rpwa::modelIntensity::setFinalStateMasses(const std::vector<double>& masses)
{
	_finalStateMasses = masses;
	if (_finalStateParticles.size() != 0 && _finalStateParticles.size() != masses.size()) {
		printErr << "number of masses given does not match" << std::endl;
		return false;
	}
	return true;
}


bool
rpwa::modelIntensity::addAmplitude(const std::complex<double> transitionAmplitude,
                                   rpwa::isobarAmplitudePtr   amplitude,
                                   const int                  reflectivity)
{
	_amplitudesInitialized = false;
	if (reflectivity != 1 && reflectivity != -1) {
		printErr << "reflectivity != +-1 (" << reflectivity << "). Aborting..." << std::endl;
		return false;
	}
	if (_initialStateParticles.size() == 0 && _finalStateParticles.size() == 0) {
		const std::vector<particlePtr>& finalStateparticles = amplitude->decayTopology()->fsParticles();
		std::vector<double> finalStateMasses;
		std::vector<std::string> finalStateParticles;
		if (finalStateparticles.size() == 0) {
			printErr << "not final state particles found" << std::endl;
			return false;
		}
		for (size_t part = 0; part < finalStateparticles.size(); ++part) {
			finalStateMasses.push_back(finalStateparticles[part]->mass());
			finalStateParticles.push_back(finalStateparticles[part]->name());
		}
		const std::vector<particlePtr>& inParticles = amplitude->decayTopology()->productionVertex()->inParticles();
		double beamMass = 0.;
		std::vector<std::string> initialStateNames;
		if (inParticles.size() == 0) {
			printErr << "no initial state particles found" << std::endl;
			return false;
		}
		for (size_t part = 0; part < inParticles.size(); ++part) {
			const std::string partName = inParticles[part]->name();
			if (partName == _target) {
				if (!setMtarget(inParticles[part]->mass())) {
					printErr << "could not set target mass" << std::endl;
				}
				continue;
			}
			beamMass = inParticles[part]->mass();
			initialStateNames.push_back(partName);
		}
		if (!setParticles(initialStateNames, finalStateParticles)) {
			printErr << "could not set particles" << std::endl;
			return false;
		}
		if (!setFinalStateMasses(finalStateMasses)) {
			printErr << "could not set final state masses" << std::endl;
			return false;
		}
		if (!setMbeam(beamMass)) {
			printErr << "could not set mBeam" << std::endl;
			return false;
		}
	}
	_allIntsLoaded = false;
	_waveNames.push_back(waveDescription::waveNameFromTopology(*(amplitude->decayTopology())));
	_reflectivities.push_back(reflectivity);
	_transitionAmplitudes.push_back(transitionAmplitude);
	_amplitudes.push_back(amplitude);
	return true;
}


bool
rpwa::modelIntensity::loadIntegrals(const rpwa::ampIntegralMatrix& integralMatrix)
{
	std::vector<double> integrals;
	const unsigned long nmbEvents = integralMatrix.nmbEvents();
	for (size_t wn = 0; wn < _waveNames.size(); ++wn) {
		if (!integralMatrix.containsWave(_waveNames[wn])) {
			std::cout << "wave '" << wn << "' not in the integral matrix. Aborting..." << std::endl;
			return false;
		}
		const std::complex<double> integral = integralMatrix.element(_waveNames[wn], _waveNames[wn]);
		integrals.push_back(std::sqrt(integral.real()*nmbEvents));
	}
	if (integrals.size() != _amplitudes.size()) {
		printErr << "integral.size() != _amplitudes.size(). Aborting..." << std::endl;
		return false;
	}
	_integrals = integrals;
	_allIntsLoaded = true;
	return true;
}


bool
rpwa::modelIntensity::setProdKinMomenta(const std::vector<TVector3>& prodKinMomenta)
{
	if (_initialStateParticles.size() != 0 && _initialStateParticles.size() != prodKinMomenta.size()) {
		printErr << "prodKinMomenta.size() = " << prodKinMomenta.size() << " != _initialStateParticles.size() = " << _initialStateParticles.size() << ". Aborting..." << std::endl;
		return false;
	}
	_prodKinMomenta = prodKinMomenta;
	return true;
}


bool
rpwa::modelIntensity::initAmplitudes(const bool fromXdecay)
{
	std::vector<std::string> initialStateParticles = _initialStateParticles;
	for (size_t amp = 0; amp < _amplitudes.size(); ++amp) {
		if (fromXdecay) {
			_amplitudes[amp]->setDecayTopology(rpwa::createIsobarDecayTopology(_amplitudes[amp]->decayTopology()->subDecayConsistent(_amplitudes[amp]->decayTopology()->XIsobarDecayVertex())));
			initialStateParticles = std::vector<std::string>(1, "X-");
		}
		_amplitudes[amp]->init();
		if (!_amplitudes[amp]->decayTopology()->initKinematicsData(initialStateParticles, _finalStateParticles)) {
			printErr << "could not initKinematicsData()" << std::endl;
			return false;
		}
	}
	_amplitudesInitialized = true;
	return true;
}


void
rpwa::modelIntensity::print() const
{
	std::cout << "status of modelIntensity" << std::endl;
	std::cout << "  waveNames:" << std::endl;
	for (size_t i = 0; i < _waveNames.size(); ++i) {
		std::cout << "   - " << _waveNames[i] << std::endl;
	}
	std::cout << "  initial particles" << std::endl;
	for (size_t i = 0; i < _initialStateParticles.size(); ++i) {
		std::cout << "   - " << _initialStateParticles[i] << std::endl;
	}
	std::cout << "  final particles" << std::endl;
	for (size_t i = 0; i < _finalStateParticles.size(); ++i) {
		std::cout << "   - " << _finalStateParticles[i] << std::endl;
	}
}


const std::string&
rpwa::modelIntensity::finalStateParticleName(const size_t particleIndex) const
{
	if (particleIndex >= nFinalState()) {
		printErr << "particle index " << particleIndex << " exceeds number of final state particles " << nFinalState() << std::endl;
		throw;
	}
	return _finalStateParticles[particleIndex];
}
