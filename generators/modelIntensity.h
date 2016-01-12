#ifndef MODELINTENSITY_H
#define MODELINTENSITY_H


#include<complex>
#include<vector>

#include<boost/shared_ptr.hpp>

#include<TVector3.h>

#include"isobarAmplitude.h"


namespace rpwa {


	class ampIntegralMatrix;
	class modelIntensity;
	typedef boost::shared_ptr<modelIntensity> modelIntensityPtr;


	class modelIntensity {

	public:

		modelIntensity();

		double getIntensity(const std::vector<TVector3>& prodKinMomenta,
		                    const std::vector<TVector3>& decayKinMomenta) const;
		double getIntensity(const std::vector<TVector3>& decayKinMomenta) const { return getIntensity(_prodKinMomenta, decayKinMomenta); }

		std::vector<std::complex<double> > getAmplitudes(const std::vector<TVector3>& prodKinMomenta,
		                                                 const std::vector<TVector3>& decayKinMomenta) const;
		std::vector<std::complex<double> > getAmplitudes(const std::vector<TVector3>& decayKinMomenta) const { return getAmplitudes(_prodKinMomenta, decayKinMomenta); }

		bool addAmplitude(const std::complex<double> transitionAmplitude,
		                  isobarAmplitudePtr amplitude,
		                  const int reflectivity = 1);
		bool loadIntegrals(const rpwa::ampIntegralMatrix& integrals);
		bool setProdKinMomenta(const std::vector<TVector3>& prodKinMomenta);
		bool initAmplitudes(const bool fromXdecay = true);

		bool setTarget(const std::string& target) { _target = target; return true; }

		double mBeam() const { return _mBeam; }
		double mTarget() const { return _mTarget; }

		size_t nFinalState() const { return _finalStateParticles.size(); }
		const std::string& finalStateParticleName(const size_t particleIndex) const;
		const std::string& tagetParticleName() const { return _target; }

		const std::vector<double>& finalStateMasses() const { return _finalStateMasses; }
		void print() const;

		bool setMbeam(const double mass) { _mBeam = mass; return true; }
		bool setMtarget(const double mass) { _mTarget = mass; return true; }
		bool setFinalStateMasses(const std::vector<double>& masses);
		bool setParticles(const std::vector<std::string>& initialState,
		                  const std::vector<std::string>& finalState);

		const std::vector<std::string>& initialStateParticles() const { return _initialStateParticles; }
		const std::vector<std::string>& finalStateParticles() const { return _finalStateParticles; }

	private:

		bool                               _allIntsLoaded;
		bool                               _amplitudesInitialized;
		double                             _mBeam;
		double                             _mTarget;
		std::string                        _target;
		std::vector<std::string>           _waveNames;
		std::vector<std::string>           _initialStateParticles;
		std::vector<std::string>           _finalStateParticles;
		std::vector<double>                _finalStateMasses;

		std::vector<TVector3>              _prodKinMomenta;
		std::vector<int>                   _reflectivities;
		std::vector<double>                _integrals;
		std::vector<std::complex<double> > _transitionAmplitudes;
		std::vector<isobarAmplitudePtr>    _amplitudes;

	};


} // namespace rpwa


#endif // MODELINTENSITY_H
