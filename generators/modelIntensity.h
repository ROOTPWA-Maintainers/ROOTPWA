#ifndef MODELINTENSITY_H
#define MODELINTENSITY_H


#include<complex>
#include<vector>

#include<boost/shared_ptr.hpp>

#include<TVector3.h>

#include"isobarAmplitude.h"
#include"fitResult.h"


namespace rpwa {


	class ampIntegralMatrix;

	class modelIntensity;
	typedef boost::shared_ptr<modelIntensity> modelIntensityPtr;


	class modelIntensity {

	public:

		modelIntensity(fitResultPtr fitResult);

		bool addDecayAmplitude     (isobarAmplitudePtr             decayAmplitude);
		bool loadPhaseSpaceIntegral(const rpwa::ampIntegralMatrix& integralMatrix);

		// initialize decay amplitudes starting from X decay
		bool initDecayAmplitudes(const std::vector<std::string>& decayKinParticleNames);

		bool initDecayAmplitudes(const std::vector<std::string>& prodKinParticleNames,
		                         const std::vector<std::string>& decayKinParticleNames);

		// get intensity of all waves except flat wave

		// get intensity if decay amplitudes have been initialized to start from X decay
		double getIntensity(const std::vector<TVector3>& decayKinMomenta) const;

		double getIntensity(const std::vector<TVector3>& prodKinMomenta,
		                    const std::vector<TVector3>& decayKinMomenta) const;

		// get intensity for set of waves

		// get intensity if decay amplitudes have been initialized to start from X decay
		double getIntensity(const std::vector<unsigned int>& waveIndices,
		                    const std::vector<TVector3>&     decayKinMomenta) const;

		double getIntensity(const std::vector<unsigned int>& waveIndices,
		                    const std::vector<TVector3>&     prodKinMomenta,
		                    const std::vector<TVector3>&     decayKinMomenta) const;

		std::ostream& print(std::ostream& out = std::cout) const;
		friend std::ostream& operator << (std::ostream&         out,
		                                  const modelIntensity& model) { return model.print(out); }

	private:

		template<typename T>
		bool initDecayAmplitudes(T&                              prodKinParticleNames,
		                         const std::vector<std::string>& decayKinParticleNames,
		                         const bool                      fromXDecay);

		std::vector<std::complex<double> > getDecayAmplitudes(const std::vector<TVector3>& prodKinMomenta,
		                                                      const std::vector<TVector3>& decayKinMomenta) const;

		fitResultPtr                       _fitResult;
		std::vector<unsigned int>          _waveIndicesWithoutFlat;

		bool                               _decayAmplitudesInitialized;
		std::vector<isobarAmplitudePtr>    _decayAmplitudes;
		std::vector<int>                   _refls;
		std::set<int>                      _allRefls;

		bool                               _phaseSpaceIntegralsLoaded;
		std::vector<double>                _phaseSpaceIntegrals;

		bool                               _decayAmplitudesFromXDecay;

	};


	inline
	double
	modelIntensity::getIntensity(const std::vector<TVector3>& decayKinMomenta) const
	{
		return getIntensity(_waveIndicesWithoutFlat, decayKinMomenta);
	}


	inline
	double
	modelIntensity::getIntensity(const std::vector<TVector3>& prodKinMomenta,
	                             const std::vector<TVector3>& decayKinMomenta) const
	{
		return getIntensity(_waveIndicesWithoutFlat, prodKinMomenta, decayKinMomenta);
	}


	inline
	double
	modelIntensity::getIntensity(const std::vector<unsigned int>& waveIndices,
	                             const std::vector<TVector3>&     decayKinMomenta) const
	{
		if (not _decayAmplitudesFromXDecay) {
			printErr << "decay amplitudes are not starting from X decay, but no production kinematics provided. Aborting..." << std::endl;
			throw;
		}

		return getIntensity(waveIndices, std::vector<TVector3>(1), decayKinMomenta);
	}


} // namespace rpwa


#endif // MODELINTENSITY_H
