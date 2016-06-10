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

		bool addAmplitude(isobarAmplitudePtr             amplitude);
		bool addIntegral (const rpwa::ampIntegralMatrix& integralMatrix);

		// initialize amplitudes starting from X decay
		bool initAmplitudes(const std::vector<std::string>& decayKinParticleNames);

		bool initAmplitudes(const std::vector<std::string>& prodKinParticleNames,
		                    const std::vector<std::string>& decayKinParticleNames,
		                    const bool                      fromXDecay = false);

		// get intensity if amplitudes have been initialized to start from X decay
		double getIntensity(const std::vector<TVector3>& decayKinMomenta) const;

		double getIntensity(const std::vector<TVector3>& prodKinMomenta,
		                    const std::vector<TVector3>& decayKinMomenta) const;

		std::ostream& print(std::ostream& out = std::cout) const;
		friend std::ostream& operator << (std::ostream&         out,
		                                  const modelIntensity& model) { return model.print(out); }

	private:

		std::vector<std::complex<double> > getAmplitudes(const std::vector<TVector3>& prodKinMomenta,
		                                                 const std::vector<TVector3>& decayKinMomenta) const;

		fitResultPtr                       _fitResult;

		bool                               _amplitudesInitialized;
		std::vector<isobarAmplitudePtr>    _amplitudes;
		std::vector<int>                   _refls;
		std::set<int>                      _allRefls;

		bool                               _integralsLoaded;
		std::vector<double>                _integrals;

		bool                               _amplitudesFromXDecay;

	};


} // namespace rpwa


#endif // MODELINTENSITY_H
