#ifndef IMPORTANCESAMPLER_H
#define IMPORTANCESAMPLER_H

#include<string>
#include<vector>

#include<BAT/BCModel.h>

#include"eventFileWriter.h"
#include"generator.h"
#include"modelIntensity.h"

class TFile;


namespace rpwa {


	class nBodyPhaseSpaceKinematics;


	class importanceSampler;
	typedef boost::shared_ptr<importanceSampler> importanceSamplerPtr;


	class importanceSampler : public BCModel {

	public:

		importanceSampler(rpwa::modelIntensityPtr model,
		                  rpwa::beamAndVertexGeneratorPtr beamAndVertexGenerator,
		                  rpwa::massAndTPrimePickerPtr massAndTPrimePicker,
		                  const rpwa::Beam& beam,
		                  const rpwa::Target& target,
		                  const rpwa::FinalState& finalState);

		// overload BAT methods
		double LogAPrioriProbability(const std::vector<double>& parameters);
		double LogLikelihood(const std::vector<double>& parameters);
		void CalculateObservables(const std::vector<double>& parameters);

		bool initializeFileWriter(TFile*             outFile,
		                          const std::string& userString         = "importanceSampledEvents",
		                          const bool         storeMassAndTPrime = true,
		                          const std::string& massVariableName   = "mass",
		                          const std::string& tPrimeVariableName = "tPrime");
		bool finalizeFileWriter();

		void setPhaseSpaceOnly(const bool input = true) { _phaseSpaceOnly = input; }
		void setMassPrior     (TF1*       prior = 0   ) { _massPrior      = prior; }

	private:

		boost::tuples::tuple<bool, double, TLorentzVector, TLorentzVector> getProductionKinematics(const double mass) const;

		bool initializeNBodyPhaseSpace(rpwa::nBodyPhaseSpaceKinematics& nBodyPhaseSpace,
		                               const std::vector<double>&       parameters,
		                               const bool                       angles = true) const;

		rpwa::modelIntensityPtr         _model;

		rpwa::beamAndVertexGeneratorPtr _beamAndVertexGenerator;
		rpwa::massAndTPrimePickerPtr    _massAndTPrimePicker;
		const rpwa::Beam                _beam;
		const rpwa::Target              _target;
		const rpwa::FinalState          _finalState;

		const size_t                    _nPart;
		std::vector<double>             _masses;
		double                          _mSum;
		bool                            _phaseSpaceOnly;
		TF1*                            _massPrior;

		rpwa::eventFileWriter           _fileWriter;
		bool                            _storeMassAndTPrime;


		// function call statistics (copied from pwaLikelihood)
	public:
		unsigned int nCalls() const { return _funcCallInfo[LOGLIKELIHOOD].nmbCalls; }
		void resetFuncInfo();
		std::ostream& printFuncInfo(std::ostream& out = std::cout) const;
	private:
		enum functionCallEnum {
			LOGAPRIORIPROBABILITY = 0,
			LOGLIKELIHOOD         = 1,
			CALCULATEOBSERVABLES  = 2,
			NMB_FUNCTIONCALLENUM  = 3
		};
		struct functionCallInfo {
			unsigned int nmbCalls;   // number of times function was called
			double       totalTime;  // total execution time of function
		};
		mutable functionCallInfo        _funcCallInfo[NMB_FUNCTIONCALLENUM];


	};


} //namespace rpwa


#endif // IMPORTANCESAMPLER_H
