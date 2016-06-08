#ifndef IMPORTANCESAMPLER_H
#define IMPORTANCESAMPLER_H

#include<string>
#include<vector>

#include<TF1.h>

#include<BAT/BCModel.h>

#include"eventFileWriter.h"
#include"generator.h"
#include"modelIntensity.h"

class TFile;


namespace rpwa {


	class nBodyPhaseSpaceKinematics;


	class importanceSampler : public BCModel {

	public:

		importanceSampler(const double mMin,
		                  const double mMax,
		                  rpwa::modelIntensityPtr model);

		// overload BAT methods
		double LogAPrioriProbability(const std::vector<double>& parameters);
		double LogLikelihood(const std::vector<double>& parameters);
		void CalculateObservables(const std::vector<double>& parameters);

		void setPhaseSpaceOnly(const bool input = true) { _phaseSpaceOnly = input; }

		bool initializeFileWriter(TFile* outFile);
		bool finalizeFileWriter();

		bool initializeProductionGenerator(rpwa::Beam& beam,
		                                   rpwa::Target& target,
		                                   rpwa::beamAndVertexGeneratorPtr generator,
		                                   rpwa::massAndTPrimePickerPtr picker);

		std::pair<std::pair<bool,double>, std::vector<TLorentzVector> > getProductionKinematics(const double mass);

		size_t nCalls() const { return _nCalls; }

	private:

		bool initializeNBodyPhaseSpace(rpwa::nBodyPhaseSpaceKinematics& nBodyPhaseSpace,
		                               const std::vector<double>&       parameters,
		                               const bool                       angles = true) const;

		bool                    _phaseSpaceOnly;
		bool                    _productionGeneratorInitialized;
		size_t                  _nPart;
		static size_t           _nCalls;
		double                  _mMin;
		double                  _mMax;
		std::vector<double>     _masses;
		double                  _mSum;
		rpwa::modelIntensityPtr _model;
		rpwa::eventFileWriter   _fileWriter;

		rpwa::Beam                      _beam;
		rpwa::Target                    _target;
		rpwa::beamAndVertexGeneratorPtr _beamAndVertexGenerator;
		rpwa::massAndTPrimePickerPtr    _pickerFunction;

	};


} //namespace rpwa


#endif // IMPORTANCESAMPLER_H
