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


	class importanceSampler : public BCModel {

	public:

		importanceSampler(const double mMin,
		                  const double mMax,
		                  rpwa::modelIntensityPtr model);

		bool setPhaseSpaceOnly(const bool input = true) { _phaseSpaceOnly = input; return true; }

		bool setMassPrior(TF1* massPrior) { _massPrior = massPrior; return true; }

		bool initializeFileWriter(TFile* outFile);
		bool finalizeFileWriter();

		double LogLikelihood(const std::vector<double>& parameters);
		double LogAPrioriProbability(const std::vector<double>& parameters);

		bool initializeProductionGenerator(rpwa::Beam& beam,
		                                   rpwa::Target& target,
		                                   rpwa::beamAndVertexGeneratorPtr generator,
		                                   rpwa::massAndTPrimePickerPtr picker);

		std::pair<std::pair<bool,double>, std::vector<TLorentzVector> > getProductionKinematics(const double mass);

		size_t nCalls() const { return _nCalls; }

		static std::pair<bool, TLorentzRotation> getInverseGJTransform(const TLorentzVector& beamLv, const TLorentzVector& XLv);

	private:

		std::vector<TLorentzVector> getFinalStateMomenta(const std::vector<double>& parameters) const;

		void CalculateObservables(const std::vector<double>& parameters);

		bool                    _phaseSpaceOnly;
		bool                    _productionGeneratorInitialized;
		bool                    _zeroBinWidth;
		size_t                  _nPart;
		static size_t           _nCalls;
		double                  _mMin;
		double                  _mMax;
		std::vector<double>     _masses;
		TF1*                    _massPrior;
		rpwa::modelIntensityPtr _model;
		rpwa::eventFileWriter   _fileWriter;

		rpwa::Beam                      _beam;
		rpwa::Target                    _target;
		rpwa::beamAndVertexGeneratorPtr _beamAndVertexGenerator;
		rpwa::massAndTPrimePickerPtr    _pickerFunction;

	};


} //namespace rpwa


#endif // IMPORTANCESAMPLER_H
