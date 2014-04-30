//-----------------------------------------------------------
//
// Description:
//    (BW) Component of mass dependent fit
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef MASSDEPFITMODEL_HH
#define MASSDEPFITMODEL_HH

#include <complex>
#include <limits>
#include <vector>

class TF1;

namespace rpwa {

	namespace massDepFit {

		class component;

		class model {

		public:

			model();
			~model();

			void add(rpwa::massDepFit::component* comp);

			TF1* getFsmdFunction() const { return _fsmdFunction; }
			void setFsmdFunction(TF1* fsmdFunction);

			bool init(const std::vector<std::string>& waveNames,
			          const std::vector<double>& massBinCenters,
			          const std::string& anchorWaveName,
			          const std::string& anchorComponentName);

			size_t getNrParameters() const { return _nrParameters; }
			void getParameters(double* par) const;
			void setParameters(const double* par);

			size_t getNrComponents() const { return _components.size(); }
			const rpwa::massDepFit::component* getComponent(size_t idx) const { return _components[idx]; }

			unsigned int getFsmdNrParameters() const { return _fsmdFreeParameterIndices.size(); }
			double getFsmdParameter(const size_t idx) const;
			void getFsmdParameterLimits(const size_t idx, double& lower, double& upper) const;

			const std::vector<std::pair<size_t, size_t> >& getComponentChannel(size_t idx) const { return _waveComponentChannel[idx]; }

			size_t getAnchorWave() const { return _idxAnchorWave; }

			std::complex<double> productionAmplitude(const size_t idxWave,
			                                         const size_t idxBin,
			                                         const double mass,
			                                         const size_t idxMass = std::numeric_limits<size_t>::max()) const;
			double intensity(const size_t idxWave,
			                 const size_t idxBin,
			                 const double mass,
			                 const size_t idxMass = std::numeric_limits<size_t>::max()) const;
			double phaseAbsolute(const size_t idxWave,
			                     const size_t idxBin,
			                     const double mass,
			                     const size_t idxMass = std::numeric_limits<size_t>::max()) const;
			std::complex<double> spinDensityMatrix(const size_t idxWave,
			                                       const size_t jdxWave,
			                                       const size_t idxBin,
			                                       const double mass,
			                                       const size_t idxMass = std::numeric_limits<size_t>::max()) const;
			double phase(const size_t idxWave,
			             const size_t jdxWave,
			             const size_t idxBin,
			             const double mass,
			             const size_t idxMass = std::numeric_limits<size_t>::max()) const;

			double calcFsmd(const double mass,
			                const size_t idxMass = std::numeric_limits<size_t>::max()) const;

			std::ostream& print(std::ostream& out) const;

		private:

			bool initMapping(const std::string& anchorWaveName,
			                 const std::string& anchorComponentName);
			bool initFsmd(const std::vector<double>& massBinCenters);

			std::vector<std::string> _waveNames;

			size_t _nrParameters;

			std::vector<rpwa::massDepFit::component*> _components;

			size_t _idxAnchorWave;
			size_t _idxAnchorComponent;
			size_t _idxAnchorChannel;

			TF1* _fsmdFunction;
			bool _fsmdFixed;
			std::vector<double> _fsmdValues;
			std::vector<size_t> _fsmdFreeParameterIndices;

			std::vector<std::vector<std::pair<size_t, size_t> > > _waveComponentChannel;

		};

	} // end namespace massDepFit

} // end namespace rpwa


inline
std::ostream&
operator<< (std::ostream& out, const rpwa::massDepFit::model& fitModel) {
	return fitModel.print(out);
}


#endif // MASSDEPFITMODEL_HH
