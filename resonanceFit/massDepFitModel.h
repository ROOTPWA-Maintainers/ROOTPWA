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
#include <iostream>
#include <limits>
#include <vector>

namespace rpwa {

	namespace massDepFit {

		class component;
		class fsmd;
		class parameters;

		class model {

		public:

			model();
			~model();

			void add(rpwa::massDepFit::component* comp);

			bool init(const std::vector<std::string>& waveNames,
			          const std::string& anchorWaveName,
			          const std::string& anchorComponentName);

			size_t getNrParameters() const { return _nrParameters; }
			void importParameters(const double* par, rpwa::massDepFit::parameters& parameters) const;

			size_t getNrComponents() const { return _components.size(); }
			const rpwa::massDepFit::component* getComponent(size_t idx) const { return _components[idx]; }

			const rpwa::massDepFit::fsmd* getFsmd() const { return _fsmd; }
			void setFsmd(rpwa::massDepFit::fsmd* fsmd);

			size_t getMaxChannelsInComponent() const { return _maxChannelsInComponent; }
			size_t getMaxParametersInComponent() const { return _maxParametersInComponent; }

			bool useBranchings() const { return _useBranchings; }
			void useBranchings(const bool val) { _useBranchings = val; }

			size_t getAnchorWave() const { return _idxAnchorWave; }

			const std::vector<std::pair<size_t, size_t> >& getComponentChannel(size_t idx) const { return _waveComponentChannel[idx]; }

			std::complex<double> productionAmplitude(const rpwa::massDepFit::parameters& fitParameters,
			                                         const size_t idxWave,
			                                         const size_t idxBin,
			                                         const double mass,
			                                         const size_t idxMass = std::numeric_limits<size_t>::max()) const;
			double intensity(const rpwa::massDepFit::parameters& fitParameters,
			                 const size_t idxWave,
			                 const size_t idxBin,
			                 const double mass,
			                 const size_t idxMass = std::numeric_limits<size_t>::max()) const;
			double phaseAbsolute(const rpwa::massDepFit::parameters& fitParameters,
			                     const size_t idxWave,
			                     const size_t idxBin,
			                     const double mass,
			                     const size_t idxMass = std::numeric_limits<size_t>::max()) const;
			std::complex<double> spinDensityMatrix(const rpwa::massDepFit::parameters& fitParameters,
			                                       const size_t idxWave,
			                                       const size_t jdxWave,
			                                       const size_t idxBin,
			                                       const double mass,
			                                       const size_t idxMass = std::numeric_limits<size_t>::max()) const;
			double phase(const rpwa::massDepFit::parameters& fitParameters,
			             const size_t idxWave,
			             const size_t jdxWave,
			             const size_t idxBin,
			             const double mass,
			             const size_t idxMass = std::numeric_limits<size_t>::max()) const;

			std::ostream& print(std::ostream& out = std::cout) const;

		private:

			bool initMapping(const std::string& anchorWaveName,
			                 const std::string& anchorComponentName);

			std::vector<std::string> _waveNames;

			size_t _nrParameters;

			std::vector<rpwa::massDepFit::component*> _components;

			rpwa::massDepFit::fsmd* _fsmd;

			size_t _maxChannelsInComponent;
			size_t _maxParametersInComponent;

			bool _useBranchings;

			size_t _idxAnchorWave;
			size_t _idxAnchorComponent;
			size_t _idxAnchorChannel;

			std::vector<std::vector<std::pair<size_t, size_t> > > _waveComponentChannel;

		};

	} // end namespace massDepFit

} // end namespace rpwa


inline
std::ostream&
operator<< (std::ostream& out, const rpwa::massDepFit::model& fitModel)
{
	return fitModel.print(out);
}


#endif // MASSDEPFITMODEL_HH
