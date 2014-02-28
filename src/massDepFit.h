#ifndef MASSDEPFIT_HH
#define MASSDEPFIT_HH

#include <map>

#include <boost/multi_array.hpp>

#include <Rtypes.h>

namespace libconfig {
	class Setting;
}
namespace ROOT {
	namespace Math {
		class Minimizer;
	}
}
class TFile;
class TTree;

namespace rpwa {

	class fitResult;

	namespace massDepFit {

		class model;

		class massDepFit {

		public:

			massDepFit();
			~massDepFit() {};

// FIXME: make private
			bool prepareMassLimits();

// FIXME: make private
			bool readConfigInput(const libconfig::Setting* configInput);
// FIXME: make private
			bool readConfigModel(const libconfig::Setting* configRoot,
			                     rpwa::massDepFit::model& fitModel);

// FIXME: make private
			bool updateConfig(libconfig::Setting* configRoot,
			                  const rpwa::massDepFit::model& fitModel,
			                  const ROOT::Math::Minimizer* minimizer,
			                  const double chi2,
			                  const int ndf,
			                  const double chi2red) const;

// FIXME: make private
			bool readInFiles(const std::string& valTreeName   = "pwa",
			                 const std::string& valBranchName = "fitResult_v2");

// FIXME: make private
			bool readSystematicsFiles(const std::string& valTreeName   = "pwa",
			                          const std::string& valBranchName = "fitResult_v2");

// FIXME: make private
			bool createPlots(const rpwa::massDepFit::model& fitModel,
			                 TFile* outFile,
			                 const bool rangePlotting) const;

// FIXME: get rid
			const std::vector<double>& getMassBinCenters() const { return _massBinCenters; }
// FIXME: get rid
			const std::vector<std::string>& getWaveNames() const { return _waveNames; }
// FIXME: get rid
			const std::string& getAnchorWaveName() const { return _anchorWaveName; }
// FIXME: get rid
			const std::string& getAnchorComponentName() const { return _anchorComponentName; }
// FIXME: get rid
			const boost::multi_array<std::pair<size_t, size_t>, 2>& getWavePairMassBinLimits() const { return _wavePairMassBinLimits; }
// FIXME: get rid
			const boost::multi_array<std::complex<double>, 3>& getInSpinDensityMatrices() const { return _inSpinDensityMatrices; }
// FIXME: get rid
			const boost::multi_array<double, 5>& getInSpinDensityCovarianceMatrices() const { return _inSpinDensityCovarianceMatrices; }

			static void setDebug(bool debug) { _debug = debug; }

		private:

			bool readConfigInputFitResults(const libconfig::Setting* configInputFitResults);
			bool readConfigInputWaves(const libconfig::Setting* configInputWaves);
			bool readConfigInputSystematics(const libconfig::Setting* configInputSystematics);

			bool readConfigModelAnchorWave(const libconfig::Setting* configAnchorWave);
			bool readConfigModelComponents(const libconfig::Setting* configComponents,
			                               rpwa::massDepFit::model& fitModel) const;
			bool readConfigModelFsmd(const libconfig::Setting* configFsmd,
			                         rpwa::massDepFit::model& fitModel) const;

			bool updateConfigModel(const libconfig::Setting* configModel,
			                       const rpwa::massDepFit::model& fitModel,
			                       const ROOT::Math::Minimizer* minimizer) const;
			bool updateConfigModelComponents(const libconfig::Setting* configComponents,
			                                 const rpwa::massDepFit::model& fitModel,
			                                 const ROOT::Math::Minimizer* minimizer) const;
			bool updateConfigModelFsmd(const libconfig::Setting* configFsmd,
			                           const rpwa::massDepFit::model& fitModel,
			                           const ROOT::Math::Minimizer* minimizer) const;

			bool readInFile(const std::string& fileName,
			                const std::vector<std::string>& overwritePhaseSpace,
			                const std::string& valTreeName   = "pwa",
			                const std::string& valBranchName = "fitResult_v2");

			bool readSystematicsFile(const size_t idxSystematics,
			                         const std::string& valTreeName   = "pwa",
			                         const std::string& valBranchName = "fitResult_v2");

			bool checkFitResultMassBins(TTree* tree,
			                            rpwa::fitResult* fit,
			                            std::vector<Long64_t>& mapping) const;
			bool readFitResultMassBins(TTree* tree,
			                           rpwa::fitResult* fit);
			bool readFitResultMatrices(TTree* tree,
			                           rpwa::fitResult* fit,
			                           const std::vector<Long64_t>& mapping,
			                           boost::multi_array<std::complex<double>, 3>& spinDensityMatrices,
			                           boost::multi_array<double, 5>& spinDensityCovarianceMatrices,
			                           boost::multi_array<double, 3>& intensities,
			                           boost::multi_array<double, 4>& phases) const;
			bool readFitResultIntegrals(TTree* tree,
			                            rpwa::fitResult* fit,
			                            const std::vector<Long64_t>& mapping,
			                            boost::multi_array<double, 2>& phaseSpaceIntegrals) const;
			bool readPhaseSpaceIntegralMatrices(const std::vector<std::string>& overwritePhaseSpace,
			                                    boost::multi_array<double, 2>& phaseSpaceIntegrals) const;

			bool createPlotsWave(const rpwa::massDepFit::model& fitModel,
			                     TFile* outFile,
			                     const bool rangePlotting,
			                     const size_t idxWave) const;
			bool createPlotsWavePair(const rpwa::massDepFit::model& fitModel,
			                         TFile* outFile,
			                         const bool rangePlotting,
			                         const size_t idxWave,
			                         const size_t jdxWave) const;

			std::vector<std::string> _inFileName;
			std::vector<std::vector<std::string> > _inOverwritePhaseSpace;

			bool _sysPlotting;
			std::vector<std::string> _sysFileNames;

			double _massMax;
			double _massMin;
			double _massStep;
			std::vector<double> _massBinCenters;

			std::vector<std::string> _waveNames;
			std::map<std::string, size_t> _waveIndices;
			std::vector<std::pair<double, double> > _waveMassLimits;
			std::vector<std::pair<size_t, size_t> > _waveMassBinLimits;

			boost::multi_array<std::pair<size_t, size_t>, 2> _wavePairMassBinLimits;

			std::string _anchorWaveName;
			std::string _anchorComponentName;

			boost::multi_array<std::complex<double>, 3> _inSpinDensityMatrices;
			boost::multi_array<double, 5> _inSpinDensityCovarianceMatrices;
			boost::multi_array<double, 2> _inPhaseSpaceIntegrals;

			boost::multi_array<double, 3> _inIntensities;
			boost::multi_array<double, 4> _inPhases;

			boost::multi_array<std::complex<double>, 4> _sysSpinDensityMatrices;
			boost::multi_array<double, 6> _sysSpinDensityCovarianceMatrices;

			boost::multi_array<double, 4> _sysIntensities;
			boost::multi_array<double, 5> _sysPhases;

			size_t _nrBins;
			size_t _nrMassBins;
			size_t _nrSystematics;
			size_t _nrWaves;

			static bool _debug;

		};

	} // end namespace massDepFit

} // end namespace rpwa

#endif
