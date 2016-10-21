///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010-2012 Sebastian Neubert (TUM)
//    Copyright 2014-2016 Sebastian Uhl (TUM)
//
//    This file is part of ROOTPWA
//
//    ROOTPWA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ROOTPWA is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ROOTPWA.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      master class of the resonance fit
//      - take care of reading the configuration file and setting up the
//        fit model including the components accordingly
//      - take care of reading the fit results from the partial-wave fit
//      - create plots of the final result
//      - update the configuration file with the results from the resonance
//        fit
//
//-------------------------------------------------------------------------


#ifndef MASSDEPFIT_HH
#define MASSDEPFIT_HH

#include <map>

#include <boost/multi_array.hpp>

#include <TMatrixT.h>

#include "forward.h"

namespace YAML {
	class Emitter;
	class Node;
}
class TFile;
class TTree;

namespace rpwa {

	class fitResult;

	namespace resonanceFit {

		class cache;
		class parameters;

		class massDepFit {

		public:

			massDepFit();
			~massDepFit() {}

			bool readConfig(const YAML::Node& configRoot,
			                rpwa::resonanceFit::dataConstPtr& fitData,
			                const rpwa::resonanceFit::modelPtr& fitModel,
			                rpwa::resonanceFit::parameters& fitParameters,
			                rpwa::resonanceFit::parameters& fitParametersError,
			                std::map<std::string, double>& fitQuality,
			                const bool useBranchings,
			                const std::string& valTreeName   = "pwa",
			                const std::string& valBranchName = "fitResult_v2");

			bool init(const rpwa::resonanceFit::dataConstPtr& fitData,
			          const rpwa::resonanceFit::modelPtr& fitModel,
			          const rpwa::resonanceFit::functionPtr& fitFunction);

			bool writeConfig(std::ostream& output,
			                 const rpwa::resonanceFit::modelConstPtr& fitModel,
			                 const rpwa::resonanceFit::parameters& fitParameters,
			                 const rpwa::resonanceFit::parameters& fitParametersError,
			                 const std::map<std::string, double>& fitQuality) const;

// FIXME: make private
			bool createPlots(const rpwa::resonanceFit::dataConstPtr& fitData,
			                 const rpwa::resonanceFit::modelConstPtr& fitModel,
			                 const rpwa::resonanceFit::parameters& fitParameters,
			                 rpwa::resonanceFit::cache& cache,
			                 TFile* outFile,
			                 const bool rangePlotting,
			                 const size_t extraBinning) const;

// FIXME: get rid
			const std::vector<std::string>& getFreeParameters() const { return _freeParameters; }

			size_t getNrBins() const { return _nrBins; }
			size_t getMaxMassBins() const { return _maxMassBins; }
			size_t getNrWaves() const { return _nrWaves; }

			static void setDebug(bool debug) { _debug = debug; }

		private:

			bool prepareMassLimits(const std::vector<size_t>& nrMassBins,
			                       const boost::multi_array<double, 2>& massBinCenters);
			bool prepareMassLimit(const size_t nrMassBins,
			                      const boost::multi_array<double, 1>& massBinCenters,
			                      const size_t idxBin);

			bool readConfigFitquality(const YAML::Node& configFitquality,
			                          std::map<std::string, double>& fitQuality) const;

			bool readConfigInput(const YAML::Node& configInput);
			bool readConfigInputFitResults(const YAML::Node& configInputFitResults);
			bool readConfigInputFitResultSystematics(const YAML::Node& configInputFitResultSystematics);
			bool readConfigInputWaves(const YAML::Node& configInputWaves);
			bool readConfigInputFreeParameters(const YAML::Node& configInputFreeParameters);

			bool readConfigModel(const YAML::Node& configModel,
			                     const rpwa::resonanceFit::modelPtr& fitModel,
			                     rpwa::resonanceFit::parameters& fitParameters,
			                     rpwa::resonanceFit::parameters& fitParametersError,
			                     const std::vector<size_t>& nrMassBins,
			                     const boost::multi_array<double, 2>& massBinCenters,
			                     const bool useBranchings);
			bool readConfigModelAnchorWave(const YAML::Node& configAnchorWave);
			bool readConfigModelComponents(const YAML::Node& configComponents,
			                               const rpwa::resonanceFit::modelPtr& fitModel,
			                               rpwa::resonanceFit::parameters& fitParameters,
			                               rpwa::resonanceFit::parameters& fitParametersError,
			                               const std::vector<size_t>& nrMassBins,
			                               const boost::multi_array<double, 2>& massBinCenters,
			                               const bool useBranchings) const;
			bool readConfigModelFsmd(const YAML::Node& configFsmd,
			                         const rpwa::resonanceFit::modelPtr& fitModel,
			                         rpwa::resonanceFit::parameters& fitParameters,
			                         rpwa::resonanceFit::parameters& fitParametersError) const;

			bool writeConfigFitquality(YAML::Emitter& yamlOutput,
			                           const std::map<std::string, double>& fitQuality) const;

			bool writeConfigInput(YAML::Emitter& yamlOutput) const;
			bool writeConfigInputFitResults(YAML::Emitter& yamlOutput) const;
			bool writeConfigInputFitResultSystematics(YAML::Emitter& yamlOutput,
			                                          const size_t idxBin) const;
			bool writeConfigInputWaves(YAML::Emitter& yamlOutput) const;
			bool writeConfigInputFreeParameters(YAML::Emitter& yamlOutput) const;

			bool writeConfigModel(YAML::Emitter& yamlOutput,
			                      const rpwa::resonanceFit::modelConstPtr& fitModel,
			                      const rpwa::resonanceFit::parameters& fitParameters,
			                      const rpwa::resonanceFit::parameters& fitParametersError) const;
			bool writeConfigModelAnchorWave(YAML::Emitter& yamlOutput) const;
			bool writeConfigModelComponents(YAML::Emitter& yamlOutput,
			                                const rpwa::resonanceFit::modelConstPtr& fitModel,
			                                const rpwa::resonanceFit::parameters& fitParameters,
			                                const rpwa::resonanceFit::parameters& fitParametersError) const;
			bool writeConfigModelFsmd(YAML::Emitter& yamlOutput,
			                          const rpwa::resonanceFit::modelConstPtr& fitModel,
			                          const rpwa::resonanceFit::parameters& fitParameters,
			                          const rpwa::resonanceFit::parameters& fitParametersError) const;

			bool readInFiles(std::vector<size_t>& nrMassBins,
			                 boost::multi_array<double, 2>& massBinCenters,
			                 const std::string& valTreeName   = "pwa",
			                 const std::string& valBranchName = "fitResult_v2");
			bool readInFile(const size_t idxBin,
			                std::vector<size_t>& nrMassBins,
			                boost::multi_array<double, 2>& massBinCenters,
			                const std::string& valTreeName   = "pwa",
			                const std::string& valBranchName = "fitResult_v2");

			bool readSystematicsFiles(const size_t idxBin,
			                          const size_t nrMassBins,
			                          const boost::multi_array<double, 1>& massBinCenters,
			                          const std::string& valTreeName   = "pwa",
			                          const std::string& valBranchName = "fitResult_v2");
			bool readSystematicsFile(const size_t idxBin,
			                         const size_t idxSystematics,
			                         const size_t nrMassBins,
			                         const boost::multi_array<double, 1>& massBinCenters,
			                         const std::string& valTreeName   = "pwa",
			                         const std::string& valBranchName = "fitResult_v2");

			bool checkFitResultMassBins(TTree* tree,
			                            rpwa::fitResult* fit,
			                            const size_t nrMassBins,
			                            const boost::multi_array<double, 1>& massBinCenters,
			                            std::vector<Long64_t>& mapping) const;
			bool readFitResultMassBins(TTree* tree,
			                           rpwa::fitResult* fit,
			                           size_t& nrMassBins,
			                           boost::multi_array<double, 1>& massBinCenters) const;
			bool readFitResultMatrices(TTree* tree,
			                           rpwa::fitResult* fit,
			                           const std::vector<Long64_t>& mapping,
			                           const double rescaleErrors,
			                           std::vector<std::string>& waveNames,
			                           boost::multi_array<std::complex<double>, 2>& productionAmplitudes,
			                           boost::multi_array<TMatrixT<double>, 1>& productionAmplitudesCovariance,
			                           boost::multi_array<std::complex<double>, 3>& spinDensityMatrices,
			                           boost::multi_array<TMatrixT<double>, 1>& spinDensityCovarianceMatrices,
			                           boost::multi_array<std::pair<double, double>, 2>& plottingIntensities,
			                           boost::multi_array<std::pair<double, double>, 3>& plottingSpinDensityMatrixElementsReal,
			                           boost::multi_array<std::pair<double, double>, 3>& plottingSpinDensityMatrixElementsImag,
			                           boost::multi_array<std::pair<double, double>, 3>& plottingPhases) const;
			bool readFitResultIntegrals(TTree* tree,
			                            rpwa::fitResult* fit,
			                            const std::vector<Long64_t>& mapping,
			                            const std::vector<std::string>& waveNames,
			                            boost::multi_array<double, 2>& phaseSpaceIntegrals) const;

			bool createPlotsWave(const rpwa::resonanceFit::dataConstPtr& fitData,
			                     const rpwa::resonanceFit::modelConstPtr& fitModel,
			                     const rpwa::resonanceFit::parameters& fitParameters,
			                     rpwa::resonanceFit::cache& cache,
			                     TDirectory* outDirectory,
			                     const bool rangePlotting,
			                     const size_t extraBinning,
			                     const size_t idxWave,
			                     const size_t idxBin) const;
			bool createPlotsWaveSum(const rpwa::resonanceFit::dataConstPtr& fitData,
			                        const rpwa::resonanceFit::modelConstPtr& fitModel,
			                        const rpwa::resonanceFit::parameters& fitParameters,
			                        rpwa::resonanceFit::cache& cache,
			                        TDirectory* outDirectory,
			                        const bool rangePlotting,
			                        const size_t extraBinning,
			                        const size_t idxWave) const;
			bool createPlotsWavePair(const rpwa::resonanceFit::dataConstPtr& fitData,
			                         const rpwa::resonanceFit::modelConstPtr& fitModel,
			                         const rpwa::resonanceFit::parameters& fitParameters,
			                         rpwa::resonanceFit::cache& cache,
			                         TDirectory* outDirectory,
			                         const bool rangePlotting,
			                         const size_t extraBinning,
			                         const size_t idxWave,
			                         const size_t jdxWave,
			                         const size_t idxBin) const;
			bool createPlotsFsmd(const rpwa::resonanceFit::dataConstPtr& fitData,
			                     const rpwa::resonanceFit::modelConstPtr& fitModel,
			                     const rpwa::resonanceFit::parameters& fitParameters,
			                     rpwa::resonanceFit::cache& cache,
			                     TDirectory* outDirectory,
			                     const bool rangePlotting,
			                     const size_t extraBinning,
			                     const size_t idxBin) const;

			std::vector<std::string> _inFileName;

			std::vector<size_t> _nrSystematics;
			std::vector<std::vector<std::string> > _sysFileNames;

			std::vector<double> _rescaleErrors;
			std::vector<double> _tPrimeMeans;

			bool _sameMassBinning;

			std::vector<std::string> _waveNames;
			std::vector<std::vector<std::string> > _waveNameAlternatives;
			std::map<std::string, size_t> _waveIndices;
			std::map<std::string, std:: vector<size_t> > _waveBins;
			std::vector<std::pair<double, double> > _waveMassLimits;
			boost::multi_array<std::pair<size_t, size_t>, 2> _waveMassBinLimits;
			boost::multi_array<std::pair<size_t, size_t>, 3> _wavePairMassBinLimits;

			std::vector<std::string> _freeParameters;

			std::string _anchorWaveName;
			std::string _anchorComponentName;

			boost::multi_array<std::complex<double>, 3> _inProductionAmplitudes;
			boost::multi_array<TMatrixT<double>, 2> _inProductionAmplitudesCovariance;
			boost::multi_array<std::complex<double>, 4> _inSpinDensityMatrices;
			boost::multi_array<TMatrixT<double>, 2> _inSpinDensityCovarianceMatrices;
			boost::multi_array<double, 3> _inPhaseSpaceIntegrals;

			boost::multi_array<std::pair<double, double>, 3> _plottingIntensities;
			boost::multi_array<std::pair<double, double>, 4> _plottingSpinDensityMatrixElementsReal;
			boost::multi_array<std::pair<double, double>, 4> _plottingSpinDensityMatrixElementsImag;
			boost::multi_array<std::pair<double, double>, 4> _plottingPhases;

			boost::multi_array<std::pair<double, double>, 3> _sysPlottingIntensities;
			boost::multi_array<std::pair<double, double>, 4> _sysPlottingSpinDensityMatrixElementsReal;
			boost::multi_array<std::pair<double, double>, 4> _sysPlottingSpinDensityMatrixElementsImag;
			boost::multi_array<std::pair<double, double>, 4> _sysPlottingPhases;

			size_t _nrBins;
			size_t _maxMassBins;
			size_t _nrWaves;

			static bool _debug;

		};

	} // end namespace resonanceFit

} // end namespace rpwa

#endif // MASSDEPFIT_HH
