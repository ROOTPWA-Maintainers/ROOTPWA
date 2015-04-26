///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010-2012 Sebastian Neubert (TUM)
//    Copyright 2014,2015 Sebastian Uhl (TUM)
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

#include <Rtypes.h>

namespace libconfig {
	class Setting;
}
class TFile;
class TTree;

namespace rpwa {

	class fitResult;

	namespace massDepFit {

		class cache;
		class function;
		class model;
		class parameters;

		class massDepFit {

		public:

			massDepFit();
			~massDepFit() {};

			bool readConfig(const libconfig::Setting* configRoot,
			                rpwa::massDepFit::model& fitModel,
			                rpwa::massDepFit::parameters& fitParameters,
			                rpwa::massDepFit::parameters& fitParametersError,
			                double& chi2,
			                unsigned int& ndf,
			                const std::string& valTreeName   = "pwa",
			                const std::string& valBranchName = "fitResult_v2");

			bool init(rpwa::massDepFit::model& fitModel,
			          rpwa::massDepFit::function& fitFunction);

			bool updateConfig(libconfig::Setting* configRoot,
			                  const rpwa::massDepFit::model& fitModel,
			                  const rpwa::massDepFit::parameters& fitParameters,
			                  const rpwa::massDepFit::parameters& fitParametersError,
			                  const double chi2,
			                  const unsigned int ndf) const;

// FIXME: make private
			bool createPlots(const rpwa::massDepFit::model& fitModel,
			                 const rpwa::massDepFit::parameters& fitParameters,
			                 rpwa::massDepFit::cache& cache,
			                 TFile* outFile,
			                 const bool rangePlotting) const;

// FIXME: get rid
			const std::vector<std::string>& getFreeParameters() const { return _freeParameters; }

			size_t getNrBins() const { return _nrBins; }
			size_t getNrMassBins() const { return _nrMassBins; }
			size_t getNrWaves() const { return _nrWaves; }

			static void setDebug(bool debug) { _debug = debug; }

		private:

			bool prepareMassLimits();

			bool readConfigFitquality(const libconfig::Setting* configFitquality,
			                          double& chi2,
			                          unsigned int& ndf) const;

			bool readConfigInput(const libconfig::Setting* configInput);
			bool readConfigInputFitResults(const libconfig::Setting* configInputFitResults);
			bool readConfigInputWaves(const libconfig::Setting* configInputWaves);
			bool readConfigInputSystematics(const libconfig::Setting* configInputSystematics);
			bool readConfigInputFreeParameters(const libconfig::Setting* configInputFreeParameters);

			bool readConfigModel(const libconfig::Setting* configModel,
			                     rpwa::massDepFit::model& fitModel,
			                     rpwa::massDepFit::parameters& fitParameters,
			                     rpwa::massDepFit::parameters& fitParametersError);
			bool readConfigModelAnchorWave(const libconfig::Setting* configAnchorWave);
			bool readConfigModelComponents(const libconfig::Setting* configComponents,
			                               rpwa::massDepFit::model& fitModel,
			                               rpwa::massDepFit::parameters& fitParameters,
			                               rpwa::massDepFit::parameters& fitParametersError) const;
			bool readConfigModelFsmd(const libconfig::Setting* configFsmd,
			                         rpwa::massDepFit::model& fitModel,
			                         rpwa::massDepFit::parameters& fitParameters,
			                         rpwa::massDepFit::parameters& fitParametersError) const;

			bool updateConfigFitquality(libconfig::Setting* configFitquality,
			                            const double chi2,
			                            const unsigned int ndf) const;

			bool updateConfigModel(const libconfig::Setting* configModel,
			                       const rpwa::massDepFit::model& fitModel,
			                       const rpwa::massDepFit::parameters& fitParameters,
			                       const rpwa::massDepFit::parameters& fitParametersError) const;
			bool updateConfigModelComponents(const libconfig::Setting* configComponents,
			                                 const rpwa::massDepFit::model& fitModel,
			                                 const rpwa::massDepFit::parameters& fitParameters,
			                                 const rpwa::massDepFit::parameters& fitParametersError) const;
			bool updateConfigModelFsmd(const libconfig::Setting* configFsmd,
			                           const rpwa::massDepFit::model& fitModel,
			                           const rpwa::massDepFit::parameters& fitParameters,
			                           const rpwa::massDepFit::parameters& fitParametersError) const;

			bool readInFiles(const std::string& valTreeName   = "pwa",
			                 const std::string& valBranchName = "fitResult_v2");
			bool readInFileFirst(const std::string& valTreeName   = "pwa",
			                     const std::string& valBranchName = "fitResult_v2");
			bool readInFile(const size_t idxBin,
			                const std::string& valTreeName   = "pwa",
			                const std::string& valBranchName = "fitResult_v2");

			bool readSystematicsFiles(const std::string& valTreeName   = "pwa",
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
			                           boost::multi_array<std::complex<double>, 2>& productionAmplitudes,
			                           boost::multi_array<double, 5>& productionAmplitudesCovariance,
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
			                     const rpwa::massDepFit::parameters& fitParameters,
			                     rpwa::massDepFit::cache& cache,
			                     TDirectory* outDirectory,
			                     const bool rangePlotting,
			                     const size_t idxWave,
			                     const size_t idxBin) const;
			bool createPlotsWavePair(const rpwa::massDepFit::model& fitModel,
			                         const rpwa::massDepFit::parameters& fitParameters,
			                         rpwa::massDepFit::cache& cache,
			                         TDirectory* outDirectory,
			                         const bool rangePlotting,
			                         const size_t idxWave,
			                         const size_t jdxWave,
			                         const size_t idxBin) const;

			std::vector<std::string> _inFileName;
			std::vector<std::vector<std::string> > _inOverwritePhaseSpace;

			bool _sysPlotting;
			std::vector<std::string> _sysFileNames;

			std::vector<double> _tPrimeMeans;

			double _massMax;
			double _massMin;
			double _massStep;
			std::vector<double> _massBinCenters;

			std::vector<std::string> _waveNames;
			std::map<std::string, size_t> _waveIndices;
			std::vector<std::pair<double, double> > _waveMassLimits;
			std::vector<std::pair<size_t, size_t> > _waveMassBinLimits;

			boost::multi_array<std::pair<size_t, size_t>, 2> _wavePairMassBinLimits;

			std::vector<std::string> _freeParameters;

			std::string _anchorWaveName;
			std::string _anchorComponentName;

			boost::multi_array<std::complex<double>, 3> _inProductionAmplitudes;
			boost::multi_array<double, 6> _inProductionAmplitudesCovariance;
			boost::multi_array<std::complex<double>, 4> _inSpinDensityMatrices;
			boost::multi_array<double, 6> _inSpinDensityCovarianceMatrices;
			boost::multi_array<double, 3> _inPhaseSpaceIntegrals;

			boost::multi_array<double, 4> _inIntensities;
			boost::multi_array<double, 5> _inPhases;

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
