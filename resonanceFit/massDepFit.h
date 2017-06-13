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
#include "function.h"
#include "information.h"

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

		void createPlots(const rpwa::resonanceFit::informationConstPtr& fitInformation,
		                 const rpwa::resonanceFit::dataConstPtr& fitData,
		                 const rpwa::resonanceFit::modelConstPtr& fitModel,
		                 const rpwa::resonanceFit::parameters& fitParameters,
		                 rpwa::resonanceFit::cache& cache,
		                 TFile* outFile,
		                 const bool rangePlotting,
		                 const size_t extraBinning);

		void setDebug(const bool debug);

		class massDepFit {

		public:

			massDepFit();
			~massDepFit() {}

			bool readConfig(const YAML::Node& configRoot,
			                rpwa::resonanceFit::informationConstPtr& fitInformation,
			                rpwa::resonanceFit::dataConstPtr& fitData,
			                const rpwa::resonanceFit::modelPtr& fitModel,
			                rpwa::resonanceFit::parameters& fitParameters,
			                rpwa::resonanceFit::parameters& fitParametersError,
			                std::map<std::string, double>& fitQuality,
			                std::vector<std::string>& freeParameters,
			                const bool useBranchings,
			                const rpwa::resonanceFit::function::useCovarianceMatrix useCovariance,
			                const std::string& valTreeName   = "pwa",
			                const std::string& valBranchName = "fitResult_v2");

			bool init(const rpwa::resonanceFit::informationConstPtr& fitInformation,
			          const rpwa::resonanceFit::dataConstPtr& fitData,
			          const rpwa::resonanceFit::modelPtr& fitModel,
			          const rpwa::resonanceFit::functionPtr& fitFunction);

			bool writeConfig(std::ostream& output,
			                 const rpwa::resonanceFit::informationConstPtr& fitInformation,
			                 const rpwa::resonanceFit::modelConstPtr& fitModel,
			                 const rpwa::resonanceFit::parameters& fitParameters,
			                 const rpwa::resonanceFit::parameters& fitParametersError,
			                 const std::map<std::string, double>& fitQuality,
			                 const std::vector<std::string>& freeParameters) const;

			static void setDebug(bool debug);

		private:

			bool readConfigModel(const YAML::Node& configModel,
			                     const rpwa::resonanceFit::informationConstPtr& fitInformation,
			                     const rpwa::resonanceFit::modelPtr& fitModel,
			                     rpwa::resonanceFit::parameters& fitParameters,
			                     rpwa::resonanceFit::parameters& fitParametersError,
			                     const boost::multi_array<std::string, 2>& waveNames,
			                     const std::vector<size_t>& nrMassBins,
			                     const boost::multi_array<double, 2>& massBinCenters,
			                     const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                     const bool useBranchings);
			bool readConfigModelAnchorWave(const YAML::Node& configAnchorWave);
			bool readConfigModelComponents(const YAML::Node& configComponents,
			                               const rpwa::resonanceFit::informationConstPtr& fitInformation,
			                               const rpwa::resonanceFit::modelPtr& fitModel,
			                               rpwa::resonanceFit::parameters& fitParameters,
			                               rpwa::resonanceFit::parameters& fitParametersError,
			                               const boost::multi_array<std::string, 2>& waveNames,
			                               const std::vector<size_t>& nrMassBins,
			                               const boost::multi_array<double, 2>& massBinCenters,
			                               const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                               const bool useBranchings) const;

			bool writeConfigModel(YAML::Emitter& yamlOutput,
			                      const rpwa::resonanceFit::modelConstPtr& fitModel,
			                      const rpwa::resonanceFit::parameters& fitParameters,
			                      const rpwa::resonanceFit::parameters& fitParametersError) const;
			bool writeConfigModelAnchorWave(YAML::Emitter& yamlOutput) const;
			bool writeConfigModelComponents(YAML::Emitter& yamlOutput,
			                                const rpwa::resonanceFit::modelConstPtr& fitModel,
			                                const rpwa::resonanceFit::parameters& fitParameters,
			                                const rpwa::resonanceFit::parameters& fitParametersError) const;

			std::string _anchorWaveName;
			std::string _anchorComponentName;

			static bool _debug;

		};

	} // end namespace resonanceFit

} // end namespace rpwa

#endif // MASSDEPFIT_HH
