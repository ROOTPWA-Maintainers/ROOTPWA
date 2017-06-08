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
#include <string>
#include <vector>

#include "forward.h"
#include "function.h"

class TDirectory;

namespace rpwa {

	class fitResult;

	namespace resonanceFit {

		class cache;
		class parameters;

		void read(const std::string& configFileName,
		          rpwa::resonanceFit::inputConstPtr& fitInput,
		          rpwa::resonanceFit::dataConstPtr& fitData,
		          rpwa::resonanceFit::modelConstPtr& fitModel,
		          rpwa::resonanceFit::parameters& fitParameters,
		          rpwa::resonanceFit::parameters& fitParametersError,
		          std::map<std::string, double>& fitQuality,
		          std::vector<std::string>& freeParameters,
		          const bool useBranchings,
		          const rpwa::resonanceFit::function::useCovarianceMatrix useCovariance,
		          const std::string& valTreeName   = "pwa",
		          const std::string& valBranchName = "fitResult_v2");

		void read(const std::string& configFileName,
		          const double maxMassBinCenter,
		          rpwa::resonanceFit::inputConstPtr& fitInput,
		          rpwa::resonanceFit::modelConstPtr& fitModel,
		          rpwa::resonanceFit::parameters& fitParameters,
		          rpwa::resonanceFit::parameters& fitParametersError,
		          std::map<std::string, double>& fitQuality,
		          std::vector<std::string>& freeParameters,
		          const bool useBranchings);

		void writeConfig(const std::string& configFileName,
		                 const rpwa::resonanceFit::inputConstPtr& fitInput,
		                 const rpwa::resonanceFit::modelConstPtr& fitModel,
		                 const rpwa::resonanceFit::parameters& fitParameters,
		                 const rpwa::resonanceFit::parameters& fitParametersError,
		                 const std::map<std::string, double>& fitQuality,
		                 const std::vector<std::string>& freeParameters);

		void createPlots(const rpwa::resonanceFit::inputConstPtr& fitInput,
		                 const rpwa::resonanceFit::dataConstPtr& fitData,
		                 const rpwa::resonanceFit::modelConstPtr& fitModel,
		                 const rpwa::resonanceFit::parameters& fitParameters,
		                 rpwa::resonanceFit::cache& cache,
		                 TDirectory* mainDirectory,
		                 const bool rangePlotting,
		                 const size_t extraBinning);

		bool debug();
		void setDebug(const bool debug);

	} // end namespace resonanceFit

} // end namespace rpwa

#endif // MASSDEPFIT_HH
