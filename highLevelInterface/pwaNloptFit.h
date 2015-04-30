#include <TTree.h>

#include <ampIntegralMatrix.h>
#include <fitResult.h>

namespace rpwa {

	namespace hli {

		rpwa::fitResultPtr pwaNloptFit(std::map<std::string, TTree*>&  ampTrees,
		                               const rpwa::ampIntegralMatrix&  normMatrix,
		                               rpwa::ampIntegralMatrix&        accMatrix,
		                               const std::vector<std::string>& waveNames,
		                               const std::vector<double>&      waveThresholds,
		                               const double                    massBinMin,
		                               const double                    massBinMax,
		                               const unsigned int              seed,
		                               const bool                      cauchy,
		                               const std::string               startValFileName,
		                               const unsigned int              accEventsOverride,
		                               const unsigned int              rank,
		                               const bool                      checkHessian,
		                               const bool                      saveSpace,
		                               const bool                      verbose);
	}

}
