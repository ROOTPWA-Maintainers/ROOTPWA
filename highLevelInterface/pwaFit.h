#include <TTree.h>

#include <ampIntegralMatrix.h>
#include <fitResult.h>

namespace rpwa {

	namespace hli {

		rpwa::fitResultPtr pwaFit(std::map<std::string, TTree*>&  ampTrees,
		                          const rpwa::ampIntegralMatrix&  normMatrix,
		                          rpwa::ampIntegralMatrix&        accMatrix,
		                          const std::vector<std::string>& waveNames,
		                          const std::vector<double>&      waveThresholds,
		                          const double                    massBinMin,
		                          const double                    massBinMax,
		                          const unsigned int              seed,
		                          const std::string               startValFileName,
		                          const unsigned int              accEventsOverride,
		                          const unsigned int              rank,
		                          const bool                      runHesse,
		                          const bool                      verbose);
	}

}
