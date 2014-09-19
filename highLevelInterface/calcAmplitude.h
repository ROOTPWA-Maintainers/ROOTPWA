
#include <complex>

#include <eventMetadata.h>
#include <isobarAmplitude.h>


namespace rpwa {

	namespace hli {

		std::vector<std::complex<double> > calcAmplitude(rpwa::eventMetadata&            eventMeta,
		                                                 const rpwa::isobarAmplitudePtr& amplitude,
		                                                 const long int                  maxNmbEvents            = -1,
		                                                 const bool                      printProgress           = true,
		                                                 const std::string&              treePerfStatOutFileName = "",         // root file name for tree performance result
		                                                 const long int                  treeCacheSize           = 25000000);

	}

}
