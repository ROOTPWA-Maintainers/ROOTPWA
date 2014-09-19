#ifndef HLI_CALCAMPLITUDE_H
#define HLI_CALCAMPLITUDE_H

#include <complex>

#include <eventMetadata.h>
#include <isobarAmplitude.h>


namespace rpwa {

	namespace hli {

		std::vector<std::complex<double> > calcAmplitude(const rpwa::eventMetadata&      eventMeta,
		                                                 const rpwa::isobarAmplitudePtr& amplitude,
		                                                 const long int                  maxNmbEvents            = -1,
		                                                 const bool                      printProgress           = true,
		                                                 const std::string&              treePerfStatOutFileName = "",         // root file name for tree performance result
		                                                 const long int                  treeCacheSize           = 25000000);

	}

}

#endif // HLI_CALCAMPLITUDE_H
