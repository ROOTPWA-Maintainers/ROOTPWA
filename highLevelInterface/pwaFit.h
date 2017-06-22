#ifndef HLI_PWAFIT_H
#define HLI_PWAFIT_H

#include <fitResult.h>
#include <pwaLikelihood.h>

namespace rpwa {

	namespace hli {

		rpwa::fitResultPtr pwaFit(const rpwa::pwaLikelihood<std::complex<double> >& L,
		                          const rpwa::multibinBoundariesType&               multibinBoundaries = rpwa::multibinBoundariesType(),
		                          const unsigned int                                seed = 0,
		                          const std::string&                                startValFileName = "",
		                          const bool                                        checkHessian = false,
		                          const bool                                        saveSpace = false,
		                          const bool                                        verbose = false);

	}

}

#endif // HLI_PWAFIT_H
