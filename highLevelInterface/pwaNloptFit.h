#ifndef HLI_PWANLOPTFIT_H
#define HLI_PWANLOPTFIT_H

#include <fitResult.h>
#include <pwaLikelihood.h>

namespace rpwa {

	namespace hli {

		rpwa::fitResultPtr pwaNloptFit(const rpwa::pwaLikelihood<std::complex<double> >& L,
		                               const binningMapType&                             binningMap = {},
		                               const unsigned int                                seed = 0,
		                               const std::string&                                startValFileName = "",
		                               const bool                                        checkHessian = false,
		                               const bool                                        saveSpace = false,
		                               const bool                                        verbose = false);

	}

}

#endif // HLI_PWANLOPTFIT_H
