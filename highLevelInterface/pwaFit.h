#include <pwaLikelihood.h>
#include <fitResult.h>

namespace rpwa {

	namespace hli {

		rpwa::fitResultPtr pwaFit(const rpwa::pwaLikelihood<std::complex<double> >& L,
		                          const unsigned int              seed,
		                          const double                    massBinMin,
		                          const double                    massBinMax,
		                          const std::string               startValFileName,
		                          const bool                      checkHessian,
		                          const bool                      verbose);
	}

}
