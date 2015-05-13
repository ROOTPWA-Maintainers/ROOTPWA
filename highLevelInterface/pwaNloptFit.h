#include <TTree.h>

#include <pwaLikelihood.h>
#include <fitResult.h>

namespace rpwa {

	namespace hli {

		rpwa::fitResultPtr pwaNloptFit(rpwa::pwaLikelihood<std::complex<double> >& L,
		                               const unsigned int              seed,
		                               const bool                      cauchy,
				                       const double                    massBinMin,
				                       const double                    massBinMax,
		                               const std::string               startValFileName,
		                               const bool                      checkHessian,
		                               const bool                      saveSpace,
		                               const bool                      verbose);
	}

}
