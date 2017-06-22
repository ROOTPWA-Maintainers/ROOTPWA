#include "pwaNloptFit_py.h"

#include <boost/python.hpp>

#include "pwaNloptFit.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;

namespace{

	rpwa::fitResultPtr pwaNloptFit_pwaNloptFit(const rpwa::pwaLikelihood<std::complex<double> >& L,
	                                           const bp::dict& multibinBoundaries,
	                                           const unsigned int seed = 0,
	                                           const std::string& startValFileName = "",
	                                           const bool checkHessian = false,
	                                           const bool saveSpace = false,
	                                           const bool verbose = false)
	                                           {
		rpwa::multibinBoundariesType cMultibinBoundaries;
		const bp::list keys = multibinBoundaries.keys();
		for (int i = 0; i < bp::len(keys); ++i) {
			std::string binningVar = bp::extract<std::string>(keys[i]);
			double lowerBound = bp::extract<double>(multibinBoundaries[binningVar][0]);
			double upperBound = bp::extract<double>(multibinBoundaries[binningVar][1]);
			cMultibinBoundaries[binningVar] = std::pair<double, double>(lowerBound, upperBound);
		}
		return rpwa::hli::pwaNloptFit(L, cMultibinBoundaries, seed, startValFileName, checkHessian, saveSpace, verbose);
	}
}


void rpwa::py::exportPwaNloptFit()
{

	bp::def(
		"pwaNloptFit"
		, &pwaNloptFit_pwaNloptFit
		, (bp::arg("likelihood"),
		   bp::arg("multibinBoundaries") = bp::dict(),
		   bp::arg("seed") = 0,
		   bp::arg("startValFileName") = "",
		   bp::arg("checkHessian") = false,
		   bp::arg("saveSpace") = false,
		   bp::arg("verbose") = false)
	);

}
