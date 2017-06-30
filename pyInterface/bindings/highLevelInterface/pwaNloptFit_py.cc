#include "pwaNloptFit_py.h"

#include <boost/python.hpp>

#include "pwaNloptFit.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;

namespace{

	rpwa::fitResultPtr pwaNloptFit_pwaNloptFit(const rpwa::pwaLikelihood<std::complex<double> >& L,
	                                           const bp::dict& pyMultibinBoundaries = bp::dict(),
	                                           const unsigned int seed = 0,
	                                           const std::string& startValFileName = "",
	                                           const bool checkHessian = false,
	                                           const bool saveSpace = false,
	                                           const bool verbose = false)
	{
		const rpwa::multibinBoundariesType multibinBoundaries = rpwa::py::convertMultibinBoundariesFromPy(pyMultibinBoundaries);
		return rpwa::hli::pwaNloptFit(L, multibinBoundaries, seed, startValFileName, checkHessian, saveSpace, verbose);
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
