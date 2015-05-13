
#include "pwaFit_py.h"

#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;

namespace {

	rpwa::fitResultPtr pwaFit(PyObject*          pyLikelihood,
	                          const int          seed = 0,
	                          const double       massBinLower = 0,
	                          const double       massBinUpper = 0,
	                          const std::string  startValFileName = "",
			                  const bool         checkHessian=false,
	                          const bool         verbose = false)
	{
		rpwa::pwaLikelihood<std::complex<double> > likelihood = bp::extract<rpwa::pwaLikelihood<std::complex<double> > >(pyLikelihood);
		return rpwa::hli::pwaFit(likelihood,
		                         seed,
		                         massBinLower,
		                         massBinUpper,
		                         startValFileName,
								 checkHessian,
		                         verbose);
	}

}

void rpwa::py::exportPwaFit()
{

	bp::def(
		"pwaFit"
		, &::pwaFit
		, (bp::arg("likelihood"),
		   bp::arg("seed") = 0,
		   bp::arg("massBinLower") = 0,
		   bp::arg("massBinUpper") = 0,
		   bp::arg("startValFileName") = "",
		   bp::arg("checkHessian") = false,
		   bp::arg("verbose") = false)
	);

}
