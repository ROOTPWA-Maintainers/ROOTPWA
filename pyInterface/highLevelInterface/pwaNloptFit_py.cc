
#include "pwaNloptFit_py.h"

#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;

namespace {

	rpwa::fitResultPtr pwaNloptFit(PyObject*          pyLikelihood,
	                               const int          seed = 0,
	                               const bool         cauchy = false,
	                               const double       massBinLower = 0,
	                               const double       massBinUpper = 0,
	                               const std::string  startValFileName = "",
	                               const bool         checkHessian=false,
	                               const bool         saveSpace=false,
	                               const bool         verbose = false)
	{
		rpwa::pwaLikelihood<std::complex<double> > likelihood = bp::extract<rpwa::pwaLikelihood<std::complex<double> > >(pyLikelihood);
		return rpwa::hli::pwaNloptFit(likelihood,
		                              seed,
		                              cauchy,
		                              massBinLower,
		                              massBinUpper,
		                              startValFileName,
		                              checkHessian,
		                              saveSpace,
		                              verbose);
	}
}

void rpwa::py::exportPwaNloptFit()
{

	bp::def(
		"pwaNloptFit"
		, &::pwaNloptFit
		, (bp::arg("likelihood"),
		   bp::arg("seed") = 0,
		   bp::arg("cauchy") = false,
		   bp::arg("massBinLower") = 0,
		   bp::arg("massBinUpper") = 0,
		   bp::arg("startValFileName") = "",
		   bp::arg("checkHessian") = false,
		   bp::arg("saveSpace") = false,
		   bp::arg("verbose") = false)
	);

}
