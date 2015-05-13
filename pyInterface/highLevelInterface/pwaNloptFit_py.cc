
#include "pwaNloptFit_py.h"

#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;

void rpwa::py::exportPwaNloptFit()
{

	bp::def(
		"pwaNloptFit"
		, &rpwa::hli::pwaNloptFit
		, (bp::arg("likelihood"),
		   bp::arg("massBinMin") = 0,
		   bp::arg("massBinMax") = 0,
		   bp::arg("seed") = 0,
		   bp::arg("startValFileName") = "",
		   bp::arg("checkHessian") = false,
		   bp::arg("saveSpace") = false,
		   bp::arg("verbose") = false)
	);

}
