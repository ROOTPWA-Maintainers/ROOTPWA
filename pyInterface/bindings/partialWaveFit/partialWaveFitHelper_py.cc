#include "partialWaveFitHelper_py.h"

#include <boost/python.hpp>

#include "partialWaveFitHelper.h"

namespace bp = boost::python;


void rpwa::py::exportPartialWaveFitHelper() {

	// the functions from the 'partialWaveFitHelper' namespace should not
	// be available from Python at the 'pyRootPwa.core.function' level,
	// but should also be included into their own module like
	// 'pyRootPwa.core.partialWaveFitHelper.function'. So below a submodule
	// 'partialWaveFitHelper' is added and the functions are added at that
	// scope.

	bp::scope current;
	std::string submoduleName(bp::extract<std::string>(current.attr("__name__")));
	submoduleName.append(".partialWaveFitHelper");

	bp::object submodule(bp::borrowed(PyImport_AddModule(submoduleName.c_str())));
	current.attr("partialWaveFitHelper") = submodule;

	{

		// switch the scope to the submodule
		bp::scope submoduleScope = submodule;

		bp::def(
			"getReflectivity"
			, &rpwa::partialWaveFitHelper::getReflectivity
			, (bp::arg("name"))
		);

	}

}
