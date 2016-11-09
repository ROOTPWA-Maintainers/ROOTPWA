#include "parameters_py.h"

#include <boost/python.hpp>

#include <resonanceFit/parameters.h>

namespace bp = boost::python;


void rpwa::py::resonanceFit::exportParameters() {

	// the classes and functions from the 'resonanceFit' namespace should
	// not be available from Python at the 'pyRootPwa.core.[...]' level,
	// but should also be included into their own module like
	// 'pyRootPwa.core.resonanceFit.[...]'. So below a submodule
	// 'resonanceFit' is added and the classes and functions are added at
	// that scope.

	bp::scope current;
	std::string submoduleName(bp::extract<std::string>(current.attr("__name__")));
	submoduleName.append(".resonanceFit");

	bp::object submodule(bp::borrowed(PyImport_AddModule(submoduleName.c_str())));
	current.attr("resonanceFit") = submodule;

	// switch the scope to the submodule
	bp::scope submoduleScope = submodule;

	bp::class_<rpwa::resonanceFit::parameters>
		(
			"parameters"
		)

		.def(
			bp::init<size_t, size_t, size_t, size_t>(bp::args("maxComponent",
			                                                  "maxChannels",
			                                                  "maxParameters",
			                                                  "maxBins"))
		)

		.def(
			bp::self_ns::str(bp::self)
		)

		.def(
			"getBranching"
			, &rpwa::resonanceFit::parameters::getBranching
			, (bp::arg("idxComponent"),
			   bp::arg("idxChannel"))
		)

		.def(
			"getCoupling"
			, &rpwa::resonanceFit::parameters::getCoupling
			, (bp::arg("idxComponent"),
			   bp::arg("idxChannel"),
			   bp::arg("idxBin"))
		)

		.def(
			"getParameter"
			, &rpwa::resonanceFit::parameters::getParameter
			, (bp::arg("idxComponent"),
			   bp::arg("idxParameter"))
		)

		.def(
			"setBranching"
			, &rpwa::resonanceFit::parameters::setBranching
			, (bp::arg("idxComponent"),
			   bp::arg("idxChannel"),
			   bp::arg("branching"))
		)

		.def(
			"setCoupling"
			, &rpwa::resonanceFit::parameters::setCoupling
			, (bp::arg("idxComponent"),
			   bp::arg("idxChannel"),
			   bp::arg("idxBin"),
			   bp::arg("coupling"))
		)

		.def(
			"setParameter"
			, &rpwa::resonanceFit::parameters::setParameter
			, (bp::arg("idxComponent"),
			   bp::arg("idxParameter"),
			   bp::arg("parameter"))
		)

		.def(
			"resize"
			, &rpwa::resonanceFit::parameters::resize
			, (bp::arg("maxComponents"),
			   bp::arg("maxChannels"),
			   bp::arg("maxParameters"),
			   bp::arg("maxBins"))
		)

		;

}
