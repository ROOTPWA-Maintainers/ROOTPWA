#include "cache_py.h"

#include <boost/python.hpp>

#include <resonanceFit/cache.h>

namespace bp = boost::python;


void rpwa::py::resonanceFit::exportCache() {

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

	bp::class_<rpwa::resonanceFit::cache>
		(
			"cache"
			, bp::init<size_t, size_t, size_t, size_t, size_t>(bp::args("maxWaves",
			                                                            "maxComponents",
			                                                            "maxChannels",
			                                                            "maxBins",
			                                                            "maxMassBins"))
		)

		.def(
			bp::self_ns::str(bp::self)
		)

		.def(
			"getCoupling"
			, &rpwa::resonanceFit::cache::getCoupling
			, (bp::arg("idxComponent"),
			   bp::arg("idxChannel"),
			   bp::arg("idxBin"),
			   bp::arg("idxMass"))
		)

		.def(
			"getComponent"
			, &rpwa::resonanceFit::cache::getComponent
			, (bp::arg("idxComponent"),
			   bp::arg("idxBin"),
			   bp::arg("idxMass"))
		)

		.def(
			"getProdAmp"
			, &rpwa::resonanceFit::cache::getProdAmp
			, (bp::arg("idxWave"),
			   bp::arg("idxBin"),
			   bp::arg("idxMass"))
		)

		.def(
			"setCoupling"
			, &rpwa::resonanceFit::cache::setCoupling
			, (bp::arg("idxComponent"),
			   bp::arg("idxChannel"),
			   bp::arg("idxBin"),
			   bp::arg("idxMass"),
			   bp::arg("coupling"))
		)

		.def(
			"setComponent"
			, &rpwa::resonanceFit::cache::setComponent
			, (bp::arg("idxComponent"),
			   bp::arg("idxBin"),
			   bp::arg("idxMass"),
			   bp::arg("component"))
		)

		.def(
			"setProdAmp"
			, &rpwa::resonanceFit::cache::setProdAmp
			, (bp::arg("idxWave"),
			   bp::arg("idxBin"),
			   bp::arg("idxMass"),
			   bp::arg("prodAmp"))
		)

		;

}
