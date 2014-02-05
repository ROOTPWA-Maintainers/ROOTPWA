
#include "generatorManager_py.h"

#include<generator.h>

namespace bp = boost::python;

void rpwa::py::exportGeneratorManager() {

	bp::class_<rpwa::generatorManager>("generatorManager")
		.def(bp::self_ns::str(bp::self))
		.def("event", &rpwa::generatorManager::event)
		.def(
			"getGenerator"
			, &rpwa::generatorManager::getGenerator
			, bp::return_internal_reference<1>())
		.def("readReactionFile", &rpwa::generatorManager::readReactionFile)
		.def("initializeGenerator", &rpwa::generatorManager::initializeGenerator)
		.def("overrideMassRange", &rpwa::generatorManager::overrideMassRange)
		.def("overrideBeamFile", &rpwa::generatorManager::overrideBeamFile)
		.def(
			"readBeamfileSequentially"
			, &rpwa::generatorManager::readBeamfileSequentially
			, bp::arg("readBeamfileSequentially")=true
		)
		.def("randomizeBeamfileStartingPosition", &rpwa::generatorManager::randomizeBeamfileStartingPosition)
		.add_static_property("debugGeneratorManager", &rpwa::generatorManager::debug, &rpwa::generatorManager::setDebug);

}
