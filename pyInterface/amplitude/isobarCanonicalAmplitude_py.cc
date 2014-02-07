#include "isobarCanonicalAmplitude_py.h"

namespace bp = boost::python;

void rpwa::py::exportIsobarCanonicalAmplitude() {

	bp::class_<rpwa::isobarCanonicalAmplitude, bp::bases<rpwa::isobarAmplitude> >("isobarCanonicalAmplitude")

		.def(bp::init<rpwa::isobarDecayTopologyPtr>())

		.def(bp::self_ns::str(bp::self))

		.def("name", &rpwa::isobarCanonicalAmplitude::name)

		.add_static_property("debugIsobarCanonicalAmplitude", &rpwa::isobarCanonicalAmplitude::debug, &rpwa::isobarCanonicalAmplitude::setDebug);

	bp::register_ptr_to_python<rpwa::isobarCanonicalAmplitudePtr>();

}

