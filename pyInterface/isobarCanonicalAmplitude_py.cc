#include "isobarCanonicalAmplitude_py.h"

namespace bp = boost::python;

namespace {

	struct isobarCanonicalAmplitudeWrapper : public rpwa::isobarCanonicalAmplitude,
	                                                bp::wrapper<rpwa::isobarCanonicalAmplitude>
	{

		isobarCanonicalAmplitudeWrapper()
			: rpwa::isobarCanonicalAmplitude(),
			  bp::wrapper<rpwa::isobarCanonicalAmplitude>() { };

		isobarCanonicalAmplitudeWrapper(const rpwa::isobarDecayTopologyPtr& decay)
			: rpwa::isobarCanonicalAmplitude(decay),
			  bp::wrapper<rpwa::isobarCanonicalAmplitude>() { };

		isobarCanonicalAmplitudeWrapper(const rpwa::isobarCanonicalAmplitude& amp)
			: rpwa::isobarCanonicalAmplitude(amp),
			  bp::wrapper<rpwa::isobarCanonicalAmplitude>() { };


	};

}

void rpwa::py::exportIsobarCanonicalAmplitude() {

	bp::class_<isobarCanonicalAmplitudeWrapper, bp::bases<rpwa::isobarAmplitude> >("isobarCanonicalAmplitude")

		.def(bp::init<rpwa::isobarDecayTopologyPtr>())

		.def(bp::self_ns::str(bp::self))

		.def("name", &isobarCanonicalAmplitudeWrapper::name)

		.add_static_property("debugIsobarCanonicalAmplitude", &isobarCanonicalAmplitudeWrapper::debug, &isobarCanonicalAmplitudeWrapper::setDebug);

	bp::register_ptr_to_python<rpwa::isobarCanonicalAmplitudePtr>();

};

