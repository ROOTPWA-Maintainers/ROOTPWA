#include "isobarHelicityAmplitude_py.h"

namespace bp = boost::python;

namespace {

	struct isobarHelicityAmplitudeWrapper : public rpwa::isobarHelicityAmplitude,
	                                                bp::wrapper<rpwa::isobarHelicityAmplitude>
	{

		isobarHelicityAmplitudeWrapper()
			: rpwa::isobarHelicityAmplitude(),
			  bp::wrapper<rpwa::isobarHelicityAmplitude>() { };

		isobarHelicityAmplitudeWrapper(const rpwa::isobarDecayTopologyPtr& decay)
			: rpwa::isobarHelicityAmplitude(decay),
			  bp::wrapper<rpwa::isobarHelicityAmplitude>() { };

		isobarHelicityAmplitudeWrapper(const rpwa::isobarHelicityAmplitude& amp)
			: rpwa::isobarHelicityAmplitude(amp),
			  bp::wrapper<rpwa::isobarHelicityAmplitude>() { };


	};

}

void rpwa::py::exportIsobarHelicityAmplitude() {

	bp::class_<isobarHelicityAmplitudeWrapper, bp::bases<rpwa::isobarAmplitude> >("isobarHelicityAmplitude")

		.def(bp::init<rpwa::isobarDecayTopologyPtr>())

		.def(bp::self_ns::str(bp::self))

		.def("name", &isobarHelicityAmplitudeWrapper::name)

		.add_static_property("debugIsobarHelicityAmplitude", &isobarHelicityAmplitudeWrapper::debug, &isobarHelicityAmplitudeWrapper::setDebug);

	bp::register_ptr_to_python<rpwa::isobarHelicityAmplitudePtr>();

};

