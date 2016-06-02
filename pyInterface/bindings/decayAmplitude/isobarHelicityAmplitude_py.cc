#include "isobarHelicityAmplitude_py.h"

#include <boost/python.hpp>

#include "isobarHelicityAmplitude.h"
#include "rootConverters_py.h"

namespace bp = boost::python;


namespace {

	PyObject* isobarHelicityAmplitude_hfTransform(PyObject* pyDaughterLv) {
		TLorentzVector* daughterLv = rpwa::py::convertFromPy<TLorentzVector*>(pyDaughterLv);
		if(daughterLv == NULL) {
			printErr<<"Got invalid input when executing rpwa::isobarHelicityAmplitude::hfTransform()."<<std::endl;
			return bp::object().ptr();
		}
		return rpwa::py::convertToPy<TLorentzRotation>(rpwa::isobarHelicityAmplitude::hfTransform(*daughterLv));
	}

}

void rpwa::py::exportIsobarHelicityAmplitude() {

	bp::class_<rpwa::isobarHelicityAmplitude, bp::bases<rpwa::isobarAmplitude> >("isobarHelicityAmplitude")

		.def(bp::init<rpwa::isobarDecayTopologyPtr>())

		.def(bp::self_ns::str(bp::self))

		.def("name", &rpwa::isobarHelicityAmplitude::name)

		.def("hfTransform", &isobarHelicityAmplitude_hfTransform)
		.staticmethod("hfTransform")

		.add_static_property("debugIsobarHelicityAmplitude", &rpwa::isobarHelicityAmplitude::debug, &rpwa::isobarHelicityAmplitude::setDebug);

	bp::register_ptr_to_python<rpwa::isobarHelicityAmplitudePtr>();

}
