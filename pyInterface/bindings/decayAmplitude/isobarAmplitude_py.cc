#include "isobarAmplitude_py.h"

#include <boost/python.hpp>

#include "isobarAmplitude.h"
#include "rootConverters_py.h"

namespace bp = boost::python;


namespace {

	struct isobarAmplitudeWrapper : public rpwa::isobarAmplitude,
	                                       bp::wrapper<rpwa::isobarAmplitude>
	{

		void init() {
			if(bp::override init = this->get_override("init")) {
				init();
			}
			rpwa::isobarAmplitude::init();
		}

		void default_init() {
			rpwa::isobarAmplitude::init();
		}

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::isobarAmplitude::name();
		}

		std::string default_name() const {
			return rpwa::isobarAmplitude::name();
		}

		std::string printParameters__() const {
			std::stringstream sstr;
			if(bp::override printParameters = this->get_override("printParameters")) {
				printParameters(sstr);
			} else {
				rpwa::isobarAmplitude::printParameters(sstr);
			}
			return sstr.str();
		}

		std::string default_printParameters__() const {
			std::stringstream sstr;
			rpwa::isobarAmplitude::printParameters(sstr);
			return sstr.str();
		}

	};

	rpwa::isobarDecayTopology& isobarAmplitude_decayTopology(const rpwa::isobarAmplitude& self) {
		return *(self.decayTopology());
	}

	PyObject* isobarAmplitude_gjTransform(PyObject* pyBeamLv, PyObject* pyXLv) {
		TLorentzVector* beamLv = rpwa::py::convertFromPy<TLorentzVector*>(pyBeamLv);
		TLorentzVector* XLv = rpwa::py::convertFromPy<TLorentzVector*>(pyXLv);
		if((beamLv == NULL) || (XLv == NULL)) {
			printErr<<"Got invalid input when executing rpwa::isobarAmplitude::gjTransform()."<<std::endl;
			return bp::object().ptr();
		}
		return rpwa::py::convertToPy<TLorentzRotation>(rpwa::isobarAmplitude::gjTransform(*beamLv, *XLv));
	}

	std::string isobarAmplitude_printParameters(const rpwa::isobarAmplitude& self) {
		std::stringstream sstr;
		self.printParameters(sstr);
		return sstr.str();
	}

}

void rpwa::py::exportIsobarAmplitude() {

	bp::class_<isobarAmplitudeWrapper, boost::noncopyable>("isobarAmplitude", bp::no_init)

		.def(bp::self_ns::str(bp::self))

		.def(
			"decayTopology"
			, &isobarAmplitude_decayTopology
			, bp::return_internal_reference<>()
		)

		.def("setDecayTopology", &rpwa::isobarAmplitude::setDecayTopology)

		.def("init", &isobarAmplitudeWrapper::init, &isobarAmplitudeWrapper::default_init)
		.def("init", &rpwa::isobarAmplitude::init)

		.add_property("reflectivityBasis", &rpwa::isobarAmplitude::reflectivityBasis, &rpwa::isobarAmplitude::enableReflectivityBasis)
		.add_property("boseSymmetrization", &rpwa::isobarAmplitude::boseSymmetrization, &rpwa::isobarAmplitude::enableBoseSymmetrization)
		.add_property("isospinSymmetrization", &rpwa::isobarAmplitude::isospinSymmetrization, &rpwa::isobarAmplitude::enableIsospinSymmetrization)
		.add_property("doSpaceInversion", &rpwa::isobarAmplitude::doSpaceInversion, &rpwa::isobarAmplitude::enableSpaceInversion)
		.add_property("doReflection", &rpwa::isobarAmplitude::doReflection, &rpwa::isobarAmplitude::enableReflection)

		.def("gjTransform", &isobarAmplitude_gjTransform)
		.staticmethod("gjTransform")

		.def("amplitude", &rpwa::isobarAmplitude::amplitude)

		.def("__call__", &rpwa::isobarAmplitude::operator())

		.def("name", &isobarAmplitudeWrapper::name, &isobarAmplitudeWrapper::default_name)
		.def("name", &isobarAmplitude::name)
		.def("printParameters", &isobarAmplitudeWrapper::printParameters__, &isobarAmplitudeWrapper::default_printParameters__)
		.def("printParameters", &isobarAmplitude_printParameters)

		.add_static_property("debugIsobarAmplitude", &rpwa::isobarAmplitude::debug, &rpwa::isobarAmplitude::setDebug);

	bp::register_ptr_to_python<rpwa::isobarAmplitudePtr>();

}
