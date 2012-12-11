#include "isobarAmplitude_py.h"

#include "rootConverters_py.h"

namespace bp = boost::python;

namespace {

	struct isobarAmplitudeWrapper : public rpwa::isobarAmplitude,
	                                       bp::wrapper<rpwa::isobarAmplitude>
	{

		static PyObject* gjTransform__(PyObject* pyBeamLv, PyObject* pyXLv) {
			TLorentzVector* beamLv = rpwa::py::convertFromPy<TLorentzVector*>(pyBeamLv);
			TLorentzVector* XLv = rpwa::py::convertFromPy<TLorentzVector*>(pyXLv);
			if((beamLv == NULL) || (XLv == NULL)) {
				printErr<<"Got invalid input when executing rpwa::isobarAmplitude::gjTransform()."<<std::endl;
				return bp::object().ptr();
			}
			return rpwa::py::convertToPy<TLorentzRotation>(rpwa::isobarAmplitude::gjTransform(*beamLv, *XLv));
		};

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::isobarAmplitude::name();
		};

		std::string default_name() const {
			return rpwa::isobarAmplitude::name();
		};

		std::string printParameters__() const {
			std::stringstream sstr;
			if(bp::override printParameters = this->get_override("printParameters")) {
				printParameters(sstr);
			} else {
				rpwa::isobarAmplitude::printParameters(sstr);
			}
			return sstr.str();
		};

		std::string default_printParameters__() const {
			std::stringstream sstr;
			rpwa::isobarAmplitude::printParameters(sstr);
			return sstr.str();
		};

	};

}

void rpwa::py::exportIsobarAmplitude() {

	bp::class_<isobarAmplitudeWrapper, boost::noncopyable>("isobarAmplitude", bp::no_init)

		.def(bp::self_ns::str(bp::self))

		.def(
			"decayTopology"
			, &isobarAmplitudeWrapper::decayTopology
			, bp::return_value_policy<bp::copy_const_reference>()
		)

		.def("setDecayTopology", &isobarAmplitudeWrapper::setDecayTopology)

		.add_property("reflectivityBasis", &isobarAmplitudeWrapper::reflectivityBasis, &isobarAmplitudeWrapper::enableReflectivityBasis)
		.add_property("boseSymmetrization", &isobarAmplitudeWrapper::boseSymmetrization, &isobarAmplitudeWrapper::enableBoseSymmetrization)
		.add_property("doSpaceInversion", &isobarAmplitudeWrapper::doSpaceInversion, &isobarAmplitudeWrapper::enableSpaceInversion)
		.add_property("doReflection", &isobarAmplitudeWrapper::doReflection, &isobarAmplitudeWrapper::enableReflection)

		.def("gjTransform", &isobarAmplitudeWrapper::gjTransform__)

		.def("amplitude", &isobarAmplitudeWrapper::amplitude)

		.def("__call__", &isobarAmplitudeWrapper::operator())

		.def("name", &isobarAmplitudeWrapper::name, &isobarAmplitudeWrapper::default_name)
		.def("printParameters__", &isobarAmplitudeWrapper::printParameters__, &isobarAmplitudeWrapper::default_printParameters__)

		.add_static_property("debugIsobarAmplitude", &isobarAmplitudeWrapper::debug, &isobarAmplitudeWrapper::setDebug);

	bp::register_ptr_to_python<rpwa::isobarAmplitudePtr>();

};

