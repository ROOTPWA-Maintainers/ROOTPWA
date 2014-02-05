
#include "generatorParameters_py.h"

#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;

namespace {

	void Target_setPosition(rpwa::Target& self, PyObject* pyPosition) {
		TVector3* position = rpwa::py::convertFromPy<TVector3*>(pyPosition);
		if(not position) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for position when executing rpwa::Target::setPosition()");
			bp::throw_error_already_set();
		}
		self.position = *position;
	}

	PyObject* Target_getPosition(const rpwa::Target& self) {
		return rpwa::py::convertToPy<TVector3>(self.position);
	}

	void FinalState_setParticles(rpwa::FinalState& self, const bp::object& pyParticles) {
		std::vector<rpwa::particle> particles;
		if(not rpwa::py::convertBPObjectToVector<rpwa::particle>(pyParticles, particles)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for particles when executing rpwa::FinalState::setParticles()");
			bp::throw_error_already_set();
		}
		self.particles = particles;
	}

	bp::list FinalState_getParticles(const rpwa::FinalState& self) {
		return bp::list(self.particles);
	}

}

void rpwa::py::exportGeneratorParameters() {

	bp::class_<rpwa::Beam>("Beam")
		.def(bp::self_ns::str(bp::self))
		.def_readwrite("particle", &Beam::particle)
		.def_readwrite("momentum", &Beam::momentum)
		.def_readwrite("momentumSigma", &Beam::momentumSigma)
		.def_readwrite("DxDz", &Beam::DxDz)
		.def_readwrite("DxDzSigma", &Beam::DxDzSigma)
		.def_readwrite("DyDz", &Beam::DyDz)
		.def_readwrite("DyDzSigma", &Beam::DyDzSigma);

	bp::class_<rpwa::Target>("Target")
		.def(bp::self_ns::str(bp::self))
		.def_readwrite("targetParticle", &Target::targetParticle)
		.def_readwrite("recoilParticle", &Target::recoilParticle)
		.add_property("position", &Target_getPosition, &Target_setPosition)
		.def_readwrite("length", &Target::length)
		.def_readwrite("radius", &Target::radius)
		.def_readwrite("interactionLength", &Target::interactionLength);

	bp::class_<rpwa::FinalState>("FinalState")
		.def(bp::self_ns::str(bp::self))
		.add_property("particles", &FinalState_getParticles, &FinalState_setParticles);

}
