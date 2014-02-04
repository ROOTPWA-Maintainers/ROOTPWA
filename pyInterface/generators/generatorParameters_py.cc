
#include "generatorParameters_py.h"

#include "rootConverters_py.h"

namespace bp = boost::python;

namespace {

	void Target_setPosition(rpwa::Target& self, PyObject* pyPosition) {
		TVector3* position = rpwa::py::convertFromPy<TVector3*>(pyPosition);
		if(not position) {
			printErr<<"Got invalid input when executing rpwa::Target::setPosition()."<<std::endl;
			throw;
		}
		self.position = *position;
	}

	PyObject* Target_getPosition(const rpwa::Target& self) {
		return rpwa::py::convertToPy<TVector3>(self.position);
	}

	void FinalState_setParticles(rpwa::FinalState& self, const bp::object& pyParticles) {
		bp::list pyListParticles = bp::extract<bp::list>(pyParticles);
		std::vector<rpwa::particle> particles(bp::len(pyListParticles));
		for(int i = 0; i < bp::len(pyListParticles); ++i) {
			particles[i] = bp::extract<rpwa::particle>(pyListParticles[i]);
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
