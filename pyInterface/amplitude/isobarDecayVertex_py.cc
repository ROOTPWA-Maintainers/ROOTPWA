#include "isobarDecayVertex_py.h"

#include<particle.h>

#include "rootConverters_py.h"

namespace bp = boost::python;

namespace {

	struct isobarDecayVertexWrapper : public rpwa::isobarDecayVertex,
	                                         bp::wrapper<rpwa::isobarDecayVertex>
	{

		isobarDecayVertexWrapper(const rpwa::particlePtr&       parent,
		                         const rpwa::particlePtr&       daughter1,
		                         const rpwa::particlePtr&       daughter2,
		                         const unsigned int       L       = 0,
		                         const unsigned int       S       = 0,
		                         const rpwa::massDependencePtr& massDep = rpwa::massDependencePtr())
			: rpwa::isobarDecayVertex(parent, daughter1, daughter2, L, S, massDep),
			  bp::wrapper<rpwa::isobarDecayVertex>() { }

		isobarDecayVertexWrapper(const isobarDecayVertex& vert)
			: rpwa::isobarDecayVertex(vert),
			  bp::wrapper<rpwa::isobarDecayVertex>() { }

		bool addInParticle(const rpwa::particlePtr& part) {
			if(bp::override addInParticle = this->get_override("addInParticle")) {
				return addInParticle(part);
			}
			return rpwa::isobarDecayVertex::addInParticle(part);
		}

		bool default_addInParticle(const rpwa::particlePtr& part) {
			return rpwa::isobarDecayVertex::addInParticle(part);
		}


		bool addOutParticle(const rpwa::particlePtr& part) {
			if(bp::override addOutParticle = this->get_override("addOutParticle")) {
				return addOutParticle(part);
			}
			return rpwa::isobarDecayVertex::addOutParticle(part);
		}

		bool default_addOutParticle(const rpwa::particlePtr& part) {
			return rpwa::isobarDecayVertex::addOutParticle(part);
		}

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::isobarDecayVertex::name();
		}

		std::string default_name() const {
			return rpwa::isobarDecayVertex::name();
		}

	};

	rpwa::particlePtr isobarDecayVertex_parent(rpwa::isobarDecayVertex& self) {
		return self.parent();
	}

	rpwa::particlePtr isobarDecayVertex_daughter1(rpwa::isobarDecayVertex& self) {
		return self.daughter1();
	}

	rpwa::particlePtr isobarDecayVertex_daughter2(rpwa::isobarDecayVertex& self) {
		return self.daughter2();
	}

	PyObject* isobarDecayVertex_calcParentLzVec(rpwa::isobarDecayVertex& self) {
		return rpwa::py::convertToPy<TLorentzVector>(self.calcParentLzVec());
	}

	const rpwa::massDependence& isobarDecayVertex_massDependence(const rpwa::isobarDecayVertex& self) {
		return *(self.massDependence());
	}

}

void rpwa::py::exportIsobarDecayVertex() {

	bp::class_<isobarDecayVertexWrapper, bp::bases<rpwa::interactionVertex> >("isobarDecayVertex", bp::no_init)

		.def(bp::init<rpwa::particlePtr, rpwa::particlePtr, rpwa::particlePtr, bp::optional<unsigned int, unsigned int, rpwa::massDependencePtr> >())

		.def(bp::self_ns::str(bp::self))

		.def(
			"clone"
			, &rpwa::isobarDecayVertex::clone
			, (bp::arg("cloneInParticles")=false,
			   bp::arg("cloneOutParticles")=false)
		)

		.def("addInParticle", &isobarDecayVertexWrapper::addInParticle, &isobarDecayVertexWrapper::default_addInParticle)
		.def("addInParticle", &rpwa::isobarDecayVertex::addInParticle)
		.def("addOutParticle", &isobarDecayVertexWrapper::addOutParticle, &isobarDecayVertexWrapper::default_addOutParticle)
		.def("addInParticle", &rpwa::isobarDecayVertex::addOutParticle)

		.def("parent", &isobarDecayVertex_parent)
		.def("daughter1", &isobarDecayVertex_daughter1)
		.def("daughter2", &isobarDecayVertex_daughter2)

		.def("calcParentLzVec", &isobarDecayVertex_calcParentLzVec)
		.def("calcParentCharge", &rpwa::isobarDecayVertex::calcParentCharge)
		.def("calcParentBaryonNmb", &rpwa::isobarDecayVertex::calcParentBaryonNmb)

		.add_property("L", &rpwa::isobarDecayVertex::L, &rpwa::isobarDecayVertex::setL)
		.add_property("S", &rpwa::isobarDecayVertex::S, &rpwa::isobarDecayVertex::setS)

		.def("massDepAmplitude", &rpwa::isobarDecayVertex::massDepAmplitude)
		.def("massDependence", &isobarDecayVertex_massDependence, bp::return_internal_reference<>())
		.def("setMassDependence", &rpwa::isobarDecayVertex::setMassDependence)

		.def("checkConsistency", &rpwa::isobarDecayVertex::checkConsistency)

		.def("name", &isobarDecayVertexWrapper::name, &isobarDecayVertexWrapper::default_name)
		.def("name", &rpwa::isobarDecayVertex::name)

		.add_static_property("debugIsobarDecayVertex", &rpwa::isobarDecayVertex::debug, &rpwa::isobarDecayVertex::setDebug);

	bp::register_ptr_to_python<rpwa::isobarDecayVertexPtr>();

}
