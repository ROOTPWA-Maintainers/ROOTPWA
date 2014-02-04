#include "interactionVertex_py.h"

#include<TLorentzRotation.h>

#include "rootConverters_py.h"
#include<particle.h>

namespace bp = boost::python;

namespace {

	struct interactionVertexWrapper : public rpwa::interactionVertex,
	                                         bp::wrapper<rpwa::interactionVertex>
	{

		interactionVertexWrapper()
			: rpwa::interactionVertex(),
			  bp::wrapper<rpwa::interactionVertex>() { }

		interactionVertexWrapper(const rpwa::interactionVertex& vert)
			: rpwa::interactionVertex(vert),
			  bp::wrapper<rpwa::interactionVertex>() { }

		void clear() {
			if(bp::override clear = this->get_override("clear")) {
				clear();
				return;
			}
			return rpwa::interactionVertex::clear();
		}

		void default_clear() {
			return rpwa::interactionVertex::clear();
		}

		bool addInParticle(const rpwa::particlePtr& part) {
			if(bp::override addInParticle = this->get_override("addInParticle")) {
				return addInParticle(part);
			}
			return rpwa::interactionVertex::addInParticle(part);
		}

		bool default_addInParticle(const rpwa::particlePtr& part) {
			return rpwa::interactionVertex::addInParticle(part);
		}

		bool addOutParticle(const rpwa::particlePtr& part) {
			if(bp::override addOutParticle = this->get_override("addOutParticle")) {
				return addOutParticle(part);
			}
			return rpwa::interactionVertex::addOutParticle(part);
		}

		bool default_addOutParticle(const rpwa::particlePtr& part) {
			return rpwa::interactionVertex::addOutParticle(part);
		}

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::interactionVertex::name();
		}

		std::string default_name() const {
			return rpwa::interactionVertex::name();
		}

	};

	void interactionVertex_transformOutParticles(rpwa::interactionVertex& self, PyObject* pyRotation) {
		TLorentzRotation* L = rpwa::py::convertFromPy<TLorentzRotation*>(pyRotation);
		if(L == NULL) {
			printErr<<"Got invalid input when executing rpwa::interactionVertex()."<<std::endl;
		} else {
			return self.transformOutParticles(*L);
		}
	}

	bp::list interactionVertex_inParticles(const rpwa::interactionVertex& self) {
		std::vector<rpwa::particlePtr> retVec = self.inParticles();
		bp::object iter = bp::iterator<std::vector<rpwa::particlePtr> >()(retVec);
		return bp::list(iter);
	}

	bp::list interactionVertex_outParticles(const rpwa::interactionVertex& self) {
		std::vector<rpwa::particlePtr> retVec = self.outParticles();
		bp::object iter = bp::iterator<std::vector<rpwa::particlePtr> >()(retVec);
		return bp::list(iter);
	}

}

void rpwa::py::exportInteractionVertex() {

	bp::class_<interactionVertexWrapper>("interactionVertex")

		.def(bp::init<interactionVertexWrapper&>())

		.def(bp::self_ns::str(bp::self))

		.def(
			"clone"
			, &rpwa::interactionVertex::clone
			, (bp::arg("cloneInParticles")=false,
			   bp::arg("cloneOutParticles")=false)
		)

		.def("clear", &interactionVertexWrapper::clear, &interactionVertexWrapper::default_clear)
		.def("clear", &rpwa::interactionVertex::clear)

		.def("addInParticle", &interactionVertexWrapper::addInParticle, &interactionVertexWrapper::default_addInParticle)
		.def("addInParticle", &rpwa::interactionVertex::addInParticle)
		.def("addOutParticle", &interactionVertexWrapper::addOutParticle, &interactionVertexWrapper::default_addOutParticle)
		.def("addOutParticle", &rpwa::interactionVertex::addOutParticle)

		.def("transformOutParticles", &interactionVertex_transformOutParticles)

		.def("inParticles", &interactionVertex_inParticles)
		.def("outParticles", &interactionVertex_outParticles)

		.add_property("nmbInParticles", &rpwa::interactionVertex::nmbInParticles)
		.add_property("nmbOutParticles", &rpwa::interactionVertex::nmbOutParticles)

		.def("name", &interactionVertexWrapper::name, &interactionVertexWrapper::default_name)
		.def("name", &rpwa::interactionVertex::name)

		.add_static_property("debugInteractionVertex", &rpwa::interactionVertex::debug, &rpwa::interactionVertex::setDebug);

	bp::register_ptr_to_python<rpwa::interactionVertexPtr>();

}
