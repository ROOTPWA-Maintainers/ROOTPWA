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
			  bp::wrapper<rpwa::interactionVertex>() { };

		interactionVertexWrapper(const rpwa::interactionVertex& vert)
			: rpwa::interactionVertex(vert),
			  bp::wrapper<rpwa::interactionVertex>() { };

		void clear() {
			if(bp::override clear = this->get_override("clear")) {
				clear();
				return;
			}
			return rpwa::interactionVertex::clear();
		};

		void default_clear() {
			return rpwa::interactionVertex::clear();
		};

		bool addInParticle(const rpwa::particlePtr& part) {
			if(bp::override addInParticle = this->get_override("addInParticle")) {
				return addInParticle(part);
			}
			return rpwa::interactionVertex::addInParticle(part);
		};

		bool default_addInParticle(const rpwa::particlePtr& part) {
			return rpwa::interactionVertex::addInParticle(part);
		};

		bool addOutParticle(const rpwa::particlePtr& part) {
			if(bp::override addOutParticle = this->get_override("addOutParticle")) {
				return addOutParticle(part);
			}
			return rpwa::interactionVertex::addOutParticle(part);
		};

		bool default_addOutParticle(const rpwa::particlePtr& part) {
			return rpwa::interactionVertex::addOutParticle(part);
		};

		void transformOutParticles(PyObject* pyRotation) {
			TLorentzRotation* L = rpwa::py::convertFromPy<TLorentzRotation*>(pyRotation);
			if(L == NULL) {
				printErr<<"Got invalid input when executing rpwa::interactionVertex()."<<std::endl;
			} else {
				return rpwa::interactionVertex::transformOutParticles(*L);
			}
		};

		bp::list inParticles() {
			std::vector<rpwa::particlePtr> retVec = rpwa::interactionVertex::inParticles();
			bp::object iter = bp::iterator<std::vector<rpwa::particlePtr> >()(retVec);
			return bp::list(iter);
		};

		bp::list outParticles() {
			std::vector<rpwa::particlePtr> retVec = rpwa::interactionVertex::outParticles();
			bp::object iter = bp::iterator<std::vector<rpwa::particlePtr> >()(retVec);
			return bp::list(iter);
		};

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::interactionVertex::name();
		};

		std::string default_name() const {
			return rpwa::interactionVertex::name();
		};

	};

}

void rpwa::py::exportInteractionVertex() {

	bp::class_<interactionVertexWrapper>("interactionVertex")

		.def(bp::init<interactionVertexWrapper&>())

		.def(bp::self_ns::str(bp::self))

		.def(
			"clone"
			, &interactionVertexWrapper::clone
			, (bp::arg("cloneInParticles")=false,
			   bp::arg("cloneOutParticles")=false)
		)

		.def("clear", &interactionVertexWrapper::clear, &interactionVertexWrapper::default_clear)

		.def("addInParticle", &interactionVertexWrapper::addInParticle, &interactionVertexWrapper::default_addInParticle)
		.def("addOutParticle", &interactionVertexWrapper::addOutParticle, &interactionVertexWrapper::default_addOutParticle)

		.def("transformOutParticles", &interactionVertexWrapper::transformOutParticles)

		.def("inParticles", &interactionVertexWrapper::inParticles)
		.def("outParticles", &interactionVertexWrapper::outParticles)

		.add_property("nmbInParticles", &interactionVertexWrapper::nmbInParticles)
		.add_property("nmbOutParticles", &interactionVertexWrapper::nmbOutParticles)

		.def("name", &interactionVertexWrapper::name, &interactionVertexWrapper::default_name)

		.add_static_property("debugInteractionVertex", &interactionVertexWrapper::debug, &interactionVertexWrapper::setDebug);

	bp::register_ptr_to_python<rpwa::interactionVertexPtr>();

};
