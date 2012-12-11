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
			  bp::wrapper<rpwa::isobarDecayVertex>() { };

		isobarDecayVertexWrapper(const isobarDecayVertex& vert)
			: rpwa::isobarDecayVertex(vert),
			  bp::wrapper<rpwa::isobarDecayVertex>() { };

		bool addInParticle(const rpwa::particlePtr& part) {
			if(bp::override addInParticle = this->get_override("addInParticle")) {
				return addInParticle(part);
			}
			return rpwa::isobarDecayVertex::addInParticle(part);
		};

		bool default_addInParticle(const rpwa::particlePtr& part) {
			return rpwa::isobarDecayVertex::addInParticle(part);
		};


		bool addOutParticle(const rpwa::particlePtr& part) {
			if(bp::override addOutParticle = this->get_override("addOutParticle")) {
				return addOutParticle(part);
			}
			return rpwa::isobarDecayVertex::addOutParticle(part);
		};

		bool default_addOutParticle(const rpwa::particlePtr& part) {
			return rpwa::isobarDecayVertex::addOutParticle(part);
		};

		rpwa::particlePtr& parent__() {
			return rpwa::isobarDecayVertex::parent();
		};

		rpwa::particlePtr& daughter1__() {
			return rpwa::isobarDecayVertex::daughter1();
		};

		rpwa::particlePtr& daughter2__() {
			return rpwa::isobarDecayVertex::daughter2();
		};

		PyObject* calcParentLzVec() {
			return rpwa::py::convertToPy<TLorentzVector>(rpwa::isobarDecayVertex::calcParentLzVec());
		};

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::isobarDecayVertex::name();
		};

		std::string default_name() const {
			return rpwa::isobarDecayVertex::name();
		};

	};

}

void rpwa::py::exportIsobarDecayVertex() {

	bp::class_<isobarDecayVertexWrapper, bp::bases<rpwa::interactionVertex> >("isobarDecayVertex", bp::no_init)

		.def(bp::init<rpwa::particlePtr, rpwa::particlePtr, rpwa::particlePtr, bp::optional<unsigned int, unsigned int, rpwa::massDependencePtr> >())

		.def(bp::self_ns::str(bp::self))

		.def(
			"clone"
			, &isobarDecayVertexWrapper::clone
			, (bp::arg("cloneInParticles")=false,
			   bp::arg("cloneOutParticles")=false)
		)

		.def("addInParticle", &isobarDecayVertexWrapper::addInParticle, &isobarDecayVertexWrapper::default_addInParticle)
		.def("addOutParticle", &isobarDecayVertexWrapper::addOutParticle, &isobarDecayVertexWrapper::default_addOutParticle)

		.def(
			"parent"
			, &isobarDecayVertexWrapper::parent__
			, bp::return_value_policy<bp::copy_non_const_reference>()
		)
		.def(
			"daughter1"
			, &isobarDecayVertexWrapper::daughter1__
			, bp::return_value_policy<bp::copy_non_const_reference>()
		)
		.def(
			"daughter2"
			, &isobarDecayVertexWrapper::daughter2__
			, bp::return_value_policy<bp::copy_non_const_reference>()
		)

		.def("calcParentLzVec", &isobarDecayVertexWrapper::calcParentLzVec)
		.def("calcParentCharge", &isobarDecayVertexWrapper::calcParentCharge)
		.def("calcParentBaryonNmb", &isobarDecayVertexWrapper::calcParentBaryonNmb)

		.add_property("L", &isobarDecayVertexWrapper::L, &isobarDecayVertexWrapper::setL)
		.add_property("S", &isobarDecayVertexWrapper::S, &isobarDecayVertexWrapper::setS)

		.def("massDepAmplitude", &isobarDecayVertexWrapper::massDepAmplitude)
		.def("massDependence", &isobarDecayVertexWrapper::massDependence, bp::return_value_policy<bp::return_by_value>())
		.def("setMassDependence", &isobarDecayVertexWrapper::setMassDependence)

		.def("checkConsistency", &isobarDecayVertexWrapper::checkConsistency)

		.def("name", &isobarDecayVertexWrapper::name, &isobarDecayVertexWrapper::default_name)

		.add_static_property("debugIsobarDecayVertex", &isobarDecayVertexWrapper::debug, &isobarDecayVertex::setDebug);

	bp::register_ptr_to_python<rpwa::isobarDecayVertexPtr>();

};
