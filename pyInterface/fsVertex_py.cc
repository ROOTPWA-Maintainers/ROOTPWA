#include "fsVertex_py.h"

#include<particle.h>

namespace bp = boost::python;

namespace {

	struct fsVertexWrapper : public rpwa::fsVertex,
	                                bp::wrapper<rpwa::fsVertex>
	{
		fsVertexWrapper(const rpwa::particlePtr& fsParticle)
			: rpwa::fsVertex(fsParticle),
			  bp::wrapper<rpwa::fsVertex>() { };

		fsVertexWrapper(const fsVertex& vert)
			: rpwa::fsVertex(vert),
			  bp::wrapper<rpwa::fsVertex>() { };

		bool addInParticle(const rpwa::particlePtr& part) {
			if(bp::override addInParticle = this->get_override("addInParticle")) {
				return addInParticle(part);
			}
			return rpwa::fsVertex::addInParticle(part);
		};

		bool default_addInParticle(const rpwa::particlePtr& part) {
			return rpwa::fsVertex::addInParticle(part);
		};


		bool addOutParticle(const rpwa::particlePtr& part) {
			if(bp::override addOutParticle = this->get_override("addOutParticle")) {
				return addOutParticle(part);
			}
			return rpwa::fsVertex::addOutParticle(part);
		};

		bool default_addOutParticle(const rpwa::particlePtr& part) {
			return rpwa::fsVertex::addOutParticle(part);
		};

		rpwa::particlePtr fsParticle__() {
			return rpwa::fsVertex::fsParticle();
		};

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::fsVertex::name();
		};

		std::string default_name() const {
			return rpwa::fsVertex::name();
		};

	};

}

void rpwa::py::exportFsVertex() {

	bp::class_<fsVertexWrapper, bp::bases<interactionVertex> >("fsVertex", bp::no_init)

		.def(bp::init<rpwa::particlePtr>())
		.def(bp::init<rpwa::fsVertex&>())

		.def(bp::self_ns::str(bp::self))

		.def(
			"clone"
			, &fsVertexWrapper::clone
			, (bp::arg("cloneInParticles")=false,
			   bp::arg("cloneOutParticles")=false)
		)

		.def("addInParticle", &fsVertexWrapper::addInParticle, &fsVertexWrapper::default_addInParticle)
		.def("addOutParticle", &fsVertexWrapper::addOutParticle, &fsVertexWrapper::default_addOutParticle)

		.def("fsParticle", &fsVertexWrapper::fsParticle__)

		.def("name", &fsVertexWrapper::name, &fsVertexWrapper::default_name)

		.add_static_property("debugFsVertex", &fsVertexWrapper::debug, &fsVertexWrapper::setDebug);

	bp::register_ptr_to_python<rpwa::fsVertexPtr>();

};
