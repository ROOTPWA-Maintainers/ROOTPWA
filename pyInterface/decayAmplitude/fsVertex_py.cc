#include "fsVertex_py.h"

#include<particle.h>

namespace bp = boost::python;

namespace {

	struct fsVertexWrapper : public rpwa::fsVertex,
	                                bp::wrapper<rpwa::fsVertex>
	{
		fsVertexWrapper(const rpwa::particlePtr& fsParticle)
			: rpwa::fsVertex(fsParticle),
			  bp::wrapper<rpwa::fsVertex>() { }

		fsVertexWrapper(const rpwa::fsVertex& vert)
			: rpwa::fsVertex(vert),
			  bp::wrapper<rpwa::fsVertex>() { }

		bool addInParticle(const rpwa::particlePtr& part) {
			if(bp::override addInParticle = this->get_override("addInParticle")) {
				return addInParticle(part);
			}
			return rpwa::fsVertex::addInParticle(part);
		}

		bool default_addInParticle(const rpwa::particlePtr& part) {
			return rpwa::fsVertex::addInParticle(part);
		}


		bool addOutParticle(const rpwa::particlePtr& part) {
			if(bp::override addOutParticle = this->get_override("addOutParticle")) {
				return addOutParticle(part);
			}
			return rpwa::fsVertex::addOutParticle(part);
		}

		bool default_addOutParticle(const rpwa::particlePtr& part) {
			return rpwa::fsVertex::addOutParticle(part);
		}

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::fsVertex::name();
		}

		std::string default_name() const {
			return rpwa::fsVertex::name();
		}

	};

	rpwa::particlePtr fsVertex_fsParticle(const rpwa::fsVertex& self) {
		return self.fsParticle();
	}

}

void rpwa::py::exportFsVertex() {

	bp::class_<fsVertexWrapper, bp::bases<interactionVertex> >("fsVertex", bp::no_init)

		.def(bp::init<rpwa::particlePtr>())
		.def(bp::init<const rpwa::fsVertex&>())

		.def(bp::self_ns::str(bp::self))

		.def(
			"clone"
			, &rpwa::fsVertex::clone
			, (bp::arg("cloneInParticles")=false,
			   bp::arg("cloneOutParticles")=false)
		)

		.def("addInParticle", &fsVertexWrapper::addInParticle, &fsVertexWrapper::default_addInParticle)
		.def("addInParticle", &rpwa::fsVertex::addInParticle)
		.def("addOutParticle", &fsVertexWrapper::addOutParticle, &fsVertexWrapper::default_addOutParticle)
		.def("addOutParticle", &rpwa::fsVertex::addOutParticle)

		.def("fsParticle", &fsVertex_fsParticle)

		.def("name", &fsVertexWrapper::name, &fsVertexWrapper::default_name)
		.def("name", &rpwa::fsVertex::name)

		.add_static_property("debugFsVertex", &rpwa::fsVertex::debug, &rpwa::fsVertex::setDebug);

	bp::register_ptr_to_python<rpwa::fsVertexPtr>();

}
