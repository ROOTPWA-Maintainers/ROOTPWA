#include "productionVertex_py.h"

#include<TClonesArray.h>

namespace bp = boost::python;

namespace {

	struct productionVertexWrapper : public rpwa::productionVertex,
                                            bp::wrapper<rpwa::productionVertex>
	{

		productionVertexWrapper()
			: rpwa::productionVertex(),
			  bp::wrapper<rpwa::productionVertex>() { };

		const TLorentzVector& referenceLzVec() const {
			return this->get_override("referenceLzVec")();
		};

		const rpwa::particlePtr& XParticle() const {
			return this->get_override("XParticle")();
		};

		std::complex<double> productionAmp() const {
			if(bp::override productionAmp = this->get_override("productionAmp")) {
				return productionAmp();
			}
			return rpwa::productionVertex::productionAmp();
		};

		std::complex<double> default_productionAmp() const {
			return rpwa::productionVertex::productionAmp();
		};

		void setXFlavorQN() {
			this->get_override("setXFlavorQN")();
		};

		bool initKinematicsData(const TClonesArray& names) {
			return this->get_override("initKinematicsData")(names);
		};

		bool readKinematicsData(const TClonesArray& momenta) {
			return this->get_override("readKinematicsData")(momenta);
		};

		bool revertMomenta() {
			return this->get_override("revertKinematicsData")();
		};

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::productionVertex::name();
		};

		std::string default_name() const {
			return rpwa::productionVertex::name();
		};

	};

}

void rpwa::py::exportProductionVertex() {

	bp::class_<productionVertexWrapper, bp::bases<rpwa::interactionVertex>, boost::noncopyable>("productionVertex", bp::no_init)

		.def(bp::self_ns::str(bp::self))

		.def(
			"referenceLzVec"
			, bp::pure_virtual(&rpwa::productionVertex::referenceLzVec)
			, bp::return_value_policy<bp::copy_const_reference>()
		)

		.def(
			"XParticle"
			, bp::pure_virtual(&rpwa::productionVertex::XParticle)
			, bp::return_value_policy<bp::copy_const_reference>()
		)

		.def("productionAmp", &productionVertexWrapper::productionAmp, &productionVertexWrapper::default_productionAmp)
		.def("productionAmp", &rpwa::productionVertex::productionAmp)
		.def("setXFlavorQN", bp::pure_virtual(&rpwa::productionVertex::setXFlavorQN))

		.def("initKinematicsData", bp::pure_virtual(&rpwa::productionVertex::initKinematicsData))
		.def("readKinematicsData", bp::pure_virtual(&rpwa::productionVertex::readKinematicsData))

		.def("revertMomenta", bp::pure_virtual(&rpwa::productionVertex::revertMomenta))

		.def("name", &productionVertexWrapper::name, &productionVertexWrapper::default_name)
		.def("name", &rpwa::productionVertex::name)

		.add_static_property("debugProductionVertex", &rpwa::productionVertex::debug, &rpwa::productionVertex::setDebug);

	bp::register_ptr_to_python<rpwa::productionVertexPtr>();

};
