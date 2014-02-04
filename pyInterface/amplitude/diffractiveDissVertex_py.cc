#include "diffractiveDissVertex_py.h"

#include<TClonesArray.h>

#include "rootConverters_py.h"

namespace bp = boost::python;

namespace {

	struct diffractiveDissVertexWrapper : public rpwa::diffractiveDissVertex,
	                                             bp::wrapper<rpwa::diffractiveDissVertex>
	{

		diffractiveDissVertexWrapper(const rpwa::particlePtr& beam,
		                             const rpwa::particlePtr& target,
		                             const rpwa::particlePtr& XParticle,
		                             const rpwa::particlePtr& recoil = rpwa::particlePtr())
			: rpwa::diffractiveDissVertex(beam, target, XParticle, recoil),
			  bp::wrapper<rpwa::diffractiveDissVertex>() { }

		diffractiveDissVertexWrapper(const rpwa::diffractiveDissVertex& vert)
			: rpwa::diffractiveDissVertex(vert),
			  bp::wrapper<rpwa::diffractiveDissVertex>() { }

		bool addInParticle(const rpwa::particlePtr& part) {
			if(bp::override addInParticle = this->get_override("addInParticle")) {
				return addInParticle(part);
			}
			return rpwa::diffractiveDissVertex::addInParticle(part);
		}

		bool default_addInParticle(const rpwa::particlePtr& part) {
			return rpwa::diffractiveDissVertex::addInParticle(part);
		}

		bool addOutParticle(const rpwa::particlePtr& part) {
			if(bp::override addOutParticle = this->get_override("addOutParticle")) {
				return addOutParticle(part);
			}
			return rpwa::diffractiveDissVertex::addOutParticle(part);
		}

		bool default_addOutParticle(const rpwa::particlePtr& part) {
			return rpwa::diffractiveDissVertex::addOutParticle(part);
		}

		PyObject* referenceLzVec__() const {
			if(bp::override referenceLzVec = this->get_override("referenceLzVec")) {
				return rpwa::py::convertToPy<TLorentzVector>(referenceLzVec());
			}
			return rpwa::py::convertToPy<TLorentzVector>(rpwa::diffractiveDissVertex::referenceLzVec());
		}

		PyObject* default_referenceLzVec__() const {
			return rpwa::py::convertToPy<TLorentzVector>(rpwa::diffractiveDissVertex::referenceLzVec());
		}

		const rpwa::particlePtr& XParticle() const {
			if(bp::override XParticle = this->get_override("XParticle")) {
				return XParticle();
			}
			return rpwa::diffractiveDissVertex::XParticle();
		}

		const rpwa::particlePtr& default_XParticle() const {
			return rpwa::diffractiveDissVertex::XParticle();
		}

		void setXFlavorQN() {
			if(bp::override setXFlavorQN = this->get_override("XFlavorQN")) {
				setXFlavorQN();
			}
			rpwa::diffractiveDissVertex::setXFlavorQN();
		}

		void default_setXFlavorQN() {
			rpwa::diffractiveDissVertex::setXFlavorQN();
		}

		bool initKinematicsData__(PyObject* pyProdKinPartNames) {
			TClonesArray* prodKinPartNames = rpwa::py::convertFromPy<TClonesArray*>(pyProdKinPartNames);
			if(prodKinPartNames == NULL) {
				printErr<<"Got invalid input when executing rpwa::diffractiveDissVertex::initKinematicsData()."<<std::endl;
				return false;
			}
			if(bp::override initKinematicsData = this->get_override("initKinematicsData")) {
				return initKinematicsData(*prodKinPartNames);
			}
			return rpwa::diffractiveDissVertex::initKinematicsData(*prodKinPartNames);
		}

		bool default_initKinematicsData__(PyObject* pyProdKinPartNames) {
			TClonesArray* prodKinPartNames = rpwa::py::convertFromPy<TClonesArray*>(pyProdKinPartNames);
			if(prodKinPartNames == NULL) {
				printErr<<"Got invalid input when executing rpwa::diffractiveDissVertex::initKinematicsData()."<<std::endl;
				return false;
			}
			return rpwa::diffractiveDissVertex::initKinematicsData(*prodKinPartNames);
		}

		bool readKinematicsData__(PyObject* pyProdKinPartNames) {
			TClonesArray* prodKinMomenta = rpwa::py::convertFromPy<TClonesArray*>(pyProdKinPartNames);
			if(prodKinMomenta == NULL) {
				printErr<<"Got invalid input when executing rpwa::diffractiveDissVertex::readKinematicsData()."<<std::endl;
				return false;
			}
			if(bp::override readKinematicsData = this->get_override("readKinematicsData")) {
				return readKinematicsData(*prodKinMomenta);
			}
			return rpwa::diffractiveDissVertex::readKinematicsData(*prodKinMomenta);
		}

		bool default_readKinematicsData__(PyObject* pyProdKinPartNames) {
			TClonesArray* prodKinMomenta = rpwa::py::convertFromPy<TClonesArray*>(pyProdKinPartNames);
			if(prodKinMomenta == NULL) {
				printErr<<"Got invalid input when executing rpwa::diffractiveDissVertex::readKinematicsData()."<<std::endl;
				return false;
			}
			return rpwa::diffractiveDissVertex::readKinematicsData(*prodKinMomenta);
		}

		bool revertMomenta() {
			if(bp::override revertMomenta = this->get_override("revertMomenta")) {
				return revertMomenta();
			}
			return rpwa::diffractiveDissVertex::revertMomenta();
		}

		bool default_revertMomenta() {
			return rpwa::diffractiveDissVertex::revertMomenta();
		}


		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::diffractiveDissVertex::name();
		}

		std::string default_name() const {
			return rpwa::diffractiveDissVertex::name();
		}

	};

	PyObject* diffractiveDissVertex_referenceLzVec(const rpwa::diffractiveDissVertex& self) {
		return rpwa::py::convertToPy<TLorentzVector>(self.referenceLzVec());
	}

	bool diffractiveDissVertex_initKinematicsData(rpwa::diffractiveDissVertex& self, PyObject* pyProdKinPartNames) {
		TClonesArray* prodKinPartNames = rpwa::py::convertFromPy<TClonesArray*>(pyProdKinPartNames);
		if(prodKinPartNames == NULL) {
			printErr<<"Got invalid input when executing rpwa::diffractiveDissVertex::initKinematicsData()."<<std::endl;
			return false;
		}
		return self.initKinematicsData(*prodKinPartNames);
	}

	bool diffractiveDissVertex_readKinematicsData(rpwa::diffractiveDissVertex& self, PyObject* pyProdKinPartNames) {
		TClonesArray* prodKinMomenta = rpwa::py::convertFromPy<TClonesArray*>(pyProdKinPartNames);
		if(prodKinMomenta == NULL) {
			printErr<<"Got invalid input when executing rpwa::diffractiveDissVertex::readKinematicsData()."<<std::endl;
			return false;
		}
		return self.readKinematicsData(*prodKinMomenta);
	}

}

void rpwa::py::exportDiffractiveDissVertex() {

	bp::class_<diffractiveDissVertexWrapper, bp::bases<rpwa::productionVertex> >("diffractiveDissVertex", bp::no_init)

		.def(bp::init<rpwa::particlePtr, rpwa::particlePtr, rpwa::particlePtr, bp::optional<rpwa::particlePtr> >())
		.def(bp::init<rpwa::diffractiveDissVertex&>())

		.def(bp::self_ns::str(bp::self))

		.def(
			"clone"
			, &rpwa::diffractiveDissVertex::clone
			, (bp::arg("cloneInParticles")=false,
			   bp::arg("cloneOutParticles")=false)
		)

		.def("addInParticle", &diffractiveDissVertexWrapper::addInParticle, &diffractiveDissVertexWrapper::default_addInParticle)
		.def("addInParticle", &rpwa::diffractiveDissVertex::addInParticle)
		.def("addOutParticle", &diffractiveDissVertexWrapper::addOutParticle, &diffractiveDissVertexWrapper::default_addOutParticle)
		.def("addOutParticle", &rpwa::diffractiveDissVertex::addOutParticle)

		.def("referenceLzVec", &diffractiveDissVertexWrapper::referenceLzVec__, &diffractiveDissVertexWrapper::default_referenceLzVec__)
		.def("referenceLzVec", &diffractiveDissVertex_referenceLzVec)
		.def(
			"XParticle"
			, &diffractiveDissVertexWrapper::XParticle
			, &diffractiveDissVertexWrapper::default_XParticle
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def(
			"XParticle"
			, &rpwa::diffractiveDissVertex::XParticle
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def("setXFlavorQN", &diffractiveDissVertexWrapper::setXFlavorQN, &diffractiveDissVertexWrapper::default_setXFlavorQN)
		.def("setXFlavorQN", &rpwa::diffractiveDissVertex::setXFlavorQN)

		.def(
			"beam"
			, &rpwa::diffractiveDissVertex::beam
			, bp::return_value_policy<bp::copy_const_reference>()
		)

		.def(
			"target"
			, &rpwa::diffractiveDissVertex::target
			, bp::return_value_policy<bp::copy_const_reference>()
		)

		.def(
			"recoil"
			, &rpwa::diffractiveDissVertex::recoil
			, bp::return_value_policy<bp::copy_const_reference>()
		)

		.def("initKinematicsData", &diffractiveDissVertexWrapper::initKinematicsData__, &diffractiveDissVertexWrapper::default_initKinematicsData__)
		.def("initKinematicsData", &diffractiveDissVertex_initKinematicsData)
		.def("readKinematicsData", &diffractiveDissVertexWrapper::readKinematicsData__, &diffractiveDissVertexWrapper::default_readKinematicsData__)
		.def("readKinematicsData", &diffractiveDissVertex_readKinematicsData)

		.def("revertMomenta", &diffractiveDissVertexWrapper::revertMomenta, &diffractiveDissVertexWrapper::default_revertMomenta)
		.def("revertMomenta", &rpwa::diffractiveDissVertex::revertMomenta)
		.def("name", &diffractiveDissVertexWrapper::name, &diffractiveDissVertexWrapper::default_name)
		.def("name", &rpwa::diffractiveDissVertex::name)

		.add_static_property("debugDiffractiveDissVertex", &rpwa::diffractiveDissVertex::debug, &rpwa::diffractiveDissVertex::setDebug);

	bp::register_ptr_to_python<rpwa::diffractiveDissVertexPtr>();

}
