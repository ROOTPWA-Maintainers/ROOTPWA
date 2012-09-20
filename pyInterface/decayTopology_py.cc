#include "decayTopology_py.h"

#include<TLorentzRotation.h>

#include<productionVertex.h>

#include "rootConverters_py.h"

namespace bp = boost::python;

namespace {

	struct decayTopologyWrapper : public rpwa::decayTopology,
	                                     bp::wrapper<rpwa::decayTopology>
	{

		decayTopologyWrapper()
			: rpwa::decayTopology(),
			  bp::wrapper<rpwa::decayTopology>() { };

		decayTopologyWrapper(const rpwa::productionVertexPtr& productionVertex,
		                     const bp::object&                pyDecayVertices,
		                     const bp::object&                pyFsParticles)
			: rpwa::decayTopology(),
			  bp::wrapper<rpwa::decayTopology>()
		{
			bp::list pyListDecayVertices = bp::extract<bp::list>(pyDecayVertices);
			bp::list pyListFsParticles = bp::extract<bp::list>(pyFsParticles);

			std::vector<rpwa::interactionVertexPtr> decayVertices(bp::len(pyListDecayVertices), rpwa::interactionVertexPtr());
			std::vector<rpwa::particlePtr> fsParticles(bp::len(pyListFsParticles), rpwa::particlePtr());

			for(int i = 0; i < bp::len(pyListDecayVertices); ++i) {
				decayVertices[i] = bp::extract<rpwa::interactionVertexPtr>(pyListDecayVertices[i]);
			}

			for(int i = 0; i < bp::len(pyListFsParticles); ++i) {
				fsParticles[i] = bp::extract<rpwa::particlePtr>(pyListFsParticles[i]);
			}

			rpwa::decayTopology::constructDecay(productionVertex, decayVertices, fsParticles);
		};


		decayTopologyWrapper(const rpwa::decayTopology& topo)
			: rpwa::decayTopology(topo),
			  bp::wrapper<rpwa::decayTopology>() { };

		void clear() {
			if(bp::override clear = this->get_override("clear")) {
				clear();
			}
			rpwa::decayTopology::clear();
		};

		void default_clear() {
			rpwa::decayTopology::clear();
		};

		bp::dict nmbIndistFsParticles__() const {
			return bp::dict(rpwa::decayTopology::nmbIndistFsParticles());
		};

		bp::list fsParticles__() const {
			return bp::list(rpwa::decayTopology::fsParticles());
		};

		bp::list decayVertices__() const {
			return bp::list(rpwa::decayTopology::decayVertices());
		};

		void transformFsParticles__(PyObject* pyL) {
			TLorentzRotation* L = rpwa::py::convertFromPy<TLorentzRotation*>(pyL);
			if(L == NULL) {
				printErr<<"Got invalid input when executing rpwa::diffractiveDissVertex::transformFsParticles()."<<std::endl;
			} else {
				rpwa::decayTopology::transformFsParticles(*L);
			}
		};

		bool initKinematicsData__(PyObject* pyProdKinParticles, PyObject* pyDecayKinParticles) {
			TClonesArray* prodKinParticles = rpwa::py::convertFromPy<TClonesArray*>(pyProdKinParticles);
			TClonesArray* decayKinParticles = rpwa::py::convertFromPy<TClonesArray*>(pyDecayKinParticles);
			if((prodKinParticles == NULL) || (decayKinParticles == NULL)) {
				printErr<<"Got invalid input when executing rpwa::diffractiveDissVertex::initKinematicsData()."<<std::endl;
				return false;
			}
			return rpwa::decayTopology::initKinematicsData(*prodKinParticles, *decayKinParticles);
		};

		bool readKinematicsData__(PyObject* pyProdKinMomenta, PyObject* pyDecayKinMomenta) {
			TClonesArray* prodKinMomenta = rpwa::py::convertFromPy<TClonesArray*>(pyProdKinMomenta);
			TClonesArray* decayKinMomenta = rpwa::py::convertFromPy<TClonesArray*>(pyDecayKinMomenta);
			if((prodKinMomenta == NULL) || (decayKinMomenta == NULL)) {
				printErr<<"Got invalid input when executing rpwa::diffractiveDissVertex::readKinematicsData()."<<std::endl;
				return false;
			}
			return rpwa::decayTopology::readKinematicsData(*prodKinMomenta, *decayKinMomenta);
		};

		bool revertMomenta__1() {
			return rpwa::decayTopology::revertMomenta();
		};

		bool revertMomenta__2(PyObject* pyIndexMap) {
			bp::list pyListIndexMap = bp::extract<bp::list>(pyIndexMap);
			std::vector<unsigned int> indexMap(bp::len(pyListIndexMap), 0);
			for(int i = 0; i < bp::len(pyListIndexMap); ++i) {
				indexMap[i] = bp::extract<unsigned int>(pyListIndexMap[i]);
			}
			return rpwa::decayTopology::revertMomenta(indexMap);
		};

	};

}

void rpwa::py::exportDecayTopology() {

	bp::class_<decayTopologyWrapper>("decayTopology")
		
		.def(bp::init<rpwa::productionVertexPtr, bp::object&, bp::object&>())
		.def(bp::init<rpwa::decayTopology>())

		.def(bp::self_ns::str(bp::self))

		.def(
			"clone"
			, &decayTopologyWrapper::clone
			, (bp::arg("cloneFsParticles")=false,
			   bp::arg("cloneProdKinematics")=false)
		)

		.def("clear", &decayTopologyWrapper::clear, &decayTopologyWrapper::default_clear)

		.def("nmbDecayVertices", &decayTopologyWrapper::nmbDecayVertices)
		.def("nmbFsParticles", &decayTopologyWrapper::nmbFsParticles)

		.def("nmbIndistFsParticles", &decayTopologyWrapper::nmbIndistFsParticles__)

		.def("fsParticlesIntrinsicParity", &decayTopologyWrapper::fsParticlesIntrinsicParity)
		.def("spaceInvEigenValue", &decayTopologyWrapper::spaceInvEigenValue)
		.def("reflectionEigenValue", &decayTopologyWrapper::reflectionEigenValue)

		.def("fsParticles", &decayTopologyWrapper::fsParticles__)
		.def("decayVertices", &decayTopologyWrapper::decayVertices__)

		.def(
			"XParticle"
			, &decayTopologyWrapper::XParticle
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def(
			"productionVertex"
			, &decayTopologyWrapper::productionVertex
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def(
			"XDecayVertex"
			, &decayTopologyWrapper::XDecayVertex
			, bp::return_value_policy<bp::copy_const_reference>()
		)

		.def("transformFsParticles", &decayTopologyWrapper::transformFsParticles__)

		.def("isProductionVertex", &decayTopologyWrapper::isProductionVertex)
		.def("isDecayVertex", &decayTopologyWrapper::isDecayVertex)
		.def("isFsVertex", &decayTopologyWrapper::isFsVertex)
		.def("isFsParticle", &decayTopologyWrapper::isFsParticle)
		.def("fsParticlesIndex", &decayTopologyWrapper::fsParticlesIndex)

		.def("checkTopology", &decayTopologyWrapper::checkTopology)
		.def("checkConsistency", &decayTopologyWrapper::checkConsistency)

// This one is missing because it returns something defined in decayGraph.hpp, which is currently omitted.
/*		.def(
			"subDecay"
			, &decayTopologyWrapper::subDecay
			, (bp::arg("startNd"),
			   bp::arg("linkToMotherTopo")=false)
		)*/

		.def("addDecay", &decayTopologyWrapper::addDecay)
		.def("setProductionVertex", &decayTopologyWrapper::setProductionVertex)

		.def("initKinematicsData", &decayTopologyWrapper::initKinematicsData__)
		.def("readKinematicsData", &decayTopologyWrapper::readKinematicsData__)

		.def("fillKinematicsDataCache", &decayTopologyWrapper::fillKinematicsDataCache)

		.def("revertMomenta", &decayTopologyWrapper::revertMomenta__1)
		.def("revertMomenta", &decayTopologyWrapper::revertMomenta__2)

		.add_static_property("debugDecayTopology", &decayTopologyWrapper::debug, &decayTopologyWrapper::setDebug);

	bp::register_ptr_to_python<rpwa::decayTopologyPtr>();

};
