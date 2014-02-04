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
			  bp::wrapper<rpwa::decayTopology>() { }

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
		}


		decayTopologyWrapper(const rpwa::decayTopology& topo)
			: rpwa::decayTopology(topo),
			  bp::wrapper<rpwa::decayTopology>() { }

		void clear() {
			if(bp::override clear = this->get_override("clear")) {
				clear();
			}
			rpwa::decayTopology::clear();
		}

		void default_clear() {
			rpwa::decayTopology::clear();
		}

	};

	const rpwa::interactionVertexPtr& decayGraph_toVertex(const rpwa::decayTopology& self,
	                                                      const rpwa::particlePtr& vertex)
	{
		return self.toVertex(vertex);
	}

	const rpwa::interactionVertexPtr& decayGraph_fromVertex(const rpwa::decayTopology& self,
		                                                    const rpwa::particlePtr& vertex)
	{
			return self.fromVertex(vertex);
	}

	bp::dict decayTopology_nmbIndistFsParticles(const rpwa::decayTopology& self) {
		bp::dict retval;
		std::map<std::string, unsigned int> nmbIndistFsParts = self.nmbIndistFsParticles();
		for(std::map<std::string, unsigned int>::const_iterator it = nmbIndistFsParts.begin(); it != nmbIndistFsParts.end(); ++it) {
			retval[it->first] = it->second;
		}
		return retval;
	}

	bp::list decayTopology_fsParticles(const rpwa::decayTopology& self) {
		return bp::list(self.fsParticles());
	}

	bp::list decayTopology_decayVertices(const rpwa::decayTopology& self) {
		return bp::list(self.decayVertices());
	}

	void decayTopology_transformFsParticles(rpwa::decayTopology& self, PyObject* pyL) {
		TLorentzRotation* L = rpwa::py::convertFromPy<TLorentzRotation*>(pyL);
		if(L == NULL) {
			printErr<<"Got invalid input when executing rpwa::diffractiveDissVertex::transformFsParticles()."<<std::endl;
		} else {
			self.transformFsParticles(*L);
		}
	}

	bool decayTopology_initKinematicsData(rpwa::decayTopology& self, PyObject* pyProdKinParticles, PyObject* pyDecayKinParticles) {
		TClonesArray* prodKinParticles = rpwa::py::convertFromPy<TClonesArray*>(pyProdKinParticles);
		TClonesArray* decayKinParticles = rpwa::py::convertFromPy<TClonesArray*>(pyDecayKinParticles);
		if((prodKinParticles == NULL) || (decayKinParticles == NULL)) {
			printErr<<"Got invalid input when executing rpwa::diffractiveDissVertex::initKinematicsData()."<<std::endl;
			return false;
		}
		return self.initKinematicsData(*prodKinParticles, *decayKinParticles);
	}

	bool decayTopology_readKinematicsData(rpwa::decayTopology& self, PyObject* pyProdKinMomenta, PyObject* pyDecayKinMomenta) {
		TClonesArray* prodKinMomenta = rpwa::py::convertFromPy<TClonesArray*>(pyProdKinMomenta);
		TClonesArray* decayKinMomenta = rpwa::py::convertFromPy<TClonesArray*>(pyDecayKinMomenta);
		if((prodKinMomenta == NULL) || (decayKinMomenta == NULL)) {
			printErr<<"Got invalid input when executing rpwa::diffractiveDissVertex::readKinematicsData()."<<std::endl;
			return false;
		}
		return self.readKinematicsData(*prodKinMomenta, *decayKinMomenta);
	}

	bool decayTopology_revertMomenta1(rpwa::decayTopology& self) {
		return self.revertMomenta();
	}

	bool decayTopology_revertMomenta2(rpwa::decayTopology& self, PyObject* pyIndexMap) {
		bp::list pyListIndexMap = bp::extract<bp::list>(pyIndexMap);
		std::vector<unsigned int> indexMap(bp::len(pyListIndexMap), 0);
		for(int i = 0; i < bp::len(pyListIndexMap); ++i) {
			indexMap[i] = bp::extract<unsigned int>(pyListIndexMap[i]);
		}
		return self.revertMomenta(indexMap);
	}

}

void rpwa::py::exportDecayTopology() {

	bp::class_<decayTopologyWrapper>("decayTopology")

		.def(bp::init<const rpwa::productionVertexPtr, const bp::object&, const bp::object&>())
		.def(bp::init<const rpwa::decayTopology&>())

		.def(bp::self_ns::str(bp::self))

		// Functions from decayGraph.hpp
		.def("toVertex", &decayGraph_toVertex, bp::return_value_policy<bp::copy_const_reference>())
		.def("fromVertex", &decayGraph_fromVertex, bp::return_value_policy<bp::copy_const_reference>())

		.def(
			"clone"
			, &rpwa::decayTopology::clone
			, (bp::arg("cloneFsParticles")=false,
			   bp::arg("cloneProdKinematics")=false)
		)

		.def("clear", &decayTopologyWrapper::clear, &decayTopologyWrapper::default_clear)
		.def("clear", &rpwa::decayTopology::clear)

		.def("nmbDecayVertices", &rpwa::decayTopology::nmbDecayVertices)
		.def("nmbFsParticles", &rpwa::decayTopology::nmbFsParticles)

		.def("nmbIndistFsParticles", &decayTopology_nmbIndistFsParticles)

		.def("fsParticlesIntrinsicParity", &rpwa::decayTopology::fsParticlesIntrinsicParity)
		.def("spaceInvEigenValue", &rpwa::decayTopology::spaceInvEigenValue)
		.def("reflectionEigenValue", &rpwa::decayTopology::reflectionEigenValue)

		.def("fsParticles", &decayTopology_fsParticles)
		.def("decayVertices", &decayTopology_decayVertices)

		.def(
			"XParticle"
			, &rpwa::decayTopology::XParticle
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def(
			"productionVertex"
			, &rpwa::decayTopology::productionVertex
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def(
			"XDecayVertex"
			, &rpwa::decayTopology::XDecayVertex
			, bp::return_value_policy<bp::copy_const_reference>()
		)

		.def("transformFsParticles", &decayTopology_transformFsParticles)

		.def("isProductionVertex", &rpwa::decayTopology::isProductionVertex)
		.def("isDecayVertex", &rpwa::decayTopology::isDecayVertex)
		.def("decayVertexIndex", &rpwa::decayTopology::decayVertexIndex)
		.def("isFsVertex", &rpwa::decayTopology::isFsVertex)
		.def("isFsParticle", &rpwa::decayTopology::isFsParticle)
		.def("fsParticlesIndex", &rpwa::decayTopology::fsParticlesIndex)

		.def("checkTopology", &rpwa::decayTopology::checkTopology)
		.def("checkConsistency", &rpwa::decayTopology::checkConsistency)

// This one is missing because it returns something defined in decayGraph.hpp, which is currently omitted.
/*		.def(
			"subDecay"
			, &rpwa::decayTopology::subDecay
			, (bp::arg("startNd"),
			   bp::arg("linkToMotherTopo")=false)
		)*/

		.def("addDecay", &rpwa::decayTopology::addDecay)
		.def("setProductionVertex", &rpwa::decayTopology::setProductionVertex)

		.def("initKinematicsData", &decayTopology_initKinematicsData)
		.def("readKinematicsData", &decayTopology_readKinematicsData)

		.def("fillKinematicsDataCache", &rpwa::decayTopology::fillKinematicsDataCache)

		.def("revertMomenta", &decayTopology_revertMomenta1)
		.def("revertMomenta", &decayTopology_revertMomenta2)

		.add_static_property("debugDecayTopology", &rpwa::decayTopology::debug, &rpwa::decayTopology::setDebug);

	bp::register_ptr_to_python<rpwa::decayTopologyPtr>();

}
