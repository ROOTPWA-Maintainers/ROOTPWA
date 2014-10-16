#include "decayTopology_py.h"

#include<TLorentzRotation.h>
#include<TVector3.h>

#include<productionVertex.h>

#include "rootConverters_py.h"
#include "stlContainers_py.h"

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
			std::vector<rpwa::interactionVertexPtr> decayVertices;
			if(not rpwa::py::convertBPObjectToVector<rpwa::interactionVertexPtr>(pyDecayVertices, decayVertices)) {
				PyErr_SetString(PyExc_TypeError, "Got invalid input for decayVertices when executing rpwa::decayTopology::decayTopology()");
				bp::throw_error_already_set();
			}
			std::vector<rpwa::particlePtr> fsParticles;
			if(not rpwa::py::convertBPObjectToVector<rpwa::particlePtr>(pyFsParticles, fsParticles)) {
				printErr<<"Got invalid input for fsParticles when executing "
				        <<"rpwa::decayTopology::decayTopology(). Aborting..."<<std::endl;
				PyErr_SetString(PyExc_TypeError, "Got invalid input for fsParticles when executing rpwa::decayTopology::decayTopology()");
				bp::throw_error_already_set();
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
		if((not prodKinParticles) or (not decayKinParticles)) {
			bp::extract<bp::list> getProdKinParticlesList(pyProdKinParticles);
			if(not getProdKinParticlesList.check()) {
				printErr<<"Got invalid input for prodKinParticles when executing rpwa::decayTopology::initKinematicsData()"<<std::endl;
				return false;
			}
			bp::list prodKinParticlesList = getProdKinParticlesList();
			bp::extract<bp::list> getDecayKinParticlesList(pyDecayKinParticles);
			if(not getProdKinParticlesList.check()) {
				printErr<<"Got invalid input for prodKinParticles when executing rpwa::decayTopology::initKinematicsData()"<<std::endl;
				return false;
			}
			bp::list decayKinParticlesList = getDecayKinParticlesList();
			std::vector<std::string> productionKinematicsParticleNames;
			if(not rpwa::py::convertBPObjectToVector<std::string>(prodKinParticlesList, productionKinematicsParticleNames))
			{
				PyErr_SetString(PyExc_TypeError, "Got invalid input for prodKinParticles when executing rpwa::decayTopology::initKinematicsData()");
				bp::throw_error_already_set();
			}

			std::vector<std::string> decayKinematicsParticleNames;
			if(not rpwa::py::convertBPObjectToVector<std::string>(decayKinParticlesList, decayKinematicsParticleNames))
			{
				PyErr_SetString(PyExc_TypeError, "Got invalid input for decayKinParticles when executing rpwa::decayTopology::initKinematicsData()");
				bp::throw_error_already_set();
			}
			return self.initKinematicsData(productionKinematicsParticleNames, decayKinematicsParticleNames);
		}
		return self.initKinematicsData(*prodKinParticles, *decayKinParticles);
	}

	bool decayTopology_readKinematicsData(rpwa::decayTopology& self, PyObject* pyProdKinMomenta, PyObject* pyDecayKinMomenta) {

		TClonesArray* prodKinMomenta = rpwa::py::convertFromPy<TClonesArray*>(pyProdKinMomenta);
		TClonesArray* decayKinMomenta = rpwa::py::convertFromPy<TClonesArray*>(pyDecayKinMomenta);
		if((not prodKinMomenta) or (not decayKinMomenta)) {

			bp::extract<bp::list> getProdKinMomentaList(pyProdKinMomenta);
			if(not getProdKinMomentaList.check()) {
				printErr<<"Got invalid input for prodKinMomenta when executing rpwa::decayTopology::readKinematicsData()"<<std::endl;
				return false;
			}
			bp::list prodKinMomentaList = getProdKinMomentaList();
			bp::extract<bp::list> getDecayKinMomentaList(pyDecayKinMomenta);
			if(not getDecayKinMomentaList.check()) {
				printErr<<"Got invalid input for decayKinMomenta when executing rpwa::decayTopology::readKinematicsData()"<<std::endl;
				return false;
			}
			bp::list decayKinMomentaList = getDecayKinMomentaList();
			std::vector<TVector3> productionKinematicsMomenta(len(prodKinMomentaList));
			for(unsigned int i = 0; i < len(prodKinMomentaList); ++i) {
				bp::object item = bp::extract<bp::object>(prodKinMomentaList[i]);
				productionKinematicsMomenta[i] = *rpwa::py::convertFromPy<TVector3*>(item.ptr());
			}
			std::vector<TVector3> decayKinematicsMomenta(len(decayKinMomentaList));
			for(unsigned int i = 0; i < len(decayKinMomentaList); ++i) {
				bp::object item = bp::extract<bp::object>(decayKinMomentaList[i]);
				decayKinematicsMomenta[i] = *rpwa::py::convertFromPy<TVector3*>(item.ptr());
			}
			return self.readKinematicsData(productionKinematicsMomenta, decayKinematicsMomenta);
		}
		return self.readKinematicsData(*prodKinMomenta, *decayKinMomenta);
	}

	bool decayTopology_revertMomenta1(rpwa::decayTopology& self) {
		return self.revertMomenta();
	}

	bool decayTopology_revertMomenta2(rpwa::decayTopology& self, const bp::object& pyIndexMap) {
		std::vector<unsigned int> indexMap;
		if(not rpwa::py::convertBPObjectToVector<unsigned int>(pyIndexMap, indexMap)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for indexMap when executing rpwa::decayTopology::revertMomenta()");
			bp::throw_error_already_set();
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
