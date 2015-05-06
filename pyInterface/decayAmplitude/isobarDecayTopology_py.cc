#include "isobarDecayTopology_py.h"

#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;

namespace {

	struct isobarDecayTopologyWrapper : public rpwa::isobarDecayTopology,
	                                           bp::wrapper<rpwa::isobarDecayTopology>
	{

		isobarDecayTopologyWrapper()
			: rpwa::isobarDecayTopology(),
			  bp::wrapper<rpwa::isobarDecayTopology>() { }

		isobarDecayTopologyWrapper(const rpwa::productionVertexPtr& productionVertex,
		                           const bp::object&                pyIsobarDecayVertices,
		                           const bp::object&                pyFsParticles,
		                           const bool                       performTopologyCheck = true)
			: rpwa::isobarDecayTopology(),
			  bp::wrapper<rpwa::isobarDecayTopology>()
		{

			// Translate the fsParticles
			std::vector<rpwa::particlePtr> fsParticles;
			if(not rpwa::py::convertBPObjectToVector<rpwa::particlePtr>(pyFsParticles, fsParticles)) {
				PyErr_SetString(PyExc_TypeError, "Got invalid input for fsParticles when executing rpwa::isobarDecayTopology::isobarDecayTopology()");
				bp::throw_error_already_set();
			}

			std::vector<rpwa::isobarDecayVertexPtr> isobarDecayVerticesIDV;
			if(rpwa::py::convertBPObjectToVector<rpwa::isobarDecayVertexPtr>(pyIsobarDecayVertices, isobarDecayVerticesIDV)) {
				rpwa::isobarDecayTopology::constructDecay(productionVertex, isobarDecayVerticesIDV, fsParticles, performTopologyCheck);
				return;
			}

			std::vector<rpwa::interactionVertexPtr> isobarDecayVerticesIV;
			if(rpwa::py::convertBPObjectToVector<rpwa::interactionVertexPtr>(pyIsobarDecayVertices, isobarDecayVerticesIV)) {
				rpwa::isobarDecayTopology::constructDecay(productionVertex, isobarDecayVerticesIV, fsParticles, performTopologyCheck);
				return;
			}

			PyErr_SetString(PyExc_TypeError, "Got invalid input for decayVertices when executing rpwa::isobarDecayTopology::isobarDecayTopology()");
			bp::throw_error_already_set();

		}

		isobarDecayTopologyWrapper(const rpwa::isobarDecayTopology& topo)
			: rpwa::isobarDecayTopology(topo),
			  bp::wrapper<rpwa::isobarDecayTopology>() { }

		isobarDecayTopologyWrapper(const rpwa::decayTopology& topo)
			: rpwa::isobarDecayTopology(topo),
			  bp::wrapper<rpwa::isobarDecayTopology>() { }

		void clear() {
			if(bp::override clear = this->get_override("clear")) {
				clear();
			}
			rpwa::isobarDecayTopology::clear();
		}

		void default_clear() {
			rpwa::isobarDecayTopology::clear();
		}

	};

	bp::list isobarDecayTopology_isobarDecayVertices(const rpwa::isobarDecayTopology& self) {
		return bp::list(self.isobarDecayVertices());
	}

	rpwa::isobarDecayTopology isobarDecayTopology_subDecay(rpwa::isobarDecayTopology& self,
	                                                       const rpwa::isobarDecayVertexPtr& startVert,
	                                                       const bool linkToParentTopo = false) {
		return self.subDecay(startVert, linkToParentTopo);
	}

	rpwa::isobarDecayTopology isobarDecayTopology_joinDaughterDecays(const rpwa::isobarDecayVertexPtr& parentVertex,
	                                                                 const rpwa::isobarDecayTopology&  daughter1Decay,
	                                                                 const rpwa::isobarDecayTopology&  daughter2Decay)
	{
		return rpwa::isobarDecayTopology::joinDaughterDecays(parentVertex, daughter1Decay, daughter2Decay);
	}

	PyObject* isobarDecayTopology_calcIsobarLzVec(rpwa::isobarDecayTopology& self) {
		return rpwa::py::convertToPy<TLorentzVector>(self.calcIsobarLzVec());
	}

	std::string isobarDecayTopology_writeGraphViz1(rpwa::isobarDecayTopology& self) {
		std::stringstream sstr;
		self.writeGraphViz(sstr);
		return sstr.str();
	}

	bool isobarDecayTopology_writeGraphViz2(rpwa::isobarDecayTopology& self, const std::string& outFileName) {
		return self.writeGraphViz(outFileName);
	}

	bp::list isobarDecayTopology_getIsospinSymmetrization(rpwa::isobarDecayTopology& self) {
		std::vector<rpwa::symTermMap> symTermMap = self.getIsospinSymmetrization();
		bp::list retval;
		for(unsigned int i = 0; i < symTermMap.size(); ++i) {
			bp::dict dict;
			dict["factor"] = symTermMap[i].factor;
			dict["fsPartPermMap"] = bp::list(symTermMap[i].fsPartPermMap);
			retval.append(dict);
		}
		return retval;
	}

	bp::list isobarDecayTopology_getBoseSymmetrization(const rpwa::isobarDecayTopology& self) {
		std::vector<rpwa::symTermMap> symTermMap = self.getBoseSymmetrization();
		bp::list retval;
		for(unsigned int i = 0; i < symTermMap.size(); ++i) {
			bp::dict dict;
			dict["factor"] = symTermMap[i].factor;
			dict["fsPartPermMap"] = bp::list(symTermMap[i].fsPartPermMap);
			retval.append(dict);
		}
		return retval;
	}

	bool isobarDecayTopology_isobarIsAffectedByPermutation(const rpwa::isobarDecayTopology& self,
	                                                       const rpwa::isobarDecayVertexPtr& vertex,
	                                                       const bp::list& pyPermutation)
	{
		std::vector<unsigned int> permutation;
		for(int i = 0; i < bp::len(pyPermutation); ++i) {
			permutation.push_back(bp::extract<unsigned int>(pyPermutation[i]));
		}
		return self.isobarIsAffectedByPermutation(vertex, permutation);
	}

	bool isobarDecayTopology_daughtersAreAffectedByPermutation(const rpwa::isobarDecayTopology& self,
		                                                       const rpwa::isobarDecayVertexPtr& vertex,
		                                                       const bp::list& pyPermutation)
	{
		std::vector<unsigned int> permutation;
		for(int i = 0; i < bp::len(pyPermutation); ++i) {
			permutation.push_back(bp::extract<unsigned int>(pyPermutation[i]));
		}
		return self.daughtersAreAffectedByPermutation(vertex, permutation);
	}

	bp::list isobarDecayTopology_getFsPartIndicesConnectedToVertex(const rpwa::isobarDecayTopology& self,
	                                                               const rpwa::isobarDecayVertexPtr& vertex)
	{
		return bp::list(self.getFsPartIndicesConnectedToVertex(vertex));
	}

	bp::list isobarDecayTopology_findIsobarBoseSymVertices(const rpwa::isobarDecayTopology& self) {
		return bp::list(self.findIsobarBoseSymVertices());
	}

}

void rpwa::py::exportIsobarDecayTopology() {

	bp::class_<isobarDecayTopologyWrapper, bp::bases<rpwa::decayTopology> >("isobarDecayTopology")

		.def(bp::init<const rpwa::productionVertexPtr, const bp::object&, const bp::object&, bp::optional<bool> >())
		.def(bp::init<const rpwa::isobarDecayTopology&>())
		.def(bp::init<const rpwa::decayTopology&>())

		.def(bp::self_ns::str(bp::self))

		.def(
			"clone"
			, &rpwa::isobarDecayTopology::clone
			, (bp::arg("cloneFsParticles")=false,
			   bp::arg("cloneProdKinematics")=false)
		)

		.def("clear", &isobarDecayTopologyWrapper::clear, &isobarDecayTopologyWrapper::default_clear)
		.def("clear", &rpwa::isobarDecayTopology::clear)

		.def("isobarDecayVertices", &isobarDecayTopology_isobarDecayVertices)
		.def(
			"XIsobarDecayVertex"
			, &rpwa::isobarDecayTopology::XIsobarDecayVertex
			, bp::return_value_policy<bp::copy_const_reference>()
		)

		.def("checkTopology", &rpwa::isobarDecayTopology::checkTopology)
		.def("checkConsistency", &rpwa::isobarDecayTopology::checkConsistency)

// This one is missing because it returns something defined in decayGraph.hpp, which is currently omitted.
/*		.def(
			"subDecay"
			, &rpwa::isobarDecayTopology::subDecay
			, (bp::arg("startNd"),
			   bp::arg("linkToParentTopo")=false)
		)*/

		.def(
			"subDecay"
			, &isobarDecayTopology_subDecay
			, (bp::arg("startVert"),
			   bp::arg("linkToParentTopo")=false)
		)

		.def("subDecayConsistent", &rpwa::isobarDecayTopology::subDecayConsistent)

		.def("addDecay", &rpwa::isobarDecayTopology::addDecay)

		.def("joinDaughterDecays", &isobarDecayTopology_joinDaughterDecays)
		.staticmethod("joinDaughterDecays")

		.def("calcIsobarLzVec", &isobarDecayTopology_calcIsobarLzVec)

		.def("writeGraphViz", &isobarDecayTopology_writeGraphViz1)
		.def("writeGraphViz", &isobarDecayTopology_writeGraphViz2)

		.def("calcIsobarCharges", &rpwa::isobarDecayTopology::calcIsobarCharges, (bp::arg("quiet")=false))
		.def("calcIsobarBaryonNmbs", &rpwa::isobarDecayTopology::calcIsobarBaryonNmbs)

		.def(
			"getIsospinClebschGordanProduct"
			, &rpwa::isobarDecayTopology::getIsospinClebschGordanProduct
			, (bp::arg("vertex")=rpwa::isobarDecayVertexPtr())
		)

		.def("getIsospinSymmetrization", &isobarDecayTopology_getIsospinSymmetrization)
		.def("getBoseSymmetrization", &isobarDecayTopology_getBoseSymmetrization)
		.def("isobarIsAffectedByPermutation", &isobarDecayTopology_isobarIsAffectedByPermutation)
		.def("daughtersAreAffectedByPermutation", &isobarDecayTopology_daughtersAreAffectedByPermutation)
		.def("getFsPartIndicesConnectedToVertex", &isobarDecayTopology_getFsPartIndicesConnectedToVertex)
		.def("findIsobarBoseSymVertices", &isobarDecayTopology_findIsobarBoseSymVertices)

		.add_static_property("debugIsobarDecayTopology", &rpwa::isobarDecayTopology::debug, &isobarDecayTopology::setDebug);

	bp::register_ptr_to_python<rpwa::isobarDecayTopologyPtr>();

}
