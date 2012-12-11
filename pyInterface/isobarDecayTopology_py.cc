#include "isobarDecayTopology_py.h"

#include "rootConverters_py.h"

namespace bp = boost::python;

namespace {

	struct isobarDecayTopologyWrapper : public rpwa::isobarDecayTopology,
	                                           bp::wrapper<rpwa::isobarDecayTopology>
	{

		isobarDecayTopologyWrapper()
			: rpwa::isobarDecayTopology(),
			  bp::wrapper<rpwa::isobarDecayTopology>() { };

		isobarDecayTopologyWrapper(const rpwa::productionVertexPtr& productionVertex,
		                           const bp::object&                pyIsobarDecayVertices,
		                           const bp::object&                pyFsParticles)
			: rpwa::isobarDecayTopology(),
			  bp::wrapper<rpwa::isobarDecayTopology>()
		{

			// Translate the fsParticles
			bp::list pyListFsParticles = bp::extract<bp::list>(pyFsParticles);
			std::vector<rpwa::particlePtr> fsParticles(bp::len(pyListFsParticles), rpwa::particlePtr());
			for(int i = 0; i < bp::len(pyListFsParticles); ++i) {
				fsParticles[i] = bp::extract<rpwa::particlePtr>(pyListFsParticles[i]);
			}

			// Translate the isobarDecayVertices
			bp::list pyListIsobarDecayVertices = bp::extract<bp::list>(pyIsobarDecayVertices);
			if(bp::len(pyListIsobarDecayVertices) == 0) {
				rpwa::isobarDecayTopology::constructDecay(productionVertex, std::vector<rpwa::isobarDecayVertexPtr>(), fsParticles);
				return;
			}

			bp::extract<rpwa::isobarDecayVertexPtr> get_iDVP(pyListIsobarDecayVertices[0]);
			if(get_iDVP.check()) {
				std::vector<rpwa::isobarDecayVertexPtr> isobarDecayVerticesIDV(bp::len(pyListIsobarDecayVertices), rpwa::isobarDecayVertexPtr());
				for(int i = 0; i < bp::len(pyListIsobarDecayVertices); ++i) {
					isobarDecayVerticesIDV[i] = bp::extract<rpwa::isobarDecayVertexPtr>(pyListIsobarDecayVertices[i]);
				}
				rpwa::isobarDecayTopology::constructDecay(productionVertex, isobarDecayVerticesIDV, fsParticles);
				return;
			}

			bp::extract<rpwa::interactionVertexPtr> get_iVP(pyListIsobarDecayVertices[0]);
			if(get_iVP.check()) {
				std::vector<rpwa::interactionVertexPtr> isobarDecayVerticesIV(bp::len(pyListIsobarDecayVertices), rpwa::isobarDecayVertexPtr());
				for(int i = 0; i < bp::len(pyListIsobarDecayVertices); ++i) {
					isobarDecayVerticesIV[i] = bp::extract<rpwa::interactionVertexPtr>(pyListIsobarDecayVertices[i]);
				}
				rpwa::isobarDecayTopology::constructDecay(productionVertex, isobarDecayVerticesIV, fsParticles);
				return;
			}

			printErr<<"Got invalid input when executing rpwa::isobarDecayTopology::isobarDecayTopology()."<<std::endl;

		};

		isobarDecayTopologyWrapper(const rpwa::isobarDecayTopology& topo)
			: rpwa::isobarDecayTopology(topo),
			  bp::wrapper<rpwa::isobarDecayTopology>() { };

		isobarDecayTopologyWrapper(const rpwa::decayTopology& topo)
			: rpwa::isobarDecayTopology(topo),
			  bp::wrapper<rpwa::isobarDecayTopology>() { };

		void clear() {
			if(bp::override clear = this->get_override("clear")) {
				clear();
			}
			rpwa::isobarDecayTopology::clear();
		};

		void default_clear() {
			rpwa::isobarDecayTopology::clear();
		};

	};

	bp::list isobarDecayTopology_isobarDecayVertices(const rpwa::isobarDecayTopology& self) {
		return bp::list(self.isobarDecayVertices());
	};

	rpwa::isobarDecayTopology isobarDecayTopology_joinDaughterDecays(const rpwa::isobarDecayVertexPtr& parentVertex,
	                                                                 const rpwa::isobarDecayTopology&  daughter1Decay,
	                                                                 const rpwa::isobarDecayTopology&  daughter2Decay)
	{
		return rpwa::isobarDecayTopology::joinDaughterDecays(parentVertex, daughter1Decay, daughter2Decay);
	};

	const PyObject* isobarDecayTopology_calcIsobarLzVec(rpwa::isobarDecayTopology& self) {
		return rpwa::py::convertToPy<TLorentzVector>(self.calcIsobarLzVec());
	};

	std::string isobarDecayTopology_writeGraphViz1(rpwa::isobarDecayTopology& self) {
		std::stringstream sstr;
		self.writeGraphViz(sstr);
		return sstr.str();
	};

	bool isobarDecayTopology_writeGraphViz2(rpwa::isobarDecayTopology& self, const std::string& outFileName) {
		return self.writeGraphViz(outFileName);
	};

}

void rpwa::py::exportIsobarDecayTopology() {

	bp::class_<isobarDecayTopologyWrapper, bp::bases<rpwa::decayTopology> >("isobarDecayTopology")

		.def(bp::init<rpwa::productionVertexPtr, bp::object&, bp::object&>())
		.def(bp::init<rpwa::isobarDecayTopology&>())
		.def(bp::init<rpwa::decayTopology&>())

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
			, &rpwa::decayTopology::subDecay
			, (bp::arg("startNd"),
			   bp::arg("linkToMotherTopo")=false)
		)*/

		.def("addDecay", &rpwa::isobarDecayTopology::addDecay)

		.def("joinDaughterDecays", &isobarDecayTopology_joinDaughterDecays)
		.staticmethod("joinDaughterDecays")

		.def(
			"calcIsobarLzVec"
			, &isobarDecayTopology_calcIsobarLzVec
			, bp::return_value_policy<bp::manage_new_object>()
		)

		.def("writeGraphViz", &isobarDecayTopology_writeGraphViz1)
		.def("writeGraphViz", &isobarDecayTopology_writeGraphViz2)


		.def("calcIsobarCharges", &rpwa::isobarDecayTopology::calcIsobarCharges)
		.def("calcIsobarBaryonNmbs", &rpwa::isobarDecayTopology::calcIsobarBaryonNmbs)

		.add_static_property("debugIsobarDecayTopology", &rpwa::isobarDecayTopology::debug, &isobarDecayTopology::setDebug);

	bp::register_ptr_to_python<rpwa::isobarDecayTopologyPtr>();

};
