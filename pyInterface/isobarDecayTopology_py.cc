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

		bp::list isobarDecayVertices__() const {
			return bp::list(rpwa::isobarDecayTopology::isobarDecayVertices());
		};

		static rpwa::isobarDecayTopology joinDaughterDecays__(const rpwa::isobarDecayVertexPtr& parentVertex,
		                                                      const rpwa::isobarDecayTopology&  daughter1Decay,
		                                                      const rpwa::isobarDecayTopology&  daughter2Decay)
		{
			return rpwa::isobarDecayTopology::joinDaughterDecays(parentVertex, daughter1Decay, daughter2Decay);
		};

		const PyObject* calcIsobarLzVec__() {
			return rpwa::py::convertToPy<TLorentzVector>(rpwa::isobarDecayTopology::calcIsobarLzVec());
		};

		std::string writeGraphViz__1() {
			std::stringstream sstr;
			rpwa::isobarDecayTopology::writeGraphViz(sstr);
			return sstr.str();
		};

		bool writeGraphViz__2(const std::string& outFileName) {
			return rpwa::isobarDecayTopology::writeGraphViz(outFileName);
		};

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
			, &isobarDecayTopologyWrapper::clone
			, (bp::arg("cloneFsParticles")=false,
			   bp::arg("cloneProdKinematics")=false)
		)

		.def("clear", &isobarDecayTopologyWrapper::clear, &isobarDecayTopologyWrapper::default_clear)

		.def("isobarDecayVertices", &isobarDecayTopologyWrapper::isobarDecayVertices__)
		.def(
			"XIsobarDecayVertex"
			, &isobarDecayTopologyWrapper::XIsobarDecayVertex
			, bp::return_value_policy<bp::copy_const_reference>()
		)

		.def("checkTopology", &isobarDecayTopologyWrapper::checkTopology)
		.def("checkConsistency", &isobarDecayTopologyWrapper::checkConsistency)

// This one is missing because it returns something defined in decayGraph.hpp, which is currently omitted.
/*		.def(
			"subDecay"
			, &decayTopologyWrapper::subDecay
			, (bp::arg("startNd"),
			   bp::arg("linkToMotherTopo")=false)
		)*/

		.def("addDecay", &isobarDecayTopologyWrapper::addDecay)

		.def("joinDaughterDecays", &isobarDecayTopologyWrapper::joinDaughterDecays__)
		.staticmethod("joinDaughterDecays")

		.def(
			"calcIsobarLzVec"
			, &isobarDecayTopologyWrapper::calcIsobarLzVec__
			, bp::return_value_policy<bp::manage_new_object>()
		)

		.def("writeGraphViz", &isobarDecayTopologyWrapper::writeGraphViz__1)
		.def("writeGraphViz", &isobarDecayTopologyWrapper::writeGraphViz__2)


		.def("calcIsobarCharges", &isobarDecayTopologyWrapper::calcIsobarCharges)
		.def("calcIsobarBaryonNmbs", &isobarDecayTopologyWrapper::calcIsobarBaryonNmbs)

		.add_static_property("debugIsobarDecayTopology", &isobarDecayTopologyWrapper::debug, &isobarDecayTopologyWrapper::setDebug);

	bp::register_ptr_to_python<rpwa::isobarDecayTopologyPtr>();

};
