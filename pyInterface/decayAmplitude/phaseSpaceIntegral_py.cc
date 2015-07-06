#include "phaseSpaceIntegral_py.h"

#include <boost/python.hpp>

#include "phaseSpaceIntegral.h"

namespace bp = boost::python;


namespace {

	bp::tuple integralTableContainer_getSubWaveNameFromVertex(const rpwa::isobarDecayVertex& vertex) {
			rpwa::isobarDecayVertexPtr vertexPtr = rpwa::isobarDecayVertexPtr();
			rpwa::isobarDecayTopologyPtr subDecay = rpwa::isobarDecayTopologyPtr();
			const std::string name = rpwa::integralTableContainer::getSubWaveNameFromVertex(vertex, vertexPtr, subDecay);
			return bp::make_tuple(name, vertexPtr, subDecay);
	}

}

void rpwa::py::exportPhaseSpaceIntegral() {

	bp::class_<phaseSpaceIntegral, boost::noncopyable>("phaseSpaceIntegral", bp::no_init)

		.add_static_property(
			"instance"
			, bp::make_function(&rpwa::phaseSpaceIntegral::instance,
			                    bp::return_value_policy<bp::reference_existing_object>())
		)

		.def("__call__", &rpwa::phaseSpaceIntegral::operator())
		.def("removeVertex", &rpwa::phaseSpaceIntegral::removeVertex);

	bp::class_<integralTableContainer>("integralTableContainer")

		.def(bp::init<const rpwa::isobarDecayVertex&>())

		.def("__call__", &rpwa::integralTableContainer::operator())

		.def("getSubWaveNameFromVertex", &integralTableContainer_getSubWaveNameFromVertex)
		.staticmethod("getSubWaveNameFromVertex");

}
