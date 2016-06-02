#include"getMassShapes_py.h"
#include"stlContainers_py.h"
#include"isobarDecayTopology.h"

namespace bp = boost::python;
namespace {
	bp::list getMassShapes_py(rpwa::isobarDecayTopologyPtr &topo,
			         double                         mass,
			         bool                           useBarrierFactors) {
		return bp::list(getMassShapes(topo, mass, useBarrierFactors));
	}
}

void rpwa::py::exportGetMassShapes(){
	bp::def(
		"getMassShapes"
		, &::getMassShapes_py
		, (bp::arg("topology"),
		   bp::arg("mass"),
		   bp::arg("useBarrierFactors") = true)
	);
}
