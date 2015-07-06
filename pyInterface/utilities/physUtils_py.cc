#include "physUtils_py.h"

#include <boost/python.hpp>

#include "physUtils.hpp"

namespace bp = boost::python;


void rpwa::py::exportPhysUtils() {

	bp::def(
		"breakupMomentumSquared"
		, &rpwa::breakupMomentumSquared
		, (bp::arg("M"),
		   bp::arg("m1"),
		   bp::arg("m2"),
		   bp::arg("allowSubThr")=false)
	);
	bp::def(
		"breakupMomentum"
		, &rpwa::breakupMomentum
		, (bp::arg("M"),
		   bp::arg("m1"),
		   bp::arg("m2"))
	);
	bp::def(
		"breakupMomentumSquared"
		, &rpwa::breakupMomentumSquared
		, (bp::arg("M"),
		   bp::arg("m1"),
		   bp::arg("m2"),
		   bp::arg("allowSubThr")=false)
	);
	bp::def(
		"barrierFactorSquared"
		, &rpwa::barrierFactorSquared
		, (bp::arg("L"),
		   bp::arg("breakupMom"),
		   bp::arg("debug") = false,
		   bp::arg("Pr") = 0.1973)
	);
	bp::def(
		"breitWigner"
		, &rpwa::breitWigner
		, (bp::arg("M"),
		   bp::arg("M0"),
		   bp::arg("Gamma0"),
		   bp::arg("L"),
		   bp::arg("q"),
		   bp::arg("q0"))
	);

}
