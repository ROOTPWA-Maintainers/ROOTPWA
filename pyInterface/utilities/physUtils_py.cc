#include "physUtils_py.h"

#include <boost/python.hpp>

#include <TLorentzVector.h>

#include "physUtils.hpp"
#include "rootConverters_py.h"

namespace bp = boost::python;


double tPrime_py(PyObject* beam, PyObject* target, PyObject* out) {
	TLorentzVector* lvBeam   = rpwa::py::convertFromPy<TLorentzVector*>(beam  );
	TLorentzVector* lvTarget = rpwa::py::convertFromPy<TLorentzVector*>(target);
	TLorentzVector* lvOut    = rpwa::py::convertFromPy<TLorentzVector*>(out   );

	return rpwa::tPrime(*lvBeam, *lvTarget, *lvOut);
}


void rpwa::py::exportPhysUtils() {

	bp::def(
		"tPrime"
		, &tPrime_py
		, (bp::arg("pBeam"),
		   bp::arg("pTarget"),
		   bp::arg("pOut"))
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
		"breakupMomentum"
		, &rpwa::breakupMomentum
		, (bp::arg("M"),
		   bp::arg("m1"),
		   bp::arg("m2"))
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
