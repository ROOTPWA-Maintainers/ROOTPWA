#include "nBodyPhaseSpaceGenerator_py.h"

#include <boost/python.hpp>

#include <TLorentzVector.h>

#include "nBodyPhaseSpaceGenerator.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;


namespace {

	double nBodyPhaseSpaceGenerator_generateDecay(rpwa::nBodyPhaseSpaceGenerator& self, PyObject* PyNBody) {
		TLorentzVector* nBody = rpwa::py::convertFromPy<TLorentzVector*>(PyNBody);
		if(not nBody) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for nBody when executing rpwa::nBodyPhaseSpace::generateDecay()");
			bp::throw_error_already_set();
		}
		return self.generateDecay(*nBody);
	}

	bool nBodyPhaseSpaceGenerator_generateDecayAccepted(rpwa::nBodyPhaseSpaceGenerator& self,
	                                                    PyObject* PyNBody,
	                                                    const double maxWeight = 0.)
	{
		TLorentzVector* nBody = rpwa::py::convertFromPy<TLorentzVector*>(PyNBody);
		if(not nBody) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for nBody when executing rpwa::nBodyPhaseSpace::generateDecayAccepted()");
			bp::throw_error_already_set();
		}
		return self.generateDecayAccepted(*nBody, maxWeight);
	}

}

void rpwa::py::exportNBodyPhaseSpaceGenerator() {

	bp::scope theScope = bp::class_<rpwa::nBodyPhaseSpaceGenerator, bp::bases<rpwa::nBodyPhaseSpaceKinematics> >("nBodyPhaseSpaceGenerator")

		.def(bp::self_ns::str(bp::self))

		.def("generateDecay", &nBodyPhaseSpaceGenerator_generateDecay)
		.def(
			"generateDecayAccepted"
			, &nBodyPhaseSpaceGenerator_generateDecayAccepted
			, (bp::arg("nBody"),
			   bp::arg("maxWeight") = 0.)
		)

		.def("pickMasses", &rpwa::nBodyPhaseSpaceGenerator::pickMasses)
		.def("pickAngles", &rpwa::nBodyPhaseSpaceGenerator::pickAngles)

		.def("setMaxWeight", &rpwa::nBodyPhaseSpaceGenerator::setMaxWeight)
		.def("maxWeight", &rpwa::nBodyPhaseSpaceGenerator::maxWeight)

		.def("estimateMaxWeight", &rpwa::nBodyPhaseSpaceGenerator::estimateMaxWeight)
		.def("eventAccepted", &rpwa::nBodyPhaseSpaceGenerator::eventAccepted);

}
