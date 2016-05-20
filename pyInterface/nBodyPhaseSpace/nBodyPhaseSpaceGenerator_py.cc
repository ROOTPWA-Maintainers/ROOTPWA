#include "nBodyPhaseSpaceGenerator_py.h"

#include <boost/python.hpp>

#include <TLorentzVector.h>

#include "nBodyPhaseSpaceGenerator.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;


namespace {

	bool nBodyPhaseSpaceGenerator_setDecay(rpwa::nBodyPhaseSpaceGenerator& self,
	                                       const bp::object& pyDaughterMasses)
	{
		std::vector<double> daughterMasses;
		if(not rpwa::py::convertBPObjectToVector<double>(pyDaughterMasses, daughterMasses)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for daughterMasses when executing rpwa::nBodyPhaseSpaceGenerator::setDecay()");
			bp::throw_error_already_set();
		}
		return self.setDecay(daughterMasses);
	}

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

	void nBodyPhaseSpaceGenerator_calcEventKinematics(rpwa::nBodyPhaseSpaceGenerator& self, PyObject* PyNBody) {
		TLorentzVector* nBody = rpwa::py::convertFromPy<TLorentzVector*>(PyNBody);
		if(not nBody) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for nBody when executing rpwa::nBodyPhaseSpace::calcEventKinematics()");
			bp::throw_error_already_set();
		}
		self.calcEventKinematics(*nBody);
	}

	PyObject* nBodyPhaseSpaceGenerator_daughter(const rpwa::nBodyPhaseSpaceGenerator& self, const int index) {
		return rpwa::py::convertToPy<TLorentzVector>(self.daughter(index));
	}

	bp::list nBodyPhaseSpaceGenerator_daughters(const rpwa::nBodyPhaseSpaceGenerator& self) {
		const std::vector<TLorentzVector>& daughters = self.daughters();
		bp::list pyDaughters;
		for(unsigned int i = 0; i < daughters.size(); ++i) {
			pyDaughters.append(bp::object(bp::handle<>(rpwa::py::convertToPy<TLorentzVector>(daughters[i]))));
		}
		return pyDaughters;
	}

}

void rpwa::py::exportNBodyPhaseSpaceGenerator() {

	bp::scope theScope = bp::class_<rpwa::nBodyPhaseSpaceGenerator>("nBodyPhaseSpaceGenerator")

		.def(bp::self_ns::str(bp::self))

		.def("setDecay", &nBodyPhaseSpaceGenerator_setDecay)

		.def("generateDecay", &nBodyPhaseSpaceGenerator_generateDecay)
		.def(
			"generateDecayAccepted"
			, &nBodyPhaseSpaceGenerator_generateDecayAccepted
			, (bp::arg("nBody"),
			   bp::arg("maxWeight") = 0.)
		)

		.def("pickMasses", &rpwa::nBodyPhaseSpaceGenerator::pickMasses)
		.def("calcWeight", &rpwa::nBodyPhaseSpaceGenerator::calcWeight)
		.def("pickAngles", &rpwa::nBodyPhaseSpaceGenerator::pickAngles)

		.def("calcEventKinematics", &nBodyPhaseSpaceGenerator_calcEventKinematics)

		.def("setKinematicsType", &rpwa::nBodyPhaseSpaceGenerator::setKinematicsType)
		.def("kinematicsType", &rpwa::nBodyPhaseSpaceGenerator::kinematicsType)

		.def("setWeightType", &rpwa::nBodyPhaseSpaceGenerator::setWeightType)
		.def("weightType", &rpwa::nBodyPhaseSpaceGenerator::weightType)

		.def("setMaxWeight", &rpwa::nBodyPhaseSpaceGenerator::setMaxWeight)
		.def("maxWeight", &rpwa::nBodyPhaseSpaceGenerator::maxWeight)
		.def("normalization", &rpwa::nBodyPhaseSpaceGenerator::normalization)
		.def("eventWeight", &rpwa::nBodyPhaseSpaceGenerator::eventWeight)
		.def("maxWeightObserved", &rpwa::nBodyPhaseSpaceGenerator::maxWeightObserved)
		.def("resetMaxWeightObserved", &rpwa::nBodyPhaseSpaceGenerator::resetMaxWeightObserved)

		.def("estimateMaxWeight", &rpwa::nBodyPhaseSpaceGenerator::estimateMaxWeight)
		.def("eventAccepted", &rpwa::nBodyPhaseSpaceGenerator::eventAccepted)

		.def("daughter", &nBodyPhaseSpaceGenerator_daughter)
		.def("daughters", &nBodyPhaseSpaceGenerator_daughters)
		.def("nmbOfDaughters", &rpwa::nBodyPhaseSpaceGenerator::nmbOfDaughters)
		.def("daughterMass", &rpwa::nBodyPhaseSpaceGenerator::daughterMass)
		.def("intermediateMass", &rpwa::nBodyPhaseSpaceGenerator::intermediateMass)
		.def("breakupMom", &rpwa::nBodyPhaseSpaceGenerator::breakupMom)
		.def("cosTheta", &rpwa::nBodyPhaseSpaceGenerator::cosTheta)
		.def("phi", &rpwa::nBodyPhaseSpaceGenerator::phi);

	bp::enum_<rpwa::nBodyPhaseSpaceGenerator::kinematicsTypeEnum>("kinematicsTypeEnum")
		.value("BLOCK", rpwa::nBodyPhaseSpaceGenerator::BLOCK)
		.value("RAUBOLD_LYNCH", rpwa::nBodyPhaseSpaceGenerator::RAUBOLD_LYNCH)
		.export_values();

	bp::enum_<rpwa::nBodyPhaseSpaceGenerator::weightTypeEnum>("weightTypeEnum")
		.value("S_U_CHUNG", rpwa::nBodyPhaseSpaceGenerator::S_U_CHUNG)
		.value("NUPHAZ", rpwa::nBodyPhaseSpaceGenerator::NUPHAZ)
		.value("GENBOD", rpwa::nBodyPhaseSpaceGenerator::GENBOD)
		.value("FLAT", rpwa::nBodyPhaseSpaceGenerator::FLAT)
		.export_values();

	theScope.attr("BLOCK") = rpwa::nBodyPhaseSpaceGenerator::BLOCK;
	theScope.attr("RAUBOLD_LYNCH") = rpwa::nBodyPhaseSpaceGenerator::RAUBOLD_LYNCH;
	theScope.attr("S_U_CHUNG") = rpwa::nBodyPhaseSpaceGenerator::S_U_CHUNG;
	theScope.attr("NUPHAZ") = rpwa::nBodyPhaseSpaceGenerator::NUPHAZ;
	theScope.attr("GENBOD") = rpwa::nBodyPhaseSpaceGenerator::GENBOD;
	theScope.attr("FLAT") = rpwa::nBodyPhaseSpaceGenerator::FLAT;

}
