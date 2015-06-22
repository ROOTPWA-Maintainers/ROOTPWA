#include "nBodyPhaseSpaceGen_py.h"

#include <boost/python.hpp>

#include <TLorentzVector.h>

#include "nBodyPhaseSpaceGen.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;


namespace {

	bool nBodyPhaseSpaceGen_setDecay(rpwa::nBodyPhaseSpaceGen& self,
	                                 const bp::object& pyDaughterMasses)
	{
		std::vector<double> daughterMasses;
		if(not rpwa::py::convertBPObjectToVector<double>(pyDaughterMasses, daughterMasses)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for daughterMasses when executing rpwa::nBodyPhaseSpaceGen::setDecay()");
			bp::throw_error_already_set();
		}
		return self.setDecay(daughterMasses);
	}

	double nBodyPhaseSpaceGen_generateDecay(rpwa::nBodyPhaseSpaceGen& self, PyObject* PyNBody) {
		TLorentzVector* nBody = rpwa::py::convertFromPy<TLorentzVector*>(PyNBody);
		if(not nBody) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for nBody when executing rpwa::nBodyPhaseSpace::generateDecay()");
			bp::throw_error_already_set();
		}
		return self.generateDecay(*nBody);
	}

	bool nBodyPhaseSpaceGen_generateDecayAccepted(rpwa::nBodyPhaseSpaceGen& self,
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

	void nBodyPhaseSpaceGen_calcEventKinematics(rpwa::nBodyPhaseSpaceGen& self, PyObject* PyNBody) {
		TLorentzVector* nBody = rpwa::py::convertFromPy<TLorentzVector*>(PyNBody);
		if(not nBody) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for nBody when executing rpwa::nBodyPhaseSpace::calcEventKinematics()");
			bp::throw_error_already_set();
		}
		self.calcEventKinematics(*nBody);
	}

	PyObject* nBodyPhaseSpaceGen_daughter(const rpwa::nBodyPhaseSpaceGen& self, const int index) {
		return rpwa::py::convertToPy<TLorentzVector>(self.daughter(index));
	}

	bp::list nBodyPhaseSpaceGen_daughters(const rpwa::nBodyPhaseSpaceGen& self) {
		const std::vector<TLorentzVector>& daughters = self.daughters();
		bp::list pyDaughters;
		for(unsigned int i = 0; i < daughters.size(); ++i) {
			pyDaughters.append(bp::object(bp::handle<>(rpwa::py::convertToPy<TLorentzVector>(daughters[i]))));
		}
		return pyDaughters;
	}

}

void rpwa::py::exportNBodyPhaseSpaceGen() {

	bp::scope theScope = bp::class_<rpwa::nBodyPhaseSpaceGen>("nBodyPhaseSpaceGen")

		.def(bp::self_ns::str(bp::self))

		.def("setDecay", &nBodyPhaseSpaceGen_setDecay)

		.def("generateDecay", &nBodyPhaseSpaceGen_generateDecay)
		.def(
			"generateDecayAccepted"
			, &nBodyPhaseSpaceGen_generateDecayAccepted
			, (bp::arg("nBody"),
			   bp::arg("maxWeight") = 0.)
		)

		.def("pickMasses", &rpwa::nBodyPhaseSpaceGen::pickMasses)
		.def("calcWeight", &rpwa::nBodyPhaseSpaceGen::calcWeight)
		.def("pickAngles", &rpwa::nBodyPhaseSpaceGen::pickAngles)

		.def("calcEventKinematics", &nBodyPhaseSpaceGen_calcEventKinematics)

		.def("setKinematicsType", &rpwa::nBodyPhaseSpaceGen::setKinematicsType)
		.def("kinematicsType", &rpwa::nBodyPhaseSpaceGen::kinematicsType)

		.def("setVerbose", &rpwa::nBodyPhaseSpaceGen::setVerbose)

		.def("setWeightType", &rpwa::nBodyPhaseSpaceGen::setWeightType)
		.def("weightType", &rpwa::nBodyPhaseSpaceGen::weightType)

		.def("setMaxWeight", &rpwa::nBodyPhaseSpaceGen::setMaxWeight)
		.def("maxWeight", &rpwa::nBodyPhaseSpaceGen::maxWeight)
		.def("normalization", &rpwa::nBodyPhaseSpaceGen::normalization)
		.def("eventWeight", &rpwa::nBodyPhaseSpaceGen::eventWeight)
		.def("maxWeightObserved", &rpwa::nBodyPhaseSpaceGen::maxWeightObserved)
		.def("resetMaxWeightObserved", &rpwa::nBodyPhaseSpaceGen::resetMaxWeightObserved)

		.def("estimateMaxWeight", &rpwa::nBodyPhaseSpaceGen::estimateMaxWeight)
		.def("eventAccepted", &rpwa::nBodyPhaseSpaceGen::eventAccepted)

		.def("daughter", &nBodyPhaseSpaceGen_daughter)
		.def("daughters", &nBodyPhaseSpaceGen_daughters)
		.def("nmbOfDaughters", &rpwa::nBodyPhaseSpaceGen::nmbOfDaughters)
		.def("daughterMass", &rpwa::nBodyPhaseSpaceGen::daughterMass)
		.def("intermediateMass", &rpwa::nBodyPhaseSpaceGen::intermediateMass)
		.def("breakupMom", &rpwa::nBodyPhaseSpaceGen::breakupMom)
		.def("cosTheta", &rpwa::nBodyPhaseSpaceGen::cosTheta)
		.def("phi", &rpwa::nBodyPhaseSpaceGen::phi);

	bp::enum_<rpwa::nBodyPhaseSpaceGen::kinematicsTypeEnum>("kinematicsTypeEnum")
		.value("BLOCK", rpwa::nBodyPhaseSpaceGen::BLOCK)
		.value("RAUBOLD_LYNCH", rpwa::nBodyPhaseSpaceGen::RAUBOLD_LYNCH)
		.export_values();

	bp::enum_<rpwa::nBodyPhaseSpaceGen::weightTypeEnum>("weightTypeEnum")
		.value("S_U_CHUNG", rpwa::nBodyPhaseSpaceGen::S_U_CHUNG)
		.value("NUPHAZ", rpwa::nBodyPhaseSpaceGen::NUPHAZ)
		.value("GENBOD", rpwa::nBodyPhaseSpaceGen::GENBOD)
		.value("FLAT", rpwa::nBodyPhaseSpaceGen::FLAT)
		.export_values();

	theScope.attr("BLOCK") = rpwa::nBodyPhaseSpaceGen::BLOCK;
	theScope.attr("RAUBOLD_LYNCH") = rpwa::nBodyPhaseSpaceGen::RAUBOLD_LYNCH;
	theScope.attr("S_U_CHUNG") = rpwa::nBodyPhaseSpaceGen::S_U_CHUNG;
	theScope.attr("NUPHAZ") = rpwa::nBodyPhaseSpaceGen::NUPHAZ;
	theScope.attr("GENBOD") = rpwa::nBodyPhaseSpaceGen::GENBOD;
	theScope.attr("FLAT") = rpwa::nBodyPhaseSpaceGen::FLAT;

}
