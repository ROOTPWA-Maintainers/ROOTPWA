#include "nBodyPhaseSpaceKinematics_py.h"

#include <boost/python.hpp>

#include <TLorentzVector.h>

#include "nBodyPhaseSpaceKinematics.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;


namespace {

	bool nBodyPhaseSpaceKinematics_setDecay(rpwa::nBodyPhaseSpaceKinematics& self,
	                                        const bp::object& pyDaughterMasses)
	{
		std::vector<double> daughterMasses;
		if(not rpwa::py::convertBPObjectToVector<double>(pyDaughterMasses, daughterMasses)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for daughterMasses when executing rpwa::nBodyPhaseSpaceKinematics::setDecay()");
			bp::throw_error_already_set();
		}
		return self.setDecay(daughterMasses);
	}

	void nBodyPhaseSpaceKinematics_calcEventKinematics(rpwa::nBodyPhaseSpaceKinematics& self, PyObject* PyNBody) {
		TLorentzVector* nBody = rpwa::py::convertFromPy<TLorentzVector*>(PyNBody);
		if(not nBody) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for nBody when executing rpwa::nBodyPhaseSpace::calcEventKinematics()");
			bp::throw_error_already_set();
		}
		self.calcEventKinematics(*nBody);
	}

	PyObject* nBodyPhaseSpaceKinematics_daughter(const rpwa::nBodyPhaseSpaceKinematics& self, const int index) {
		return rpwa::py::convertToPy<TLorentzVector>(self.daughter(index));
	}

	bp::list nBodyPhaseSpaceKinematics_daughters(const rpwa::nBodyPhaseSpaceKinematics& self) {
		const std::vector<TLorentzVector>& daughters = self.daughters();
		bp::list pyDaughters;
		for(unsigned int i = 0; i < daughters.size(); ++i) {
			pyDaughters.append(bp::object(bp::handle<>(rpwa::py::convertToPy<TLorentzVector>(daughters[i]))));
		}
		return pyDaughters;
	}

}

void rpwa::py::exportNBodyPhaseSpaceKinematics() {

	bp::scope theScope = bp::class_<rpwa::nBodyPhaseSpaceKinematics>("nBodyPhaseSpaceKinematics")

		.def(bp::self_ns::str(bp::self))

		.def("setDecay", &nBodyPhaseSpaceKinematics_setDecay)

		.def("calcBreakupMomenta", &rpwa::nBodyPhaseSpaceKinematics::calcBreakupMomenta)

		.def("calcWeight", &rpwa::nBodyPhaseSpaceKinematics::calcWeight)

		.def("calcEventKinematics", &nBodyPhaseSpaceKinematics_calcEventKinematics)

		.def("setKinematicsType", &rpwa::nBodyPhaseSpaceKinematics::setKinematicsType)
		.def("kinematicsType", &rpwa::nBodyPhaseSpaceKinematics::kinematicsType)

		.def("setWeightType", &rpwa::nBodyPhaseSpaceKinematics::setWeightType)
		.def("weightType", &rpwa::nBodyPhaseSpaceKinematics::weightType)

		.def("normalization", &rpwa::nBodyPhaseSpaceKinematics::normalization)
		.def("eventWeight", &rpwa::nBodyPhaseSpaceKinematics::eventWeight)
		.def("maxWeightObserved", &rpwa::nBodyPhaseSpaceKinematics::maxWeightObserved)
		.def("resetMaxWeightObserved", &rpwa::nBodyPhaseSpaceKinematics::resetMaxWeightObserved)

		.def("daughter", &nBodyPhaseSpaceKinematics_daughter)
		.def("daughters", &nBodyPhaseSpaceKinematics_daughters)
		.def("nmbOfDaughters", &rpwa::nBodyPhaseSpaceKinematics::nmbOfDaughters)
		.def("daughterMass", &rpwa::nBodyPhaseSpaceKinematics::daughterMass)
		.def("intermediateMass", &rpwa::nBodyPhaseSpaceKinematics::intermediateMass)
		.def("breakupMom", &rpwa::nBodyPhaseSpaceKinematics::breakupMom)
		.def("cosTheta", &rpwa::nBodyPhaseSpaceKinematics::cosTheta)
		.def("phi", &rpwa::nBodyPhaseSpaceKinematics::phi);

	bp::enum_<rpwa::nBodyPhaseSpaceKinematics::kinematicsTypeEnum>("kinematicsTypeEnum")
		.value("BLOCK", rpwa::nBodyPhaseSpaceKinematics::BLOCK)
		.value("RAUBOLD_LYNCH", rpwa::nBodyPhaseSpaceKinematics::RAUBOLD_LYNCH)
		.export_values();

	bp::enum_<rpwa::nBodyPhaseSpaceKinematics::weightTypeEnum>("weightTypeEnum")
		.value("S_U_CHUNG", rpwa::nBodyPhaseSpaceKinematics::S_U_CHUNG)
		.value("NUPHAZ", rpwa::nBodyPhaseSpaceKinematics::NUPHAZ)
		.value("GENBOD", rpwa::nBodyPhaseSpaceKinematics::GENBOD)
		.value("FLAT", rpwa::nBodyPhaseSpaceKinematics::FLAT)
		.export_values();

	theScope.attr("BLOCK") = rpwa::nBodyPhaseSpaceKinematics::BLOCK;
	theScope.attr("RAUBOLD_LYNCH") = rpwa::nBodyPhaseSpaceKinematics::RAUBOLD_LYNCH;
	theScope.attr("S_U_CHUNG") = rpwa::nBodyPhaseSpaceKinematics::S_U_CHUNG;
	theScope.attr("NUPHAZ") = rpwa::nBodyPhaseSpaceKinematics::NUPHAZ;
	theScope.attr("GENBOD") = rpwa::nBodyPhaseSpaceKinematics::GENBOD;
	theScope.attr("FLAT") = rpwa::nBodyPhaseSpaceKinematics::FLAT;

}
