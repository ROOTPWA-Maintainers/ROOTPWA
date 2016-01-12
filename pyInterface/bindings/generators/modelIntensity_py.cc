#include"modelIntensity_py.h"

#include<boost/python.hpp>

#include<TVector3.h>

#include"ampIntegralMatrix.h"
#include"modelIntensity.h"
#include"rootConverters_py.h"
#include"stlContainers_py.h"

namespace bp = boost::python;


namespace {

	bool
	modelIntensity_loadIntegrals(rpwa::modelIntensity& self,
	                             PyObject*             pyIntegralMatrix)
	{
		rpwa::ampIntegralMatrix* integralMatrix = rpwa::py::convertFromPy<rpwa::ampIntegralMatrix* >(pyIntegralMatrix);
		if(not integralMatrix) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for integralMatrix when executing modelIntensity::loadIntegrals");
			bp::throw_error_already_set();
		}
		return self.loadIntegrals(*integralMatrix);
	}


	bool
	modelIntensity_setFinalStateMasses(rpwa::modelIntensity& self,
	                                   const bp::list&       pyMasses)
	{
		std::vector<double> masses;
		if(not rpwa::py::convertBPObjectToVector<double>(pyMasses, masses)) {
			PyErr_SetString(PyExc_TypeError, "invalid masses gotten");
			bp::throw_error_already_set();
		}
		return self.setFinalStateMasses(masses);
	}


	bool
	modelIntensity_setParticles(rpwa::modelIntensity& self,
	                            const bp::list&       pyInitial,
	                            const bp::list&       pyFinal)
	{
		std::vector<std::string> initial;
		if (not rpwa::py::convertBPObjectToVector<std::string>(pyInitial, initial)) {
			PyErr_SetString(PyExc_TypeError, "invalid initial particles gotten");
			bp::throw_error_already_set();
		}

		std::vector<std::string> final;
		if (not rpwa::py::convertBPObjectToVector<std::string>(pyFinal, final)) {
			PyErr_SetString(PyExc_TypeError, "invalid final particles gotten");
			bp::throw_error_already_set();
		}

		return self.setParticles(initial, final);
	}


	bool
	modelIntensity_setProdKinMomenta(rpwa::modelIntensity& self,
	                                 PyObject*             pyProdKinMomenta)
	{
		bp::extract<bp::list> getProdKinMomentaList(pyProdKinMomenta);
		if(not getProdKinMomentaList.check()) {
			printErr<<"Got invalid input for prodKinMomenta when executing rpwa::modelIntensity::setProdKinMomenta()"<<std::endl;
			return false;
		}
		bp::list prodKinMomentaList = getProdKinMomentaList();

		std::vector<TVector3> prodKinMomenta(len(prodKinMomentaList));
		for(int i = 0; i < len(prodKinMomentaList); ++i) {
			bp::object item = bp::extract<bp::object>(prodKinMomentaList[i]);
			prodKinMomenta[i] = *rpwa::py::convertFromPy<TVector3*>(item.ptr());
		}

		return self.setProdKinMomenta(prodKinMomenta);
	}


	double
	modelIntensity_getIntensity1(rpwa::modelIntensity& self,
	                             PyObject*             pyDecayKinMomenta)
	{
		bp::extract<bp::list> getDecayKinMomentaList(pyDecayKinMomenta);
		if(not getDecayKinMomentaList.check()) {
			PyErr_SetString(PyExc_TypeError,"Got invalid input for decayKinMomenta when executing rpwa::modelIntensity::getIntensity()");
			bp::throw_error_already_set();
		}
		bp::list decayKinMomentaList = getDecayKinMomentaList();

		std::vector<TVector3> decayKinMomenta(len(decayKinMomentaList));
		for(int i = 0; i < len(decayKinMomentaList); ++i) {
			bp::object item = bp::extract<bp::object>(decayKinMomentaList[i]);
			decayKinMomenta[i] = *rpwa::py::convertFromPy<TVector3*>(item.ptr());
		}

		return self.getIntensity(decayKinMomenta);
	}


	double
	modelIntensity_getIntensity2(rpwa::modelIntensity& self,
	                             PyObject*             pyProdKinMomenta,
	                             PyObject*             pyDecayKinMomenta)
	{
		bp::extract<bp::list> getProdKinMomentaList(pyProdKinMomenta);
		if(not getProdKinMomentaList.check()) {
			PyErr_SetString(PyExc_TypeError,"Got invalid input for prodKinMomenta when executing rpwa::modelIntensity::getIntensity()");
			bp::throw_error_already_set();
		}
		bp::list prodKinMomentaList = getProdKinMomentaList();

		std::vector<TVector3> prodKinMomenta(len(prodKinMomentaList));
		for (int i = 0; i < len(prodKinMomentaList); ++i) {
			bp::object item = bp::extract<bp::object>(prodKinMomentaList[i]);
			prodKinMomenta[i] = *rpwa::py::convertFromPy<TVector3*>(item.ptr());
		}

		bp::extract<bp::list> getDecayKinMomentaList(pyDecayKinMomenta);
		if(not getDecayKinMomentaList.check()) {
			PyErr_SetString(PyExc_TypeError,"Got invalid input for decayKinMomenta when executing rpwa::modelIntensity::getIntensity()");
			bp::throw_error_already_set();
		}
		bp::list decayKinMomentaList = getDecayKinMomentaList();

		std::vector<TVector3> decayKinMomenta(len(decayKinMomentaList));
		for(int i = 0; i < len(decayKinMomentaList); ++i) {
			bp::object item = bp::extract<bp::object>(decayKinMomentaList[i]);
			decayKinMomenta[i] = *rpwa::py::convertFromPy<TVector3*>(item.ptr());
		}

		return self.getIntensity(prodKinMomenta, decayKinMomenta);
	}

}


void rpwa::py::exportModelIntensity() {

	bp::class_<rpwa::modelIntensity>("modelIntensity")

		.def(
			"getIntensity"
			, &::modelIntensity_getIntensity1
			, (bp::arg("decayKinMomenta"))
		)
		.def(
			"getIntensity"
			, &::modelIntensity_getIntensity2
			, (bp::arg("prodKinMomenta"),
			   bp::arg("decayKinMomenta"))
		)

		.def(
			"loadIntegrals"
			, &::modelIntensity_loadIntegrals
			, (bp::arg("integralMatrix"))
		)
		.def(
			"loadIntegrals"
			, &rpwa::modelIntensity::loadIntegrals
			, (bp::arg("integralMatrix"))
		)

		.def(
			"addAmplitude"
			, &rpwa::modelIntensity::addAmplitude
			, (bp::arg("transitionAmplitude"),
			   bp::arg("amplitude"),
			   bp::arg("reflectivity") = 1)
		)
		.def(
			"initAmplitudes"
			, &rpwa::modelIntensity::initAmplitudes
			, (bp::arg("fromXdecay") = true)
		)

		.def(
			"setProdKinMomenta"
			, &::modelIntensity_setProdKinMomenta
			, (bp::arg("prodKinMomenta"))
		)

		.def("mBeam", &rpwa::modelIntensity::mBeam)
		.def(
			"setMbeam"
			, &rpwa::modelIntensity::setMbeam
			, (bp::arg("mBeam"))
		)

		.def("mTarget", &rpwa::modelIntensity::mTarget)
		.def(
			"setMtarget"
			, &rpwa::modelIntensity::setMtarget
			, (bp::arg("mTarget"))
		)
		.def(
			"setTarget"
			, &rpwa::modelIntensity::setTarget
			, (bp::arg("target"))
		)

		.def(
			"setParticles"
			, &::modelIntensity_setParticles
			, (bp::arg("intialState"),
			   bp::arg("finalState"))
		)
		.def(
			"setFinalStateMasses"
			, &::modelIntensity_setFinalStateMasses
			, (bp::arg("masses"))
		)

		.def("Print", &rpwa::modelIntensity::print)

	;

}
