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
	modelIntensity_addIntegral(rpwa::modelIntensity& self,
	                           PyObject*             pyIntegralMatrix)
	{
		rpwa::ampIntegralMatrix* integralMatrix = rpwa::py::convertFromPy<rpwa::ampIntegralMatrix* >(pyIntegralMatrix);
		if(not integralMatrix) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for integralMatrix when executing modelIntensity::addIntegral");
			bp::throw_error_already_set();
		}
		return self.addIntegral(*integralMatrix);
	}


	double
	modelIntensity_initAmplitudes_1(rpwa::modelIntensity& self,
	                                const bp::list&       pyDecayKinParticleNames)
	{
		std::vector<std::string> decayKinParticleNames;
		if(not rpwa::py::convertBPObjectToVector<std::string>(pyDecayKinParticleNames, decayKinParticleNames)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for decayKinParticleNames when executing rpwa::modelIntensity::initAmplitudes()");
			bp::throw_error_already_set();
		}

		return self.initAmplitudes(decayKinParticleNames);
	}


	double
	modelIntensity_initAmplitudes_2(rpwa::modelIntensity& self,
	                                const bp::list&       pyProdKinParticleNames,
	                                const bp::list&       pyDecayKinParticleNames,
	                                const bool            fromXDecay)
	{
		std::vector<std::string> prodKinParticleNames;
		if(not rpwa::py::convertBPObjectToVector<std::string>(pyProdKinParticleNames, prodKinParticleNames)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for prodKinParticleNames when executing rpwa::modelIntensity::initAmplitudes()");
			bp::throw_error_already_set();
		}

		std::vector<std::string> decayKinParticleNames;
		if(not rpwa::py::convertBPObjectToVector<std::string>(pyDecayKinParticleNames, decayKinParticleNames)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for decayKinParticleNames when executing rpwa::modelIntensity::initAmplitudes()");
			bp::throw_error_already_set();
		}

		return self.initAmplitudes(prodKinParticleNames, decayKinParticleNames, fromXDecay);
	}


	double
	modelIntensity_getIntensity_1(rpwa::modelIntensity& self,
	                              const bp::list&       pyDecayKinMomenta)
	{
		std::vector<TVector3> decayKinMomenta(len(pyDecayKinMomenta));
		for(unsigned int i = 0; i < len(pyDecayKinMomenta); ++i) {
			bp::object item = bp::extract<bp::object>(pyDecayKinMomenta[i]);
			decayKinMomenta[i] = *rpwa::py::convertFromPy<TVector3*>(item.ptr());
		}

		return self.getIntensity(decayKinMomenta);
	}


	double
	modelIntensity_getIntensity_2(rpwa::modelIntensity& self,
	                              const bp::list&       pyList,
	                              const bp::list&       pyDecayKinMomenta)
	{
		bool                      indices = false;
		std::vector<unsigned int> waveIndices   (len(pyList));
		bool                      momenta = false;
		std::vector<TVector3>     prodKinMomenta(len(pyList));
		for(unsigned int i = 0; i < len(pyList); ++i) {
			bp::extract<unsigned int> itemU(pyList[i]);
			if (itemU.check()) {
				indices = true;
				if (momenta) {
					PyErr_SetString(PyExc_TypeError, "First argument of modelIntensity::getIntensity is a mixed list of 'unsigned int' and 'TVector3'.");
					bp::throw_error_already_set();
				}
				waveIndices[i] = itemU();
				continue;
			}
			bp::extract<bp::object> itemV(pyList[i]);
			if (itemV.check()) {
				TVector3* ptrV = rpwa::py::convertFromPy<TVector3*>(itemV().ptr());
				if (ptrV) {
					momenta = true;
					if (indices) {
						PyErr_SetString(PyExc_TypeError, "First argument of modelIntensity::getIntensity is a mixed list of 'unsigned int' and 'TVector3'.");
						bp::throw_error_already_set();
					}
					prodKinMomenta[i] = *ptrV;
					continue;
				}
			}
			PyErr_SetString(PyExc_TypeError, "First argument of modelIntensity::getIntensity is not a list of either 'unsigned int' or 'TVector3'.");
			bp::throw_error_already_set();
		}

		std::vector<TVector3> decayKinMomenta(len(pyDecayKinMomenta));
		for(unsigned int i = 0; i < len(pyDecayKinMomenta); ++i) {
			bp::object item = bp::extract<bp::object>(pyDecayKinMomenta[i]);
			decayKinMomenta[i] = *rpwa::py::convertFromPy<TVector3*>(item.ptr());
		}

		if (indices)
			return self.getIntensity(waveIndices, decayKinMomenta);
		if (momenta)
			return self.getIntensity(prodKinMomenta, decayKinMomenta);

		PyErr_SetString(PyExc_TypeError, "Could not determine type of first argument of modelIntensity::getIntensity.");
		bp::throw_error_already_set();
		return 0;
	}


	double
	modelIntensity_getIntensity_3(rpwa::modelIntensity& self,
	                              const bp::list&       pyWaveIndices,
	                              const bp::list&       pyProdKinMomenta,
	                              const bp::list&       pyDecayKinMomenta)
	{
		std::vector<unsigned int> waveIndices;
		if (not rpwa::py::convertBPObjectToVector<unsigned int>(pyWaveIndices, waveIndices)) {
			PyErr_SetString(PyExc_TypeError, "Cannot convert first argument of modelIntensity::getIntensity to a vector of 'unsigned int'.");
			bp::throw_error_already_set();
		}

		std::vector<TVector3> prodKinMomenta(len(pyProdKinMomenta));
		for(unsigned int i = 0; i < len(pyProdKinMomenta); ++i) {
			bp::object item = bp::extract<bp::object>(pyProdKinMomenta[i]);
			prodKinMomenta[i] = *rpwa::py::convertFromPy<TVector3*>(item.ptr());
		}

		std::vector<TVector3> decayKinMomenta(len(pyDecayKinMomenta));
		for(unsigned int i = 0; i < len(pyDecayKinMomenta); ++i) {
			bp::object item = bp::extract<bp::object>(pyDecayKinMomenta[i]);
			decayKinMomenta[i] = *rpwa::py::convertFromPy<TVector3*>(item.ptr());
		}

		return self.getIntensity(waveIndices, prodKinMomenta, decayKinMomenta);
	}

}


void rpwa::py::exportModelIntensity() {

	bp::class_<rpwa::modelIntensity>("modelIntensity", bp::init<fitResultPtr>())

		.def(bp::self_ns::str(bp::self))

		.def(
			"addAmplitude"
			, &rpwa::modelIntensity::addAmplitude
			, (bp::arg("amplitude"))
		)

		.def(
			"addIntegral"
			, &::modelIntensity_addIntegral
			, (bp::arg("integralMatrix"))
		)
		.def(
			"addIntegral"
			, &rpwa::modelIntensity::addIntegral
			, (bp::arg("integralMatrix"))
		)

		.def(
			"initAmplitudes"
			, &modelIntensity_initAmplitudes_1
			, (bp::arg("decayKinParticleNames"))
		)
		.def(
			"initAmplitudes"
			, &modelIntensity_initAmplitudes_2
			, (bp::arg("prodKinParticleNames"),
			   bp::arg("decayKinParticleNames"),
			   bp::arg("fromXdecay") = false)
		)

		.def(
			"getIntensity"
			, &modelIntensity_getIntensity_1
			, (bp::arg("decayKinMomenta"))
		)
		.def(
			"getIntensity"
			, &modelIntensity_getIntensity_2
			, (bp::arg("prodKinMomenta"),
			   bp::arg("decayKinMomenta"))
		)
		.def(
			"getIntensity"
			, &modelIntensity_getIntensity_2
			, (bp::arg("waveIndices"),
			   bp::arg("decayKinMomenta"))
		)
		.def(
			"getIntensity"
			, &modelIntensity_getIntensity_3
			, (bp::arg("waveIndices"),
			   bp::arg("prodKinMomenta"),
			   bp::arg("decayKinMomenta"))
		)

	;

	bp::register_ptr_to_python<rpwa::modelIntensityPtr>();

}
