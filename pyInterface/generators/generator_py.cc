#include "generator_py.h"

#include <boost/python.hpp>

#include "generator.h"
#include "generatorPickerFunctions.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;


namespace {

	struct generatorWrapper : public rpwa::generator,
	                                 bp::wrapper<rpwa::generator>
	{

		generatorWrapper()
			: rpwa::generator(),
			  bp::wrapper<rpwa::generator>() { }

		const rpwa::particle& getGeneratedBeam() const
		{
			if(bp::override getGeneratedBeam = this->get_override("getGeneratedBeam")) {
				return getGeneratedBeam();
			}
			return rpwa::generator::getGeneratedBeam();
		}

		const rpwa::particle& default_getGeneratedBeam() const
		{
			return rpwa::generator::getGeneratedBeam();
		}

		const rpwa::particle& getGeneratedRecoil() const
		{
			if(bp::override getGeneratedRecoil = this->get_override("getGeneratedRecoil")) {
				return getGeneratedRecoil();
			}
			return rpwa::generator::getGeneratedRecoil();
		}

		const rpwa::particle& default_getGeneratedRecoil() const
		{
			return rpwa::generator::getGeneratedRecoil();
		}

		bp::list getGeneratedFinalState__() const
		{
			if(bp::override getGeneratedFinalState = this->get_override("getGeneratedFinalState")) {
				return bp::list(getGeneratedFinalState());
			}
			return bp::list(rpwa::generator::getGeneratedFinalState());
		}

		bp::list default_getGeneratedFinalState__() const
		{
			return bp::list(rpwa::generator::getGeneratedFinalState());
		}

		PyObject* getGeneratedVertex__() const
		{
			TVector3 retval;
			if(bp::override getGeneratedVertex = this->get_override("getGeneratedVertex")) {
				retval = getGeneratedVertex();
			} else {
				retval = rpwa::generator::getGeneratedVertex();
			}
			return rpwa::py::convertToPy<TVector3>(retval);
		}

		PyObject* default_getGeneratedVertex__() const
		{
			return rpwa::py::convertToPy<TVector3>(rpwa::generator::getGeneratedVertex());
		}

		void setBeam(const rpwa::Beam& beam) {
			if(bp::override setBeam = this->get_override("setBeam")) {
				setBeam(beam);
			} else {
				rpwa::generator::setBeam(beam);
			}
		}

		void default_setBeam(const rpwa::Beam& beam) {
			rpwa::generator::setBeam(beam);
		}

		void setTarget(const rpwa::Target& target) {
			if(bp::override setTarget = this->get_override("setTarget")) {
				setTarget(target);
			} else {
				rpwa::generator::setTarget(target);
			}
		}

		void default_setTarget(const rpwa::Target& target) {
			rpwa::generator::setTarget(target);
		}

		void setTPrimeAndMassPicker(const rpwa::massAndTPrimePickerPtr& pickerFunction) {
			if(bp::override setTPrimeAndMassPicker = this->get_override("setTPrimeAndMassPicker")) {
				setTPrimeAndMassPicker(pickerFunction);
			} else {
				rpwa::generator::setTPrimeAndMassPicker(pickerFunction);
			}
		}

		void default_setTPrimeAndMassPicker(const rpwa::massAndTPrimePickerPtr& pickerFunction) {
			rpwa::generator::setTPrimeAndMassPicker(pickerFunction);
		}

		const rpwa::massAndTPrimePickerPtr& getTPrimeAndMassPicker() {
			if(bp::override getTPrimeAndMassPicker = this->get_override("getTPrimeAndMassPicker")) {
				return getTPrimeAndMassPicker();
			}
			return rpwa::generator::getTPrimeAndMassPicker();
		}

		const rpwa::massAndTPrimePickerPtr& default_getTPrimeAndMassPicker() {
			return rpwa::generator::getTPrimeAndMassPicker();
		}

		void setPrimaryVertexGenerator(const rpwa::beamAndVertexGeneratorPtr& beamAndVertexGenerator) {
			if(bp::override setPrimaryVertexGenerator = this->get_override("setPrimaryVertexGenerator")) {
				setPrimaryVertexGenerator(beamAndVertexGenerator);
			} else {
				rpwa::generator::setPrimaryVertexGenerator(beamAndVertexGenerator);
			}
		}

		void default_setPrimaryVertexGenerator(const rpwa::beamAndVertexGeneratorPtr& beamAndVertexGenerator) {
			rpwa::generator::setPrimaryVertexGenerator(beamAndVertexGenerator);
		}

		void setDecayProducts__(bp::object& pyParticles) {
			std::vector<rpwa::particle> particles;
			if(not rpwa::py::convertBPObjectToVector<rpwa::particle>(pyParticles, particles)) {
				PyErr_SetString(PyExc_TypeError, "Got invalid input for particles when executing rpwa::generator::setDecayProducts()");
				bp::throw_error_already_set();
			}
			if(bp::override setDecayProducts = this->get_override("setDecayProducts")) {
				setDecayProducts(particles);
			} else {
				rpwa::generator::setDecayProducts(particles);
			}
		}

		void default_setDecayProducts__(bp::object& pyParticles) {
			std::vector<rpwa::particle> particles;
			if(not rpwa::py::convertBPObjectToVector<rpwa::particle>(pyParticles, particles)) {
				PyErr_SetString(PyExc_TypeError, "Got invalid input for particles when executing rpwa::generator::setDecayProducts()");
				bp::throw_error_already_set();
			}
			rpwa::generator::setDecayProducts(particles);
		}

		void addDecayProduct(const rpwa::particle& particle) {
			if(bp::override addDecayProduct = this->get_override("addDecayProduct")) {
				addDecayProduct(particle);
			} else {
				rpwa::generator::addDecayProduct(particle);
			}
		}

		void default_addDecayProduct(const rpwa::particle& particle) {
			rpwa::generator::addDecayProduct(particle);
		}

	};

	bp::list generator_getGeneratedFinalState(const rpwa::generator& self) {
		return bp::list(self.getGeneratedFinalState());
	}

	PyObject* generator_getGeneratedVertex(const rpwa::generator& self) {
		return rpwa::py::convertToPy<TVector3>(self.getGeneratedVertex());
	}

	void generator_setDecayProducts(rpwa::generator& self, bp::object& pyParticles) {
		std::vector<rpwa::particle> particles;
		if(not rpwa::py::convertBPObjectToVector<rpwa::particle>(pyParticles, particles)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for particles when executing rpwa::generator::setDecayProducts()");
			bp::throw_error_already_set();
		}
		self.setDecayProducts(particles);
	}

	bp::str generator_convertEventToAscii(const rpwa::particle& beam,
	                                      const bp::object& pyFinalState)
	{
		std::vector<rpwa::particle> finalState;
		if(not rpwa::py::convertBPObjectToVector<rpwa::particle>(pyFinalState, finalState)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for finalState when executing rpwa::generator::convertEventToAscii()");
			bp::throw_error_already_set();
		}
		std::stringstream sstr;
		rpwa::generator::convertEventToAscii(sstr, beam, finalState);
		return bp::str(sstr.str());
	}

	bp::str generator_convertEventToComgeant(const rpwa::particle& beam,
	                                         const rpwa::particle& recoil,
	                                         PyObject* pyVertex,
	                                         const bp::object& pyFinalState,
	                                         bool writeBinary = false)
	{
		TVector3* vertex = rpwa::py::convertFromPy<TVector3*>(pyVertex);
		if(not vertex) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input when executing rpwa::generator::convertEventToAscii()");
			bp::throw_error_already_set();
		}
		std::vector<rpwa::particle> finalState;
		if(not rpwa::py::convertBPObjectToVector<rpwa::particle>(pyFinalState, finalState)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for finalState when executing rpwa::generator::convertEventToComgeant()");
			bp::throw_error_already_set();
		}
		std::stringstream sstr;
		rpwa::generator::convertEventToComgeant(sstr, beam, recoil, *vertex, finalState, writeBinary);
		return bp::str(sstr.str());
	}

}

void rpwa::py::exportGenerator() {

	bp::class_<generatorWrapper, boost::noncopyable>("generator", bp::no_init)
		.def("event", bp::pure_virtual(&rpwa::generator::event))
		.def(
			"getGeneratedBeam"
			, &generatorWrapper::getGeneratedBeam
			, &generatorWrapper::default_getGeneratedBeam
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def("getGeneratedBeam"
			, &rpwa::generator::getGeneratedBeam
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def(
			"getGeneratedRecoil"
			, &generatorWrapper::getGeneratedRecoil
			, &generatorWrapper::default_getGeneratedRecoil
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def("getGeneratedRecoil"
			, &rpwa::generator::getGeneratedRecoil
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def("getGeneratedFinalState"
			, &generatorWrapper::getGeneratedFinalState__
			, &generatorWrapper::default_getGeneratedFinalState__
		)
		.def("getGeneratedFinalState", &generator_getGeneratedFinalState)
		.def("getGeneratedVertex"
			, &generatorWrapper::getGeneratedVertex__
			, &generatorWrapper::default_getGeneratedVertex__
		)
		.def("getGeneratedVertex", &generator_getGeneratedVertex)
		.def("setBeam", &generatorWrapper::setBeam, &generatorWrapper::default_setBeam)
		.def("setBeam", &rpwa::generator::setBeam)
		.def("setTarget", &generatorWrapper::setTarget, &generatorWrapper::default_setTarget)
		.def("setTarget", &rpwa::generator::setTarget)
		.def(
			"setTPrimeAndMassPicker"
			, &generatorWrapper::setTPrimeAndMassPicker
			, &generatorWrapper::default_setTPrimeAndMassPicker
		)
		.def("setTPrimeAndMassPicker", &rpwa::generator::setTPrimeAndMassPicker)
		.def("getTPrimeAndMassPicker"
			, &generatorWrapper::getTPrimeAndMassPicker
			, &generatorWrapper::default_getTPrimeAndMassPicker
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def("getTPrimeAndMassPicker"
			, &rpwa::generator::getTPrimeAndMassPicker
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def(
			"setPrimaryVertexGenerator"
			, &generatorWrapper::setPrimaryVertexGenerator
			, &generatorWrapper::default_setPrimaryVertexGenerator
		)
		.def("setPrimaryVertexGenerator", &rpwa::generator::setPrimaryVertexGenerator)
		.def("setDecayProducts", &generatorWrapper::setDecayProducts__, &generatorWrapper::default_setDecayProducts__)
		.def("setDecayProducts", &generator_setDecayProducts)
		.def("addDecayProduct", &generatorWrapper::addDecayProduct, &generatorWrapper::default_addDecayProduct)
		.def("addDecayProduct", &rpwa::generator::addDecayProduct)
		.def("convertEventToAscii", &generator_convertEventToAscii)
		.staticmethod("convertEventToAscii")
		.def(
			"convertEventToComgeant"
			, &generator_convertEventToComgeant
			, (bp::arg("beam"),
			   bp::arg("recoil"),
			   bp::arg("vertex"),
			   bp::arg("finalState"),
			   bp::arg("writeBinary")=false)
		)
		.staticmethod("convertEventToComgeant");

}
