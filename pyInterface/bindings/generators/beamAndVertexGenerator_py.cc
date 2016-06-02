#include "beamAndVertexGenerator_py.h"

#include <boost/python.hpp>

#include "beamAndVertexGenerator.h"
#include "generatorParameters.hpp"
#include "rootConverters_py.h"

namespace bp = boost::python;


namespace {

	struct beamAndVertexGeneratorWrapper : rpwa::beamAndVertexGenerator,
	                                       bp::wrapper<rpwa::beamAndVertexGenerator>
	{

		beamAndVertexGeneratorWrapper()
			: rpwa::beamAndVertexGenerator(),
			  bp::wrapper<rpwa::beamAndVertexGenerator>() { }

		beamAndVertexGeneratorWrapper(const rpwa::beamAndVertexGenerator& gen)
			: rpwa::beamAndVertexGenerator(gen),
			  bp::wrapper<rpwa::beamAndVertexGenerator>() { }

		bool loadBeamFile(const std::string& beamFileName) {
			if(bp::override loadBeamFile = this->get_override("loadBeamFile")) {
				return loadBeamFile(beamFileName);
			}
			return rpwa::beamAndVertexGenerator::loadBeamFile(beamFileName);
		}

		bool default_loadBeamFile(const std::string& beamFileName) {
			return rpwa::beamAndVertexGenerator::loadBeamFile(beamFileName);
		}

		void setBeamfileSequentialReading(bool sequentialReading = true)
		{
			if(bp::override setBeamfileSequentialReading = this->get_override("setBeamfileSequentialReading")) {
				setBeamfileSequentialReading(sequentialReading);
			} else {
				rpwa::beamAndVertexGenerator::setBeamfileSequentialReading(sequentialReading);
			}
		}

		void default_setBeamfileSequentialReading(bool sequentialReading = true)
		{
			rpwa::beamAndVertexGenerator::setBeamfileSequentialReading(sequentialReading);
		}

		void randomizeBeamfileStartingPosition()
		{
			if(bp::override randomizeBeamfileStartingPosition = this->get_override("randomizeBeamfileStartingPosition")) {
				randomizeBeamfileStartingPosition();
			} else {
				rpwa::beamAndVertexGenerator::randomizeBeamfileStartingPosition();
			}
		}

		void default_randomizeBeamfileStartingPosition()
		{
			rpwa::beamAndVertexGenerator::randomizeBeamfileStartingPosition();
		}

		bool check() const {
			if(bp::override check = this->get_override("check")) {
				return check();
			}
			return rpwa::beamAndVertexGenerator::check();
		}

		bool default_check() const {
			return rpwa::beamAndVertexGenerator::check();
		}

		bool event(const rpwa::Target& target, const rpwa::Beam& beam) {
			if(bp::override event = this->get_override("event")) {
				return event(target, beam);
			}
			return rpwa::beamAndVertexGenerator::event(target, beam);
		}

		bool default_event(const rpwa::Target& target, const rpwa::Beam& beam) {
			return rpwa::beamAndVertexGenerator::event(target, beam);
		}

		PyObject* getVertex__() const {
			TVector3 retval;
			if(bp::override getVertex = this->get_override("getVertex")) {
				retval = getVertex();
			} else {
				retval = rpwa::beamAndVertexGenerator::getVertex();
			}
			return rpwa::py::convertToPy<TVector3>(retval);
		}

		PyObject* default_getVertex__() const {
			return rpwa::py::convertToPy<TVector3>(rpwa::beamAndVertexGenerator::getVertex());
		}

		PyObject* getBeam__() const {
			TLorentzVector retval;
			if(bp::override getBeam = this->get_override("getBeam")) {
				retval = getBeam();
			} else {
				retval = rpwa::beamAndVertexGenerator::getBeam();
			}
			return rpwa::py::convertToPy<TLorentzVector>(retval);
		}

		PyObject* default_getBeam__() const {
			return rpwa::py::convertToPy<TLorentzVector>(rpwa::beamAndVertexGenerator::getBeam());
		}

		void setSigmaScalingFactor(const double& scalingFactor)
		{
			if(bp::override setSigmaScalingFactor = this->get_override("setSigmaScalingFactor")) {
				setSigmaScalingFactor(scalingFactor);
			} else {
				rpwa::beamAndVertexGenerator::setSigmaScalingFactor(scalingFactor);
			}
		}

		void default_setSigmaScalingFactor(const double& scalingFactor)
		{
			rpwa::beamAndVertexGenerator::setSigmaScalingFactor(scalingFactor);
		}

	};

	PyObject* beamAndVertexGenerator_getVertex(const rpwa::beamAndVertexGenerator& self) {
		return rpwa::py::convertToPy<TVector3>(self.getVertex());
	}

	PyObject* beamAndVertexGenerator_getBeam(const rpwa::beamAndVertexGenerator& self) {
		return rpwa::py::convertToPy<TLorentzVector>(self.getBeam());
	}

}

void rpwa::py::exportBeamAndVertexGenerator() {

	bp::class_<beamAndVertexGeneratorWrapper>("beamAndVertexGenerator")
		.def(bp::self_ns::str(bp::self))
		.def(
			"loadBeamFile"
			, &beamAndVertexGeneratorWrapper::loadBeamFile
			, &beamAndVertexGeneratorWrapper::default_loadBeamFile
		)
		.def("loadBeamFile", &rpwa::beamAndVertexGenerator::loadBeamFile)
		.def(
			"setBeamfileSequentialReading"
			, &beamAndVertexGeneratorWrapper::setBeamfileSequentialReading
			, &beamAndVertexGeneratorWrapper::default_setBeamfileSequentialReading
			, bp::arg("sequentialReading")=true
		)
		.def("setBeamfileSequentialReading", &rpwa::beamAndVertexGenerator::setBeamfileSequentialReading)
		.def(
			"randomizeBeamfileStartingPosition"
			, &beamAndVertexGeneratorWrapper::randomizeBeamfileStartingPosition
			, &beamAndVertexGeneratorWrapper::default_randomizeBeamfileStartingPosition
		)
		.def("randomizeBeamfileStartingPosition", &rpwa::beamAndVertexGenerator::randomizeBeamfileStartingPosition)
		.def("check", &beamAndVertexGeneratorWrapper::check, &beamAndVertexGeneratorWrapper::default_check)
		.def("check", &rpwa::beamAndVertexGenerator::check)
		.def("event", &beamAndVertexGeneratorWrapper::event, &beamAndVertexGeneratorWrapper::default_event)
		.def("event", &rpwa::beamAndVertexGenerator::event)
		.def(
			"getVertex"
			, &beamAndVertexGeneratorWrapper::getVertex__
			, &beamAndVertexGeneratorWrapper::default_getVertex__
		)
		.def("getVertex", &beamAndVertexGenerator_getVertex)
		.def(
			"getBeam"
			, &beamAndVertexGeneratorWrapper::getBeam__
			, &beamAndVertexGeneratorWrapper::default_getBeam__
		)
		.def("getBeam", &beamAndVertexGenerator_getBeam)
		.def(
			"setSigmaScalingFactor"
			, &beamAndVertexGeneratorWrapper::setSigmaScalingFactor
			, &beamAndVertexGeneratorWrapper::default_setSigmaScalingFactor
		)
		.def("setSigmaScalingFactor", &rpwa::beamAndVertexGenerator::setSigmaScalingFactor);


}
