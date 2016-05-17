#include "ampIntegralMatrix_py.h"

#include <boost/python.hpp>

#include <TDirectory.h>

#include "ampIntegralMatrix.h"
#include "amplitudeMetadata.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;


namespace {

/*
	const bp::list ampIntegralMatrix_waveDescriptions(const rpwa::ampIntegralMatrix& self) {
		return bp::list(self());
	}
*/
	const rpwa::waveDescription& ampIntegralMatrix_waveDesc1(const rpwa::ampIntegralMatrix& self,
	                                                         const unsigned int waveIndex)
	{
		return *(self.waveDesc(waveIndex));
	}

	const rpwa::waveDescription& ampIntegralMatrix_waveDesc2(const rpwa::ampIntegralMatrix& self,
	                                                         const std::string& waveName) {
		return *(self.waveDesc(waveName));
	}

	std::complex<double> ampIntegralMatrix_element1(const rpwa::ampIntegralMatrix& self,
	                                                const unsigned int waveIndexI,
	                                                const unsigned int waveIndexJ)
	{
		return self.element(waveIndexI, waveIndexJ);
	}

	std::complex<double> ampIntegralMatrix_element2(const rpwa::ampIntegralMatrix& self,
	                                                const std::string& waveNameI,
	                                                const std::string& waveNameJ)
	{
		return self.element(waveNameI, waveNameJ);
	}

	bool ampIntegralMatrix_integrate(rpwa::ampIntegralMatrix& self,
	                                 const bp::object& pyAmplitudeMetadata,
	                                 const unsigned long maxNmbEvents,
	                                 const std::string& weightFileName,
	                                 const rpwa::eventMetadata* eventMeta,
	                                 const bp::dict& pyOtfBin)
	{
		std::vector<const rpwa::amplitudeMetadata*> amplitudeMeta;
		if(not rpwa::py::convertBPObjectToVector<const rpwa::amplitudeMetadata*>(pyAmplitudeMetadata, amplitudeMeta)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for amplitudeMetadata when executing rpwa::ampIntegralMatrix::integrate()");
			bp::throw_error_already_set();
		}
		std::map<std::string, std::pair<double, double> > otfBin;
		if(bp::len(pyOtfBin.keys()) > 0) {
			bp::list keys = pyOtfBin.keys();
			for(unsigned int i = 0; i < len(keys); ++i) {
				otfBin[bp::extract<std::string>(keys[i])] =
					std::pair<double, double>(bp::extract<double>(pyOtfBin[pyOtfBin.keys()[i]][0]),
					                          bp::extract<double>(pyOtfBin[pyOtfBin.keys()[i]][1]));
			}
		}
		return self.integrate(amplitudeMeta, maxNmbEvents, weightFileName, eventMeta, otfBin);
	}

	bool ampIntegralMatrix_setWaveNames(rpwa::ampIntegralMatrix& self,
	                                    const bp::object&        waveNamesPy)
	{
		std::vector<std::string> waveNames;
		if(not rpwa::py::convertBPObjectToVector<std::string>(waveNamesPy, waveNames)) {
			printErr << "conversion of waveNames failed. Abort..." << std::endl;
			return false;
		}
		return self.setWaveNames(waveNames);
	}

	bool ampIntegralMatrix_addEvent(rpwa::ampIntegralMatrix& self,
	                                const bp::dict&          pyAmplitudes) {
		std::map<std::string, std::complex<double> > amplitudes;
		const bp::list keys = pyAmplitudes.keys();
		for (int i = 0; i < bp::len(keys); ++i) {
			std::string waveName = bp::extract<std::string>(keys[i]);
			std::complex<double> ampl = bp::extract<std::complex<double> >(pyAmplitudes[waveName]);
			amplitudes[waveName] = ampl;
		}
		return self.addEvent(amplitudes);
	}

	bool ampIntegralMatrix_writeAscii(const rpwa::ampIntegralMatrix& self, const std::string& outFileName) {
		return self.writeAscii(outFileName);
	}

	bool ampIntegralMatrix_readAscii(rpwa::ampIntegralMatrix& self, const std::string& inFileName) {
		return self.readAscii(inFileName);
	}

	int ampIntegralMatrix_Write(const rpwa::ampIntegralMatrix& self, const char* name = 0) {
		return self.Write(name);
	}

}

void rpwa::py::exportAmpIntegralMatrix() {

	bp::class_<rpwa::ampIntegralMatrix>("ampIntegralMatrix")

		.def(bp::init<const rpwa::ampIntegralMatrix&>())

		.def(bp::self_ns::str(bp::self))
		.def(bp::self == bp::self)
		.def(bp::self != bp::self)
		.def(bp::self + bp::self)
		.def(bp::self - bp::self)
		.def(bp::self * double())
		.def(bp::self / double())
		.def(double() * bp::self)
		.def(bp::self += bp::self)
		.def(bp::self -= bp::self)
		.def(bp::self *= double())
		.def(bp::self /= double())

		.def("clear", &rpwa::ampIntegralMatrix::clear)
		.def("nmbWaves", &rpwa::ampIntegralMatrix::nmbWaves)
		.def("nmbEvents", &rpwa::ampIntegralMatrix::nmbEvents)
		.def("setNmbEvents", &rpwa::ampIntegralMatrix::setNmbEvents)
		.def("containsWave", &rpwa::ampIntegralMatrix::containsWave)
		.def("waveIndex", &rpwa::ampIntegralMatrix::waveIndex)
		.def(
			"waveName"
			, &rpwa::ampIntegralMatrix::waveName
			, bp::return_value_policy<bp::return_by_value>()
		)
//		Disabled because of missing == operator in rpwa::waveDescription
//		See also http://stackoverflow.com/questions/10680691/why-do-i-need-comparison-operators-in-boost-python-vector-indexing-suite
//		.def("waveDescriptions", &ampIntegralMatrix_waveDescriptions)
		.def(
			"waveDesc"
			, &ampIntegralMatrix_waveDesc1
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def(
			"waveDesc"
			, &ampIntegralMatrix_waveDesc2
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def("allWavesHaveDesc", &rpwa::ampIntegralMatrix::allWavesHaveDesc)

//		Commenting this until it is decided how the boost::multi_array should be handled in python
//		.def("matrix", &rpwa::ampIntegralMatrix::matrix)

		.def("element", &ampIntegralMatrix_element1)
		.def("element", &ampIntegralMatrix_element2)

		.def("integrate"
		     , &ampIntegralMatrix_integrate
		     , (bp::arg("amplitudeMetadata"),
		        bp::arg("maxNmbEvents")=0,
		        bp::arg("weightFileName")="",
		        bp::arg("eventMeta")=bp::object(),
		        bp::arg("otfBin")=bp::dict())
		)
		.def("setWaveNames"
		     , &ampIntegralMatrix_setWaveNames
		     , bp::arg("waveNames"))

		.def("addEvent"
		     , &::ampIntegralMatrix_addEvent
		     , bp::arg("waveNameAmplitudeMap"))

		.def("renormalize", &rpwa::ampIntegralMatrix::renormalize)
		.def("writeAscii", &ampIntegralMatrix_writeAscii)
		.def("readAscii", &ampIntegralMatrix_readAscii)

		.def("Write", &ampIntegralMatrix_Write, bp::arg("name")=0)
		.def("getFromTDirectory", &rpwa::py::getFromTDirectory<rpwa::ampIntegralMatrix>, bp::return_value_policy<bp::manage_new_object>())
		.staticmethod("getFromTDirectory")

		.add_static_property("debugAmpIntegralMatrix", &rpwa::ampIntegralMatrix::debug, &rpwa::ampIntegralMatrix::setDebug)
		.def_readonly("integralObjectName", &rpwa::ampIntegralMatrix::integralObjectName);

}
